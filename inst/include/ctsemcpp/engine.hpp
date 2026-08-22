#ifndef CTSEMCPP_ENGINE_HPP
#define CTSEMCPP_ENGINE_HPP

// The multi-subject objective: value, adjoint gradient, and an L-BFGS driver.
//
// Each subject is filtered forward once (with tracing on when a gradient is
// wanted), then its tape is replayed backwards and its contribution accumulated
// into the shared gradient. The cost of a gradient is therefore independent of
// the number of free parameters, which is the whole point of having an adjoint:
// forward-mode needs one pass per chunk of parameters, and `indvarying` models
// inflate that count without bound.

#include <RcppEigen.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "filter.hpp"
#include "model.hpp"
#include "reverse.hpp"

namespace ctsemcpp {

struct CppObjective {
  CppModel model;
  std::vector<SubjectData> subjects;

  // Owned copies of the data, so the objective is safe to hold across R calls
  // independently of the R objects it was built from.
  std::vector<double> y, tdpreds, tipreds, times;
  std::vector<int> starts, stops;

  FilterWorkspace fws;
  AdjointWorkspace aws;
  AdjointTape tape;

  // Prior specification, as index/scale pairs over the raw parameter vector.
  // Deliberately data rather than logic: which raw parameter carries which
  // prior is ctsem semantics the R side already knows, and all the engine has
  // to do is evaluate a normal log-density and its derivative.
  std::vector<int> priorIndex;      // 0-based
  std::vector<double> priorScale;
  double priorWeight = 1.0;

  int nvalues() const { return model.nvalues; }

  // Matches the generated Stan model term for term: `normal_lpdf(x/scale|0,1)`,
  // i.e. -x^2/(2 scale^2) - log(2pi)/2. The missing -log(scale) is missing in
  // Stan too -- it takes the density of the scaled quantity -- and while that
  // is a constant, dropping it here would leave this engine's log probability a
  // constant away from Stan's, which the parity tests would report as a
  // mismatch.
  double logPrior(const double* values) const {
    double total = 0.0;
    for (std::size_t k = 0; k < priorIndex.size(); ++k) {
      const double scaled = values[priorIndex[k]] / priorScale[k];
      total += -0.5 * scaled * scaled - 0.5 * kLog2Pi;
    }
    return priorWeight * total;
  }

  void addLogPriorGradient(const double* values, double* gradient_out,
                           double weight = 1.0) const {
    const double scaled_weight = priorWeight * weight;
    for (std::size_t k = 0; k < priorIndex.size(); ++k) {
      const int idx = priorIndex[k];
      const double scale = priorScale[k];
      gradient_out[idx] -= scaled_weight * values[idx] / (scale * scale);
    }
  }

  void prepare() {
    fws.resize(model);
    aws.resize(model);
  }

  double value(const double* values) {
    double total = 0.0;
    // The discretization caches are not invalidated between subjects or
    // between evaluations: every entry is guarded by an exact comparison
    // against the input matrices that produced it, so a parameter change
    // misses and a genuinely identical input hits. Sharing one workspace
    // across subjects therefore also shares the linear-model exponential and
    // Lyapunov solutions, which are the same at every row of every subject.
    for (SubjectData& s : subjects) {
      const double ll = filterSubject(model, fws, values, s, nullptr);
      if (!std::isfinite(ll)) return std::nan("");
      total += ll;
    }
    return total + logPrior(values);
  }

  // Returns the value; `gradient` is overwritten. An invalid trial point
  // propagates as a non-finite value *and* a non-finite gradient, matching the
  // Julia backend and what an optimizer's finiteness guards expect. There is
  // deliberately no silent fallback to a different gradient method.
  double gradient(const double* values, double* gradient_out) {
    const int p = model.nvalues;
    std::fill(gradient_out, gradient_out + p, 0.0);
    // A previous call that bailed out on an invalid trial point can leave a
    // queued Frechet direction behind; it belongs to that call's cotangent.
    aws.frechetPending = false;
    aws.jaxBarDeferred.setZero();
    const bool shared = model.parameterLayerShareable();
    if (shared) std::fill(aws.theta_bar.begin(), aws.theta_bar.end(), 0.0);

    double total = 0.0;
    const std::vector<double>* lastSubjectValues = nullptr;
    const double* lastTipreds = nullptr;

    for (SubjectData& s : subjects) {
      tape.reset();
      const double ll = filterSubject(model, fws, values, s, &tape);
      if (!std::isfinite(ll)) {
        std::fill(gradient_out, gradient_out + p, std::nan(""));
        return ll;
      }
      total += ll;
      if (!shared) std::fill(aws.theta_bar.begin(), aws.theta_bar.end(), 0.0);
      reverseTape(model, aws, tape);
      lastSubjectValues = &tape.subject_values;
      lastTipreds = s.tipreds;
      if (!shared) parameterLayer(model, aws, tape.subject_values, s.tipreds, gradient_out);
    }
    if (shared && lastSubjectValues != nullptr) {
      parameterLayer(model, aws, *lastSubjectValues, lastTipreds, gradient_out);
    }
    // The prior is a closed-form function of the raw parameters alone, so it is
    // added once rather than per subject.
    total += logPrior(values);
    addLogPriorGradient(values, gradient_out);

    // The deferred matrix-exponential Frechet contribution, pushed through the
    // parameter layer in one extra pass. It only ever touches JAx cells, and
    // `frechetDeferrable()` guarantees every subject shares the same parameter
    // layer for those, so one pass covers all of them.
    if (aws.deferFrechet && lastSubjectValues != nullptr) {
      flushFrechet(model, aws);
      std::fill(aws.theta_bar.begin(), aws.theta_bar.end(), 0.0);
      const int n = model.nlatent;
      for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
          aws.theta_bar[model.offJAx + j * n + i] = aws.jaxBarDeferred(i, j);
      parameterLayer(model, aws, *lastSubjectValues, lastTipreds, gradient_out);
    }
    return total;
  }

  // Per-subject gradient contributions -- the score matrix `scorecalc()`
  // produces for the Stan backend, and what the OPG, sandwich and
  // score-bootstrap uncertainty methods consume. `scores` is nsubjects x
  // nvalues, row-major.
  //
  // The summed gradient takes two shortcuts that have to be switched off for
  // the rows to be individually correct: the shared parameter layer, which
  // unwinds the transform layer once for all subjects, and the deferred
  // Frechet contribution, which is batched across them. Both are flags.
  double subjectGradients(const double* values, double* scores) {
    const int p = model.nvalues;
    const int nsubject = static_cast<int>(subjects.size());
    std::fill(scores, scores + static_cast<std::size_t>(nsubject) * p, 0.0);
    const bool defer = aws.deferFrechet;
    aws.deferFrechet = false;
    double total = 0.0;
    for (int s = 0; s < nsubject; ++s) {
      aws.frechetPending = false;
      tape.reset();
      const double ll = filterSubject(model, fws, values, subjects[s], &tape);
      if (!std::isfinite(ll)) {
        aws.deferFrechet = defer;
        std::fill(scores, scores + static_cast<std::size_t>(nsubject) * p,
                  std::nan(""));
        return ll;
      }
      total += ll;
      std::fill(aws.theta_bar.begin(), aws.theta_bar.end(), 0.0);
      reverseTape(model, aws, tape);
      parameterLayer(model, aws, tape.subject_values, subjects[s].tipreds,
                     scores + static_cast<std::size_t>(s) * p);
    }
    aws.deferFrechet = defer;

    // Each subject's row carries 1/nsubjects of the prior, matching what
    // `scorecalc()` does on the Stan side (it sets priormod = 1/nsubjects
    // before taking per-subject gradients). The rows then still sum to the
    // full posterior gradient.
    if (!priorIndex.empty() && nsubject > 0) {
      const double share = 1.0 / nsubject;
      for (int s = 0; s < nsubject; ++s) {
        addLogPriorGradient(values, scores + static_cast<std::size_t>(s) * p, share);
      }
    }
    return total + logPrior(values);
  }
};

// ---------------------------------------------------------------------------
// L-BFGS
// ---------------------------------------------------------------------------
//
// A compact two-loop-recursion L-BFGS with a backtracking Armijo line search,
// driving the negative log-likelihood. The Julia backend calls Optim.LBFGS();
// this is the same algorithm with a simpler line search, kept in C++ so a fit
// does not pay an R round-trip per iteration. Like the Julia driver it has no
// multi-start or restart robustness -- see the note in HANDOFF.md.

struct LbfgsResult {
  std::vector<double> minimizer;
  double value = 0.0;
  std::vector<double> gradient;
  int iterations = 0;
  bool converged = false;
};

template <typename Fn>
inline LbfgsResult lbfgs(Fn&& fg, const std::vector<double>& start, int maxiter,
                         double gtol, int history = 10) {
  const int p = static_cast<int>(start.size());
  std::vector<double> x = start, g(p), xNew(p), gNew(p), dir(p);
  std::vector<std::vector<double>> sHist, yHist;
  std::vector<double> rhoHist;

  double f = fg(x.data(), g.data());
  LbfgsResult out;
  out.minimizer = x;
  out.value = f;
  out.gradient = g;

  auto infnorm = [](const std::vector<double>& v) {
    double best = 0.0;
    for (double e : v) best = std::max(best, std::fabs(e));
    return best;
  };

  if (!std::isfinite(f)) return out;

  for (int iter = 0; iter < maxiter; ++iter) {
    if (infnorm(g) < gtol) { out.converged = true; break; }

    // two-loop recursion
    dir = g;
    const int mHist = static_cast<int>(sHist.size());
    std::vector<double> alpha(mHist, 0.0);
    for (int i = mHist - 1; i >= 0; --i) {
      double dot = 0.0;
      for (int j = 0; j < p; ++j) dot += sHist[i][j] * dir[j];
      alpha[i] = rhoHist[i] * dot;
      for (int j = 0; j < p; ++j) dir[j] -= alpha[i] * yHist[i][j];
    }
    if (mHist > 0) {
      double ys = 0.0, yy = 0.0;
      for (int j = 0; j < p; ++j) {
        ys += sHist[mHist - 1][j] * yHist[mHist - 1][j];
        yy += yHist[mHist - 1][j] * yHist[mHist - 1][j];
      }
      const double scale = (yy > 0.0) ? ys / yy : 1.0;
      for (int j = 0; j < p; ++j) dir[j] *= scale;
    }
    for (int i = 0; i < mHist; ++i) {
      double dot = 0.0;
      for (int j = 0; j < p; ++j) dot += yHist[i][j] * dir[j];
      const double beta = rhoHist[i] * dot;
      for (int j = 0; j < p; ++j) dir[j] += sHist[i][j] * (alpha[i] - beta);
    }
    for (int j = 0; j < p; ++j) dir[j] = -dir[j];

    double slope = 0.0;
    for (int j = 0; j < p; ++j) slope += dir[j] * g[j];
    if (slope >= 0.0) {  // a bad curvature estimate; restart from steepest descent
      for (int j = 0; j < p; ++j) dir[j] = -g[j];
      slope = 0.0;
      for (int j = 0; j < p; ++j) slope += dir[j] * g[j];
      sHist.clear(); yHist.clear(); rhoHist.clear();
    }

    // backtracking line search with the Armijo condition; an invalid trial
    // point (non-finite objective) is treated as a failed step, not accepted.
    double step = 1.0;
    bool ok = false;
    double fNew = f;
    for (int ls = 0; ls < 40; ++ls) {
      for (int j = 0; j < p; ++j) xNew[j] = x[j] + step * dir[j];
      fNew = fg(xNew.data(), gNew.data());
      if (std::isfinite(fNew) && fNew <= f + 1e-4 * step * slope) { ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;

    std::vector<double> s(p), y(p);
    double ys = 0.0;
    for (int j = 0; j < p; ++j) {
      s[j] = xNew[j] - x[j];
      y[j] = gNew[j] - g[j];
      ys += s[j] * y[j];
    }
    if (ys > 1e-12) {
      sHist.push_back(s);
      yHist.push_back(y);
      rhoHist.push_back(1.0 / ys);
      if (static_cast<int>(sHist.size()) > history) {
        sHist.erase(sHist.begin());
        yHist.erase(yHist.begin());
        rhoHist.erase(rhoHist.begin());
      }
    }
    x = xNew;
    g = gNew;
    f = fNew;
    out.iterations = iter + 1;
    if (infnorm(g) < gtol) { out.converged = true; break; }
  }

  out.minimizer = x;
  out.value = f;
  out.gradient = g;
  return out;
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_ENGINE_HPP
