#ifndef CTSEMCPP_SUMMARY_HPP
#define CTSEMCPP_SUMMARY_HPP

// Model-implied parameter matrices for a given raw parameter vector.
//
// This is the primitive the whole summary/plot architecture rests on, and it
// deliberately runs *inside the engine* rather than re-deriving the transforms
// in R. The engine already materializes every model matrix from the raw vector
// on its way to a log likelihood; asking it for that same materialization is
// the only way to guarantee that what a summary reports is what the likelihood
// used. A parallel R implementation of the transforms is precisely the kind of
// duplication that has historically let ctsem's backends drift apart.
//
// Two things here are choices rather than consequences, and both are surfaced
// to the caller rather than hidden:
//
//   * State-dependent cells (the predict/td/update transform groups) have no
//     single value -- they are functions of the latent state. They are
//     evaluated at a caller-supplied state, defaulting to T0MEANS, with
//     tdpreds fixed at zero. For a linear model this is exact, because no cell
//     depends on the state; for a nonlinear one it is a conditional value, and
//     stateDependentPositions() reports which cells those are so a summary can
//     say so rather than pretending a single number describes them.
//
//   * The derived matrices (DIFFUSIONcov, MANIFESTcov, T0cov, asymDIFFUSIONcov,
//     asymCINT) follow the generated Stan code's definitions, including its
//     restriction of the diffusion-related ones to the states that carry their
//     own diffusion.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "engine.hpp"

namespace ctsemcpp {

// The flat layout of everything parameterMatrices() returns: the engine's own
// model matrices in their table order, then the derived ones.
struct SummaryLayout {
  std::vector<std::string> name;
  std::vector<int> nrow, ncol, offset;
  int size = 0;
  int nbase = 0;  // how many entries came from model.layouts

  void add(const std::string& nm, int r, int c) {
    name.push_back(nm);
    nrow.push_back(r);
    ncol.push_back(c);
    offset.push_back(size);
    size += r * c;
  }
};

inline SummaryLayout summaryLayout(const CppModel& model) {
  SummaryLayout out;
  for (const MatrixLayout& L : model.layouts) out.add(L.name, L.nrow, L.ncol);
  out.nbase = static_cast<int>(model.layouts.size());
  const int n = model.nlatent;
  const int m = model.nmanifest;
  out.add("DIFFUSIONcov", n, n);
  out.add("MANIFESTcov", m, m);
  out.add("T0cov", n, n);
  out.add("asymDIFFUSIONcov", n, n);
  out.add("asymCINT", n, 1);
  return out;
}

// Flat all_params positions written by a state-dependent transform group, so a
// caller can mark those cells as conditional on the state they were evaluated
// at rather than reporting them as if they were constants.
inline std::vector<int> stateDependentPositions(const CppModel& model) {
  std::vector<int> out;
  const TransformGroup* groups[3] = {&model.predict, &model.td, &model.update};
  for (int g = 0; g < 3; ++g) {
    out.insert(out.end(), groups[g]->pos.begin(), groups[g]->pos.end());
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Pack whatever is currently in `ws.all_params` into the flat summary layout,
// deriving the covariance and asymptotic matrices from it.
//
// Split out from parameterMatrices() because there are two ways to arrive at a
// materialized parameter vector -- transform it from a raw vector, or take the
// one a subject's filter pass ended with -- and only the first half differs.
// `out` must have room for layout.size doubles.
inline void packMatrices(const CppModel& model, FilterWorkspace& ws,
                         const SummaryLayout& layout, double* out) {
  const int n = model.nlatent;
  const int m = model.nmanifest;

  for (int i = 0; i < layout.nbase; ++i) {
    const MatrixLayout& L = model.layouts[static_cast<std::size_t>(i)];
    const int count = L.nrow * L.ncol;
    std::copy(ws.all_params.begin() + L.offset, ws.all_params.begin() + L.offset + count,
              out + layout.offset[static_cast<std::size_t>(i)]);
  }

  Eigen::Map<const MatrixXd> DIFFUSION(ws.all_params.data() + model.offDIFFUSION, n, n);
  Eigen::Map<const MatrixXd> MANIFESTVAR(ws.all_params.data() + model.offMANIFESTVAR, m, m);
  Eigen::Map<const MatrixXd> T0VAR(ws.all_params.data() + model.offT0VAR, n, n);
  Eigen::Map<const MatrixXd> DRIFT(ws.all_params.data() + model.offDRIFT, n, n);
  Eigen::Map<const MatrixXd> CINT(ws.all_params.data() + model.offCINT, n, 1);

  MatrixXd diffusioncov(n, n), manifestcov(m, m), t0cov(n, n);
  sdcovsqrt2cov(DIFFUSION, n, diffusioncov, ws.covO, ws.covB);
  sdcovsqrt2cov(MANIFESTVAR, m, manifestcov, ws.covO, ws.covB);
  sdcovsqrt2cov(T0VAR, n, t0cov, ws.covO, ws.covB);

  // Stan restricts DIFFUSIONcov to the states carrying their own diffusion; the
  // augmented carrier states an intoverpop model adds are structurally zero
  // there, and writing that explicitly keeps it visible.
  const std::vector<int>& dyn = model.diffusionStates;
  const int k = static_cast<int>(dyn.size());
  MatrixXd restricted = MatrixXd::Zero(n, n);
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) restricted(dyn[i], dyn[j]) = diffusioncov(dyn[i], dyn[j]);

  MatrixXd asymdiffusion = MatrixXd::Zero(n, n);
  MatrixXd asymcint = MatrixXd::Zero(n, 1);
  if (k > 0) {
    MatrixXd Ad(k, k), Qd(k, k), X;
    VectorXd cintd(k);
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < k; ++i) Ad(i, j) = DRIFT(dyn[i], dyn[j]);
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < k; ++i) Qd(i, j) = restricted(dyn[i], dyn[j]);
    for (int i = 0; i < k; ++i) cintd(i) = CINT(dyn[i], 0);

    // An unstable or singular drift has no asymptotic form. That is a property
    // of the parameter value, not an error in the caller, so it comes back NaN
    // rather than thrown: a posterior sample that happens to be non-stationary
    // should not abort the summary of the other 199.
    //
    // Discrete time solves the discrete Lyapunov equation X = A X A' + Q via
    // (I - A (x) A) vec(X) = vec(Q), and the asymptotic intercept becomes
    // (I - A)^-1 c. Both are Stan's own definitions for a discrete model.
    bool ok = true;
    try {
      if (model.continuousTime) {
        lyapSolve(Ad, Qd, X, ws.lyapws);
      } else {
        const int kk = k * k;
        MatrixXd system = MatrixXd::Identity(kk, kk);
        for (int b = 0; b < k; ++b)
          for (int a = 0; a < k; ++a)
            for (int j = 0; j < k; ++j)
              for (int i = 0; i < k; ++i) {
                system(i + k * j, a + k * b) -= Ad(i, a) * Ad(j, b);
              }
        Eigen::Map<const VectorXd> q(Qd.data(), kk);
        Eigen::FullPivLU<MatrixXd> lyaplu(system);
        if (!lyaplu.isInvertible()) throw std::runtime_error("singular");
        VectorXd solved = lyaplu.solve(q);
        X = Eigen::Map<const MatrixXd>(solved.data(), k, k);
      }
    } catch (...) {
      ok = false;
    }
    if (ok && X.allFinite()) {
      for (int j = 0; j < k; ++j)
        for (int i = 0; i < k; ++i) asymdiffusion(dyn[i], dyn[j]) = X(i, j);
    } else {
      asymdiffusion.setConstant(std::nan(""));
    }

    MatrixXd intercept = -Ad;
    if (!model.continuousTime) {
      for (int i = 0; i < k; ++i) intercept(i, i) += 1.0;
    }
    Eigen::FullPivLU<MatrixXd> lu(intercept);
    if (lu.isInvertible()) {
      VectorXd solved = lu.solve(cintd);
      for (int i = 0; i < k; ++i) asymcint(dyn[i], 0) = solved(i);
    } else {
      asymcint.setConstant(std::nan(""));
    }
  }

  const int nnames = static_cast<int>(layout.name.size());
  for (int i = layout.nbase; i < nnames; ++i) {
    const std::string& nm = layout.name[static_cast<std::size_t>(i)];
    const MatrixXd* source = nullptr;
    if (nm == "DIFFUSIONcov") source = &restricted;
    else if (nm == "MANIFESTcov") source = &manifestcov;
    else if (nm == "T0cov") source = &t0cov;
    else if (nm == "asymDIFFUSIONcov") source = &asymdiffusion;
    else if (nm == "asymCINT") source = &asymcint;
    if (source == nullptr) continue;
    std::copy(source->data(), source->data() + source->size(),
              out + layout.offset[static_cast<std::size_t>(i)]);
  }
}

// `tipreds` must have model.ntipred entries (all zero gives the population
// values); `stateIn` may be null, in which case T0MEANS is used.
inline void parameterMatrices(const CppModel& model, FilterWorkspace& ws,
                              const SummaryLayout& layout, const double* values,
                              const double* tipreds, const double* stateIn, double time,
                              double dt, double* out) {
  const int n = model.nlatent;

  materializeParameters(model, ws, values, tipreds);

  if (stateIn != nullptr) {
    for (int i = 0; i < n; ++i) ws.state(i) = stateIn[i];
  } else {
    for (int i = 0; i < n; ++i) ws.state(i) = ws.all_params[model.offT0MEANS + i];
  }

  std::fill(ws.tdrow.begin(), ws.tdrow.end(), 0.0);
  ExprContext base;
  base.cells = ws.all_params.data();
  base.ncells = model.nall;
  base.state = ws.state.data();
  base.nstate = n;
  base.tdpreds = ws.tdrow.data();
  base.ntdpred = model.ntdpred;
  base.tipreds = tipreds;
  base.ntipred = model.ntipred;
  base.time = time;
  base.dt = dt;
  // All three groups, in filter order, so every Jacobian block is populated.
  applyGroup(model, ws, model.predict, base);
  applyGroup(model, ws, model.td, base);
  applyGroup(model, ws, model.update, base);

  packMatrices(model, ws, layout, out);
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_SUMMARY_HPP
