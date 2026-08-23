#ifndef CTSEMCPP_REVERSE_HPP
#define CTSEMCPP_REVERSE_HPP

// The reverse-mode pass: walk one subject's tape backwards.
//
// Conventions, matching the Julia adjoint this is ported from:
//
//   * `x_bar` and `P_bar` are the running cotangents of "the state and
//     covariance at this point in the forward pass". Every reverse step
//     consumes them as the cotangent of its outputs and leaves them as the
//     cotangent of its inputs. Identity copies in the primal (P_update =
//     P_predict between substeps, and on a fully missing row) therefore need
//     no reverse code at all.
//   * `theta_bar` is the cotangent of `all_params`, in exactly the primal's
//     flat layout, so a model matrix's cotangent is a contiguous block.
//   * `Theta_bar` is the running cotangent of the manifest covariance, kept
//     separate from `theta_bar` because the point at which the forward pass
//     builds it from MANIFESTVAR is a tape entry of its own.
//   * Covariance cotangents are symmetrised before use. The primal maintains
//     symmetry explicitly, and those operations are the identity on symmetric
//     perturbations, which is what makes that valid. Several steps below then
//     *rely* on that exact symmetry to avoid recomputing a product that would
//     otherwise be formed twice; those places are marked.
//   * The +1e-10 ridges are affine, so they are invisible to the reverse pass
//     except through the *value* at which downstream derivatives are
//     evaluated -- hence the reverse recomputes `Pr = P + 1e-10 I` rather than
//     reusing the unridged `P`.
//
// Every matrix below lives in `ReverseScratch`, resized once per model and
// reused for every row of every subject of every gradient evaluation, and every
// product is written through `.noalias()`. Both matter: the reverse pass forms
// roughly twenty small dense products per row, and at that size Eigen's own
// temporaries and their allocation are a large fraction of the arithmetic.

#include <RcppEigen.h>

#include <vector>

#include "filter.hpp"
#include "linalg.hpp"
#include "model.hpp"
#include "tape.hpp"

namespace ctsemcpp {

// Preallocated working storage for one reverse step. Sized for the largest
// possible observed subset (all manifest variables); rows with fewer observed
// variables use leading blocks of the same buffers.
struct ReverseScratch {
  MatrixXd Pr, M, Ps, Mbar, PbarNew, tmpNN, tmpNN2, Abar, JAxbar, QcBar, diffusionBar;
  MatrixXd PHt, G, PHtBar, Gbar, tmpNM;
  MatrixXd S, Sinv, Sbar, Sbar0, Rbar, tmpMM, eye;
  MatrixXd Hbar, Lbar;
  VectorXd yt, alpha, alphabar, beta, ytbar, xbarNew, dINTbar;
  MatrixXd Ad, JAxd, Qb, Xbar, Adbar, JAxdbar, Qcdbar, tmpKK, scaled, frechet;
  VectorXd sbar, affinebar, dINTbarDyn;
  Eigen::LLT<MatrixXd> llt;
  Eigen::PartialPivLU<MatrixXd> lu;

  void resize(int n, int m, int k) {
    Pr.resize(n, n); M.resize(n, n); Ps.resize(n, n); Mbar.resize(n, n);
    PbarNew.resize(n, n); tmpNN.resize(n, n); tmpNN2.resize(n, n);
    Abar.resize(n, n); JAxbar.resize(n, n); QcBar.resize(n, n); diffusionBar.resize(n, n);
    PHt.resize(n, m); G.resize(n, m); PHtBar.resize(n, m); Gbar.resize(n, m); tmpNM.resize(n, m);
    S.resize(m, m); Sinv.resize(m, m); Sbar.resize(m, m); Sbar0.resize(m, m);
    Rbar.resize(m, m); tmpMM.resize(m, m);
    Hbar.resize(m, n); Lbar.resize(m, n);
    eye = MatrixXd::Identity(m, m);
    yt.resize(m); alpha.resize(m); alphabar.resize(m); beta.resize(m); ytbar.resize(m);
    xbarNew.resize(n); dINTbar.resize(n);
    Ad.resize(k, k); JAxd.resize(k, k); Qb.resize(k, k); Xbar.resize(k, k);
    Adbar.resize(k, k); JAxdbar.resize(k, k); Qcdbar.resize(k, k); tmpKK.resize(k, k);
    scaled.resize(n, n);
    sbar.resize(k); affinebar.resize(k); dINTbarDyn.resize(k);
  }
};

struct AdjointWorkspace {
  std::vector<double> theta_bar;          // cotangent of all_params
  std::vector<double> subject_values_bar;
  std::vector<double> group_scratch;      // replay buffer for transform groups
  VectorXd x_bar;
  MatrixXd P_bar, Theta_bar;
  LyapWorkspace lyapws;
  FrechetWorkspace frechetws;
  ReverseScratch sc;

  // Deferred matrix-exponential Frechet derivative; see `flushFrechet`.
  bool frechetPending = false;
  bool deferFrechet = false;
  double frechetDt = 0.0;
  MatrixXd frechetA;        // the `JAx * dt` the pending directions belong to
  MatrixXd frechetAccum;    // sum of the directions still to be pushed through
  MatrixXd jaxBarDeferred;  // Frechet contributions held back for one final pass
  std::vector<double> exprWork, exprAdj;
  std::vector<double> groupSaved, groupWritten;

  void resize(const CppModel& model) {
    theta_bar.assign(model.nall, 0.0);
    subject_values_bar.assign(model.nvalues, 0.0);
    // NaN, not zero: only a group's relevant cells are ever written here, so a
    // transform reading anything that was not recorded poisons the gradient
    // visibly instead of quietly using a stale value.
    group_scratch.assign(model.nall, std::numeric_limits<double>::quiet_NaN());
    x_bar.setZero(model.nlatent);
    P_bar.setZero(model.nlatent, model.nlatent);
    Theta_bar.setZero(model.nmanifest, model.nmanifest);
    sc.resize(model.nlatent, model.nmanifest,
              static_cast<int>(model.diffusionStates.size()));
    frechetPending = false;
    deferFrechet = model.frechetDeferrable();
    frechetA.setZero(model.nlatent, model.nlatent);
    frechetAccum.setZero(model.nlatent, model.nlatent);
    jaxBarDeferred.setZero(model.nlatent, model.nlatent);
  }
};

namespace detail {

template <typename Derived>
inline void symmetrizeInto(MatrixXd& out, const Eigen::MatrixBase<Derived>& A) {
  const int n = static_cast<int>(A.rows());
  out.resize(n, n);
  for (int j = 0; j < n; ++j) {
    out(j, j) = A(j, j);
    for (int i = j + 1; i < n; ++i) {
      const double v = 0.5 * (A(i, j) + A(j, i));
      out(i, j) = v;
      out(j, i) = v;
    }
  }
}

}  // namespace detail

// Push any deferred matrix-exponential Frechet derivative into the JAx
// cotangent.
//
// `A = exp(JAx dt)` is the single most expensive step in the reverse pass: its
// adjoint is a matrix exponential of *twice* the state dimension, which at 20
// latents measured 540 us against ~30 us for every other step in one prediction
// substep -- roughly two thirds of the whole reverse pass.
//
// But the Frechet derivative is **linear in its direction argument**, and a
// balanced panel design hands it the same `JAx * dt` at every substep of every
// subject. So instead of `sum_steps dt * L(A', Abar_step)`, accumulate the
// directions and evaluate `dt * L(A', sum_steps Abar_step)` once. That is exact,
// not an approximation, and it turns an O(rows) count of block exponentials into
// O(distinct (JAx, dt) pairs) -- one, for a linear model on an equally spaced
// panel.
//
// The guard is an exact comparison against the actual `JAx * dt` that produced
// the pending directions, so nothing here assumes linearity: a state-dependent
// model changes JAx per row, misses every time, and simply pays the comparison,
// while irregular observation times batch within each distinct interval.
//
// Correctness depends on nothing consuming the JAx cotangent while a flush is
// outstanding. Two things can: a state-dependent transform group that *writes*
// a JAx cell (whose reverse zeroes that cell's cotangent), and the parameter
// layer. Both flush first.
inline void flushFrechet(const CppModel& model, AdjointWorkspace& aws) {
  if (!aws.frechetPending) return;
  aws.frechetPending = false;
  const int n = model.nlatent;
  expmFrechetAdjoint(aws.frechetA, aws.frechetAccum, aws.sc.frechet, aws.frechetws);
  if (aws.deferFrechet) {
    // Held back so the batch can span subjects too; `gradient()` pushes this
    // through the parameter layer once at the end.
    aws.jaxBarDeferred.noalias() += aws.frechetDt * aws.sc.frechet;
    return;
  }
  double* theta = aws.theta_bar.data();
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      theta[model.offJAx + j * n + i] += aws.frechetDt * aws.sc.frechet(i, j);
}

// Undo one prediction substep.
//
// Forward, on the dynamic index subset D of size k:
//     A          = exp(JAx dt)                (also the discrete drift)
//     affine[i]  = CINT[Di] + sum_j (DRIFT[Di,j] - JAx[Di,j]) x[j]
//     s[i]       = -affine[i] + sum_q A[Di,Dq] affine[q]
//     dINT[D]    = JAx[D,D] \ s
//     X          = lyap(JAx[D,D], Qc[D,D])
//     dDIFF[D,D] = X - A[D,D] X A[D,D]'
//     x+         = A x + dINT
//     P+         = A (P + eps I) A' + dDIFF
inline void reversePredict(const CppModel& model, AdjointWorkspace& aws,
                           const PredictRecord& rec) {
  const int n = model.nlatent;
  const std::vector<int>& dyn = model.diffusionStates;
  const int k = static_cast<int>(dyn.size());
  ReverseScratch& sc = aws.sc;

  const MatrixXd& A = rec.A;
  const MatrixXd& JAx = rec.JAx;
  const VectorXd& x = rec.state_in;

  double* thetaAll = aws.theta_bar.data();

  // The discrete-time reverse pass is the continuous one with its three hard
  // pieces removed rather than a second implementation of it: A is JAx (no
  // Frechet derivative of an exponential), dDIFFUSION is the diffusion
  // covariance itself (no Lyapunov pullback), and dINT is the affine offset
  // (no linear solve). What remains is the same mean/covariance recursion,
  // written out here because sharing it with the continuous branch would mean
  // branching inside every step of it.
  if (!model.continuousTime) {
    sc.Abar.noalias() = aws.x_bar * x.transpose();
    sc.dINTbar = aws.x_bar;
    sc.xbarNew.noalias() = A.transpose() * aws.x_bar;

    detail::symmetrizeInto(sc.Ps, aws.P_bar);
    sc.Pr = rec.P_in;
    for (int i = 0; i < n; ++i) sc.Pr(i, i) += kRidge;
    sc.tmpNN.noalias() = sc.Ps * A;
    sc.Abar.noalias() += 2.0 * (sc.tmpNN * sc.Pr);
    sc.tmpNN2.noalias() = A.transpose() * sc.Ps;
    sc.PbarNew.noalias() = sc.tmpNN2 * A;

    // dDIFFUSION[D,D] = Qc[D,D], so its cotangent passes straight through.
    for (int j = 0; j < k; ++j) {
      for (int i = 0; i < k; ++i) {
        sc.Qcdbar(i, j) = 0.5 * (sc.Ps(dyn[i], dyn[j]) + sc.Ps(dyn[j], dyn[i]));
      }
    }

    // dINT[i] = CINT[i] + sum_j (DRIFT[i,j] - JAx[i,j]) x[j], over every state.
    sc.JAxbar.setZero();
    for (int i = 0; i < n; ++i) thetaAll[model.offCINT + i] += sc.dINTbar(i);
    for (int j = 0; j < n; ++j) {
      double* driftCol = thetaAll + model.offDRIFT + j * n;
      double acc = 0.0;
      for (int i = 0; i < n; ++i) {
        const double contribution = sc.dINTbar(i) * x(j);
        driftCol[i] += contribution;
        sc.JAxbar(i, j) -= contribution;
        acc += (rec.DRIFT(i, j) - JAx(i, j)) * sc.dINTbar(i);
      }
      sc.xbarNew(j) += acc;
    }
    // A *is* JAx here, so its cotangent simply adds.
    sc.JAxbar += sc.Abar;
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) thetaAll[model.offJAx + j * n + i] += sc.JAxbar(i, j);
    }

    sc.QcBar.setZero();
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < k; ++i) sc.QcBar(dyn[i], dyn[j]) = sc.Qcdbar(i, j);
    sc.diffusionBar.setZero();
    sdcovsqrt2covPullback(sc.diffusionBar, rec.DIFFUSION, sc.QcBar, n);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        thetaAll[model.offDIFFUSION + j * n + i] += sc.diffusionBar(i, j);
      }
    }

    aws.x_bar = sc.xbarNew;
    aws.P_bar = sc.PbarNew;
    return;
  }

  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) { sc.Ad(i, j) = A(dyn[i], dyn[j]); sc.JAxd(i, j) = JAx(dyn[i], dyn[j]); }

  // mean: x+ = A x + dINT
  sc.Abar.noalias() = aws.x_bar * x.transpose();
  sc.dINTbar = aws.x_bar;
  sc.xbarNew.noalias() = A.transpose() * aws.x_bar;

  // covariance: P+ = A (P + eps I) A' + dDIFF
  detail::symmetrizeInto(sc.Ps, aws.P_bar);
  sc.Pr = rec.P_in;
  for (int i = 0; i < n; ++i) sc.Pr(i, i) += kRidge;
  // d(A Pt A')/dA contracted with Pbar is Pbar A Pt' + Pbar' A Pt; both are
  // symmetric here, so that collapses to twice one term.
  sc.tmpNN.noalias() = sc.Ps * A;
  sc.Abar.noalias() += 2.0 * (sc.tmpNN * sc.Pr);
  sc.tmpNN2.noalias() = A.transpose() * sc.Ps;
  sc.PbarNew.noalias() = sc.tmpNN2 * A;

  // dDIFF[D,D] = X - Ad X Ad'; the cotangent of dDIFF is the symmetrised Pbar.
  const MatrixXd& X = rec.Xlyap;
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) sc.Qb(i, j) = 0.5 * (sc.Ps(dyn[i], dyn[j]) + sc.Ps(dyn[j], dyn[i]));
  sc.tmpKK.noalias() = sc.Ad.transpose() * sc.Qb;
  sc.Xbar = sc.Qb;
  sc.Xbar.noalias() -= sc.tmpKK * sc.Ad;
  // Qb is symmetric by construction, so `Qb' Ad` is the same product as
  // `Qb Ad` and the two terms of the Ad cotangent share it.
  sc.tmpKK.noalias() = sc.Qb * sc.Ad;
  sc.Adbar.noalias() = sc.tmpKK * X.transpose();
  sc.Adbar.noalias() += sc.tmpKK * X;
  sc.Adbar *= -1.0;

  // X = lyap(JAx[D,D], Qc[D,D])
  lyapPullback(sc.JAxd, X, sc.Xbar, sc.JAxdbar, sc.Qcdbar, aws.lyapws);

  // dINT[D] = JAx[D,D] \ s -- a linear solve
  for (int i = 0; i < k; ++i) sc.dINTbarDyn(i) = sc.dINTbar(dyn[i]);
  sc.tmpKK = sc.JAxd.transpose();
  sc.lu.compute(sc.tmpKK);
  sc.sbar = sc.lu.solve(sc.dINTbarDyn);
  sc.JAxdbar.noalias() -= sc.sbar * rec.dINTdyn.transpose();

  // s = -affine + Ad affine
  sc.affinebar.noalias() = sc.Ad.transpose() * sc.sbar;
  sc.affinebar -= sc.sbar;
  sc.Adbar.noalias() += sc.sbar * rec.affine.transpose();

  // affine[i] = CINT[Di] + sum_j (DRIFT[Di,j] - JAx[Di,j]) x[j]
  double* theta = aws.theta_bar.data();
  sc.JAxbar.setZero();
  for (int i = 0; i < k; ++i) theta[model.offCINT + dyn[i]] += sc.affinebar(i);
  for (int j = 0; j < n; ++j) {
    double* driftCol = theta + model.offDRIFT + j * n;
    double acc = 0.0;
    for (int i = 0; i < k; ++i) {
      const double contribution = sc.affinebar(i) * x(j);
      driftCol[dyn[i]] += contribution;
      sc.JAxbar(dyn[i], j) -= contribution;
      acc += (rec.DRIFT(dyn[i], j) - JAx(dyn[i], j)) * sc.affinebar(i);
    }
    sc.xbarNew(j) += acc;
  }

  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < k; ++i) {
      sc.Abar(dyn[i], dyn[j]) += sc.Adbar(i, j);
      sc.JAxbar(dyn[i], dyn[j]) += sc.JAxdbar(i, j);
    }
  }

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) theta[model.offJAx + j * n + i] += sc.JAxbar(i, j);

  // A = exp(JAx * dt). The direction is queued rather than pushed through
  // immediately, so substeps sharing a `JAx * dt` cost one block exponential
  // between them instead of one each -- see `flushFrechet`.
  sc.scaled.noalias() = rec.dt * JAx;
  if (aws.frechetPending && aws.frechetDt == rec.dt &&
      detail::blocksIdentical(aws.frechetA, sc.scaled)) {
    aws.frechetAccum += sc.Abar;
  } else {
    flushFrechet(model, aws);
    aws.frechetA = sc.scaled;
    aws.frechetDt = rec.dt;
    aws.frechetAccum = sc.Abar;
    aws.frechetPending = true;
  }

  // Qc = sdcovsqrt2cov(DIFFUSION); only the dynamic block was consumed.
  sc.QcBar.setZero();
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) sc.QcBar(dyn[i], dyn[j]) = sc.Qcdbar(i, j);
  sc.diffusionBar.setZero();
  sdcovsqrt2covPullback(sc.diffusionBar, rec.DIFFUSION, sc.QcBar, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) theta[model.offDIFFUSION + j * n + i] += sc.diffusionBar(i, j);

  aws.x_bar = sc.xbarNew;
  aws.P_bar = sc.PbarNew;
}

// Undo one TD-predictor impulse: x+ = x + TDPREDEFFECT td, P+ = Jtd P Jtd'.
inline void reverseTd(const CppModel& model, AdjointWorkspace& aws, const TdRecord& rec) {
  const int n = model.nlatent;
  const int ntd = static_cast<int>(rec.tdpreds.size());
  ReverseScratch& sc = aws.sc;
  double* theta = aws.theta_bar.data();
  for (int j = 0; j < ntd; ++j) {
    double* col = theta + model.offTDPREDEFFECT + j * n;
    const double td = rec.tdpreds(j);
    for (int i = 0; i < n; ++i) col[i] += aws.x_bar(i) * td;
  }

  detail::symmetrizeInto(sc.Ps, aws.P_bar);
  const MatrixXd& Jtd = rec.Jtd;
  // Ps is exactly symmetric, so Ps' * Jtd == Ps * Jtd and the two terms of the
  // Jtd cotangent share one product.
  sc.tmpNN.noalias() = sc.Ps * Jtd;
  sc.Mbar.noalias() = sc.tmpNN * rec.P_in.transpose();
  sc.Mbar.noalias() += sc.tmpNN * rec.P_in;
  for (int j = 0; j < n; ++j) {
    double* col = theta + model.offJtd + j * n;
    for (int i = 0; i < n; ++i) col[i] += sc.Mbar(i, j);
  }
  sc.tmpNN2.noalias() = Jtd.transpose() * sc.Ps;
  aws.P_bar.noalias() = sc.tmpNN2 * Jtd;
  // The mean update is a pure translation, so x_bar passes through unchanged.
}

// Undo one measurement update, including its log-likelihood contribution.
//
// Forward (on the observed subset, with eps the Stan-matching ridge):
//     Pr  = P + eps I
//     PHt = Pr H'
//     S   = sym(H PHt + R) + eps I
//     yt  = y - (Lambda x + mu)
//     a   = S^-1 yt
//     x+  = x + PHt a
//     G   = PHt S^-1;  M = I - G H
//     P+  = M P M' + G R G'          (Joseph form, on the unridged P)
//     ll  = -0.5 (m log 2pi + logdet S + yt' a)
//
// The seed for `ll` is 1: the reverse pass differentiates the summed
// log-likelihood, so every row contributes with unit weight.
inline void reverseUpdate(const CppModel& model, AdjointWorkspace& aws,
                          const UpdateRecord& rec) {
  const int n = model.nlatent;
  const int m = static_cast<int>(rec.observed.size());
  ReverseScratch& sc = aws.sc;
  const MatrixXd& H = rec.H;
  const MatrixXd& L = rec.Lambda;
  const MatrixXd& R = rec.R;
  const VectorXd& x = rec.state_in;

  auto PHt = sc.PHt.leftCols(m);
  auto G = sc.G.leftCols(m);
  auto PHtBar = sc.PHtBar.leftCols(m);
  auto Gbar = sc.Gbar.leftCols(m);
  auto tmpNM = sc.tmpNM.leftCols(m);
  auto S = sc.S.topLeftCorner(m, m);
  auto Sinv = sc.Sinv.topLeftCorner(m, m);
  auto Sbar = sc.Sbar.topLeftCorner(m, m);
  auto Rbar = sc.Rbar.topLeftCorner(m, m);
  auto tmpMM = sc.tmpMM.topLeftCorner(m, m);
  auto Hbar = sc.Hbar.topRows(m);
  auto Lbar = sc.Lbar.topRows(m);
  auto yt = sc.yt.head(m);
  auto alpha = sc.alpha.head(m);
  auto alphabar = sc.alphabar.head(m);
  auto beta = sc.beta.head(m);
  auto ytbar = sc.ytbar.head(m);

  sc.Pr = rec.P_in;
  for (int i = 0; i < n; ++i) sc.Pr(i, i) += kRidge;
  PHt.noalias() = sc.Pr * H.transpose();
  sc.tmpMM.topLeftCorner(m, m).noalias() = H * PHt;
  sc.tmpMM.topLeftCorner(m, m) += R;
  S = 0.5 * (sc.tmpMM.topLeftCorner(m, m) + sc.tmpMM.topLeftCorner(m, m).transpose());
  for (int i = 0; i < m; ++i) S(i, i) += kRidge;

  // S is symmetric positive definite here (the primal already Cholesky-
  // factorized it successfully), so an LLT inverse is both cheaper and better
  // conditioned than a general LU one.
  sc.llt.compute(S);
  Sinv = sc.llt.solve(sc.eye.topLeftCorner(m, m));

  yt.noalias() = L * x;
  for (int i = 0; i < m; ++i) yt(i) = rec.y(i) - (yt(i) + rec.manifestmeans(i));
  alpha.noalias() = Sinv * yt;
  G.noalias() = PHt * Sinv;
  sc.M.setIdentity();
  sc.M.noalias() -= G * H;

  // log-likelihood contribution (unit seed)
  Sbar = -0.5 * Sinv;
  Sbar.noalias() += 0.5 * (alpha * alpha.transpose());
  ytbar = -alpha;

  // x+ = x + PHt a
  sc.xbarNew = aws.x_bar;
  PHtBar.noalias() = aws.x_bar * alpha.transpose();
  alphabar.noalias() = PHt.transpose() * aws.x_bar;

  // a = S^-1 yt
  beta.noalias() = Sinv * alphabar;
  ytbar += beta;
  Sbar.noalias() -= beta * alpha.transpose();

  // P+ = M P M' + G R G'
  detail::symmetrizeInto(sc.Ps, aws.P_bar);
  const MatrixXd& Pin = rec.P_in;
  // Ps is exactly symmetric, so `Ps' M` equals `Ps M` and the two terms of
  // each cotangent below share one product rather than forming it twice.
  sc.tmpNN.noalias() = sc.Ps * sc.M;
  sc.Mbar.noalias() = sc.tmpNN * Pin.transpose();
  sc.Mbar.noalias() += sc.tmpNN * Pin;
  sc.tmpNN2.noalias() = sc.M.transpose() * sc.Ps;
  sc.PbarNew.noalias() = sc.tmpNN2 * sc.M;
  tmpNM.noalias() = sc.Ps * G;
  Gbar.noalias() = tmpNM * R.transpose();
  Gbar.noalias() += tmpNM * R;
  Rbar.noalias() = tmpNM.transpose() * G;

  // M = I - G H
  Gbar.noalias() -= sc.Mbar * H.transpose();
  Hbar.noalias() = G.transpose() * sc.Mbar;
  Hbar *= -1.0;

  // G = PHt S^-1
  PHtBar.noalias() += Gbar * Sinv;
  tmpMM.noalias() = G.transpose() * Gbar;
  Sbar.noalias() -= tmpMM * Sinv;

  // S = sym(H PHt + R) + eps I
  auto Sbar0 = sc.Sbar0.topLeftCorner(m, m);
  Sbar0 = 0.5 * (Sbar + Sbar.transpose());
  Hbar.noalias() += Sbar0 * PHt.transpose();
  PHtBar.noalias() += H.transpose() * Sbar0;
  Rbar += Sbar0;

  // PHt = Pr H'
  sc.PbarNew.noalias() += PHtBar * H;
  Hbar.noalias() += PHtBar.transpose() * sc.Pr;

  // yt = y - (Lambda x + mu)
  Lbar.noalias() = ytbar * x.transpose();
  Lbar *= -1.0;
  sc.xbarNew.noalias() -= L.transpose() * ytbar;

  double* theta = aws.theta_bar.data();
  for (int i = 0; i < m; ++i) {
    const int oi = rec.observed[i];
    theta[model.offMANIFESTMEANS + oi] -= ytbar(i);
    for (int j = 0; j < n; ++j) {
      theta[model.offLAMBDA + j * model.nmanifest + oi] += Lbar(i, j);
      theta[model.offJy + j * model.nmanifest + oi] += Hbar(i, j);
    }
    for (int j = 0; j < m; ++j) aws.Theta_bar(oi, rec.observed[j]) += Rbar(i, j);
  }

  aws.x_bar = sc.xbarNew;
  aws.P_bar = sc.PbarNew;
}

// Reverse one recorded group of state-dependent transforms.
//
// Walks the group backwards, because within a group a later transform may read
// a cell an earlier one has already written. For each transform this takes the
// cotangent currently on the cell it wrote, zeroes it (the cell was
// *overwritten*, so the value there beforehand did not reach the likelihood
// through this path), and adds cotangent * d f / d input to every parameter
// cell and state entry the transform reads.
//
// Each transform's derivative is evaluated at the values that transform
// actually saw, which are not the group's final values: transform k runs after
// 1..k-1 have written their cells but before k+1..K have written theirs. The
// group is therefore replayed forward once from `params_before`, and the
// reverse walk unwinds those writes one at a time as it descends.
//
// The derivative itself comes from one reverse sweep over the expression's AST,
// which yields the partials with respect to every cell and state entry the
// expression reads at once. The Julia backend instead seeds a ForwardDiff dual
// once per input, so it pays `reads + state_dim` evaluations of the expression
// where this pays one forward and one reverse.
inline void reverseGroup(const CppModel& model, AdjointWorkspace& aws,
                         const GroupRecord& rec) {
  const TransformGroup& group = rec.group == 0 ? model.predict
                              : rec.group == 1 ? model.td
                                               : model.update;
  const int K = static_cast<int>(group.pos.size());
  if (K == 0) return;

  double* working = aws.group_scratch.data();
  for (std::size_t j = 0; j < group.relevant.size(); ++j) {
    working[group.relevant[j]] = rec.params_before[j];
  }

  ExprContext ctx;
  ctx.cells = working;
  ctx.ncells = model.nall;
  ctx.state = rec.state.data();
  ctx.nstate = model.nlatent;
  ctx.tdpreds = rec.tdpreds.empty() ? nullptr : rec.tdpreds.data();
  ctx.ntdpred = static_cast<int>(rec.tdpreds.size());
  ctx.time = rec.time;
  ctx.dt = rec.dt;

  aws.groupSaved.resize(K);
  aws.groupWritten.resize(K);
  for (int k = 0; k < K; ++k) aws.groupSaved[k] = working[group.pos[k]];
  for (int k = 0; k < K; ++k) {
    const double value = group.expr[k].eval(ctx, aws.exprWork);
    aws.groupWritten[k] = value;
    working[group.pos[k]] = value;
  }
  // `working` now holds post-replay values; put the last transform's own cell
  // back, since it ran before its write landed.
  working[group.pos[K - 1]] = aws.groupSaved[K - 1];

  for (int k = K - 1; k >= 0; --k) {
    const int idx = group.pos[k];
    const double cotangent = aws.theta_bar[idx];
    aws.theta_bar[idx] = 0.0;
    if (cotangent != 0.0) {
      group.expr[k].eval(ctx, aws.exprWork);
      group.expr[k].backward(ctx, aws.exprWork, cotangent, aws.exprAdj,
                             aws.theta_bar.data(), nullptr, aws.x_bar.data());
    }
    if (k > 0) working[group.pos[k - 1]] = aws.groupSaved[k - 1];
  }
}

// Walk one subject's tape backwards, accumulating into `theta_bar`.
// Deliberately stops there rather than going on to the free parameters:
// whether `theta_bar` can be shared across subjects is a property of the model,
// not of one tape. `theta_bar` is not zeroed here; the caller owns that.
inline void reverseTape(const CppModel& model, AdjointWorkspace& aws, const AdjointTape& tape) {
  const int n = model.nlatent;
  const int m = model.nmanifest;
  aws.x_bar.setZero(n);
  aws.P_bar.setZero(n, n);
  aws.Theta_bar.setZero(m, m);
  double* theta = aws.theta_bar.data();

  for (int e = static_cast<int>(tape.program.size()) - 1; e >= 0; --e) {
    const TapeKind kind = tape.program[e].first;
    const int index = tape.program[e].second;
    switch (kind) {
      case TapeKind::Update: reverseUpdate(model, aws, tape.updates[index]); break;
      case TapeKind::Td: reverseTd(model, aws, tape.tds[index]); break;
      case TapeKind::Predict: reversePredict(model, aws, tape.predicts[index]); break;
      case TapeKind::Group:
        // Only a group that writes a JAx cell can consume a queued Frechet
        // contribution (its reverse zeroes that cell's cotangent). Flushing
        // before every group instead would break the batch on ctsem's default
        // model, whose individually varying MANIFESTMEANS puts a group between
        // every pair of prediction substeps.
        if (model.groupsWriteJAx) flushFrechet(model, aws);
        reverseGroup(model, aws, tape.groups[index]);
        break;
      case TapeKind::Theta: {
        aws.sc.tmpMM.setZero(m, m);
        sdcovsqrt2covPullback(aws.sc.tmpMM, tape.thetas[index].MANIFESTVAR, aws.Theta_bar, m);
        for (int j = 0; j < m; ++j)
          for (int i = 0; i < m; ++i) theta[model.offMANIFESTVAR + j * m + i] += aws.sc.tmpMM(i, j);
        aws.Theta_bar.setZero();
        break;
      }
      case TapeKind::Init: {
        for (int i = 0; i < n; ++i) theta[model.offT0MEANS + i] += aws.x_bar(i);
        aws.x_bar.setZero();
        aws.sc.tmpNN.setZero(n, n);
        detail::symmetrizeInto(aws.sc.Ps, aws.P_bar);
        sdcovsqrt2covPullback(aws.sc.tmpNN, tape.inits[index].T0VAR, aws.sc.Ps, n);
        for (int j = 0; j < n; ++j)
          for (int i = 0; i < n; ++i) theta[model.offT0VAR + j * n + i] += aws.sc.tmpNN(i, j);
        aws.P_bar.setZero();
        break;
      }
    }
  }
}

// Unwind the parameter layer: whatever cotangent is left on `all_params` was
// put there by the regular (non-state-dependent) transforms, so push it through
// them and then through the TI-predictor effects into `values_bar`.
inline void parameterLayer(const CppModel& model, AdjointWorkspace& aws,
                           const std::vector<double>& subject_values,
                           const double* tipreds, double* values_bar) {
  if (!aws.deferFrechet) flushFrechet(model, aws);
  std::fill(aws.subject_values_bar.begin(), aws.subject_values_bar.end(), 0.0);
  ExprContext ctx;
  ctx.params = subject_values.data();
  ctx.nparams = model.nvalues;
  for (std::size_t k = 0; k < model.regPos.size(); ++k) {
    const double cotangent = aws.theta_bar[model.regPos[k]];
    if (cotangent == 0.0) continue;
    model.regExpr[k].eval(ctx, aws.exprWork);
    model.regExpr[k].backward(ctx, aws.exprWork, cotangent, aws.exprAdj, nullptr,
                              aws.subject_values_bar.data(), nullptr);
  }
  for (int i = 0; i < model.nvalues; ++i) values_bar[i] += aws.subject_values_bar[i];
  for (std::size_t i = 0; i < model.tiParameter.size(); ++i) {
    values_bar[model.tiCoefficient[i]] +=
        aws.subject_values_bar[model.tiParameter[i]] * tipreds[model.tiPredictor[i]];
  }
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_REVERSE_HPP
