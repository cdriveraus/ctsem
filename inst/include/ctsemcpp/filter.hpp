#ifndef CTSEMCPP_FILTER_HPP
#define CTSEMCPP_FILTER_HPP

// The continuous-time extended Kalman filter: one primal loop, which optionally
// records the adjoint tape as it goes.
//
// There is deliberately no second, traced copy of the loop. A parallel traced
// loop is precisely how a hand-written reverse pass silently drifts out of sync
// with the forward pass it is meant to mirror; the cost of keeping the
// recording inline is that adding a step here means adding its record.
//
// The arithmetic mirrors `ContinuousTimeSEM/src/kalman_filters.jl` and
// `discrete_time_form.jl` step for step, including Stan's `makesym()` +1e-10
// ridges at the three points Stan applies them (before propagating the
// posterior covariance, before projecting the prior covariance into observation
// space, and before the innovation Cholesky). Those ridges are not defensive
// decoration: with MANIFESTVAR fixed at zero and LAMBDA the identity the
// innovation covariance is numerically the predicted state covariance, which
// sits exactly on the PD boundary for small trial DIFFUSION values, and a port
// that omits them rejects trial points Stan accepts.

#include <RcppEigen.h>

#include <cmath>
#include <limits>
#include <vector>

#include "linalg.hpp"
#include "model.hpp"
#include "tape.hpp"

namespace ctsemcpp {

constexpr double kRidge = 1e-10;
constexpr double kLog2Pi = 1.8378770664093454835606594728112;

// Last-value cache for the two expensive, row-invariant pieces of the
// discretization. Both are guarded by a comparison against the actual input
// matrices rather than by a model-level "is this linear" flag, so a
// state-dependent model simply misses every time and pays only the comparison,
// a cache cannot survive a parameter change between optimizer iterations
// (the parameters *are* the inputs being compared), and irregular observation
// times miss the exponential cache while still hitting the Lyapunov one.
//
// The Lyapunov entry is the important one: its solve takes JAx and the
// diffusion covariance and does *not* take dt, so for any model whose drift and
// diffusion are not state-dependent it is identical at every row of every
// subject, yet sits in the middle of a routine whose whole purpose is to
// produce dt-dependent quantities.
struct DiscretizationCache {
  MatrixXd expJAx, expOut;
  double expDt = 0.0;
  bool expValid = false;
  MatrixXd lyapJAx, lyapQ, lyapOut;
  bool lyapValid = false;

  void invalidate() { expValid = false; lyapValid = false; }
};

struct SubjectData {
  // manifest variables by observation; NaN marks a missing entry
  const double* y = nullptr;
  const double* tdpreds = nullptr;
  const double* tipreds = nullptr;
  const double* times = nullptr;
  int nobs = 0;
  int subject = 1;
};

struct FilterWorkspace {
  std::vector<double> subject_values;
  std::vector<double> all_params;

  VectorXd state;
  MatrixXd P_predict, P_update;
  MatrixXd Qc, Theta;             // continuous diffusion covariance, manifest covariance
  MatrixXd covO, covB;            // scratch for sdcovsqrt2cov
  MatrixXd eJAx, dDIFFUSION;
  VectorXd dINT;

  // dynamic-block scratch
  MatrixXd JAxd, Qcd, Xlyap, AdX;
  VectorXd affine, sVec, dINTdyn;

  MatrixXd Pr, PHt, S, Kgain, Mjoseph, Ptmp;
  VectorXd innovation, alpha;
  Eigen::LLT<MatrixXd> llt;

  LyapWorkspace lyapws;
  DiscretizationCache cache;

  std::vector<double> exprWork;
  std::vector<int> observed;
  std::vector<double> tdrow;

  void resize(const CppModel& model) {
    const int n = model.nlatent;
    const int m = model.nmanifest;
    const int k = static_cast<int>(model.diffusionStates.size());
    subject_values.assign(model.nvalues, 0.0);
    // NaN rather than zero: a state-dependent cell is written by its transform
    // group and read nowhere before that, so anything reading one early gets a
    // loud NaN instead of a plausible zero.
    all_params.assign(model.nall, std::numeric_limits<double>::quiet_NaN());
    state.setZero(n);
    P_predict.setZero(n, n);
    P_update.setZero(n, n);
    Qc.setZero(n, n);
    Theta.setZero(m, m);
    eJAx.setZero(n, n);
    dDIFFUSION.setZero(n, n);
    dINT.setZero(n);
    JAxd.setZero(k, k);
    Qcd.setZero(k, k);
    Xlyap.setZero(k, k);
    AdX.setZero(k, k);
    affine.setZero(k);
    sVec.setZero(k);
    dINTdyn.setZero(k);
    Pr.setZero(n, n);
    Ptmp.setZero(n, n);
    Mjoseph.setZero(n, n);
    tdrow.assign(model.ntdpred, 0.0);
    cache.invalidate();
    lyapws.reset();
  }
};

namespace detail {

inline void symmetrize(MatrixXd& A) {
  A = (0.5 * (A + A.transpose())).eval();
}

// Exact equality, not a tolerance: a cache that reused a factorization computed
// at a nearby-but-different parameter value would be a silently wrong gradient.
template <typename A, typename B>
inline bool blocksIdentical(const A& a, const B& b) {
  return a.rows() == b.rows() && a.cols() == b.cols() && (a.array() == b.array()).all();
}

}  // namespace detail

// Materialize `subject_values` (population values plus TI-predictor effects)
// and then the full transformed parameter vector.
inline void materializeParameters(const CppModel& model, FilterWorkspace& ws,
                                  const double* values, const double* tipreds) {
  for (int i = 0; i < model.nvalues; ++i) ws.subject_values[i] = values[i];
  for (std::size_t i = 0; i < model.tiParameter.size(); ++i) {
    ws.subject_values[model.tiParameter[i]] +=
        values[model.tiCoefficient[i]] * tipreds[model.tiPredictor[i]];
  }

  ExprContext ctx;
  ctx.params = ws.subject_values.data();
  ctx.nparams = model.nvalues;
  for (std::size_t k = 0; k < model.regPos.size(); ++k) {
    ws.all_params[model.regPos[k]] = model.regExpr[k].eval(ctx, ws.exprWork);
  }
  for (std::size_t k = 0; k < model.fixedPos.size(); ++k) {
    ws.all_params[model.fixedPos[k]] = model.fixedVal[k];
  }
}

inline void applyGroup(const CppModel& model, FilterWorkspace& ws,
                       const TransformGroup& group, const ExprContext& base) {
  if (group.empty()) return;
  ExprContext ctx = base;
  ctx.cells = ws.all_params.data();
  ctx.ncells = model.nall;
  for (std::size_t k = 0; k < group.pos.size(); ++k) {
    ws.all_params[group.pos[k]] = group.expr[k].eval(ctx, ws.exprWork);
  }
}

inline void recordGroup(AdjointTape* tape, const CppModel& model, FilterWorkspace& ws,
                        const TransformGroup& group, int groupId, const ExprContext& base,
                        int row) {
  if (tape == nullptr || group.empty()) return;
  GroupRecord& rec = tape->newGroup();
  rec.group = groupId;
  rec.params_before.resize(group.relevant.size());
  for (std::size_t j = 0; j < group.relevant.size(); ++j) {
    rec.params_before[j] = ws.all_params[group.relevant[j]];
  }
  rec.state = ws.state;
  rec.tdpreds.assign(base.tdpreds, base.tdpreds + base.ntdpred);
  rec.time = base.time;
  rec.dt = base.dt;
  rec.row = row;
}

// Continuous -> discrete for one interval.
//
// Stan propagates the local affine EKF model, not the raw DRIFT matrix:
// f(x) = DRIFT x + CINT, J = JAx, c = f(x) - J x. The Lyapunov solve and the
// discrete-intercept solve are restricted to the states with their own
// diffusion; the augmented static "carrier" states an `intoverpop` model adds
// have structurally zero diffusion and a singular Jacobian block, while their
// covariance still propagates through the full exponential.
inline void computeDiscreteTimeForm(const CppModel& model, FilterWorkspace& ws, double dt) {
  const int n = model.nlatent;
  const int k = static_cast<int>(model.diffusionStates.size());
  const std::vector<int>& dyn = model.diffusionStates;
  Eigen::Map<const MatrixXd> JAx(ws.all_params.data() + model.offJAx, n, n);
  Eigen::Map<const MatrixXd> DRIFT(ws.all_params.data() + model.offDRIFT, n, n);
  Eigen::Map<const MatrixXd> CINT(ws.all_params.data() + model.offCINT, n, 1);

  if (ws.cache.expValid && ws.cache.expDt == dt && detail::blocksIdentical(ws.cache.expJAx, JAx)) {
    ws.eJAx = ws.cache.expOut;
  } else {
    ws.Ptmp = JAx * dt;
    expm(ws.Ptmp, ws.eJAx);
    ws.cache.expJAx = JAx;
    ws.cache.expDt = dt;
    ws.cache.expOut = ws.eJAx;
    ws.cache.expValid = true;
  }

  ws.dDIFFUSION.setZero();
  ws.dINT.setZero();

  for (int i = 0; i < k; ++i) {
    const int ii = dyn[i];
    double affine = CINT(ii, 0);
    for (int j = 0; j < n; ++j) affine += (DRIFT(ii, j) - JAx(ii, j)) * ws.state(j);
    ws.affine(i) = affine;
  }
  for (int i = 0; i < k; ++i) {
    const int ii = dyn[i];
    double correction = -ws.affine(i);
    for (int q = 0; q < k; ++q) correction += ws.eJAx(ii, dyn[q]) * ws.affine(q);
    ws.sVec(i) = correction;
  }
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) ws.JAxd(i, j) = JAx(dyn[i], dyn[j]);
  ws.dINTdyn = ws.JAxd.partialPivLu().solve(ws.sVec);
  for (int i = 0; i < k; ++i) ws.dINT(dyn[i]) = ws.dINTdyn(i);

  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) ws.Qcd(i, j) = ws.Qc(dyn[i], dyn[j]);

  if (ws.cache.lyapValid && detail::blocksIdentical(ws.cache.lyapJAx, ws.JAxd) &&
      detail::blocksIdentical(ws.cache.lyapQ, ws.Qcd)) {
    ws.Xlyap = ws.cache.lyapOut;
  } else {
    lyapSolve(ws.JAxd, ws.Qcd, ws.Xlyap, ws.lyapws);
    ws.cache.lyapJAx = ws.JAxd;
    ws.cache.lyapQ = ws.Qcd;
    ws.cache.lyapOut = ws.Xlyap;
    ws.cache.lyapValid = true;
  }

  // dDIFFUSION[D,D] = X - Ad X Ad'
  for (int r = 0; r < k; ++r) {
    for (int i = 0; i < k; ++i) {
      double value = 0.0;
      for (int q = 0; q < k; ++q) value += ws.eJAx(dyn[i], dyn[q]) * ws.Xlyap(q, r);
      ws.AdX(i, r) = value;
    }
  }
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < k; ++i) {
      double value = ws.Xlyap(i, j);
      for (int r = 0; r < k; ++r) value -= ws.AdX(i, r) * ws.eJAx(dyn[j], dyn[r]);
      ws.dDIFFUSION(dyn[i], dyn[j]) = value;
    }
  }
}

inline void predictStep(const CppModel& model, FilterWorkspace& ws, double dt,
                        AdjointTape* tape) {
  const int n = model.nlatent;
  Eigen::Map<const MatrixXd> DIFFUSION(ws.all_params.data() + model.offDIFFUSION, n, n);
  sdcovsqrt2cov(DIFFUSION, n, ws.Qc, ws.covO, ws.covB);

  VectorXd stateIn;
  MatrixXd Pin;
  if (tape) { stateIn = ws.state; Pin = ws.P_update; }

  computeDiscreteTimeForm(model, ws, dt);

  ws.state = (ws.eJAx * ws.state + ws.dINT).eval();

  // Stan bakes the +1e-10 ridge into the transition itself, not just into a
  // defensive check before a Cholesky.
  for (int i = 0; i < n; ++i) ws.P_update(i, i) += kRidge;
  ws.Ptmp.noalias() = ws.eJAx * ws.P_update;
  ws.P_predict.noalias() = ws.Ptmp * ws.eJAx.transpose();
  ws.P_predict += ws.dDIFFUSION;
  detail::symmetrize(ws.P_predict);

  if (tape) {
    Eigen::Map<const MatrixXd> JAx(ws.all_params.data() + model.offJAx, n, n);
    Eigen::Map<const MatrixXd> DRIFT(ws.all_params.data() + model.offDRIFT, n, n);
    PredictRecord& rec = tape->newPredict();
    rec.state_in = stateIn;
    rec.P_in = Pin;
    rec.A = ws.eJAx;
    rec.JAx = JAx;
    rec.DRIFT = DRIFT;
    rec.DIFFUSION = DIFFUSION;
    rec.Xlyap = ws.Xlyap;
    rec.affine = ws.affine;
    rec.dINTdyn = ws.dINTdyn;
    rec.dt = dt;
  }
}

inline void applyTdImpulse(const CppModel& model, FilterWorkspace& ws,
                           const double* tdpreds, AdjointTape* tape) {
  if (model.ntdpred == 0) return;
  const int n = model.nlatent;
  Eigen::Map<const MatrixXd> TDPREDEFFECT(ws.all_params.data() + model.offTDPREDEFFECT, n,
                                          model.ntdpred);
  Eigen::Map<const MatrixXd> Jtd(ws.all_params.data() + model.offJtd, n, n);
  if (tape) {
    TdRecord& rec = tape->newTd();
    rec.P_in = ws.P_predict;
    rec.Jtd = Jtd;
    rec.tdpreds = Eigen::Map<const VectorXd>(tdpreds, model.ntdpred);
  }
  Eigen::Map<const VectorXd> td(tdpreds, model.ntdpred);
  ws.state.noalias() += TDPREDEFFECT * td;
  ws.Ptmp.noalias() = Jtd * ws.P_predict;
  ws.P_predict.noalias() = ws.Ptmp * Jtd.transpose();
}

// One measurement update restricted to the observed manifest rows. There is one
// implementation, not a separate fully-observed fast path, matching how Stan's
// own generated code always operates on the observed subset.
inline bool maskedUpdateStep(const CppModel& model, FilterWorkspace& ws,
                             const std::vector<int>& observed, const double* yrow,
                             double& loglik) {
  const int n = model.nlatent;
  const int mo = static_cast<int>(observed.size());
  Eigen::Map<const MatrixXd> LAMBDA(ws.all_params.data() + model.offLAMBDA, model.nmanifest, n);
  Eigen::Map<const MatrixXd> Jy(ws.all_params.data() + model.offJy, model.nmanifest, n);
  Eigen::Map<const MatrixXd> MANIFESTMEANS(ws.all_params.data() + model.offMANIFESTMEANS,
                                           model.nmanifest, 1);

  MatrixXd Lv(mo, n), Hv(mo, n), Rv(mo, mo);
  VectorXd yv(mo);
  for (int i = 0; i < mo; ++i) {
    const int oi = observed[i];
    for (int j = 0; j < n; ++j) { Lv(i, j) = LAMBDA(oi, j); Hv(i, j) = Jy(oi, j); }
    for (int j = 0; j < mo; ++j) Rv(i, j) = ws.Theta(oi, observed[j]);
    yv(i) = yrow[oi];
  }
  // LAMBDA (evaluated at the current state) predicts the manifest mean for the
  // innovation, while Jy (its Jacobian) propagates covariance. These coincide
  // for a fixed or linear LAMBDA but differ for a state-dependent one.
  ws.innovation.noalias() = Lv * ws.state;
  for (int i = 0; i < mo; ++i) {
    ws.innovation(i) = yv(i) - (ws.innovation(i) + MANIFESTMEANS(observed[i], 0));
  }

  ws.Pr = ws.P_predict;
  for (int i = 0; i < n; ++i) ws.Pr(i, i) += kRidge;

  ws.PHt.noalias() = ws.Pr * Hv.transpose();
  ws.S.noalias() = Hv * ws.PHt;
  ws.S += Rv;
  detail::symmetrize(ws.S);
  for (int i = 0; i < mo; ++i) ws.S(i, i) += kRidge;

  ws.llt.compute(ws.S);
  if (ws.llt.info() != Eigen::Success) return false;

  ws.alpha = ws.llt.solve(ws.innovation);
  ws.state.noalias() += ws.PHt * ws.alpha;

  // Joseph form, on the *unridged* prior covariance, matching Stan.
  ws.Kgain = ws.llt.solve(ws.PHt.transpose()).transpose();  // PHt * S^-1
  ws.Mjoseph = MatrixXd::Identity(n, n);
  ws.Mjoseph.noalias() -= ws.Kgain * Hv;
  ws.Ptmp.noalias() = ws.Mjoseph * ws.P_predict;
  ws.P_update.noalias() = ws.Ptmp * ws.Mjoseph.transpose();
  MatrixXd GR = ws.Kgain * Rv;
  ws.P_update.noalias() += GR * ws.Kgain.transpose();
  detail::symmetrize(ws.P_update);

  double logdet_half = 0.0;
  const MatrixXd& L = ws.llt.matrixLLT();
  for (int i = 0; i < mo; ++i) logdet_half += std::log(L(i, i));
  loglik = -0.5 * (mo * kLog2Pi + 2.0 * logdet_half + ws.innovation.dot(ws.alpha));
  return true;
}

// Whether a manifest entry counts as observed, matching the Julia backend's
// `!ismissing(x) && isfinite(x)` (R hands NA across as NaN).
inline bool isObserved(double x) { return std::isfinite(x); }

// One row's measurement update, including the fully-missing case.
//
// A fully missing row leaves the predicted state and covariance in place and
// contributes zero to the likelihood, matching Stan's own
// `if(si==0 || nobs_y[rowi] > 0 || dosmoother)` gate around the whole
// measurement block.
inline bool updateObserved(const CppModel& model, FilterWorkspace& ws, const double* yrow,
                           AdjointTape* tape, double& loglik) {
  const int m = model.nmanifest;
  const int n = model.nlatent;
  ws.observed.clear();
  for (int i = 0; i < m; ++i) {
    if (isObserved(yrow[i])) ws.observed.push_back(i);
  }
  if (ws.observed.empty()) {
    ws.P_update = ws.P_predict;
    loglik = 0.0;
    return true;
  }

  if (tape) {
    Eigen::Map<const MatrixXd> LAMBDA(ws.all_params.data() + model.offLAMBDA, m, n);
    Eigen::Map<const MatrixXd> Jy(ws.all_params.data() + model.offJy, m, n);
    Eigen::Map<const MatrixXd> MANIFESTMEANS(ws.all_params.data() + model.offMANIFESTMEANS, m, 1);
    const int mo = static_cast<int>(ws.observed.size());
    UpdateRecord& rec = tape->newUpdate();
    rec.observed = ws.observed;
    rec.state_in = ws.state;  // recorded before the update overwrites it
    rec.P_in = ws.P_predict;
    rec.Lambda.resize(mo, n);
    rec.H.resize(mo, n);
    rec.R.resize(mo, mo);
    rec.manifestmeans.resize(mo);
    rec.y.resize(mo);
    for (int i = 0; i < mo; ++i) {
      const int oi = ws.observed[i];
      for (int j = 0; j < n; ++j) { rec.Lambda(i, j) = LAMBDA(oi, j); rec.H(i, j) = Jy(oi, j); }
      for (int j = 0; j < mo; ++j) rec.R(i, j) = ws.Theta(oi, ws.observed[j]);
      rec.manifestmeans(i) = MANIFESTMEANS(oi, 0);
      rec.y(i) = yrow[oi];
    }
  }
  return maskedUpdateStep(model, ws, ws.observed, yrow, loglik);
}

// Evaluate one subject's log-likelihood, optionally recording the tape.
// Returns NaN for an invalid trial point (a failed innovation Cholesky), which
// poisons the total and the gradient rather than producing a finite-but-wrong
// answer the optimizer would accept.
inline double filterSubject(const CppModel& model, FilterWorkspace& ws,
                            const double* values, const SubjectData& data,
                            AdjointTape* tape) {
  const int n = model.nlatent;
  const int m = model.nmanifest;
  const int ntd = model.ntdpred;

  materializeParameters(model, ws, values, data.tipreds);
  if (tape) {
    tape->subject_values = ws.subject_values;
  }

  Eigen::Map<const MatrixXd> T0VAR(ws.all_params.data() + model.offT0VAR, n, n);
  Eigen::Map<const MatrixXd> T0MEANS(ws.all_params.data() + model.offT0MEANS, n, 1);
  sdcovsqrt2cov(T0VAR, n, ws.P_predict, ws.covO, ws.covB);
  detail::symmetrize(ws.P_predict);
  ws.state = T0MEANS.col(0);
  if (tape) tape->newInit().T0VAR = T0VAR;

  ExprContext base;
  base.cells = ws.all_params.data();
  base.ncells = model.nall;
  base.state = ws.state.data();
  base.nstate = n;
  base.tdpreds = ws.tdrow.data();
  base.ntdpred = ntd;
  base.tipreds = data.tipreds;
  base.ntipred = model.ntipred;

  auto loadTdRow = [&](int col) {
    for (int j = 0; j < ntd; ++j) ws.tdrow[j] = data.tdpreds[static_cast<std::size_t>(col) * ntd + j];
  };

  // Every row follows one contract: prediction, TD impulse, measurement. The
  // first row has no prediction interval but can still carry an impulse.
  loadTdRow(0);
  base.time = data.times[0];
  base.dt = 0.0;
  base.state = ws.state.data();
  recordGroup(tape, model, ws, model.td, 1, base, 1);
  applyGroup(model, ws, model.td, base);
  applyTdImpulse(model, ws, ws.tdrow.data(), tape);
  base.state = ws.state.data();
  recordGroup(tape, model, ws, model.update, 2, base, 1);
  applyGroup(model, ws, model.update, base);

  // Unlike the Julia backend, the manifest covariance for the first row is
  // built *after* the update group rather than before it. For a MANIFESTVAR
  // that is not state-dependent the two are identical; for one that is, the
  // Julia order reads the cell before its transform has ever written it.
  {
    Eigen::Map<const MatrixXd> MANIFESTVAR(ws.all_params.data() + model.offMANIFESTVAR, m, m);
    sdcovsqrt2cov(MANIFESTVAR, m, ws.Theta, ws.covO, ws.covB);
    if (tape) tape->newTheta().MANIFESTVAR = MANIFESTVAR;
  }

  double ll = 0.0;
  double rowll = 0.0;
  if (!updateObserved(model, ws, data.y, tape, rowll)) return std::nan("");
  ll += rowll;

  double prev = data.times[0];
  for (int t = 1; t < data.nobs; ++t) {
    const double now = data.times[t];
    const double dt = now - prev;
    loadTdRow(t);

    // Match Stan's nonlinear integration contract: re-materialize the local
    // affine model at each bounded substep before predicting.
    int nsub = 1;
    if (std::isfinite(model.maxTimestep) && model.maxTimestep > 0 && dt > model.maxTimestep) {
      nsub = static_cast<int>(std::ceil(dt / model.maxTimestep));
    }
    const double sdt = dt / nsub;
    for (int s = 1; s <= nsub; ++s) {
      base.time = prev + s * sdt;
      base.dt = sdt;
      base.state = ws.state.data();
      recordGroup(tape, model, ws, model.predict, 0, base, t + 1);
      applyGroup(model, ws, model.predict, base);
      predictStep(model, ws, sdt, tape);
      ws.P_update = ws.P_predict;  // the next bounded step starts from this one
    }

    base.time = now;
    base.dt = dt;
    base.state = ws.state.data();
    recordGroup(tape, model, ws, model.td, 1, base, t + 1);
    applyGroup(model, ws, model.td, base);
    applyTdImpulse(model, ws, ws.tdrow.data(), tape);

    base.state = ws.state.data();
    recordGroup(tape, model, ws, model.update, 2, base, t + 1);
    applyGroup(model, ws, model.update, base);
    {
      Eigen::Map<const MatrixXd> MANIFESTVAR(ws.all_params.data() + model.offMANIFESTVAR, m, m);
      sdcovsqrt2cov(MANIFESTVAR, m, ws.Theta, ws.covO, ws.covB);
      if (tape) tape->newTheta().MANIFESTVAR = MANIFESTVAR;
    }

    if (!updateObserved(model, ws, data.y + static_cast<std::size_t>(t) * m, tape, rowll)) {
      return std::nan("");
    }
    ll += rowll;
    prev = now;
  }
  return ll;
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_FILTER_HPP
