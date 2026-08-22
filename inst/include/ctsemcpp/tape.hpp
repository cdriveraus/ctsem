#ifndef CTSEMCPP_TAPE_HPP
#define CTSEMCPP_TAPE_HPP

// The adjoint tape.
//
// Deliberately an ordered program of small typed records rather than a
// fixed-shape array of checkpoints, for the same reason the Julia backend's is:
// the forward pass has no canonical shape. The number of prediction substeps
// depends on `max_timestep` and the observation spacing, TD impulses only occur
// when there are TD predictors, state-dependent transform groups only exist for
// nonlinear models, and a fully missing row skips the measurement update
// entirely. Recording the order the forward pass actually took is what lets one
// reverse pass cover every case without special-casing any of them.
//
// Records are pooled and reused across subjects and across gradient
// evaluations: `reset()` rewinds the counters without freeing the Eigen
// storage, so a steady-state gradient evaluation allocates nothing here.

#include <RcppEigen.h>

#include <vector>

namespace ctsemcpp {

using Eigen::MatrixXd;
using Eigen::VectorXd;

enum class TapeKind : int { Init = 0, Theta, Group, Predict, Td, Update };

// One prediction substep: everything the reverse pass needs to undo it.
struct PredictRecord {
  VectorXd state_in;    // state entering the substep
  MatrixXd P_in;        // posterior covariance entering the substep (unridged)
  MatrixXd A;           // eJAx = exp(JAx * dt), which is also the discrete drift
  MatrixXd JAx;
  MatrixXd DRIFT;
  MatrixXd DIFFUSION;   // raw SD/correlation-sqrt parameters
  MatrixXd Xlyap;       // Lyapunov solution on the dynamic block
  VectorXd affine;      // the local affine correction, `r` in the derivation
  VectorXd dINTdyn;     // solved discrete intercept on the dynamic block
  double dt = 0.0;
};

struct TdRecord {
  MatrixXd P_in;
  MatrixXd Jtd;
  VectorXd tdpreds;
};

// One measurement update. Stores only inputs: the innovation, gain and Cholesky
// factor are cheap to recompute and doing so keeps the primal free of tracing
// code (it reuses its gain buffer as scratch partway through, so not all of its
// intermediates are live at the end anyway).
struct UpdateRecord {
  std::vector<int> observed;
  VectorXd state_in;
  MatrixXd P_in;          // prior covariance, unridged
  MatrixXd Lambda;        // LAMBDA[observed, :]
  VectorXd manifestmeans; // MANIFESTMEANS[observed]
  MatrixXd H;             // Jy[observed, :]
  MatrixXd R;             // Theta[observed, observed]
  VectorXd y;             // observed data for this row
};

// One application of a group of state-dependent transforms.
//
// `params_before` holds the parameter values from *before* the group ran, and
// both halves of that matter. Within a group the transforms are applied in
// ascending flattened-index order into the same buffer they read from, so a
// transform can read a cell a later transform in the same group is about to
// overwrite -- and when it does, it sees the previous row's value. Recording
// the pre-group state lets the reverse pass reconstruct, for each transform,
// the values that transform actually saw. Only the group's relevant cells are
// stored; everything else in the parameter vector is unreachable from the
// group's expressions.
struct GroupRecord {
  int group = 0;  // 0 = predict, 1 = td, 2 = update
  std::vector<double> params_before;
  VectorXd state;
  std::vector<double> tdpreds;
  double time = 0.0;
  double dt = 0.0;
  int row = 0;
};

struct ThetaRecord { MatrixXd MANIFESTVAR; };
struct InitRecord { MatrixXd T0VAR; };

struct AdjointTape {
  std::vector<std::pair<TapeKind, int>> program;
  std::vector<PredictRecord> predicts;
  std::vector<TdRecord> tds;
  std::vector<UpdateRecord> updates;
  std::vector<GroupRecord> groups;
  std::vector<ThetaRecord> thetas;
  std::vector<InitRecord> inits;
  std::vector<double> subject_values;

  int nPredict = 0, nTd = 0, nUpdate = 0, nGroup = 0, nTheta = 0, nInit = 0;

  void reset() {
    program.clear();
    nPredict = nTd = nUpdate = nGroup = nTheta = nInit = 0;
  }

  template <typename T>
  static T& take(std::vector<T>& pool, int& counter) {
    if (counter == static_cast<int>(pool.size())) pool.emplace_back();
    return pool[counter++];
  }

  PredictRecord& newPredict() { program.emplace_back(TapeKind::Predict, nPredict); return take(predicts, nPredict); }
  TdRecord& newTd() { program.emplace_back(TapeKind::Td, nTd); return take(tds, nTd); }
  UpdateRecord& newUpdate() { program.emplace_back(TapeKind::Update, nUpdate); return take(updates, nUpdate); }
  GroupRecord& newGroup() { program.emplace_back(TapeKind::Group, nGroup); return take(groups, nGroup); }
  ThetaRecord& newTheta() { program.emplace_back(TapeKind::Theta, nTheta); return take(thetas, nTheta); }
  InitRecord& newInit() { program.emplace_back(TapeKind::Init, nInit); return take(inits, nInit); }
};

}  // namespace ctsemcpp

#endif  // CTSEMCPP_TAPE_HPP
