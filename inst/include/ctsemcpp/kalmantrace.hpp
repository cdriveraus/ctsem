#ifndef CTSEMCPP_KALMANTRACE_HPP
#define CTSEMCPP_KALMANTRACE_HPP

// Per-row Kalman output, for ctKalman()/ctPredict() rather than for the
// likelihood.
//
// This is a plain data sink, deliberately, so that `filter.hpp` can take a
// pointer to one without depending on the smoother or on anything R-facing.
// It is the second nullable out-parameter the filter carries, alongside the
// adjoint tape, and for the same reason: there is one forward loop in this
// engine, and a second copy of it written "for prediction" is exactly how the
// prediction output and the likelihood drift apart.
//
// Everything here is indexed by *global* data row, not by position within a
// subject, so a caller allocates once for the whole dataset and each subject
// writes its own slice. Rows are grouped by kind:
//
//   Prior  -- state and observation before this row's measurement update
//   Upd    -- after it (the filtered estimate)
//   Smooth -- after the backward RTS pass over the subject
//
// The memory here is O(nrows * nlatentpop^2), which is why it is opt-in: this
// is a diagnostic path, not the one an optimizer iterates.

#include <RcppEigen.h>

#include <vector>

namespace ctsemcpp {

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct KalmanTrace {
  enum Kind { Prior = 0, Upd = 1, Smooth = 2, nKinds = 3 };

  int nlatent = 0;
  int nmanifest = 0;
  int nrows = 0;
  int nsubjects = 0;

  // [kind * nrows + row]
  std::vector<VectorXd> eta, y;
  std::vector<MatrixXd> etacov, ycov;
  std::vector<double> llrow;

  // Per row, for the backward pass: the state transition that *reached* this
  // row, and the measurement Jacobian used at it.
  std::vector<MatrixXd> transition, Jy;

  // Per subject: the materialized parameter vector as of the subject's last
  // row, and the smoothed latent state at its first row. Together these are
  // what individual-difference parameter estimates are read from -- in an
  // intoverpop model the random effects *are* augmented latent states, so a
  // subject's parameters are its smoothed t0 state.
  std::vector<std::vector<double> > subjectParams;
  std::vector<VectorXd> subjectT0;

  // Which rows a subject actually wrote, so a caller can tell a genuinely
  // filtered row from an untouched allocation.
  std::vector<int> rowSubject;

  void resize(int nlatent1, int nmanifest1, int nrows1, int nsubjects1, int nall) {
    nlatent = nlatent1;
    nmanifest = nmanifest1;
    nrows = nrows1;
    nsubjects = nsubjects1;
    const std::size_t total = static_cast<std::size_t>(nKinds) * nrows;
    eta.assign(total, VectorXd::Zero(nlatent));
    y.assign(total, VectorXd::Zero(nmanifest));
    etacov.assign(total, MatrixXd::Zero(nlatent, nlatent));
    ycov.assign(total, MatrixXd::Zero(nmanifest, nmanifest));
    llrow.assign(static_cast<std::size_t>(nrows), 0.0);
    transition.assign(static_cast<std::size_t>(nrows), MatrixXd::Identity(nlatent, nlatent));
    Jy.assign(static_cast<std::size_t>(nrows), MatrixXd::Zero(nmanifest, nlatent));
    subjectParams.assign(static_cast<std::size_t>(nsubjects),
                         std::vector<double>(static_cast<std::size_t>(nall), 0.0));
    subjectT0.assign(static_cast<std::size_t>(nsubjects), VectorXd::Zero(nlatent));
    rowSubject.assign(static_cast<std::size_t>(nrows), 0);
  }

  std::size_t at(int kind, int row) const {
    return static_cast<std::size_t>(kind) * nrows + row;
  }
};

// Posterior-predictive data generation: the same filter, drawing each row's
// observation from its own prior predictive instead of reading it.
//
// This is not a recorder -- it changes what the filter consumes, so the state
// it carries forward is conditioned on the drawn data rather than the real
// data, which is exactly what makes the result a draw from the model. It rides
// on the same pass for the same reason the recorders do: a separate "simulator"
// would have to re-derive the prior predictive at every row, and would be free
// to disagree with the filter about it.
//
// The standard normal draws come from the caller rather than from an engine
// RNG, so that a seed set in R gives the same data whichever engine ran it, and
// so that `set.seed()` means what a user expects.
//
// `base` and `out` are nmanifest-by-nrows in the same layout as the observed
// data. Missing entries are never generated: the generated dataset keeps the
// original missingness, since the point of it is comparison against the
// observations that are actually there.
struct GenerateSpec {
  const double* base = nullptr;   // standard normal draws, one per manifest per row
  double* out = nullptr;          // generated observations; untouched where missing
  double* llrow = nullptr;        // per-row log likelihood of the generated data
};

}  // namespace ctsemcpp

#endif  // CTSEMCPP_KALMANTRACE_HPP
