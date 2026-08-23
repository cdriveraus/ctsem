#ifndef CTSEMCPP_KALMAN_HPP
#define CTSEMCPP_KALMAN_HPP

// Driving the filter for prediction rather than for a likelihood.
//
// Everything here is orchestration: the arithmetic lives in `filter.hpp` (the
// forward pass and the backward RTS pass) and in `summary.hpp` (packing a
// materialized parameter vector into named matrices). That split is the point.
// ctPredict()/ctKalman() need per-row prior, filtered and smoothed estimates,
// and for a model with random effects they need each subject's own parameters;
// all of those are byproducts of the pass the likelihood already makes, and
// writing a second "prediction filter" to produce them is how the two would
// start disagreeing.

#include <algorithm>
#include <vector>

#include "summary.hpp"

namespace ctsemcpp {

// One traced pass over every subject. Returns each subject's log likelihood,
// NaN for one whose filter failed -- a failed subject leaves its smoothed rows
// untouched, so the caller has to be able to tell.
inline std::vector<double> runKalman(CppObjective& objective, const double* values,
                                     KalmanTrace& trace) {
  const CppModel& model = objective.model;
  int nrows = 0;
  for (std::size_t s = 0; s < objective.subjects.size(); ++s) {
    const SubjectData& subject = objective.subjects[s];
    nrows = std::max(nrows, subject.firstRow + subject.nobs);
  }
  trace.resize(model.nlatent, model.nmanifest, nrows,
               static_cast<int>(objective.subjects.size()), model.nall);

  // A private workspace: this must not disturb the cached one the objective
  // uses for likelihoods and gradients.
  FilterWorkspace ws;
  ws.resize(model);

  std::vector<double> loglik(objective.subjects.size(), 0.0);
  for (std::size_t s = 0; s < objective.subjects.size(); ++s) {
    loglik[s] = filterSubject(model, ws, values, objective.subjects[s], nullptr, &trace);
  }
  return loglik;
}

// One posterior-predictive dataset.
//
// `base` holds one standard normal per manifest per row; `out` receives the
// generated observations and is left untouched wherever the original data was
// missing, so the generated dataset has the same missingness as the real one.
// Returns each subject's log likelihood *of the generated data*, which is what
// a posterior predictive check compares against the observed one.
inline std::vector<double> runGenerate(CppObjective& objective, const double* values,
                                       const double* base, double* out, double* llrow) {
  const CppModel& model = objective.model;
  GenerateSpec generate;
  generate.base = base;
  generate.out = out;
  generate.llrow = llrow;

  FilterWorkspace ws;
  ws.resize(model);
  std::vector<double> loglik(objective.subjects.size(), 0.0);
  for (std::size_t s = 0; s < objective.subjects.size(); ++s) {
    loglik[s] = filterSubject(model, ws, values, objective.subjects[s], nullptr, nullptr,
                              &generate);
  }
  return loglik;
}

// One subject's model matrices: the parameter vector its filter pass ended
// with, with T0MEANS replaced by its smoothed initial state.
//
// For a model with individually varying parameters this is where the
// individual differences come out. ctsem carries such a parameter as an
// augmented latent state with no drift and no diffusion, so the smoothed t0
// estimate of that state *is* the subject's value for it, and the rest of the
// subject's matrices already reflect it because the forward pass evaluated them
// against the same carrier states at the subject's last row.
inline void subjectMatrices(const CppModel& model, FilterWorkspace& ws,
                            const SummaryLayout& layout, const KalmanTrace& trace,
                            int subject, double* out) {
  ws.all_params = trace.subjectParams[static_cast<std::size_t>(subject)];
  const VectorXd& t0 = trace.subjectT0[static_cast<std::size_t>(subject)];
  for (int i = 0; i < model.nlatent; ++i) {
    ws.all_params[static_cast<std::size_t>(model.offT0MEANS + i)] = t0(i);
  }
  packMatrices(model, ws, layout, out);
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_KALMAN_HPP
