// Rcpp entry points for the hand-written C++ ctsem backend.
//
// The engine itself is header-only under `inst/include/ctsemcpp/`, so it can be
// compiled standalone (via Rcpp::sourceCpp with that directory on the include
// path) without reinstalling the package. This file is only the R boundary.

// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>

#include "ctsemcpp/rinterface.hpp"
#include "ctsemcpp/summary.hpp"

using namespace Rcpp;

namespace {

ctsemcpp::CppObjective* fromPtr(SEXP ptr) {
  XPtr<ctsemcpp::CppObjective> handle(ptr);
  if (!handle) stop("ctsem C++ backend: the objective handle is no longer valid; rebuild it.");
  return handle.get();
}

}  // namespace

// [[Rcpp::export(.ctsemCppBuild)]]
SEXP ctsemCppBuild(List spec) {
  std::unique_ptr<ctsemcpp::CppObjective> objective = ctsemcpp::buildObjective(spec);
  XPtr<ctsemcpp::CppObjective> ptr(objective.release(), true);
  return ptr;
}

// [[Rcpp::export(.ctsemCppNpars)]]
int ctsemCppNpars(SEXP handle) { return fromPtr(handle)->nvalues(); }

// [[Rcpp::export(.ctsemCppEvaluate)]]
List ctsemCppEvaluate(SEXP handle, NumericVector pars, bool gradient = true,
                      bool contributions = false) {
  ctsemcpp::CppObjective* objective = fromPtr(handle);
  const int p = objective->nvalues();
  if (pars.size() != p) {
    stop("ctsem C++ backend: expected %d free parameters, got %d.", p,
         static_cast<int>(pars.size()));
  }
  List out;
  if (gradient) {
    NumericVector grad(p);
    const double value = objective->gradient(pars.begin(), grad.begin());
    out = List::create(_["value"] = value, _["gradient"] = grad);
  } else {
    out = List::create(_["value"] = objective->value(pars.begin()),
                       _["gradient"] = R_NilValue);
  }
  if (contributions) {
    NumericVector subject(objective->subjects.size());
    for (std::size_t s = 0; s < objective->subjects.size(); ++s) {
      subject[static_cast<R_xlen_t>(s)] =
          ctsemcpp::filterSubject(objective->model, objective->fws, pars.begin(),
                                  objective->subjects[s], nullptr);
    }
    out["subject_loglik"] = subject;
  }
  return out;
}

// [[Rcpp::export(.ctsemCppOptimize)]]
List ctsemCppOptimize(SEXP handle, NumericVector start, int maxiter = 1000,
                      double gtol = 1e-8) {
  ctsemcpp::CppObjective* objective = fromPtr(handle);
  const int p = objective->nvalues();
  if (start.size() != p) {
    stop("ctsem C++ backend: expected %d starting values, got %d.", p,
         static_cast<int>(start.size()));
  }
  std::vector<double> x0(start.begin(), start.end());
  // The optimizer minimizes, so it drives the negative log-likelihood.
  auto fg = [objective, p](const double* x, double* g) {
    const double value = objective->gradient(x, g);
    for (int i = 0; i < p; ++i) g[i] = -g[i];
    return -value;
  };
  ctsemcpp::LbfgsResult result = ctsemcpp::lbfgs(fg, x0, maxiter, gtol);

  NumericVector minimizer(p), gradient(p);
  for (int i = 0; i < p; ++i) {
    minimizer[i] = result.minimizer[static_cast<std::size_t>(i)];
    gradient[i] = -result.gradient[static_cast<std::size_t>(i)];
  }
  NumericVector subject(objective->subjects.size());
  for (std::size_t s = 0; s < objective->subjects.size(); ++s) {
    subject[static_cast<R_xlen_t>(s)] =
        ctsemcpp::filterSubject(objective->model, objective->fws, minimizer.begin(),
                                objective->subjects[s], nullptr);
  }
  return List::create(_["minimizer"] = minimizer,
                      _["maximum_loglik"] = -result.value,
                      _["gradient"] = gradient,
                      _["subject_loglik"] = subject,
                      _["iterations"] = result.iterations,
                      _["converged"] = result.converged);
}

// Per-subject gradient contributions: the score matrix the OPG, sandwich and
// score-bootstrap uncertainty methods consume. Rows are subjects, columns are
// free parameters, and the rows sum to the full gradient.
// [[Rcpp::export(.ctsemCppSubjectGradients)]]
List ctsemCppSubjectGradients(SEXP handle, NumericVector pars) {
  ctsemcpp::CppObjective* objective = fromPtr(handle);
  const int p = objective->nvalues();
  if (pars.size() != p) {
    stop("ctsem C++ backend: expected %d free parameters, got %d.", p,
         static_cast<int>(pars.size()));
  }
  const int nsubject = static_cast<int>(objective->subjects.size());
  std::vector<double> rowmajor(static_cast<std::size_t>(nsubject) * p);
  const double value = objective->subjectGradients(pars.begin(), rowmajor.data());
  NumericMatrix scores(nsubject, p);
  for (int s = 0; s < nsubject; ++s) {
    for (int j = 0; j < p; ++j) {
      scores(s, j) = rowmajor[static_cast<std::size_t>(s) * p + j];
    }
  }
  return List::create(_["value"] = value, _["scores"] = scores);
}

// Diagnostic: report the flat layout the engine derived from a parameter table,
// so an R-side test can check the C++ and Julia views of a model agree without
// evaluating anything.
// [[Rcpp::export(.ctsemCppLayout)]]
List ctsemCppLayout(SEXP handle) {
  ctsemcpp::CppObjective* objective = fromPtr(handle);
  const ctsemcpp::CppModel& model = objective->model;
  CharacterVector names(model.layouts.size());
  IntegerVector nrow(model.layouts.size()), ncol(model.layouts.size()), offset(model.layouts.size());
  for (std::size_t i = 0; i < model.layouts.size(); ++i) {
    names[static_cast<R_xlen_t>(i)] = model.layouts[i].name;
    nrow[static_cast<R_xlen_t>(i)] = model.layouts[i].nrow;
    ncol[static_cast<R_xlen_t>(i)] = model.layouts[i].ncol;
    offset[static_cast<R_xlen_t>(i)] = model.layouts[i].offset;
  }
  return List::create(_["matrix"] = names, _["nrow"] = nrow, _["ncol"] = ncol,
                      _["offset"] = offset, _["nall"] = model.nall,
                      _["nlatent"] = model.nlatent, _["nmanifest"] = model.nmanifest,
                      _["ntdpred"] = model.ntdpred, _["nvalues"] = model.nvalues,
                      _["npredict"] = static_cast<int>(model.predict.size()),
                      _["nupdate"] = static_cast<int>(model.update.size()),
                      _["ntd"] = static_cast<int>(model.td.size()),
                      _["parameter_layer_shareable"] = model.parameterLayerShareable());
}

// Layout of the matrices `.ctsemCppParMatrices` returns, plus which cells are
// state dependent. Queried once per model; the R side reshapes the flat columns
// with it.
// [[Rcpp::export(.ctsemCppSummaryLayout)]]
List ctsemCppSummaryLayout(SEXP handle) {
  ctsemcpp::CppObjective* objective = fromPtr(handle);
  const ctsemcpp::CppModel& model = objective->model;
  const ctsemcpp::SummaryLayout layout = ctsemcpp::summaryLayout(model);
  const R_xlen_t nmat = static_cast<R_xlen_t>(layout.name.size());

  CharacterVector names(nmat);
  IntegerVector nrow(nmat), ncol(nmat), offset(nmat);
  for (R_xlen_t i = 0; i < nmat; ++i) {
    names[i] = layout.name[static_cast<std::size_t>(i)];
    nrow[i] = layout.nrow[static_cast<std::size_t>(i)];
    ncol[i] = layout.ncol[static_cast<std::size_t>(i)];
    offset[i] = layout.offset[static_cast<std::size_t>(i)];
  }

  // State-dependent cells reported as (matrix, row, col) rather than as flat
  // engine offsets, because the caller thinks in model matrices and should not
  // have to reimplement the engine's addressing to find out which entries of a
  // summary are conditional on a state.
  const std::vector<int> statedep = ctsemcpp::stateDependentPositions(model);
  CharacterVector sdmat(statedep.size());
  IntegerVector sdrow(statedep.size()), sdcol(statedep.size());
  for (std::size_t s = 0; s < statedep.size(); ++s) {
    const int position = statedep[s];
    for (std::size_t i = 0; i < model.layouts.size(); ++i) {
      const ctsemcpp::MatrixLayout& L = model.layouts[i];
      if (position < L.offset || position >= L.offset + L.nrow * L.ncol) continue;
      const int local = position - L.offset;
      sdmat[static_cast<R_xlen_t>(s)] = L.name;
      sdrow[static_cast<R_xlen_t>(s)] = local % L.nrow + 1;
      sdcol[static_cast<R_xlen_t>(s)] = local / L.nrow + 1;
      break;
    }
  }

  return List::create(_["matrix"] = names, _["nrow"] = nrow, _["ncol"] = ncol,
                      _["offset"] = offset, _["size"] = layout.size,
                      _["nlatent"] = model.nlatent, _["nmanifest"] = model.nmanifest,
                      _["ntipred"] = model.ntipred,
                      _["statedep"] = List::create(_["matrix"] = sdmat, _["row"] = sdrow,
                                                   _["col"] = sdcol));
}

// Materialize every model matrix for one or many raw parameter vectors.
//
// `pars` is npar x nsamples so that a whole posterior costs one call rather
// than one call per sample; the returned matrix is (flat layout) x nsamples.
// [[Rcpp::export(.ctsemCppParMatrices)]]
NumericMatrix ctsemCppParMatrices(SEXP handle, NumericMatrix pars, SEXP tipreds, SEXP state,
                                  double time = 0.0, double dt = 0.0) {
  ctsemcpp::CppObjective* objective = fromPtr(handle);
  const ctsemcpp::CppModel& model = objective->model;
  const int p = objective->nvalues();
  if (pars.nrow() != p) {
    stop("ctsem C++ backend: expected %d free parameters per column, got %d.", p,
         static_cast<int>(pars.nrow()));
  }

  std::vector<double> tivalues(static_cast<std::size_t>(std::max(model.ntipred, 0)), 0.0);
  if (!Rf_isNull(tipreds)) {
    NumericVector supplied(tipreds);
    if (supplied.size() != model.ntipred) {
      stop("ctsem C++ backend: expected %d TI predictor values, got %d.", model.ntipred,
           static_cast<int>(supplied.size()));
    }
    for (int i = 0; i < model.ntipred; ++i) tivalues[static_cast<std::size_t>(i)] = supplied[i];
  }

  std::vector<double> statevalues;
  const double* statepointer = nullptr;
  if (!Rf_isNull(state)) {
    NumericVector supplied(state);
    if (supplied.size() != model.nlatent) {
      stop("ctsem C++ backend: expected %d latent state values, got %d.", model.nlatent,
           static_cast<int>(supplied.size()));
    }
    statevalues.assign(supplied.begin(), supplied.end());
    statepointer = statevalues.data();
  }

  const ctsemcpp::SummaryLayout layout = ctsemcpp::summaryLayout(model);
  // A private workspace: this must not disturb the cached filter workspace the
  // objective uses for likelihoods and gradients.
  ctsemcpp::FilterWorkspace ws;
  ws.resize(model);

  NumericMatrix out(layout.size, pars.ncol());
  std::vector<double> column(static_cast<std::size_t>(layout.size));
  for (R_xlen_t s = 0; s < pars.ncol(); ++s) {
    NumericMatrix::Column values = pars(_, s);
    std::vector<double> raw(values.begin(), values.end());
    ctsemcpp::parameterMatrices(model, ws, layout, raw.data(), tivalues.data(), statepointer,
                                time, dt, column.data());
    for (int i = 0; i < layout.size; ++i) out(i, s) = column[static_cast<std::size_t>(i)];
  }
  return out;
}
