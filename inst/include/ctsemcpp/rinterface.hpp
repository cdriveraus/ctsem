#ifndef CTSEMCPP_RINTERFACE_HPP
#define CTSEMCPP_RINTERFACE_HPP

// Turn the R-side model specification into a `CppObjective`.
//
// The specification is exactly the one the Julia backend consumes -- the same
// `parameter_table`, `ti_effects`, `dynamic_state_indices`, `subject_starts`,
// `times`, `manifest_data`, `tdpred_data`, `tipred_data` and `max_timestep`
// that `.ctJuliaPrepare()` builds. Nothing on the R side had to change to feed
// this engine, which is the point: the transform strings that Julia parses and
// evaluates and that Stan compiles are the same strings interpreted here.

#include <Rcpp.h>

#include <memory>
#include <string>
#include <vector>

#include "engine.hpp"

namespace ctsemcpp {

inline std::string asStringOrEmpty(SEXP element) {
  if (element == NA_STRING) return std::string();
  return Rcpp::as<std::string>(Rcpp::CharacterVector::create(element));
}

inline std::vector<int> asIntVector(SEXP x) {
  Rcpp::IntegerVector v(x);
  std::vector<int> out(v.size());
  for (R_xlen_t i = 0; i < v.size(); ++i) out[static_cast<std::size_t>(i)] = v[i];
  return out;
}

inline std::unique_ptr<CppObjective> buildObjective(const Rcpp::List& spec) {
  using namespace Rcpp;
  const List table = spec["parameter_table"];
  CharacterVector matName = table["matrix"];
  IntegerVector matRow = table["row"];
  IntegerVector matCol = table["col"];
  IntegerVector parnumber = table["parnumber"];
  NumericVector value = table["value"];
  CharacterVector transform = table["transform"];
  CharacterVector predicttransform = table["predicttransform"];
  CharacterVector updatetransform = table["updatetransform"];
  CharacterVector tdtransform = table["tdtransform"];

  const R_xlen_t nrows = matName.size();
  std::vector<ParameterTableRow> rows(static_cast<std::size_t>(nrows));
  auto cell = [](const CharacterVector& v, R_xlen_t i) -> std::string {
    if (CharacterVector::is_na(v[i])) return std::string();
    return std::string(v[i]);
  };
  for (R_xlen_t i = 0; i < nrows; ++i) {
    ParameterTableRow& r = rows[static_cast<std::size_t>(i)];
    r.matrix = std::string(matName[i]);
    r.row = matRow[i];
    r.col = matCol[i];
    r.parnumber = IntegerVector::is_na(parnumber[i]) ? 0 : parnumber[i];
    r.hasValue = !NumericVector::is_na(value[i]);
    r.value = r.hasValue ? value[i] : 0.0;
    r.transform = cell(transform, i);
    r.predicttransform = cell(predicttransform, i);
    r.updatetransform = cell(updatetransform, i);
    r.tdtransform = cell(tdtransform, i);
  }

  const List effects = spec["ti_effects"];
  std::vector<int> tiPar, tiPred, tiCoef;
  if (effects.size() > 0) {
    tiPar = asIntVector(effects["parameter"]);
    tiPred = asIntVector(effects["predictor"]);
    tiCoef = asIntVector(effects["coefficient"]);
  }
  std::vector<int> dyn = asIntVector(spec["dynamic_state_indices"]);

  NumericMatrix tipredData = spec["tipred_data"];
  const double maxTimestep = as<double>(spec["max_timestep"]);

  auto objective = std::unique_ptr<CppObjective>(new CppObjective());
  objective->model = buildModel(rows, tiPar, tiPred, tiCoef, dyn,
                                static_cast<int>(tipredData.ncol()), maxTimestep);

  NumericMatrix manifest = spec["manifest_data"];      // manifest x observations
  NumericMatrix tdpred = spec["tdpred_data"];          // tdpred x observations
  NumericVector times = spec["times"];
  IntegerVector starts = spec["subject_starts"];

  const int m = objective->model.nmanifest;
  const int ntd = objective->model.ntdpred;
  if (manifest.nrow() != m) {
    Rcpp::stop("ctsem C++ backend: manifest data has %d rows but the model has %d manifest variables.",
               static_cast<int>(manifest.nrow()), m);
  }
  if (ntd > 0 && tdpred.nrow() != ntd) {
    Rcpp::stop("ctsem C++ backend: TD predictor data has %d rows but the model expects %d.",
               static_cast<int>(tdpred.nrow()), ntd);
  }

  const R_xlen_t nobs = times.size();
  objective->y.assign(manifest.begin(), manifest.end());
  objective->times.assign(times.begin(), times.end());
  if (ntd > 0) objective->tdpreds.assign(tdpred.begin(), tdpred.end());

  const int nsubject = static_cast<int>(starts.size());
  const int ntipred = static_cast<int>(tipredData.ncol());
  objective->subjects.resize(nsubject);
  for (int s = 0; s < nsubject; ++s) {
    const int begin = starts[s] - 1;
    const int end = (s + 1 < nsubject) ? starts[s + 1] - 1 : static_cast<int>(nobs);
    if (begin < 0 || end <= begin || end > nobs) {
      Rcpp::stop("ctsem C++ backend: subject_starts are outside the observation range.");
    }
    SubjectData& data = objective->subjects[static_cast<std::size_t>(s)];
    data.y = objective->y.data() + static_cast<std::size_t>(begin) * m;
    data.tdpreds = ntd > 0 ? objective->tdpreds.data() + static_cast<std::size_t>(begin) * ntd : nullptr;
    data.times = objective->times.data() + begin;
    data.nobs = end - begin;
    data.subject = s + 1;
  }

  // tipred_data arrives as subjects x predictors, which R stores column-major,
  // so one subject's predictors are strided; repack them contiguously.
  std::vector<double> tiBySubject(static_cast<std::size_t>(nsubject) * ntipred);
  for (int s = 0; s < nsubject; ++s) {
    for (int j = 0; j < ntipred; ++j) {
      tiBySubject[static_cast<std::size_t>(s) * ntipred + j] = tipredData(s, j);
    }
  }
  objective->tipreds.swap(tiBySubject);
  for (int s = 0; s < nsubject; ++s) {
    objective->subjects[static_cast<std::size_t>(s)].tipreds =
        ntipred > 0 ? objective->tipreds.data() + static_cast<std::size_t>(s) * ntipred : nullptr;
  }

  objective->prepare();
  return objective;
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_RINTERFACE_HPP
