#ifndef CTSEMCPP_MODEL_HPP
#define CTSEMCPP_MODEL_HPP

// The runtime model specification.
//
// Built once per model from exactly the same parameter table the Julia backend
// consumes (`.ctJuliaParameterTable` / `.ctJuliaAugmentRandomEffects` in
// R/ctJuliaBackend.R): one row per model-matrix cell, carrying the cell's
// coordinates, its fixed value or free-parameter number, and its transform as
// an expression string. Nothing about the model shape is compiled in.
//
// Flat layout: matrices appear in order of first appearance in the table and
// each occupies a contiguous, column-major block of `all_params`, so an
// `Eigen::Map<const MatrixXd>` over that block is the model matrix itself with
// no copy.

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "expr.hpp"

namespace ctsemcpp {

// One group of state/TD-dependent cells, re-evaluated at every row (and every
// prediction substep). The three groups mirror the Julia backend's split, and
// the split is load-bearing: `predict` cells must be current before the
// discretization, `td` before the impulse, `update` before the measurement.
struct TransformGroup {
  std::vector<int> pos;                // flat all_params index written
  std::vector<Expr> expr;
  std::vector<int> relevant;           // union of reads and writes, sorted
  bool empty() const { return pos.empty(); }
  std::size_t size() const { return pos.size(); }
};

struct CppModel {
  // --- matrix layout -------------------------------------------------------
  std::vector<MatrixLayout> layouts;
  std::map<std::string, MatrixLayout> byName;
  int nall = 0;

  int nlatent = 0;    // rows of DRIFT (the augmented state count)
  int nmanifest = 0;  // rows of LAMBDA
  int ntdpred = 0;    // columns of TDPREDEFFECT
  int ntipred = 0;
  int nvalues = 0;    // free parameters, including TI-effect coefficients

  int offT0MEANS = -1, offT0VAR = -1, offDRIFT = -1, offDIFFUSION = -1;
  int offCINT = -1, offLAMBDA = -1, offMANIFESTVAR = -1, offMANIFESTMEANS = -1;
  int offJAx = -1, offJy = -1, offJtd = -1, offTDPREDEFFECT = -1;

  // --- regular (population) transforms -------------------------------------
  std::vector<int> regPos;
  std::vector<Expr> regExpr;

  // --- fixed values --------------------------------------------------------
  std::vector<int> fixedPos;
  std::vector<double> fixedVal;

  // --- state-dependent transform groups ------------------------------------
  TransformGroup predict, td, update;

  // --- TI predictor effects (all 0-based here) ------------------------------
  std::vector<int> tiParameter, tiPredictor, tiCoefficient;

  // --- states with their own diffusion (0-based) ----------------------------
  std::vector<int> diffusionStates;

  // Whether any state-dependent transform group writes a JAx cell. The reverse
  // pass defers the matrix-exponential Frechet derivative (see `flushFrechet`),
  // and a group whose reverse zeroes a JAx cotangent is the one thing that can
  // invalidate an outstanding one.
  bool groupsWriteJAx = false;

  double maxTimestep = 1e300;

  // Whether the parameter layer is identical for every subject, so its
  // cotangent can be accumulated across subjects and unwound once instead of
  // once per subject. TI effects break it (each subject's `subject_values`
  // differs) and so do state-dependent transforms (their reverse *zeroes* the
  // cotangent on the cells they write, which would destroy another subject's
  // accumulated value). This is a correctness constraint, not conservatism.
  bool parameterLayerShareable() const {
    return tiParameter.empty() && predict.empty() && td.empty() && update.empty();
  }

  // Whether the deferred Frechet contribution can be carried across subjects
  // and pushed through the parameter layer once at the end. It can when no
  // group writes JAx (so nothing consumes the pending cotangent mid-tape) and
  // there are no TI-predictor effects (so every subject sees the same parameter
  // layer, and one pass through it covers them all). That is a weaker condition
  // than `parameterLayerShareable`, which additionally forbids *any* group --
  // and it has to be, because ctsem's default model makes MANIFESTMEANS
  // individually varying, which puts a state-dependent group between every pair
  // of prediction substeps in an otherwise entirely linear model.
  bool frechetDeferrable() const { return tiParameter.empty() && !groupsWriteJAx; }

  const MatrixLayout& layout(const std::string& name) const {
    auto it = byName.find(name);
    if (it == byName.end()) {
      throw std::runtime_error("ctsem C++ backend: model has no matrix '" + name + "'");
    }
    return it->second;
  }

  int offsetOf(const std::string& name) const { return layout(name).offset; }
};

// One row of the R-side parameter table.
struct ParameterTableRow {
  std::string matrix;
  int row = 0;
  int col = 0;
  int parnumber = 0;        // 0 when the cell is not a free parameter
  bool hasValue = false;
  double value = 0.0;
  std::string transform;
  std::string predicttransform;
  std::string updatetransform;
  std::string tdtransform;
};

namespace detail {

inline void collectRelevant(TransformGroup& group) {
  std::set<int> relevant(group.pos.begin(), group.pos.end());
  for (const Expr& e : group.expr) {
    for (int idx : e.cell_reads()) relevant.insert(idx);
  }
  group.relevant.assign(relevant.begin(), relevant.end());
}

}  // namespace detail

inline CppModel buildModel(const std::vector<ParameterTableRow>& rows,
                           const std::vector<int>& tiParameter1,
                           const std::vector<int>& tiPredictor1,
                           const std::vector<int>& tiCoefficient1,
                           const std::vector<int>& diffusionStates1,
                           int ntipred, double maxTimestep) {
  CppModel model;

  // Matrix order follows first appearance in the table, and dimensions come
  // from the largest row/column index seen, so a table that lists cells in any
  // order still yields the same layout.
  std::vector<std::string> order;
  std::map<std::string, std::pair<int, int>> dims;
  for (const ParameterTableRow& r : rows) {
    auto it = dims.find(r.matrix);
    if (it == dims.end()) {
      order.push_back(r.matrix);
      dims[r.matrix] = {r.row, r.col};
    } else {
      it->second.first = std::max(it->second.first, r.row);
      it->second.second = std::max(it->second.second, r.col);
    }
  }
  int offset = 0;
  for (const std::string& name : order) {
    MatrixLayout layout;
    layout.name = name;
    layout.nrow = dims[name].first;
    layout.ncol = dims[name].second;
    layout.offset = offset;
    offset += layout.nrow * layout.ncol;
    model.layouts.push_back(layout);
    model.byName[name] = layout;
  }
  model.nall = offset;

  model.nlatent = model.layout("DRIFT").nrow;
  model.nmanifest = model.layout("LAMBDA").nrow;
  model.offT0MEANS = model.offsetOf("T0MEANS");
  model.offT0VAR = model.offsetOf("T0VAR");
  model.offDRIFT = model.offsetOf("DRIFT");
  model.offDIFFUSION = model.offsetOf("DIFFUSION");
  model.offCINT = model.offsetOf("CINT");
  model.offLAMBDA = model.offsetOf("LAMBDA");
  model.offMANIFESTVAR = model.offsetOf("MANIFESTVAR");
  model.offMANIFESTMEANS = model.offsetOf("MANIFESTMEANS");
  model.offJAx = model.offsetOf("JAx");
  model.offJy = model.offsetOf("Jy");
  if (model.byName.count("TDPREDEFFECT")) {
    model.offTDPREDEFFECT = model.offsetOf("TDPREDEFFECT");
    model.ntdpred = model.layout("TDPREDEFFECT").ncol;
  }
  if (model.byName.count("Jtd")) model.offJtd = model.offsetOf("Jtd");
  if (model.ntdpred > 0 && model.offJtd < 0) {
    // The TD impulse propagates the covariance through Jtd, TDPREDEFFECT's
    // Jacobian, so a model with TD predictors and no Jtd is malformed rather
    // than a case to handle.
    throw std::runtime_error(
        "ctsem C++ backend: the model has TD predictors but no Jtd matrix.");
  }
  model.ntipred = ntipred;

  ExprParser parser(&model.byName);

  int maxpar = 0;
  for (const ParameterTableRow& r : rows) {
    const MatrixLayout& layout = model.layout(r.matrix);
    const int flat = layout.flat(r.row, r.col);

    if (r.parnumber > 0) {
      maxpar = std::max(maxpar, r.parnumber);
      Expr e;
      const std::string text =
          r.transform.empty() ? ("param[" + std::to_string(r.parnumber) + "]") : r.transform;
      parser.parse(text, e);
      if (!e.cell_reads().empty() || e.reads_state()) {
        throw std::runtime_error(
            "ctsem C++ backend: population transform for " + r.matrix + " reads a model "
            "matrix cell or the state; such cells belong in a predict/update/td transform.");
      }
      model.regPos.push_back(flat);
      model.regExpr.push_back(e);
    } else if (r.hasValue) {
      model.fixedPos.push_back(flat);
      model.fixedVal.push_back(r.value);
    }

    auto addComplex = [&](const std::string& text, TransformGroup& group) {
      if (text.empty()) return;
      Expr e;
      parser.parse(text, e);
      if (!e.param_reads().empty()) {
        throw std::runtime_error(
            "ctsem C++ backend: a state-dependent expression may not reference param[]; "
            "the R side renders those as population transforms.");
      }
      group.pos.push_back(flat);
      group.expr.push_back(e);
    };
    addComplex(r.predicttransform, model.predict);
    addComplex(r.updatetransform, model.update);
    addComplex(r.tdtransform, model.td);
  }

  // Within a group the transforms are applied in ascending flattened-index
  // order, which is what makes a state-dependent DRIFT cell see this row's
  // PARS value (PARS sorts first in the parameter axis) rather than the
  // previous row's. The reverse pass depends on this order, so sort here
  // rather than relying on the table's row order.
  auto sortGroup = [](TransformGroup& group) {
    std::vector<int> idx(group.pos.size());
    for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = static_cast<int>(i);
    std::stable_sort(idx.begin(), idx.end(),
                     [&](int a, int b) { return group.pos[a] < group.pos[b]; });
    std::vector<int> pos(group.pos.size());
    std::vector<Expr> expr(group.expr.size());
    for (std::size_t i = 0; i < idx.size(); ++i) {
      pos[i] = group.pos[idx[i]];
      expr[i] = group.expr[idx[i]];
    }
    group.pos.swap(pos);
    group.expr.swap(expr);
    detail::collectRelevant(group);
  };
  sortGroup(model.predict);
  sortGroup(model.td);
  sortGroup(model.update);

  {
    const MatrixLayout& jax = model.layout("JAx");
    const int lo = jax.offset, hi = jax.offset + jax.nrow * jax.ncol;
    auto writesJAx = [&](const TransformGroup& group) {
      for (int pos : group.pos) if (pos >= lo && pos < hi) return true;
      return false;
    };
    model.groupsWriteJAx = writesJAx(model.predict) || writesJAx(model.td) ||
                           writesJAx(model.update);
  }

  for (std::size_t i = 0; i < tiParameter1.size(); ++i) {
    model.tiParameter.push_back(tiParameter1[i] - 1);
    model.tiPredictor.push_back(tiPredictor1[i] - 1);
    model.tiCoefficient.push_back(tiCoefficient1[i] - 1);
    maxpar = std::max(maxpar, tiCoefficient1[i]);
  }
  model.nvalues = maxpar;

  if (diffusionStates1.empty()) {
    for (int i = 0; i < model.nlatent; ++i) model.diffusionStates.push_back(i);
  } else {
    for (int s : diffusionStates1) {
      if (s < 1 || s > model.nlatent) {
        throw std::runtime_error("ctsem C++ backend: diffusion-state index outside the latent range");
      }
      model.diffusionStates.push_back(s - 1);
    }
  }
  model.maxTimestep = maxTimestep;
  return model;
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_MODEL_HPP
