#ifndef CTSEMCPP_EXPR_HPP
#define CTSEMCPP_EXPR_HPP

// Runtime expression interpreter for ctsem model-matrix cells.
//
// This is the piece that makes the C++ engine model-agnostic in the same way
// the Julia backend is: the R side already renders every free parameter's
// transform, and every state/TD-dependent matrix cell, as a short arithmetic
// expression string (see `.ctJuliaParameterTable` in R/ctJuliaBackend.R).
// Julia `Meta.parse`s those strings and `eval`s them into JIT-compiled
// closures; Stan pastes them into generated Stan code and compiles it. Here
// they are parsed once into a flat AST and interpreted.
//
// Two things keep the interpretive cost off the critical path:
//
//   1. The AST is a flat, topologically ordered array of nodes, so evaluation
//      is a single forward loop over contiguous memory with no pointer chasing
//      and no recursion.
//   2. Derivatives come from one reverse sweep over that same array, giving
//      the partials with respect to *every* input (all parameter cells and all
//      state entries) in one pass. The Julia backend instead seeds a
//      ForwardDiff dual once per input, so a transform reading `s` parameter
//      cells in an `n`-state model costs it `s + n` evaluations of the
//      expression where this costs one forward plus one reverse.
//
// Recognised leaves:
//   param[k]            free/raw parameter k (1-based, as R renders it)
//   state[i]            latent state i
//   PARS[i,j], DRIFT[i,j], ... any model matrix cell
//   tdpreds[i] / ctx.tdpreds[i]   time-dependent predictor i for the row
//   tipreds[i]          time-independent predictor i for the subject
//   time, dt            row time and elapsed interval
//
// Recognised operators: + - * / ^ and unary -/+.
// Recognised calls: exp log log1p log1p_exp sqrt square inv_logit logit
//                   tanh sin cos abs fabs pow

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ctsemcpp {

enum class ExprOp : int {
  Const = 0,
  ParamRef,   // idx = 0-based index into the raw parameter vector
  CellRef,    // idx = 0-based flat index into all_params
  StateRef,   // idx = 0-based state index
  TdRef,      // idx = 0-based TD predictor index
  TiRef,      // idx = 0-based TI predictor index
  TimeRef,
  DtRef,
  Add, Sub, Mul, Div, Pow,
  Neg,
  Exp, Log, Log1p, Log1pExp, Sqrt, Square, InvLogit, Logit, Tanh, Sin, Cos, Abs
};

struct ExprNode {
  ExprOp op;
  int a;        // first child (index into nodes), -1 when unused
  int b;        // second child, -1 when unused
  int idx;      // leaf index
  double val;   // constant value
};

// `log1p(exp(x))` in the numerically stable form Stan and the Julia backend
// both use; the transform strings ctsem renders call this by name.
inline double log1p_exp(double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

// Inputs an expression can read. Anything absent is a hard error at evaluation
// time rather than a silent zero: a transform reading a quantity the caller did
// not supply is a model-specification bug, not a value of zero.
struct ExprContext {
  const double* params = nullptr;      // raw (subject) parameter vector
  int nparams = 0;
  const double* cells = nullptr;       // all_params
  int ncells = 0;
  const double* state = nullptr;
  int nstate = 0;
  const double* tdpreds = nullptr;
  int ntdpred = 0;
  const double* tipreds = nullptr;
  int ntipred = 0;
  double time = 0.0;
  double dt = 0.0;
};

class Expr {
 public:
  Expr() = default;

  const std::vector<ExprNode>& nodes() const { return nodes_; }
  bool empty() const { return nodes_.empty(); }
  const std::string& source() const { return source_; }

  // Distinct leaf indices this expression reads, discovered structurally from
  // the parsed AST rather than by probing derivatives (a derivative that
  // happens to vanish at the probe point would hide a real dependency).
  const std::vector<int>& cell_reads() const { return cell_reads_; }
  const std::vector<int>& param_reads() const { return param_reads_; }
  bool reads_state() const { return reads_state_; }

  // Forward evaluation. `work` is resized to the node count and left holding
  // every node's value, so `backward` can reuse it without a second pass.
  double eval(const ExprContext& ctx, std::vector<double>& work) const;

  double eval(const ExprContext& ctx) const {
    std::vector<double> work;
    return eval(ctx, work);
  }

  // Reverse sweep: given the forward values in `work` (from `eval`) and an
  // output cotangent `seed`, accumulate partials into the supplied
  // accumulators. Any accumulator may be null, in which case that class of
  // input is skipped.
  void backward(const ExprContext& ctx, const std::vector<double>& work,
                double seed, std::vector<double>& adj,
                double* cell_bar, double* param_bar, double* state_bar) const;

 private:
  std::vector<ExprNode> nodes_;
  std::vector<int> cell_reads_;
  std::vector<int> param_reads_;
  bool reads_state_ = false;
  std::string source_;

  friend class ExprParser;
};

// Layout the parser needs to turn `DRIFT[2,1]` into a flat all_params index.
struct MatrixLayout {
  std::string name;
  int nrow = 0;
  int ncol = 0;
  int offset = 0;  // start of this matrix in the flat all_params vector

  int flat(int row1, int col1) const {  // 1-based row/col, column-major
    return offset + (col1 - 1) * nrow + (row1 - 1);
  }
};

class ExprParser {
 public:
  explicit ExprParser(const std::map<std::string, MatrixLayout>* layouts)
      : layouts_(layouts) {}

  void parse(const std::string& text, Expr& out) {
    text_ = &text;
    pos_ = 0;
    nodes_.clear();
    cells_.clear();
    params_.clear();
    reads_state_ = false;
    int root = parseExpression();
    skipSpace();
    if (pos_ != text.size()) {
      throw std::runtime_error("ctsem C++ backend: trailing input in expression '" +
                               text + "' at position " + std::to_string(pos_));
    }
    // The evaluator assumes the root is the last node, which the bottom-up
    // construction below guarantees; assert it rather than trusting it.
    if (root != static_cast<int>(nodes_.size()) - 1) {
      throw std::runtime_error("ctsem C++ backend: expression AST is not topologically ordered");
    }
    out.nodes_ = nodes_;
    out.source_ = text;
    dedupe(cells_);
    dedupe(params_);
    out.cell_reads_ = cells_;
    out.param_reads_ = params_;
    out.reads_state_ = reads_state_;
  }

 private:
  const std::map<std::string, MatrixLayout>* layouts_;
  const std::string* text_ = nullptr;
  std::size_t pos_ = 0;
  std::vector<ExprNode> nodes_;
  std::vector<int> cells_;
  std::vector<int> params_;
  bool reads_state_ = false;

  static void dedupe(std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }

  int push(ExprOp op, int a, int b, int idx, double val) {
    nodes_.push_back(ExprNode{op, a, b, idx, val});
    return static_cast<int>(nodes_.size()) - 1;
  }

  void skipSpace() {
    while (pos_ < text_->size() && std::isspace(static_cast<unsigned char>((*text_)[pos_]))) ++pos_;
  }

  bool consume(char c) {
    skipSpace();
    if (pos_ < text_->size() && (*text_)[pos_] == c) { ++pos_; return true; }
    return false;
  }

  void expect(char c) {
    if (!consume(c)) {
      throw std::runtime_error(std::string("ctsem C++ backend: expected '") + c +
                               "' in expression '" + *text_ + "'");
    }
  }

  [[noreturn]] void fail(const std::string& what) {
    throw std::runtime_error("ctsem C++ backend: " + what + " in expression '" + *text_ + "'");
  }

  int parseExpression() {
    int lhs = parseTerm();
    for (;;) {
      skipSpace();
      if (pos_ >= text_->size()) return lhs;
      char c = (*text_)[pos_];
      if (c == '+' || c == '-') {
        ++pos_;
        int rhs = parseTerm();
        lhs = push(c == '+' ? ExprOp::Add : ExprOp::Sub, lhs, rhs, -1, 0.0);
      } else {
        return lhs;
      }
    }
  }

  int parseTerm() {
    int lhs = parseUnary();
    for (;;) {
      skipSpace();
      if (pos_ >= text_->size()) return lhs;
      char c = (*text_)[pos_];
      if (c == '*' || c == '/') {
        ++pos_;
        int rhs = parseUnary();
        lhs = push(c == '*' ? ExprOp::Mul : ExprOp::Div, lhs, rhs, -1, 0.0);
      } else {
        return lhs;
      }
    }
  }

  int parseUnary() {
    skipSpace();
    if (pos_ < text_->size() && ((*text_)[pos_] == '-' || (*text_)[pos_] == '+')) {
      char c = (*text_)[pos_];
      ++pos_;
      int operand = parseUnary();
      if (c == '-') return push(ExprOp::Neg, operand, -1, -1, 0.0);
      return operand;
    }
    return parsePower();
  }

  int parsePower() {
    int base = parsePrimary();
    skipSpace();
    if (pos_ < text_->size() && (*text_)[pos_] == '^') {
      ++pos_;
      int exponent = parseUnary();  // right associative, and binds a unary minus
      return push(ExprOp::Pow, base, exponent, -1, 0.0);
    }
    return base;
  }

  // Index expressions in ctsem-rendered strings are always integer literals
  // (`DRIFT[2,1]`, `state[3]`). Anything else would be a model that indexes a
  // matrix by a computed value, which no ctsem backend supports.
  int parseIndex() {
    skipSpace();
    bool negative = false;
    if (pos_ < text_->size() && ((*text_)[pos_] == '-' || (*text_)[pos_] == '+')) {
      negative = (*text_)[pos_] == '-';
      ++pos_;
      skipSpace();
    }
    std::size_t start = pos_;
    while (pos_ < text_->size() && std::isdigit(static_cast<unsigned char>((*text_)[pos_]))) ++pos_;
    if (pos_ == start) fail("expected an integer index");
    int value = std::stoi(text_->substr(start, pos_ - start));
    return negative ? -value : value;
  }

  int parsePrimary() {
    skipSpace();
    if (pos_ >= text_->size()) fail("unexpected end of expression");
    char c = (*text_)[pos_];

    if (c == '(') {
      ++pos_;
      int inner = parseExpression();
      expect(')');
      return inner;
    }

    if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
      std::size_t start = pos_;
      while (pos_ < text_->size() &&
             (std::isdigit(static_cast<unsigned char>((*text_)[pos_])) || (*text_)[pos_] == '.')) ++pos_;
      if (pos_ < text_->size() && ((*text_)[pos_] == 'e' || (*text_)[pos_] == 'E')) {
        std::size_t save = pos_;
        ++pos_;
        if (pos_ < text_->size() && ((*text_)[pos_] == '+' || (*text_)[pos_] == '-')) ++pos_;
        if (pos_ < text_->size() && std::isdigit(static_cast<unsigned char>((*text_)[pos_]))) {
          while (pos_ < text_->size() && std::isdigit(static_cast<unsigned char>((*text_)[pos_]))) ++pos_;
        } else {
          pos_ = save;  // a trailing 'e' that is not an exponent
        }
      }
      return push(ExprOp::Const, -1, -1, -1, std::stod(text_->substr(start, pos_ - start)));
    }

    if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
      std::size_t start = pos_;
      while (pos_ < text_->size() &&
             (std::isalnum(static_cast<unsigned char>((*text_)[pos_])) ||
              (*text_)[pos_] == '_' || (*text_)[pos_] == '.')) ++pos_;
      std::string name = text_->substr(start, pos_ - start);
      // `.ctJuliaTDExpression` rewrites Stan-style `tdpreds[rowi,j]` to the
      // Julia context form; accept either spelling so the same table feeds
      // both engines unchanged.
      if (name.rfind("ctx.", 0) == 0) name = name.substr(4);
      if (name == "pars") {  // `ctx.pars.DRIFT[...]`, if ever handed to us
        if (pos_ < text_->size() && (*text_)[pos_] == '.') {
          ++pos_;
          std::size_t s2 = pos_;
          while (pos_ < text_->size() &&
                 (std::isalnum(static_cast<unsigned char>((*text_)[pos_])) || (*text_)[pos_] == '_')) ++pos_;
          name = text_->substr(s2, pos_ - s2);
        }
      }
      skipSpace();
      if (pos_ < text_->size() && (*text_)[pos_] == '(') {
        ++pos_;
        return parseCall(name);
      }
      if (pos_ < text_->size() && (*text_)[pos_] == '[') {
        ++pos_;
        return parseIndexed(name);
      }
      if (name == "time") return push(ExprOp::TimeRef, -1, -1, -1, 0.0);
      if (name == "dt") return push(ExprOp::DtRef, -1, -1, -1, 0.0);
      if (name == "pi") return push(ExprOp::Const, -1, -1, -1, 3.14159265358979323846);
      fail("unknown symbol '" + name + "'");
    }

    fail(std::string("unexpected character '") + c + "'");
  }

  int parseIndexed(const std::string& name) {
    int i = parseIndex();
    int j = 1;
    bool twod = false;
    skipSpace();
    if (pos_ < text_->size() && (*text_)[pos_] == ',') {
      ++pos_;
      j = parseIndex();
      twod = true;
    }
    expect(']');

    if (name == "param") {
      if (twod) fail("param[] takes one index");
      params_.push_back(i - 1);
      return push(ExprOp::ParamRef, -1, -1, i - 1, 0.0);
    }
    if (name == "state") {
      // ctsem sometimes renders the state as a column matrix (`state[i,1]`).
      if (twod && j != 1) fail("state[] must be indexed as state[i] or state[i,1]");
      reads_state_ = true;
      return push(ExprOp::StateRef, -1, -1, i - 1, 0.0);
    }
    if (name == "tdpreds") {
      if (twod) { i = j; }  // `tdpreds[rowi, k]` -> predictor k
      return push(ExprOp::TdRef, -1, -1, i - 1, 0.0);
    }
    if (name == "tipreds") {
      if (twod) { i = j; }
      return push(ExprOp::TiRef, -1, -1, i - 1, 0.0);
    }

    auto it = layouts_->find(name);
    if (it == layouts_->end()) fail("unknown matrix '" + name + "'");
    const MatrixLayout& layout = it->second;
    if (i < 1 || i > layout.nrow || j < 1 || j > layout.ncol) {
      fail("index out of range for matrix '" + name + "'");
    }
    int flat = layout.flat(i, j);
    cells_.push_back(flat);
    return push(ExprOp::CellRef, -1, -1, flat, 0.0);
  }

  int parseCall(const std::string& name) {
    std::vector<int> args;
    skipSpace();
    if (pos_ < text_->size() && (*text_)[pos_] == ')') {
      ++pos_;
    } else {
      for (;;) {
        args.push_back(parseExpression());
        skipSpace();
        if (consume(',')) continue;
        expect(')');
        break;
      }
    }

    auto unary = [&](ExprOp op) {
      if (args.size() != 1) fail("function '" + name + "' takes one argument");
      return push(op, args[0], -1, -1, 0.0);
    };

    if (name == "exp") return unary(ExprOp::Exp);
    if (name == "log") return unary(ExprOp::Log);
    if (name == "log1p") return unary(ExprOp::Log1p);
    if (name == "log1p_exp") return unary(ExprOp::Log1pExp);
    if (name == "sqrt") return unary(ExprOp::Sqrt);
    if (name == "square") return unary(ExprOp::Square);
    if (name == "inv_logit") return unary(ExprOp::InvLogit);
    if (name == "logit") return unary(ExprOp::Logit);
    if (name == "tanh") return unary(ExprOp::Tanh);
    if (name == "sin") return unary(ExprOp::Sin);
    if (name == "cos") return unary(ExprOp::Cos);
    if (name == "abs" || name == "fabs") return unary(ExprOp::Abs);
    if (name == "pow") {
      if (args.size() != 2) fail("pow() takes two arguments");
      return push(ExprOp::Pow, args[0], args[1], -1, 0.0);
    }
    fail("unknown function '" + name + "'");
  }
};

inline double Expr::eval(const ExprContext& ctx, std::vector<double>& work) const {
  const std::size_t n = nodes_.size();
  if (n == 0) return 0.0;
  work.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    const ExprNode& node = nodes_[k];
    double v = 0.0;
    switch (node.op) {
      case ExprOp::Const: v = node.val; break;
      case ExprOp::ParamRef:
        if (node.idx >= ctx.nparams) throw std::runtime_error("ctsem C++ backend: param index out of range");
        v = ctx.params[node.idx];
        break;
      case ExprOp::CellRef:
        if (node.idx >= ctx.ncells) throw std::runtime_error("ctsem C++ backend: matrix cell index out of range");
        v = ctx.cells[node.idx];
        break;
      case ExprOp::StateRef:
        if (node.idx >= ctx.nstate) throw std::runtime_error("ctsem C++ backend: state index out of range");
        v = ctx.state[node.idx];
        break;
      case ExprOp::TdRef:
        if (node.idx >= ctx.ntdpred) throw std::runtime_error("ctsem C++ backend: TD predictor index out of range");
        v = ctx.tdpreds[node.idx];
        break;
      case ExprOp::TiRef:
        if (node.idx >= ctx.ntipred) throw std::runtime_error("ctsem C++ backend: TI predictor index out of range");
        v = ctx.tipreds[node.idx];
        break;
      case ExprOp::TimeRef: v = ctx.time; break;
      case ExprOp::DtRef: v = ctx.dt; break;
      case ExprOp::Add: v = work[node.a] + work[node.b]; break;
      case ExprOp::Sub: v = work[node.a] - work[node.b]; break;
      case ExprOp::Mul: v = work[node.a] * work[node.b]; break;
      case ExprOp::Div: v = work[node.a] / work[node.b]; break;
      case ExprOp::Pow: v = std::pow(work[node.a], work[node.b]); break;
      case ExprOp::Neg: v = -work[node.a]; break;
      case ExprOp::Exp: v = std::exp(work[node.a]); break;
      case ExprOp::Log: v = std::log(work[node.a]); break;
      case ExprOp::Log1p: v = std::log1p(work[node.a]); break;
      case ExprOp::Log1pExp: v = log1p_exp(work[node.a]); break;
      case ExprOp::Sqrt: v = std::sqrt(work[node.a]); break;
      case ExprOp::Square: v = work[node.a] * work[node.a]; break;
      case ExprOp::InvLogit: v = 1.0 / (1.0 + std::exp(-work[node.a])); break;
      case ExprOp::Logit: v = std::log(work[node.a] / (1.0 - work[node.a])); break;
      case ExprOp::Tanh: v = std::tanh(work[node.a]); break;
      case ExprOp::Sin: v = std::sin(work[node.a]); break;
      case ExprOp::Cos: v = std::cos(work[node.a]); break;
      case ExprOp::Abs: v = std::fabs(work[node.a]); break;
    }
    work[k] = v;
  }
  return work[n - 1];
}

inline void Expr::backward(const ExprContext& ctx, const std::vector<double>& work,
                           double seed, std::vector<double>& adj,
                           double* cell_bar, double* param_bar,
                           double* state_bar) const {
  const std::size_t n = nodes_.size();
  if (n == 0 || seed == 0.0) return;
  adj.assign(n, 0.0);
  adj[n - 1] = seed;
  for (std::size_t k = n; k-- > 0;) {
    const double g = adj[k];
    if (g == 0.0) continue;
    const ExprNode& node = nodes_[k];
    switch (node.op) {
      case ExprOp::Const:
      case ExprOp::TimeRef:
      case ExprOp::DtRef:
      case ExprOp::TdRef:
      case ExprOp::TiRef:
        break;
      case ExprOp::ParamRef: if (param_bar) param_bar[node.idx] += g; break;
      case ExprOp::CellRef:  if (cell_bar) cell_bar[node.idx] += g; break;
      case ExprOp::StateRef: if (state_bar) state_bar[node.idx] += g; break;
      case ExprOp::Add: adj[node.a] += g; adj[node.b] += g; break;
      case ExprOp::Sub: adj[node.a] += g; adj[node.b] -= g; break;
      case ExprOp::Mul: adj[node.a] += g * work[node.b]; adj[node.b] += g * work[node.a]; break;
      case ExprOp::Div: {
        const double denom = work[node.b];
        adj[node.a] += g / denom;
        adj[node.b] -= g * work[node.a] / (denom * denom);
        break;
      }
      case ExprOp::Pow: {
        const double base = work[node.a];
        const double power = work[node.b];
        adj[node.a] += g * power * std::pow(base, power - 1.0);
        // Only pay for the log when the exponent is genuinely a variable; a
        // constant exponent (the common `param^3` case) would otherwise take
        // log() of a negative base and produce a NaN in a dead branch.
        if (nodes_[node.b].op != ExprOp::Const) {
          adj[node.b] += g * work[k] * std::log(base);
        }
        break;
      }
      case ExprOp::Neg: adj[node.a] -= g; break;
      case ExprOp::Exp: adj[node.a] += g * work[k]; break;
      case ExprOp::Log: adj[node.a] += g / work[node.a]; break;
      case ExprOp::Log1p: adj[node.a] += g / (1.0 + work[node.a]); break;
      case ExprOp::Log1pExp: {
        const double x = work[node.a];
        adj[node.a] += g * (x > 0.0 ? 1.0 / (1.0 + std::exp(-x)) : std::exp(x) / (1.0 + std::exp(x)));
        break;
      }
      case ExprOp::Sqrt: adj[node.a] += g * 0.5 / work[k]; break;
      case ExprOp::Square: adj[node.a] += g * 2.0 * work[node.a]; break;
      case ExprOp::InvLogit: adj[node.a] += g * work[k] * (1.0 - work[k]); break;
      case ExprOp::Logit: {
        const double x = work[node.a];
        adj[node.a] += g / (x * (1.0 - x));
        break;
      }
      case ExprOp::Tanh: adj[node.a] += g * (1.0 - work[k] * work[k]); break;
      case ExprOp::Sin: adj[node.a] += g * std::cos(work[node.a]); break;
      case ExprOp::Cos: adj[node.a] -= g * std::sin(work[node.a]); break;
      case ExprOp::Abs: adj[node.a] += g * (work[node.a] < 0.0 ? -1.0 : 1.0); break;
    }
  }
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_EXPR_HPP
