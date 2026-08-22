#ifndef CTSEMCPP_LINALG_HPP
#define CTSEMCPP_LINALG_HPP

// Matrix kernels and their hand-written reverse-mode pullbacks.
//
// This is the C++ counterpart of `ContinuousTimeSEM/src/adjoint_primitives.jl`
// and every derivation there applies unchanged; the comments below record the
// result rather than repeating the derivation. Three kernels need hand-written
// pullbacks because their forward implementations bottom out in dense
// factorizations (LU, real Schur) rather than in differentiable arithmetic:
//
//   expm      d/dA exp(A) in direction E is the Frechet derivative L(A,E), the
//             top-right block of exp([A E; 0 A]). Its adjoint is
//             Abar = L(A', Ybar): one block exponential of twice the dimension,
//             once per reverse predict step, not once per parameter.
//
//   lyap      For A X + X A' + Q = 0: Qbar = W where A' W + W A + sym(Xbar) = 0,
//             and Abar = W X' + W' X. The pullback of a Lyapunov solve is
//             another Lyapunov solve, so it costs one forward solve.
//
//   linsolve  X = A^-1 B gives Bbar = A^-T Xbar and Abar = -Bbar X'.
//
// The fourth, `sdcovsqrt2cov`, is differentiable by ordinary AD but a naive
// Jacobian is O(d^5); the row-separable structure of the correlation constraint
// brings its pullback down to O(d^3).

#include <RcppEigen.h>
#include <unsupported/Eigen/MatrixFunctions>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace ctsemcpp {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ---------------------------------------------------------------------------
// Matrix exponential
// ---------------------------------------------------------------------------

// Eigen's scaling-and-squaring Pade implementation picks its degree adaptively
// from the matrix norm, which matters: the Julia backend measured its own
// fixed-degree-13 kernel as up to 2x slower than an adaptive one on the
// small-norm `JAx * dt` matrices this filter actually forms.
inline void expm(const MatrixXd& A, MatrixXd& out) { out = A.exp(); }

inline MatrixXd expm(const MatrixXd& A) { return A.exp(); }

// L(A, E), the Frechet derivative of exp at A in direction E.
inline MatrixXd expm_frechet(const MatrixXd& A, const MatrixXd& E) {
  const int n = static_cast<int>(A.rows());
  const double scale = E.cwiseAbs().maxCoeff();
  if (!(scale > 0.0)) return MatrixXd::Zero(n, n);
  MatrixXd block = MatrixXd::Zero(2 * n, 2 * n);
  block.topLeftCorner(n, n) = A;
  block.topRightCorner(n, n) = E / scale;   // see expmFrechetAdjoint for why
  block.bottomRightCorner(n, n) = A;
  MatrixXd expblock = block.exp();
  return expblock.topRightCorner(n, n) * scale;
}

// Abar with <Ybar, d exp(A)> = <Abar, dA>: the Frechet derivative at A' in
// direction Ybar.
inline MatrixXd expm_frechet_adjoint(const MatrixXd& A, const MatrixXd& Ybar) {
  return expm_frechet(A.transpose(), Ybar);
}

// Buffered form. The reverse pass calls this once per prediction substep, so
// the 2n x 2n block and its exponential are held in a workspace rather than
// allocated per call.
struct FrechetWorkspace {
  MatrixXd block, expblock;
  // How many block exponentials have actually been formed. Cheap, and the only
  // way to confirm from outside that the batching in `flushFrechet` is hitting.
  std::size_t calls = 0;
};

inline void expmFrechetAdjoint(const Eigen::Ref<const MatrixXd>& A,
                               const Eigen::Ref<const MatrixXd>& Ybar,
                               MatrixXd& out, FrechetWorkspace& ws) {
  const int n = static_cast<int>(A.rows());
  if (ws.block.rows() != 2 * n) ws.block.setZero(2 * n, 2 * n);

  // Normalise the direction before forming the block, and undo it afterwards.
  //
  // L(A, E) is *linear* in E, so this is exact -- but it is not cosmetic. A
  // scaling-and-squaring matrix exponential chooses its squaring count from the
  // norm of the matrix it is given, and the block `[A' E; 0 A']` inherits the
  // norm of E. Cotangents arriving here routinely have norms several orders of
  // magnitude above the drift's, so the unnormalised block was being squared
  // many extra times: measured at 608 us per prediction substep on a 20-latent
  // model (66% of the whole reverse pass) against ~30 us for every other step
  // in it. Scaling E to unit max-norm makes the block's norm that of A' alone,
  // which is what actually governs how hard this exponential is.
  const double scale = Ybar.cwiseAbs().maxCoeff();
  if (!(scale > 0.0)) { out.setZero(n, n); return; }
  ++ws.calls;

  ws.block.topLeftCorner(n, n) = A.transpose();
  ws.block.topRightCorner(n, n) = Ybar / scale;
  ws.block.bottomRightCorner(n, n) = A.transpose();
  ws.block.bottomLeftCorner(n, n).setZero();
  ws.expblock = ws.block.exp();
  out = ws.expblock.topRightCorner(n, n) * scale;
}

// ---------------------------------------------------------------------------
// Continuous Lyapunov solve:  A X + X A' + Q = 0
// ---------------------------------------------------------------------------
//
// Bartels-Stewart. Eigen ships a real Schur decomposition but no Sylvester
// back-substitution (LAPACK's `trsyl`, which the Julia path calls), so the
// block back-substitution is written out here. Blocks are 1x1 or 2x2 (a 2x2
// block is a complex-conjugate eigenvalue pair), so every small system solved
// below is at most 4x4.
//
// Reusable across calls: the Schur factorization of `A` is the expensive part
// and the reverse pass solves a second system with `A'`, so the factors are
// held in a small workspace rather than recomputed.

struct LyapWorkspace {
  Eigen::RealSchur<MatrixXd> schur;
  MatrixXd S;       // quasi-upper-triangular Schur form of A
  MatrixXd U;       // orthogonal Schur vectors
  MatrixXd C;       // right-hand side in Schur coordinates
  MatrixXd tmp;
  MatrixXd lastA;
  MatrixXd Y, rhs, At;   // scratch for the back-substitution and the pullback
  bool valid = false;
  std::vector<int> blockStart;  // 1x1/2x2 block partition of S
  std::vector<int> blockSize;

  void reset() { valid = false; }
};

namespace detail {

// Split a real Schur form into its 1x1 and 2x2 diagonal blocks.
inline void schurBlocks(const MatrixXd& S, std::vector<int>& start, std::vector<int>& size) {
  const int n = static_cast<int>(S.rows());
  start.clear();
  size.clear();
  int i = 0;
  while (i < n) {
    if (i + 1 < n && S(i + 1, i) != 0.0) {
      start.push_back(i);
      size.push_back(2);
      i += 2;
    } else {
      start.push_back(i);
      size.push_back(1);
      i += 1;
    }
  }
}

// Solve the small system  Ablk * Y + Y * Bblk = R  with Ablk (a x a) and
// Bblk (b x b), a, b in {1, 2}. Column-major vec gives
// (I_b kron Ablk + Bblk' kron I_a) vec(Y) = vec(R).
inline void solveSmallSylvester(const Eigen::Ref<const MatrixXd>& Ablk,
                                const Eigen::Ref<const MatrixXd>& Bblk,
                                const Eigen::Ref<const MatrixXd>& R,
                                Eigen::Ref<MatrixXd> Y) {
  const int a = static_cast<int>(Ablk.rows());
  const int b = static_cast<int>(Bblk.rows());
  const int d = a * b;
  Eigen::Matrix<double, 4, 4> M = Eigen::Matrix<double, 4, 4>::Zero();
  Eigen::Matrix<double, 4, 1> rhs = Eigen::Matrix<double, 4, 1>::Zero();
  for (int q = 0; q < b; ++q) {
    for (int p = 0; p < a; ++p) {
      const int row = q * a + p;
      rhs(row) = R(p, q);
      for (int k = 0; k < a; ++k) M(row, q * a + k) += Ablk(p, k);
      for (int k = 0; k < b; ++k) M(row, k * a + p) += Bblk(k, q);
    }
  }
  Eigen::Matrix<double, 4, 1> sol = Eigen::Matrix<double, 4, 1>::Zero();
  sol.head(d) = M.topLeftCorner(d, d).partialPivLu().solve(rhs.head(d));
  for (int q = 0; q < b; ++q)
    for (int p = 0; p < a; ++p) Y(p, q) = sol(q * a + p);
}

// Solve S Y + Y S' = C for Y, with S quasi-upper-triangular.
//
// The trailing-block corrections are accumulated with whole-row/column products
// rather than a loop over individual blocks: for a 20x20 system the block
// partition has ~20 members, so a per-block-pair inner loop means thousands of
// 1x1 and 2x2 Eigen products whose call overhead swamps their arithmetic.
inline void solveSchurLyapunov(const MatrixXd& S, const std::vector<int>& start,
                               const std::vector<int>& size, MatrixXd& Y,
                               const MatrixXd& C) {
  const int p = static_cast<int>(start.size());
  const int n = static_cast<int>(C.rows());
  Y.setZero(n, n);
  Eigen::Matrix<double, 2, 2> R;
  Eigen::Matrix<double, 2, 2> Yij;
  for (int bi = p - 1; bi >= 0; --bi) {
    const int i0 = start[bi], ni = size[bi];
    const int itail = i0 + ni;
    for (int bj = p - 1; bj >= 0; --bj) {
      const int j0 = start[bj], nj = size[bj];
      const int jtail = j0 + nj;
      auto Rblk = R.topLeftCorner(ni, nj);
      Rblk = C.block(i0, j0, ni, nj);
      // S is zero below its block diagonal, so only rows/columns past the
      // current block contribute, and those are already solved for.
      if (itail < n) {
        Rblk.noalias() -= S.block(i0, itail, ni, n - itail) * Y.block(itail, j0, n - itail, nj);
      }
      if (jtail < n) {
        Rblk.noalias() -= Y.block(i0, jtail, ni, n - jtail) *
                          S.block(j0, jtail, nj, n - jtail).transpose();
      }
      solveSmallSylvester(S.block(i0, i0, ni, ni), S.block(j0, j0, nj, nj).transpose(),
                          Rblk, Yij.topLeftCorner(ni, nj));
      Y.block(i0, j0, ni, nj) = Yij.topLeftCorner(ni, nj);
    }
  }
}

}  // namespace detail

// Solve A X + X A' + Q = 0. `Q` must be symmetric; `X` comes back symmetric.
inline void lyapSolve(const MatrixXd& A, const MatrixXd& Q, MatrixXd& X,
                      LyapWorkspace& ws) {
  const int n = static_cast<int>(A.rows());
  if (n == 0) { X.resize(0, 0); return; }
  if (n == 1) {
    if (A(0, 0) == 0.0) throw std::runtime_error("ctsem C++ backend: singular Lyapunov system");
    X.resize(1, 1);
    X(0, 0) = -Q(0, 0) / (2.0 * A(0, 0));
    return;
  }
  if (!ws.valid || ws.lastA.rows() != n || ws.lastA != A) {
    ws.schur.compute(A, true);
    if (ws.schur.info() != Eigen::Success) {
      throw std::runtime_error("ctsem C++ backend: Schur decomposition failed in the Lyapunov solve");
    }
    ws.S = ws.schur.matrixT();
    ws.U = ws.schur.matrixU();
    detail::schurBlocks(ws.S, ws.blockStart, ws.blockSize);
    ws.lastA = A;
    ws.valid = true;
  }
  ws.tmp.noalias() = ws.U.transpose() * Q;
  ws.C.noalias() = ws.tmp * ws.U;
  ws.C *= -1.0;
  detail::solveSchurLyapunov(ws.S, ws.blockStart, ws.blockSize, ws.Y, ws.C);
  ws.tmp.noalias() = ws.U * ws.Y;
  X.noalias() = ws.tmp * ws.U.transpose();
  // The solution is symmetric in exact arithmetic; enforce it so downstream
  // Cholesky factorizations and the reverse pass see exactly that.
  ws.Y = X;
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) {
      const double v = 0.5 * (ws.Y(i, j) + ws.Y(j, i));
      X(i, j) = v;
      X(j, i) = v;
    }
}

// Reverse of X = lyap(A, Q).
//
// The incoming cotangent is symmetrised first, and it must be: `X` is always
// symmetric, so dX is symmetric and <Xbar, dX> = <sym(Xbar), dX>; feeding a
// non-symmetric right-hand side to a solver that only spans symmetric unknowns
// returns a plausible-looking wrong answer rather than an error.
inline void lyapPullback(const MatrixXd& A, const MatrixXd& X, const MatrixXd& Xbar,
                         MatrixXd& Abar, MatrixXd& Qbar, LyapWorkspace& ws) {
  const int n = static_cast<int>(A.rows());
  if (n == 0) { Abar.resize(0, 0); Qbar.resize(0, 0); return; }
  ws.rhs = 0.5 * (Xbar + Xbar.transpose());
  ws.At = A.transpose();
  lyapSolve(ws.At, ws.rhs, Qbar, ws);
  Abar.noalias() = Qbar * X.transpose();
  Abar.noalias() += Qbar.transpose() * X;
}

// ---------------------------------------------------------------------------
// sdcovsqrt2cov: SD + unconstrained correlation square root -> covariance
// ---------------------------------------------------------------------------
//
// Forward, mirroring `constrain_cor_sqrt.jl` statement for statement:
//     v_i        = row i of mat read symmetrically from its lower triangle
//     s_i        = eps + sum_{j != i} v_j
//     ss_i       = eps + sum_{j != i} v_j^2
//     r_i        = R(s_i, ss_i)
//     O[i,j]     = v_j / r_i           (j != i)
//     O[i,i]     = sqrt(1 - sum_j O[i,j]^2 + eps)
//     B          = diag(mat) * O
//     cov        = B B'

inline double sdcovRowScale(double s, double ss) {
  const double abs_s = std::fabs(s);
  const double tmp = std::sqrt(std::log1p(std::exp(2.0 * (abs_s - s - 1.0) - 4.0)));
  return std::sqrt(ss + (tmp * (abs_s / std::sqrt(ss) - 1.0) + 1.0) * tmp + 1.0);
}

// Both partials of the row scale, by hand (the Julia port takes them from one
// two-partial dual evaluation; the closed forms are short enough to write out).
inline void sdcovRowScaleDerivs(double s, double ss, double& r, double& dr_ds, double& dr_dss) {
  const double abs_s = std::fabs(s);
  const double sign_s = (s < 0.0) ? -1.0 : 1.0;
  const double z = 2.0 * (abs_s - s - 1.0) - 4.0;
  const double ez = std::exp(z);
  const double L = std::log1p(ez);
  const double tmp = std::sqrt(L);
  const double sqrt_ss = std::sqrt(ss);
  const double g = abs_s / sqrt_ss - 1.0;
  const double inner = ss + (tmp * g + 1.0) * tmp + 1.0;
  r = std::sqrt(inner);

  const double dz_ds = 2.0 * (sign_s - 1.0);
  const double dL_ds = (ez / (1.0 + ez)) * dz_ds;
  const double dtmp_ds = (tmp > 0.0) ? 0.5 * dL_ds / tmp : 0.0;
  const double dg_ds = sign_s / sqrt_ss;
  // inner = ss + tmp^2 g + tmp + 1
  const double dinner_ds = 2.0 * tmp * dtmp_ds * g + tmp * tmp * dg_ds + dtmp_ds;
  dr_ds = 0.5 * dinner_ds / r;

  const double dg_dss = -0.5 * abs_s / (ss * sqrt_ss);
  const double dinner_dss = 1.0 + tmp * tmp * dg_dss;
  dr_dss = 0.5 * dinner_dss / r;
}

// `mat` is read symmetrically from its lower triangle, matching the Julia
// `_sym_lower_get`.
inline double symLowerGet(const Eigen::Ref<const MatrixXd>& mat, int i, int j) {
  return (j >= i) ? mat(j, i) : mat(i, j);
}

inline void constrainCorSqrt(const Eigen::Ref<const MatrixXd>& mat, int d, MatrixXd& O,
                             double epsilon = 1e-5) {
  O.resize(d, d);
  for (int i = 0; i < d; ++i) {
    double s = epsilon, ss = epsilon;
    for (int j = 0; j < d; ++j) {
      if (j == i) continue;
      const double v = symLowerGet(mat, i, j);
      s += v;
      ss += v * v;
    }
    const double r = sdcovRowScale(s, ss);
    const double inv_r = 1.0 / r;
    double sq = 0.0;
    for (int j = 0; j < d; ++j) {
      if (j == i) { O(i, j) = 0.0; continue; }
      const double o = symLowerGet(mat, i, j) * inv_r;
      O(i, j) = o;
      sq += o * o;
    }
    O(i, i) = std::sqrt(1.0 - sq + epsilon);
  }
}

inline void sdcovsqrt2cov(const Eigen::Ref<const MatrixXd>& mat, int d, MatrixXd& cov,
                          MatrixXd& scratchO, MatrixXd& scratchB, double epsilon = 1e-5) {
  constrainCorSqrt(mat, d, scratchO, epsilon);
  scratchB.resize(d, d);
  for (int j = 0; j < d; ++j)
    for (int i = 0; i < d; ++i) scratchB(i, j) = mat(i, i) * scratchO(i, j);
  cov.noalias() = scratchB * scratchB.transpose();
}

// Pullback of one row of the correlation constraint. Hand-derived; O(d) rather
// than the O(d^2) a per-input dual seeding would cost.
inline void corrSqrtRowPullback(const std::vector<double>& v, const std::vector<double>& obar,
                                int i, int d, double epsilon, std::vector<double>& vbar) {
  double s = epsilon, ss = epsilon;
  for (int j = 0; j < d; ++j) {
    if (j == i) continue;
    s += v[j];
    ss += v[j] * v[j];
  }
  double r, dr_ds, dr_dss;
  sdcovRowScaleDerivs(s, ss, r, dr_ds, dr_dss);
  const double inv_r = 1.0 / r;

  double sq = 0.0;
  for (int j = 0; j < d; ++j) {
    if (j == i) continue;
    const double o = v[j] * inv_r;
    sq += o * o;
  }
  const double diagonal = std::sqrt(1.0 - sq + epsilon);
  const double sq_bar = -obar[i] * 0.5 / diagonal;

  double r_bar = 0.0;
  vbar.assign(d, 0.0);
  for (int j = 0; j < d; ++j) {
    if (j == i) continue;
    const double o = v[j] * inv_r;
    const double o_bar = obar[j] + 2.0 * sq_bar * o;
    vbar[j] = o_bar * inv_r;
    r_bar -= o_bar * o * inv_r;
  }
  const double s_bar = r_bar * dr_ds;
  const double ss_bar = r_bar * dr_dss;
  for (int j = 0; j < d; ++j) {
    if (j == i) continue;
    vbar[j] += s_bar + 2.0 * ss_bar * v[j];
  }
}

// Accumulate into `mat_bar` the cotangent of `mat` for cov = sdcovsqrt2cov(mat).
// `mat_bar` is added to, never overwritten, so repeated uses of one parameter
// matrix across rows accumulate correctly. Only the lower triangle and diagonal
// are written, which is where the free parameters live.
inline void sdcovsqrt2covPullback(MatrixXd& mat_bar, const Eigen::Ref<const MatrixXd>& mat,
                                  const Eigen::Ref<const MatrixXd>& cov_bar, int d,
                                  double epsilon = 1e-5) {
  if (d == 0) return;
  // A zero cotangent is routine (a fully missing row leaves the manifest
  // covariance untouched) and this map is O(d^3) to pull back.
  if (cov_bar.isZero(0.0)) return;

  MatrixXd O;
  constrainCorSqrt(mat, d, O, epsilon);
  MatrixXd B(d, d);
  for (int j = 0; j < d; ++j)
    for (int i = 0; i < d; ++i) B(i, j) = mat(i, i) * O(i, j);

  // Bbar = (Cbar + Cbar') B. cov is symmetric, so this is correct whether or
  // not the incoming cotangent was symmetrised.
  MatrixXd Csym = cov_bar + cov_bar.transpose();
  MatrixXd Bbar = Csym * B;

  for (int i = 0; i < d; ++i) {
    double acc = 0.0;
    for (int j = 0; j < d; ++j) acc += Bbar(i, j) * O(i, j);
    mat_bar(i, i) += acc;
  }

  std::vector<double> v(d), obar(d), vbar(d);
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) obar[j] = mat(i, i) * Bbar(i, j);
    for (int j = 0; j < d; ++j) v[j] = symLowerGet(mat, i, j);
    corrSqrtRowPullback(v, obar, i, d, epsilon, vbar);
    for (int j = 0; j < d; ++j) {
      if (j == i) continue;
      if (j < i) mat_bar(i, j) += vbar[j];
      else mat_bar(j, i) += vbar[j];
    }
  }
}

}  // namespace ctsemcpp

#endif  // CTSEMCPP_LINALG_HPP
