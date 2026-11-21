// Identification of the disc model parameters from measured disc angle histories:
// bounded Levenberg-Marquardt on the decimated residuals of all runs.
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fr {

static DiscParams with(const DiscParams& base, const std::array<double, 4>& x, bool fric) {
    DiscParams p = base;
    p.k = x[0]; p.preload = x[1]; p.mu = x[2]; p.Mc = fric ? x[3] : 0.0; p.friction = fric;
    return p;
}

Metrics metrics(const std::vector<const Run*>& runs, const DiscParams& p, int nsub) {
    Metrics M;
    double se = 0, ae = 0, be = 0; size_t n = 0;
    std::vector<double> all;
    for (auto* r : runs) {
        auto m = simulate(r->Mf, r->dt, p, nsub);
        double s = 0, a = 0, b = 0; size_t c = 0;
        std::vector<double> meas;
        for (size_t i = 0; i < m.size(); ++i) if (r->in_soil[i]) {
            double e = m[i] - r->alpha_deg[i];
            s += e * e; a += std::fabs(e); b += e; ++c; meas.push_back(r->alpha_deg[i]);
        }
        double mm = mean(meas), ss = 0;
        for (double x : meas) ss += (x - mm) * (x - mm);
        M.per[r->name] = {std::sqrt(s / c), a / c, b / c, 1 - s / ss};
        se += s; ae += a; be += b; n += c;
        all.insert(all.end(), meas.begin(), meas.end());
    }
    double mm = mean(all), ss = 0;
    for (double x : all) ss += (x - mm) * (x - mm);
    M.rmse = std::sqrt(se / n); M.mae = ae / n; M.bias = be / n; M.r2 = 1 - se / ss;
    return M;
}

// Solve (A) x = y by Gaussian elimination with partial pivoting.
static std::vector<double> solve(std::vector<std::vector<double>> A, std::vector<double> y) {
    const int n = (int)y.size();
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
        std::swap(A[c], A[piv]); std::swap(y[c], y[piv]);
        if (std::fabs(A[c][c]) < 1e-300) throw std::runtime_error("singular matrix");
        for (int r = c + 1; r < n; ++r) {
            double f = A[r][c] / A[c][c];
            for (int k = c; k < n; ++k) A[r][k] -= f * A[c][k];
            y[r] -= f * y[c];
        }
    }
    std::vector<double> x(n);
    for (int r = n - 1; r >= 0; --r) {
        double s = y[r];
        for (int k = r + 1; k < n; ++k) s -= A[r][k] * x[k];
        x[r] = s / A[r][r];
    }
    return x;
}
static std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& A) {
    const size_t n = A.size();
    std::vector<std::vector<double>> inv(n, std::vector<double>(n));
    for (size_t j = 0; j < n; ++j) {
        std::vector<double> e(n, 0.0); e[j] = 1;
        auto c = solve(A, e);
        for (size_t i = 0; i < n; ++i) inv[i][j] = c[i];
    }
    return inv;
}

FitResult fit(const std::vector<const Run*>& runs, const FitOptions& o) {
    const int np = o.friction ? 4 : 3;
    std::vector<int> fr;
    for (int i = 0; i < np; ++i) if (!o.fixed.count(i)) fr.push_back(i);
    const int nf = (int)fr.size();
    auto full = [&](const std::vector<double>& q) {
        std::array<double, 4> x = o.x0;
        if (!o.friction) x[3] = 0;
        for (int i = 0; i < nf; ++i) x[fr[i]] = q[i];
        for (auto& kv : o.fixed) x[kv.first] = kv.second;
        return x;
    };
    auto resid = [&](const std::vector<double>& q) {
        auto p = with(o.base, full(q), o.friction);
        std::vector<double> r;
        for (auto* run : runs) {
            auto m = simulate(run->Mf, run->dt, p, o.nsub);
            // residuals at every 'step'-th sample of the soil phase (as in the published fit)
            size_t c = 0;
            for (size_t i = 0; i < m.size(); ++i) if (run->in_soil[i]) {
                if (c % o.step == 0) r.push_back(m[i] - run->alpha_deg[i]);
                ++c;
            }
        }
        return r;
    };
    auto cost = [](const std::vector<double>& r) { double s = 0; for (double v : r) s += v * v; return 0.5 * s; };
    std::vector<double> q(nf), lb(nf), ub(nf), sc(nf);
    const std::array<double, 4> xs{10, 5, 0.1, 0.5};
    for (int i = 0; i < nf; ++i) { q[i] = o.x0[fr[i]]; lb[i] = o.lb[fr[i]]; ub[i] = o.ub[fr[i]]; sc[i] = xs[fr[i]]; }
    auto r = resid(q);
    double C = cost(r), lambda = 1e-3;
    std::vector<std::vector<double>> JtJ(nf, std::vector<double>(nf));
    auto normal_eq = [&](const std::vector<double>& qc, const std::vector<double>& rc,
                         std::vector<std::vector<double>>& A, std::vector<double>& g) {
        // forward-difference Jacobian, step relative to the parameter scale
        std::vector<std::vector<double>> Jc(nf);
        for (int j = 0; j < nf; ++j) {
            auto qq = qc;
            double h = 1e-4 * sc[j];
            if (qq[j] + h > ub[j]) h = -h;
            qq[j] += h;
            auto rr = resid(qq);
            Jc[j].resize(rc.size());
            for (size_t i = 0; i < rc.size(); ++i) Jc[j][i] = (rr[i] - rc[i]) / h;
        }
        g.assign(nf, 0.0);
        for (int a = 0; a < nf; ++a) {
            for (int b = 0; b < nf; ++b) {
                double s2 = 0; for (size_t i = 0; i < rc.size(); ++i) s2 += Jc[a][i] * Jc[b][i];
                A[a][b] = s2;
            }
            double s2 = 0; for (size_t i = 0; i < rc.size(); ++i) s2 += Jc[a][i] * rc[i];
            g[a] = s2;
        }
    };
    int it = 0;
    for (; it < o.max_iter; ++it) {
        std::vector<double> g;
        normal_eq(q, r, JtJ, g);
        bool improved = false;
        for (int tries = 0; tries < 12; ++tries) {
            auto A = JtJ;
            for (int a = 0; a < nf; ++a) A[a][a] += lambda * std::max(JtJ[a][a], 1e-12);
            std::vector<double> mg(nf); for (int a = 0; a < nf; ++a) mg[a] = -g[a];
            auto d = solve(A, mg);
            auto qn = q;
            for (int a = 0; a < nf; ++a) qn[a] = std::min(std::max(q[a] + d[a], lb[a]), ub[a]);
            auto rn = resid(qn);
            double Cn = cost(rn);
            if (Cn < C) {
                double rel = (C - Cn) / C;
                double step = 0; for (int a = 0; a < nf; ++a) step = std::max(step, std::fabs(qn[a] - q[a]) / sc[a]);
                q = qn; r = rn; C = Cn; lambda = std::max(lambda / 3, 1e-9); improved = true;
                if (rel < 1e-10 || step < 1e-8) it = o.max_iter;  // converged
                break;
            }
            lambda *= 4;
        }
        if (!improved) break;
    }
    FitResult R;
    R.p = with(o.base, full(q), o.friction);
    R.cost = C; R.iterations = it; R.free_idx = fr;
    // parameter correlation from the Gauss-Newton covariance (J^T J)^-1 at the solution
    try {
        std::vector<double> g0; normal_eq(q, r, JtJ, g0);
        auto Ci = inverse(JtJ);
        R.corr.assign(nf, std::vector<double>(nf));
        for (int a = 0; a < nf; ++a) for (int b = 0; b < nf; ++b) R.corr[a][b] = Ci[a][b] / std::sqrt(Ci[a][a] * Ci[b][b]);
    } catch (...) {}
    R.m = metrics(runs, R.p, 20);
    return R;
}
}  // namespace fr
