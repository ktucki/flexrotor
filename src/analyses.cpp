// Order analysis, tine kinematics, spring selection maps, mean-load model and the
// coupling selection rule.
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <stdexcept>

namespace fr {
using cd = std::complex<double>;

// ------------------------------------------------------------------ order analysis
OrderResult order_analysis(const Run& r, const DiscParams& p, int n) {
    OrderResult o;
    o.name = r.name;
    o.rpm = r.omega * 60 / (2 * PI);
    o.Omega = n * r.omega;
    // soil phase as used for identification: middle 60% of the samples with alpha > 1 deg
    size_t s0 = r.alpha_deg.size(), cnt = 0;
    for (size_t i = 0; i < r.alpha_deg.size(); ++i) if (r.in_soil[i]) { s0 = std::min(s0, i); ++cnt; }
    if (!cnt) throw std::runtime_error("order_analysis: no soil phase");
    size_t i0 = s0 + (size_t)(0.2 * cnt), i1 = s0 + (size_t)(0.8 * cnt);
    const auto& th = r.theta_rad;
    const double th0 = th[i0];
    o.nrev = (int)((th[i1] - th0) / (2 * PI));
    std::vector<char> m(th.size(), 0);
    for (size_t i = i0; i < th.size(); ++i) m[i] = th[i] >= th0 && th[i] < th0 + o.nrev * 2 * PI;
    auto coef = [&](const std::vector<double>& x, int ord) {
        double mu = 0; size_t c = 0;
        for (size_t i = 0; i < x.size(); ++i) if (m[i]) { mu += x[i]; ++c; }
        mu /= c;
        cd s = 0;
        for (size_t i = 0; i < x.size(); ++i) if (m[i]) s += (x[i] - mu) * std::exp(cd(0, -ord * th[i]));
        return 2.0 * s / double(c);
    };
    const double J = p.J, kr = k_si(p.k);
    const cd a_n = coef(r.alpha_deg, n) * PI / 180.0;
    const cd Mt_n = coef(r.Mg, n);
    const cd Mg_n = Mt_n - J * o.Omega * o.Omega * a_n;  // J a'' = Mg - Mt
    const cd C = 1.0 / (kr - J * o.Omega * o.Omega + cd(0, o.Omega * p.mu));
    o.Mt_n = std::abs(Mt_n); o.Mg_n = std::abs(Mg_n);
    o.alpha_n_deg = rad2deg(std::abs(a_n));
    o.alpha_n_lin_deg = rad2deg(std::abs(Mg_n * C));
    o.T_exp = std::abs(Mt_n / Mg_n);
    o.T_lin = std::abs((kr + cd(0, o.Omega * p.mu)) * C);
    return o;
}

// ------------------------------------------------------------------ kinematics
KinResult kinematics(const Run& r, const DiscParams& p, int z, double radius, double u) {
    KinResult k;
    k.name = r.name;
    const double om = r.omega, Tz = 2 * PI / (z * om);
    const int nz = (int)std::lround(Tz / r.dt);
    k.lambda = radius * om / u;
    k.k_crit = 1 - 1 / k.lambda;
    std::vector<double> da, f;
    for (size_t i = 0; i + nz < r.alpha_deg.size(); ++i) if (r.in_soil[i] && r.in_soil[i + nz]) {
        double d = deg2rad(r.alpha_deg[i + nz] - r.alpha_deg[i]);
        da.push_back(std::fabs(d));
        f.push_back((2 * PI / z) / (2 * PI / z - d));
    }
    if (da.empty()) throw std::runtime_error("kinematics: no soil phase");
    k.dalpha_max = rad2deg(*std::max_element(da.begin(), da.end()));
    k.dalpha_p95 = rad2deg(percentile(da, 95));
    k.f_min = *std::min_element(f.begin(), f.end());
    k.f_max = *std::max_element(f.begin(), f.end());
    k.f_p025 = percentile(f, 2.5);
    k.f_p975 = percentile(f, 97.5);
    // power dissipated in the coupling vs power absorbed by the soil
    auto a = simulate(r.Mf, r.dt, p, 20);
    std::vector<double> ar(a.size());
    for (size_t i = 0; i < a.size(); ++i) ar[i] = deg2rad(a[i]);
    double Pd = 0, Pg = 0; size_t c = 0;
    for (size_t i = 0; i < ar.size(); ++i) {
        size_t lo = i ? i - 1 : 0, hi = i + 1 < ar.size() ? i + 1 : i;
        double v = (ar[hi] - ar[lo]) / ((hi - lo) * r.dt);
        if (!r.in_soil[i]) continue;
        Pd += p.mu * v * v + p.Mc * std::fabs(std::tanh(v / p.v0)) * std::fabs(v);
        Pg += r.Mf[i] * om; ++c;
    }
    k.Pd = Pd / c; k.Pg = Pg / c;
    return k;
}

// ------------------------------------------------------------------ selection map
SelectionMap selection_map(const std::vector<const Run*>& runs, const DiscParams& p0,
                           const std::vector<double>& ks, const std::vector<double>& Ms, int nsub) {
    SelectionMap S;
    S.ks = ks; S.Ms = Ms;
    S.low.assign(Ms.size(), std::vector<double>(ks.size()));
    S.mid = S.up = S.low;
    for (size_t j = 0; j < ks.size(); ++j)
        for (size_t i = 0; i < Ms.size(); ++i) {
            DiscParams p = p0;
            p.k = ks[j];
            p.preload = rad2deg(Ms[i] / k_si(ks[j]));
            double nl = 0, nm = 0, nu = 0;
            for (auto* r : runs) {
                auto a = simulate(r->Mf, r->dt, p, nsub);
                for (size_t t = 0; t < a.size(); ++t) if (r->in_soil[t]) {
                    if (a[t] < 0.5) ++nl; else if (a[t] > p.amax - 0.5) ++nu; else ++nm;
                }
            }
            double tot = nl + nm + nu;
            S.low[i][j] = nl / tot; S.mid[i][j] = nm / tot; S.up[i][j] = nu / tot;
        }
    return S;
}

// ------------------------------------------------------------------ mean-load model
static std::vector<double> lstsq(const std::vector<std::vector<double>>& X, const std::vector<double>& y) {
    const size_t p = X[0].size();
    std::vector<std::vector<double>> A(p, std::vector<double>(p + 1, 0.0));
    for (size_t r = 0; r < X.size(); ++r)
        for (size_t a = 0; a < p; ++a) {
            for (size_t b = 0; b < p; ++b) A[a][b] += X[r][a] * X[r][b];
            A[a][p] += X[r][a] * y[r];
        }
    for (size_t c = 0; c < p; ++c) {
        size_t piv = c;
        for (size_t r = c + 1; r < p; ++r) if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
        std::swap(A[c], A[piv]);
        for (size_t r = 0; r < p; ++r) if (r != c) {
            double f = A[r][c] / A[c][c];
            for (size_t k = c; k <= p; ++k) A[r][k] -= f * A[c][k];
        }
    }
    std::vector<double> b(p);
    for (size_t i = 0; i < p; ++i) b[i] = A[i][p] / A[i][i];
    return b;
}

Characteristic characteristic(const std::vector<Variant>& v, double amin, double amax) {
    std::vector<std::vector<double>> X; std::vector<double> y;
    for (auto& x : v) if (x.alpha > amin && x.alpha < amax) { X.push_back({1.0, x.alpha}); y.push_back(x.M); }
    if (X.size() < 2) throw std::runtime_error("characteristic: too few variants");
    auto b = lstsq(X, y);
    return {b[0], b[1]};
}

Characteristic characteristic_binned(const std::vector<const Run*>& runs, double width, double lo, double hi,
                                     std::vector<Bin>* bins, size_t min_n) {
    std::map<int, std::vector<double>> g;
    for (auto* r : runs) {
        const int w = (int)std::lround(0.2 / r->dt);
        auto a = moving_average(r->alpha_deg, w), m = moving_average(r->Mg, w);
        for (size_t i = 0; i < a.size(); ++i) if (r->in_soil[i]) g[(int)std::floor(a[i] / width)].push_back(m[i]);
    }
    std::vector<std::vector<double>> X; std::vector<double> y;
    for (auto& kv : g) {
        Bin B;
        B.center = (kv.first + 0.5) * width;
        B.q25 = percentile(kv.second, 25); B.median = percentile(kv.second, 50); B.q75 = percentile(kv.second, 75);
        B.n = kv.second.size();
        if (bins) bins->push_back(B);
        if (kv.first * width >= lo - 1e-9 && (kv.first + 1) * width <= hi + 1e-9 && B.n >= min_n) { X.push_back({1.0, B.center}); y.push_back(B.median); }
    }
    if (X.size() < 2) throw std::runtime_error("characteristic_binned: too few bins");
    auto b = lstsq(X, y);
    return {b[0], b[1]};
}

double LoadModel::torque(const std::string& s, double u, double omega) const {
    auto it = offset.find(s);
    return Ma + (it == offset.end() ? 0.0 : it->second) + K * bite(u, omega);
}

static std::vector<std::string> settings_of(const std::vector<Variant>& v) {
    std::vector<std::string> s;
    for (auto& x : v) if (std::find(s.begin(), s.end(), x.setting) == s.end()) s.push_back(x.setting);
    return s;
}
static std::vector<double> row(const Variant& x, const std::vector<std::string>& S, const LoadModel& L) {
    std::vector<double> r{1.0};
    for (size_t i = 1; i < S.size(); ++i) r.push_back(x.setting == S[i] ? 1.0 : 0.0);
    r.push_back(L.bite(x.u, x.omega));
    return r;
}

LoadModel fit_load_model(const std::vector<Variant>& v, double h, double b, int z) {
    LoadModel L; L.h = h; L.b = b; L.z = z;
    L.settings = settings_of(v);
    std::vector<std::vector<double>> X; std::vector<double> y;
    for (auto& x : v) { X.push_back(row(x, L.settings, L)); y.push_back(x.M); }
    auto c = lstsq(X, y);
    L.Ma = c[0];
    L.offset[L.settings[0]] = 0;
    for (size_t i = 1; i < L.settings.size(); ++i) L.offset[L.settings[i]] = c[i];
    L.K = c.back();                     // N*m per mm of bite length
    L.e = 2 * PI * L.K * 1e3 / (z * h * b);   // J/m^3: work of cutting a unit volume of soil, e = 2*pi*K/(z*h*b)
    return L;
}

LooResult loo_chain(const std::vector<Variant>& v, const std::map<std::string, Characteristic>& ch, double h, double b, int z) {
    LooResult R;
    LoadModel L; L.h = h; L.b = b; L.z = z; L.settings = settings_of(v);
    double se = 0, seM = 0;
    std::vector<double> am;
    for (size_t j = 0; j < v.size(); ++j) {
        std::vector<std::vector<double>> X; std::vector<double> y;
        for (size_t i = 0; i < v.size(); ++i) if (i != j) { X.push_back(row(v[i], L.settings, L)); y.push_back(v[i].M); }
        auto c = lstsq(X, y);
        auto x = row(v[j], L.settings, L);
        double Mp = 0; for (size_t k = 0; k < x.size(); ++k) Mp += x[k] * c[k];
        const auto& C = ch.at(v[j].setting);
        double ap = std::min(std::max((Mp - C.M0) / C.k, 0.0), 90.0);
        R.M_pred.push_back(Mp); R.alpha_pred.push_back(ap);
        se += (ap - v[j].alpha) * (ap - v[j].alpha); seM += (Mp - v[j].M) * (Mp - v[j].M);
        am.push_back(v[j].alpha);
    }
    double m = mean(am), ss = 0; for (double a : am) ss += (a - m) * (a - m);
    R.rmse_alpha = std::sqrt(se / v.size()); R.r2_alpha = 1 - se / ss; R.rmse_M = std::sqrt(seM / v.size());
    return R;
}

// Minimum stiffness: the torque range expected over the operating range (mean load from
// the bite length, widened by z_q standard deviations of the slow fluctuations) must be
// covered by the disc rotation range.
SelectionRule selection_rule(double Ma, double K, double s_min, double s_max, double sigma, double z_q,
                             double amax_deg, double soil_scale) {
    SelectionRule R;
    R.M_lo = soil_scale * (Ma + K * s_min) - z_q * sigma;
    R.M_hi = soil_scale * (Ma + K * s_max) + z_q * sigma;
    R.k_min = (R.M_hi - R.M_lo) / amax_deg * 1000.0;
    return R;
}

// ------------------------------------------------------------------ rotor pair
PairResult pair_analysis(const Run& r, double shift, double phi0_deg, int n, int z, double spacing, double radius,
                         const std::vector<double>* alpha_sim) {
    PairResult o;
    o.name = r.name; o.shift = shift;
    const std::vector<double>& al = alpha_sim ? *alpha_sim : r.alpha_deg;
    if (al.size() != r.alpha_deg.size()) throw std::runtime_error("pair_analysis: alpha size");
    size_t s0 = r.alpha_deg.size(), cnt = 0;
    for (size_t i = 0; i < r.in_soil.size(); ++i) if (r.in_soil[i]) { s0 = std::min(s0, i); ++cnt; }
    if (!cnt) throw std::runtime_error("pair_analysis: no soil phase");
    const size_t i0 = s0 + (size_t)(0.2 * cnt), i1 = s0 + (size_t)(0.8 * cnt);
    const size_t s = (size_t)std::lround(shift / r.dt);
    if (i0 + s + 10 >= i1) throw std::runtime_error("pair_analysis: shift longer than the soil phase");
    const size_t e = i1 - s;  // rotor 1 samples i in [i0, e), rotor 2 sample j = i + s
    const auto& th = r.theta_rad;
    const double d2r = PI / 180.0;
    // correlation of the slow soil loads and relative disc angle
    {
        double m1 = 0, m2 = 0, c = 0, v1 = 0, v2 = 0; const double N = double(e - i0);
        for (size_t i = i0; i < e; ++i) { m1 += r.Mf[i]; m2 += r.Mf[i + s]; }
        m1 /= N; m2 /= N;
        for (size_t i = i0; i < e; ++i) { double a1 = r.Mf[i] - m1, a2 = r.Mf[i + s] - m2; c += a1 * a2; v1 += a1 * a1; v2 += a2 * a2; }
        o.corr = (v1 > 0 && v2 > 0) ? c / std::sqrt(v1 * v2) : 1.0;
        std::vector<double> da;
        for (size_t i = i0; i < e; ++i) da.push_back(al[i] - al[i + s]);
        o.da_mean = mean(da); o.da_sd = stdev(da);
        o.da_lo = percentile(da, 2.5); o.da_hi = percentile(da, 97.5);
        for (double d : da) o.da_max = std::max(o.da_max, std::fabs(d));
    }
    // tine-order torque of the timed pair. The tine-entry torque of each rotor is locked to its tine
    // angle (shaft angle minus disc angle), so with equal pulse amplitudes the order-n component of
    // the pair sum is A |1 + exp(i n (phi + dalpha))|, i.e. relative to the sum of the two rotors
    // rho = |cos(n (phi + dalpha) / 2)|; for rigid tines dalpha = 0. dalpha varies slowly compared
    // with the tine period, so rho is averaged over time.
    std::vector<double> da;
    for (size_t i = i0; i < e; ++i) da.push_back((al[i] - al[i + s]) * d2r);
    auto rho = [&](double phi) {
        double acc = 0;
        for (double d : da) acc += std::fabs(std::cos(0.5 * n * (phi * d2r + d)));
        return acc / double(da.size());
    };
    o.nrev = (int)((th[e - 1] - th[i0]) / (2 * PI));
    o.rho = rho(phi0_deg);
    o.rho_rigid = std::fabs(std::cos(0.5 * n * phi0_deg * d2r));
    for (double p = 0; p <= 180 + 1e-9; p += 1) {
        o.phi_deg.push_back(p);
        o.ratio.push_back(rho(p));
        o.ratio_rigid.push_back(std::fabs(std::cos(0.5 * n * p * d2r)));
    }
    // tip-to-tip distance of the two rotors (rotor 2 mirrored: counter-rotating)
    if (spacing > 0 && radius > 0 && spacing < 2 * radius) {
        auto dist = [&](double psi1, double psi2) {
            double best = 1e9;
            for (int p = 0; p < z; ++p)
                for (int q = 0; q < z; ++q) {
                    const double a1 = psi1 + 2 * PI * p / z, a2 = psi2 + 2 * PI * q / z;
                    const double dx = (spacing - radius * std::cos(a2)) - radius * std::cos(a1);
                    const double dy = radius * std::sin(a2) - radius * std::sin(a1);
                    best = std::min(best, std::hypot(dx, dy));
                }
            return best;
        };
        o.dmin = 1e9; o.dmin_rigid = 1e9;
        double rev_start = th[i0], rev_min = 1e9, rev_t0 = r.t[i0];
        for (size_t i = i0; i < e; ++i) {
            const double t1 = th[i];
            const double d = dist(t1 - al[i] * d2r, t1 + phi0_deg * d2r - al[i + s] * d2r);
            o.dmin = std::min(o.dmin, d);
            rev_min = std::min(rev_min, d);
            if (t1 >= rev_start + 2 * PI) {   // one shaft revolution completed
                o.rev_t.push_back(0.5 * (rev_t0 + r.t[i]) - r.t[i0]); o.rev_dmin.push_back(rev_min);
                rev_start = t1; rev_min = 1e9; rev_t0 = r.t[i];
            }
        }
        for (int k = 0; k < 3600; ++k) {
            const double t1 = 2 * PI * k / 3600.0;
            o.dmin_rigid = std::min(o.dmin_rigid, dist(t1, t1 + phi0_deg * d2r));
        }
    }
    return o;
}
}  // namespace fr
