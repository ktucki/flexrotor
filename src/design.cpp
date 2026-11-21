// Design path: synthetic soil load, torsion spring dimensioning after EN 13906-3 and the tine
// geometry of two adjacent rotors of a row.
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fr {

// ------------------------------------------------------------------ load fluctuation
static void soil_window(const Run& r, size_t& i0, size_t& i1) {
    size_t s0 = r.in_soil.size(), cnt = 0;
    for (size_t i = 0; i < r.in_soil.size(); ++i) if (r.in_soil[i]) { s0 = std::min(s0, i); ++cnt; }
    if (!cnt) throw std::runtime_error("no soil phase in " + r.name);
    i0 = s0 + (size_t)(0.2 * cnt); i1 = s0 + (size_t)(0.8 * cnt);
}

Fluct fluctuation(const Run& r) {
    size_t i0, i1; soil_window(r, i0, i1);
    std::vector<double> x(r.Mf.begin() + i0, r.Mf.begin() + i1);
    Fluct f; f.mean = mean(x);
    for (double& v : x) v -= f.mean;
    double c0 = 0; for (double v : x) c0 += v * v;
    f.sd = std::sqrt(c0 / x.size());
    const size_t n = x.size();
    f.tau = 0.5 * n * r.dt;
    double prev = 1.0;
    for (size_t L = 1; L < n / 2; ++L) {
        double c = 0; for (size_t i = 0; i + L < n; ++i) c += x[i] * x[i + L];
        c /= c0;
        if (c < std::exp(-1.0)) {   // linear interpolation of the 1/e crossing
            const double frac = (prev - std::exp(-1.0)) / (prev - c);
            f.tau = (L - 1 + frac) * r.dt;
            break;
        }
        prev = c;
    }
    return f;
}

// ------------------------------------------------------------------ random numbers
Rng::Rng(unsigned seed) : m_(seed) {}
double Rng::uniform() { return (double(m_()) + 0.5) / 4294967296.0; }
double Rng::normal() {
    if (has_) { has_ = false; return next_; }
    const double u1 = uniform(), u2 = uniform(), rr = std::sqrt(-2.0 * std::log(u1)), t = 2 * PI * u2;
    next_ = rr * std::sin(t); has_ = true;
    return rr * std::cos(t);
}

std::vector<double> ou_process(double sd, double tau, double dt, size_t n, Rng& g) {
    std::vector<double> x(n, 0.0);
    if (n == 0 || sd <= 0) return x;
    if (tau <= 0 || dt <= 0) throw std::runtime_error("ou_process: tau and dt must be positive");
    const double c = std::exp(-dt / tau), s = sd * std::sqrt(1 - c * c);
    x[0] = sd * g.normal();
    for (size_t i = 1; i < n; ++i) x[i] = c * x[i - 1] + s * g.normal();
    return x;
}

DiscStats disc_stats(const std::vector<double>& a, double amax, size_t i0, size_t i1) {
    if (i1 == 0 || i1 > a.size()) i1 = a.size();
    std::vector<double> w(a.begin() + i0, a.begin() + i1);
    if (w.empty()) throw std::runtime_error("disc_stats: empty window");
    DiscStats s;
    s.mean = mean(w); s.sd = stdev(w); s.p05 = percentile(w, 5); s.p95 = percentile(w, 95);
    for (double v : w) { if (v < 0.5) s.low += 1; else if (v > amax - 0.5) s.up += 1; else s.mid += 1; }
    s.low /= w.size(); s.mid /= w.size(); s.up /= w.size();
    return s;
}

// ------------------------------------------------------------------ torsion spring (EN 13906-3)
// Minimum tensile strength of patented cold drawn spring steel wire, grades SH and DH, EN 10270-1.
double wire_rm_dh(double d) {
    static const double T[][2] = {{3.0, 1840}, {3.2, 1820}, {3.4, 1790}, {3.6, 1770}, {3.8, 1750}, {4.0, 1740},
                                  {4.25, 1710}, {4.5, 1690}, {4.75, 1680}, {5.0, 1660}, {5.3, 1640}, {5.6, 1620},
                                  {6.0, 1590}, {6.3, 1570}, {6.5, 1560}, {7.0, 1540}, {7.5, 1510}, {8.0, 1490}};
    const int n = sizeof(T) / sizeof(T[0]);
    if (d < T[0][0] - 1e-9 || d > T[n - 1][0] + 1e-9) throw std::runtime_error("wire_rm_dh: diameter outside 3-8 mm; give Rm");
    for (int i = 0; i + 1 < n; ++i)
        if (d <= T[i + 1][0] + 1e-12) return T[i][1] + (T[i + 1][1] - T[i][1]) * (d - T[i][0]) / (T[i + 1][0] - T[i][0]);
    return T[n - 1][1];
}

static const double KDEG = 64.0 * 180.0 / PI;   // 3667.3: k (N*mm/deg) = E d^4 / (KDEG D n)
double spring_rate(const TorsionSpring& s) { return s.E * std::pow(s.d, 4) / (KDEG * s.D * s.n); }
double coils_for_rate(double k, double d, double D, double E) { return E * std::pow(d, 4) / (KDEG * D * k); }

double gerber_se(double Sr, double Sut) {
    const double h = Sr / 2.0;
    return h / (1.0 - (h / Sut) * (h / Sut));
}

double gerber_nf(double sa, double sm, double Se, double Sut) {
    if (sa <= 0) return 1e9;
    if (sm <= 0) return Se / sa;
    const double r = 2.0 * sm * Se / (Sut * sa);
    return 0.5 * (Sut / sm) * (Sut / sm) * (sa / Se) * (-1.0 + std::sqrt(1.0 + r * r));
}

SpringCheck check_spring(const TorsionSpring& s, double Mstar, double M_lo, double M_hi, double amax_deg, double stress_factor,
                         double M_fl, double r6, double r5) {
    SpringCheck c;
    c.k = spring_rate(s);
    c.beta = Mstar * 1e3 / c.k;
    c.M_end = Mstar + c.k * amax_deg * 1e-3;
    const double W = PI * std::pow(s.d, 3) / 32.0;          // section modulus, mm^3
    const double w = s.D / s.d;
    c.q = (w + 0.07) / (w - 0.75);                           // stress correction factor (fatigue)
    c.sigma_end = c.M_end * 1e3 / W;                         // static: without q
    c.sigma_lo = c.q * M_lo * 1e3 / W;
    c.sigma_hi = c.q * M_hi * 1e3 / W;
    const double Rm = s.Rm > 0 ? s.Rm : wire_rm_dh(s.d);
    c.sigma_allow = stress_factor * Rm;
    c.stress_ok = c.sigma_end <= c.sigma_allow;
    const double th = c.beta + amax_deg;                     // twist at the end stop, deg (winding direction)
    c.Di_end = s.D * s.n / (s.n + th / 360.0) - s.d;
    c.Do = s.D + s.d;
    c.Lk = (s.n + 1.5) * s.d;
    c.mass = 7850e-9 * (PI * s.d * s.d / 4) * (PI * s.D * s.n);
    if (r5 > 0) c.nf_range = gerber_nf((c.sigma_hi - c.sigma_lo) / 2, (c.sigma_hi + c.sigma_lo) / 2, gerber_se(r5 * Rm, Rm), Rm);
    if (M_fl > 0 && r6 > 0) {
        c.sigma_fa = c.q * M_fl * 1e3 / W;
        c.sigma_fm = c.q * (M_hi - M_fl) * 1e3 / W;
        c.nf_fl = gerber_nf(c.sigma_fa, c.sigma_fm, gerber_se(r6 * Rm, Rm), Rm);
    }
    return c;
}

std::vector<std::pair<TorsionSpring, SpringCheck>> design_springs(double k_min, double Mstar, double M_lo, double M_hi,
                                                                  double amax_deg, const SpringLimits& lim, double E,
                                                                  double stress_factor) {
    std::vector<double> ds = lim.d_list;
    if (ds.empty()) ds = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};
    std::vector<std::pair<TorsionSpring, SpringCheck>> out;
    for (double d : ds) {
        bool found = false; std::pair<TorsionSpring, SpringCheck> best;
        for (double D = lim.Di_min + d; D + d <= lim.Do_max + 1e-9; D += lim.D_step) {
            TorsionSpring s; s.d = d; s.D = D; s.E = E;
            s.n = std::floor(coils_for_rate(k_min, d, D, E) * 4.0) / 4.0;   // quarter turns, rate >= k_min
            if (s.n < 2) continue;
            auto c = check_spring(s, Mstar, M_lo, M_hi, amax_deg, stress_factor);
            if (!c.stress_ok || c.Di_end < lim.Di_min || c.Lk > lim.L_max) continue;
            if (!found || c.Do < best.second.Do - 1e-9) { best = {s, c}; found = true; }   // most compact
        }
        if (found) out.push_back(best);
    }
    return out;
}

// ------------------------------------------------------------------ rotor row
double tip_distance(const RowGeometry& g, double psi1, double psi2) {
    double best = 1e9;
    for (int p = 0; p < g.z; ++p)
        for (int q = 0; q < g.z; ++q) {
            const double a1 = psi1 + 2 * PI * p / g.z, a2 = psi2 + 2 * PI * q / g.z;
            const double dy = (g.a - g.r * std::cos(a2)) - g.r * std::cos(a1);
            const double dx = g.r * std::sin(a2) - g.r * std::sin(a1);
            best = std::min(best, std::hypot(dx, dy));
        }
    return best;
}

namespace {
// Exact squared Euclidean distance transform in 1D (Felzenszwalb & Huttenlocher).
void edt1(const double* f, double* d, int n, std::vector<int>& v, std::vector<double>& z) {
    v.assign(n, 0); z.assign(n + 1, 0);
    int k = 0; v[0] = 0; z[0] = -1e300; z[1] = 1e300;
    for (int q = 1; q < n; ++q) {
        double s;
        while (true) {
            s = ((f[q] + double(q) * q) - (f[v[k]] + double(v[k]) * v[k])) / (2.0 * q - 2.0 * v[k]);
            if (s <= z[k] && k > 0) --k; else break;
        }
        if (s <= z[k]) { v[0] = q; z[0] = -1e300; z[1] = 1e300; k = 0; continue; }
        ++k; v[k] = q; z[k] = s; z[k + 1] = 1e300;
    }
    k = 0;
    for (int q = 0; q < n; ++q) {
        while (z[k + 1] < q) ++k;
        d[q] = double(q - v[k]) * (q - v[k]) + f[v[k]];
    }
}

struct Raster {
    int nx = 0, ny = 0; double x0 = 0, res = 0.001, per = 0;
    std::vector<unsigned char> m;
    void mark(double x, double y) {
        int i = (int)std::floor((x - x0) / res);
        double yy = std::fmod(y, per); if (yy < 0) yy += per;
        int j = (int)std::floor(yy / res); if (j >= ny) j = ny - 1;
        if (i >= 0 && i < nx) m[(size_t)i * ny + j] = 1;
    }
    // largest empty circle over x in [xa, xb]; returns diameter and centre
    std::array<double, 3> largest(double xa, double xb) const {
        const double INF = 1e18;
        const int nyt = 3 * ny;                        // lateral period: tile three times
        std::vector<double> g((size_t)nx * nyt), col(std::max(nx, nyt)), out(std::max(nx, nyt));
        std::vector<int> v; std::vector<double> z;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < nyt; ++j) col[j] = m[(size_t)i * ny + (j % ny)] ? 0.0 : INF;
            edt1(col.data(), out.data(), nyt, v, z);
            for (int j = 0; j < nyt; ++j) g[(size_t)i * nyt + j] = out[j];
        }
        double best = -1; int bi = 0, bj = 0;
        const int ia = std::max(0, (int)std::ceil((xa - x0) / res)), ib = std::min(nx - 1, (int)std::floor((xb - x0) / res));
        for (int j = ny; j < 2 * ny; ++j) {
            for (int i = 0; i < nx; ++i) col[i] = g[(size_t)i * nyt + j];
            edt1(col.data(), out.data(), nx, v, z);
            for (int i = ia; i <= ib; ++i) if (out[i] > best) { best = out[i]; bi = i; bj = j - ny; }
        }
        return {x0 + (bi + 0.5) * res, (bj + 0.5) * res, 2.0 * std::sqrt(best) * res};
    }
};
}  // namespace

RowResult row_analysis(const RowGeometry& g, double omega, double u, const std::vector<double>& al1,
                       const std::vector<double>& al2, double dt, double res, double plot_len, double uncut_len) {
    if (al1.size() != al2.size() || al1.size() < 2) throw std::runtime_error("row_analysis: disc angle series");
    if (g.z < 1 || g.r <= 0 || g.a <= 0) throw std::runtime_error("row_analysis: geometry");
    const size_t N = al1.size();
    // the uncut element is evaluated over the first uncut_len metres of travel (+ the tip circle)
    const size_t Nr = std::min(N, (size_t)std::ceil((uncut_len + 2 * g.r + 0.04) / (u * dt)) + 1);
    const double T = (Nr - 1) * dt, d2r = PI / 180.0, phi = g.phi_deg * d2r;
    const double xa = g.r + 0.01, xb = u * T - g.r - 0.01;
    if (xb - xa < 0.05) throw std::runtime_error("row_analysis: record too short for the tip paths");
    RowResult R;
    // clearance
    R.clear_rigid = 1e9;
    for (int k = 0; k < 7200; ++k) { const double t1 = 2 * PI * k / 7200.0; R.clear_rigid = std::min(R.clear_rigid, tip_distance(g, t1, t1 + phi)); }
    R.clear_min = 1e9;
    double rev_min = 1e9, rev_start = 0;
    for (size_t i = 0; i < N; ++i) {
        const double th = omega * i * dt;
        const double d = tip_distance(g, th - al1[i] * d2r, th + phi - al2[i] * d2r);
        R.clear_min = std::min(R.clear_min, d); rev_min = std::min(rev_min, d);
        if (th >= rev_start + 2 * PI) { R.rev_clear.push_back(rev_min); rev_min = 1e9; rev_start = th; }
    }
    // largest uncut element
    const int sub = std::max(1, (int)std::ceil((g.r * omega + u) * dt / (0.4 * res)));
    auto sweep = [&](bool rigid, std::vector<std::array<double, 3>>* path, double p0) {
        Raster ra; ra.res = res; ra.per = 2 * g.a; ra.x0 = -g.r - 0.01;
        ra.nx = (int)std::ceil((u * T + 2 * g.r + 0.02) / res); ra.ny = (int)std::lround(ra.per / res);
        ra.m.assign((size_t)ra.nx * ra.ny, 0);
        const double p1 = p0 + plot_len;
        for (size_t i = 0; i + 1 < Nr; ++i)
            for (int s = 0; s < sub; ++s) {
                const double f = double(s) / sub, t = (i + f) * dt, th = omega * t;
                const double a1 = rigid ? 0 : (al1[i] + f * (al1[i + 1] - al1[i])) * d2r;
                const double a2 = rigid ? 0 : (al2[i] + f * (al2[i + 1] - al2[i])) * d2r;
                for (int p = 0; p < g.z; ++p) {
                    const double s1 = th - a1 + 2 * PI * p / g.z, s2 = th + phi - a2 + 2 * PI * p / g.z;
                    const double x1 = u * t + g.r * std::sin(s1), y1 = g.r * std::cos(s1);
                    const double x2 = u * t + g.r * std::sin(s2), y2 = g.a - g.r * std::cos(s2);
                    ra.mark(x1, y1); ra.mark(x2, y2);
                    if (path) {   // tine tips of rotor 1 and rotor 2 within the plotted window
                        if (x1 >= p0 && x1 <= p1) path->push_back({x1, y1, 1.0 + 0.1 * p});
                        if (x2 >= p0 && x2 <= p1) path->push_back({x2, y2, 2.0 + 0.1 * p});
                    }
                }
            }
        return ra.largest(xa, xb);
    };
    auto c = sweep(false, nullptr, 0);
    R.uncut = c[2]; R.uncut_circle = c;
    R.uncut_rigid = sweep(true, nullptr, 0)[2];
    if (plot_len > 0)   // tip paths over a window centred on the largest uncut element
        sweep(false, &R.path, std::min(std::max(c[0] - 0.5 * plot_len, xa), std::max(xa, xb - plot_len)));
    return R;
}
}  // namespace fr
