// Unit and regression tests for FlexRotor (no test framework needed).
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

using namespace fr;
static int failures = 0, checks = 0;
static void check(bool ok, const std::string& what, double got, double ref) {
    ++checks;
    std::printf("%s %-58s got %.6g ref %.6g\n", ok ? "  ok  " : "  FAIL", what.c_str(), got, ref);
    if (!ok) ++failures;
}
static void near(double got, double ref, double tol, const std::string& what) { check(std::fabs(got - ref) <= tol, what, got, ref); }

static const std::string DATA = std::string(FLEXROTOR_SOURCE_DIR) + "/examples/data/";

int main() {
    // ---- signal processing against SciPy 1.x reference values
    {
        std::vector<double> x(2000);
        for (int i = 0; i < 2000; ++i)
            x[i] = std::sin(i * 0.002 * 2 * PI * 0.7) + 0.3 * std::sin(i * 0.002 * 2 * PI * 20) + 0.001 * i;
        auto y = butter_lowpass_filtfilt(x, 1.5, 500, 4);
        near(y[0], 0.18725391083643014, 1e-7, "filtfilt Butterworth, first sample");
        near(y[1000], 1.5863558796959085, 1e-7, "filtfilt Butterworth, middle sample");
        near(y[1999], 1.1975946049264563, 1e-7, "filtfilt Butterworth, last sample");
        auto d = savgol(x, 21, 2, 1, 0.002);
        near(d[0], 23.933404719122837, 1e-6, "Savitzky-Golay derivative, edge");
        near(d[1000], 14.209850780675467, 1e-6, "Savitzky-Golay derivative, middle");
        near(d[1999], 6.433870574428943, 1e-6, "Savitzky-Golay derivative, edge");
    }
    // ---- disc model: static equilibrium and threshold
    {
        DiscParams p; p.friction = false; p.k = 120; p.preload = 90; p.mu = 0.5;
        std::vector<double> M(5000, 15.0);
        auto a = simulate(M, 0.002, p, 10);
        double ref = rad2deg(15.0 / k_si(p.k)) - p.preload;
        near(a.back(), ref, 1e-3, "equilibrium angle (M - M*)/k");
        std::vector<double> lo(3000, 0.9 * p.threshold());
        auto b = simulate(lo, 0.002, p, 10);
        near(b.back(), 0.0, 0.05, "below threshold the disc stays on the initial stop");
        std::vector<double> hi(3000, p.threshold() + k_si(p.k) * deg2rad(120));
        auto c = simulate(hi, 0.002, p, 10);
        near(c.back(), 90.0, 0.1, "above range the disc rests on the end stop");
        DiscParams q = p; q.friction = true; q.Mc = 1.5;
        auto e = simulate_play(M, q);
        check(std::fabs(e.back() - ref) <= rad2deg(q.Mc / k_si(q.k)) + 1e-9, "play model stays within the friction band", e.back(), ref);
    }
    // ---- identification recovers known parameters from synthetic data
    {
        Run r = process(read_csv(DATA + "u2024092707.csv"), [] { ProcessOptions o; o.nominal_rpm = 150; return o; }());
        DiscParams t; t.k = 140; t.preload = 85; t.mu = 0.5; t.Mc = 1.2;
        r.alpha_deg = simulate(r.Mf, r.dt, t, 20);
        for (size_t i = 0; i < r.alpha_deg.size(); ++i) r.in_soil[i] = r.alpha_deg[i] > 1.0;
        FitOptions fo; fo.x0 = {120, 95, 0.4, 1.0};
        auto F = fit({&r}, fo);
        near(F.p.threshold(), t.threshold(), 0.02 * t.threshold(), "synthetic: threshold torque M* recovered");
        near(F.p.Mc, t.Mc, 0.05, "synthetic: dry friction recovered");
        near(F.m.rmse, 0.0, 0.05, "synthetic: residual RMSE");
    }
    // ---- regression: published results (Tucki, Gruz, Mieszkalski; setting N2, 27.09.2024)
    {
        const char* f[9] = {"07", "08", "09", "02", "01", "03", "06", "05", "04"};
        const int rpm[9] = {150, 150, 150, 200, 200, 200, 250, 250, 250};
        std::vector<Run> R;
        for (int i = 0; i < 9; ++i) {
            ProcessOptions o; o.nominal_rpm = rpm[i];
            R.push_back(process(read_csv(DATA + "u20240927" + f[i] + ".csv"), o));
        }
        near(R[0].Mf[5000], 0.6468035388, 1e-6, "soil resistance torque after low-pass filter");
        std::vector<const Run*> P; for (auto& r : R) P.push_back(&r);
        DiscParams p; p.k = 131.7751691258267; p.preload = 89.54239396135536; p.mu = 0.4607769742832293; p.Mc = 1.4225520236415723;
        auto m = metrics(P, p, 20);
        near(m.rmse, 3.9279, 0.002, "published disc model: RMSE (deg)");
        near(m.r2, 0.9571, 0.001, "published disc model: R2");
        near(p.threshold(), 11.8, 0.05, "published threshold torque M* (N*m)");
        auto o2 = order_analysis(R[1], p, 2);
        check(o2.alpha_n_deg > 0.02 && o2.alpha_n_deg < 0.25, "order-2 disc response 0.03-0.23 deg", o2.alpha_n_deg, 0.13);
    }
    // ---- mean-load model and selection rule (setting N1 and N2)
    {
        std::vector<Variant> V;
        const double om[3] = {26.3, 21.3, 15.8}, u[3] = {0.325, 0.276, 0.227};
        for (int s = 0; s < 2; ++s)
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    Variant v; v.setting = s ? "N2" : "N1"; v.omega = om[i]; v.u = u[j];
                    LoadModel L; L.z = 2;
                    v.M = 13.0 - 1.5 * s + 0.19 * L.bite(v.u, v.omega);
                    V.push_back(v);
                }
        auto L = fit_load_model(V, 0.14, 0.25, 2);
        near(L.K, 0.19, 1e-9, "load model: slope recovered");
        near(L.offset["N2"], -1.5, 1e-9, "load model: setting offset recovered");
        near(L.e / 1e3, PI * 0.19 * 1e3 / (0.14 * 0.25) / 1e3, 1e-9, "specific work of cutting");
        auto sr = selection_rule(13.36 - 0.81, 0.1913, 26.5, 62.9, 1.5, 1.645, 90, 1.0);
        near(sr.k_min, (0.1913 * (62.9 - 26.5) + 2 * 1.645 * 1.5) / 90 * 1000, 1e-9, "selection rule: minimum stiffness (N*mm/deg)");
    }
    // ---- rig configuration file: the distributed file equals the built-in defaults
    {
        RigConfig c = read_config(std::string(FLEXROTOR_SOURCE_DIR) + "/examples/rig_sggw.cfg");
        RigConfig d;
        bool same = c.ch.counts_per_rev == d.ch.counts_per_rev && c.ch.torque_scale == d.ch.torque_scale &&
                    c.ch.cart_counts_per_m == d.ch.cart_counts_per_m && c.loss.c0 == d.loss.c0 &&
                    c.loss.c1 == d.loss.c1 && c.loss.c2 == d.loss.c2 && c.disc.J == d.disc.J &&
                    c.disc.amax == d.disc.amax && c.z == d.z && c.depth == d.depth && c.width == d.width &&
                    c.spacing == d.spacing && c.phase_deg == d.phase_deg && c.spring.d == d.spring.d && c.spring.D == d.spring.D &&
                    c.spring.n == d.spring.n && c.spring.E == d.spring.E && c.spring.Rm == d.spring.Rm && c.stress_factor == d.stress_factor && c.fatigue_r6 == d.fatigue_r6 && c.fatigue_r5 == d.fatigue_r5 &&
                    c.limits.Do_max == d.limits.Do_max && c.limits.Di_min == d.limits.Di_min && c.limits.L_max == d.limits.L_max &&
                    c.load_corr == d.load_corr && c.radius == d.radius;
        check(same, "rig_sggw.cfg reproduces the built-in defaults", same, 1);
        bool threw = false;
        try { read_config(std::string(FLEXROTOR_SOURCE_DIR) + "/tests/bad.cfg"); } catch (const std::exception&) { threw = true; }
        check(threw, "unknown configuration key is rejected", threw, 1);
    }

    // ---- rotor pair: identical loads reduce to the rigid timed pair; tip geometry
    {
        Run r = process(read_csv(DATA + "u2024092707.csv"), [] { ProcessOptions o; o.nominal_rpm = 150; return o; }());
        auto q = pair_analysis(r, 0.0, 90.0, 2, 2, 0.20, 0.11);
        near(q.da_max, 0.0, 1e-12, "pair, zero shift: relative disc angle");
        near(q.rho, q.rho_rigid, 1e-12, "pair, zero shift: tine-order ratio equals rigid");
        near(q.rho_rigid, 0.0, 1e-12, "pair, 90 deg offset, two tines: rigid ratio cos(90 deg)");
        auto q0 = pair_analysis(r, 0.0, 0.0, 2, 2, 0.20, 0.11);
        near(q0.dmin_rigid, 0.0, 2e-4, "pair, no offset, overlapping circles: tips meet");
        // 90 deg offset, rigid: closed form min |tip1 - tip2| over a revolution (brute force, fine grid)
        double best = 1e9;
        for (int k = 0; k < 360000; ++k) {
            double a = 2 * PI * k / 360000.0;
            for (int p = 0; p < 2; ++p) for (int m = 0; m < 2; ++m) {
                double a1 = a + PI * p, a2 = a + PI / 2 + PI * m;
                best = std::min(best, std::hypot(0.20 - 0.11 * std::cos(a2) - 0.11 * std::cos(a1), 0.11 * std::sin(a2) - 0.11 * std::sin(a1)));
            }
        }
        near(q.dmin_rigid, best, 1e-4, "pair, 90 deg offset: rigid tip distance");
    }

    // ---- design path: torsion spring after EN 13906-3, prototype springs of Gruz (2026), Table 14
    {
        const double T14[4][4] = {{5.5, 46.2, 3.5, 319.10}, {5.5, 46.2, 4.5, 248.19}, {5.0, 46.5, 4.5, 168.60}, {5.0, 46.5, 5.5, 137.95}};
        double worst = 0;
        for (auto& t : T14) { TorsionSpring s; s.d = t[0]; s.D = t[1]; s.n = t[2]; s.E = 207000; worst = std::max(worst, std::fabs(spring_rate(s) / t[3] - 1)); }
        near(worst, 0.0, 2e-3, "spring rate E d^4/(3667 D n) vs design table (rel.)");
        near(wire_rm_dh(5.0), 1660, 1e-9, "EN 10270-1 SH/DH minimum Rm, d = 5 mm");
        TorsionSpring s; s.d = 5; s.D = 46.5; s.n = 5.5;
        auto c = check_spring(s, 11.8, 17.0, 24.0, 90.0);
        near(c.sigma_end, c.M_end * 1e3 * 32 / (PI * 125.0), 1e-9, "static bending stress 32 M / (pi d^3)");
        near(c.Di_end, 46.5 * 5.5 / (5.5 + (c.beta + 90) / 360) - 5, 1e-9, "inner diameter at the end-stop twist");
        SpringLimits lim; lim.Do_max = 80; lim.L_max = 120;
        auto D = design_springs(134, 15.3, 15.3, 27.4, 90.0, lim);
        bool ok = !D.empty();
        for (auto& x : D) {
            TorsionSpring t = x.first; t.n += 0.25;
            ok = ok && x.second.k >= 134 && x.second.stress_ok && x.second.Di_end >= 37.7 && x.second.Lk <= 120 && x.second.Do <= 80 + 1e-9 && spring_rate(t) < 134;
        }
        check(ok, "designed springs: rate, stress, diameters, length, fewest coils", ok, 1);
        const bool none = design_springs(134, 15.3, 15.3, 27.4, 90.0, SpringLimits()).empty();
        check(none, "no 90 deg spring within the 40 mm space of the prototype", none, 1);
        // Gerber criterion: n sa/Se + (n sm/Sut)^2 = 1; endurance of a repeated stress Sr: amplitude = mean = Sr/2 on the parabola
        const double Sut = 1560, Se = gerber_se(0.5 * Sut, Sut), sa = 120, sm = 1000, nf = gerber_nf(sa, sm, Se, Sut);
        near(nf * sa / Se + (nf * sm / Sut) * (nf * sm / Sut), 1.0, 1e-12, "Gerber safety factor on the failure parabola");
        near((0.25 * Sut) / gerber_se(0.5 * Sut, Sut) + 0.0625, 1.0, 1e-12, "Gerber endurance from the repeated-stress strength");
        TorsionSpring e; e.d = 6.5; e.D = 63.2; e.n = 4.5;
        auto ce = check_spring(e, 15.3, 15.3, 27.4, 36.0, 0.7, 1.645 * 1.5);
        check(ce.nf_fl > 1 && ce.nf_range > 1 && c.nf_range < 1, "fatigue: designed spring passes, prototype spring fails", ce.nf_fl, 1);
    }
    // ---- design path: synthetic load (Ornstein-Uhlenbeck) and its estimation from a record
    {
        Rng g(7);
        auto x = ou_process(1.0, 0.5, 0.002, 1000000, g);
        const double m = mean(x), sd = stdev(x);
        double c0 = 0, c1 = 0; for (size_t i = 0; i < x.size(); ++i) { c0 += (x[i] - m) * (x[i] - m); if (i + 250 < x.size()) c1 += (x[i] - m) * (x[i + 250] - m); }
        near(sd, 1.0, 0.03, "OU process: standard deviation");
        near(c1 / c0, std::exp(-1.0), 0.03, "OU process: autocorrelation at lag tau");
        Run r; r.dt = 0.002; r.name = "syn";
        Rng h(3); r.Mf = ou_process(1.2, 0.4, 0.002, 60000, h);
        for (double& v : r.Mf) v += 20;
        r.in_soil.assign(r.Mf.size(), 1);
        auto f = fluctuation(r);
        near(f.mean, 20.0, 0.15, "fluctuation: mean of a synthetic record");
        near(f.sd, 1.2, 0.12, "fluctuation: SD of a synthetic record");
        near(f.tau, 0.4, 0.08, "fluctuation: correlation time of a synthetic record");
    }
    // ---- design path: rotor row geometry
    {
        RowGeometry G; G.r = 0.11; G.z = 2;
        std::vector<double> z0((size_t)(4.0 / 0.002), 0.0);
        G.a = 0.27; G.phi_deg = 0;   // no overlap: both tips at the inner position at the same time
        auto R = row_analysis(G, 250 * PI / 30, 0.2, z0, z0, 0.002);
        near(R.clear_rigid, 0.27 - 0.22, 1e-6, "row, no offset: clearance a - 2r");
        near(R.uncut_rigid, 0.27 - 0.22, 0.003, "row, a > 2r: untilled strip a - 2r (largest empty circle)");
        G.a = 0.217; G.phi_deg = 90;
        auto Q = row_analysis(G, 200 * PI / 30, 0.28, z0, z0, 0.002);
        near(Q.clear_min, Q.clear_rigid, 1e-4, "row, rigid discs: clearance over time equals the rigid value");
        Run r = process(read_csv(DATA + "u2024092707.csv"), [] { ProcessOptions o; o.nominal_rpm = 150; return o; }());
        auto P = pair_analysis(r, 0.0, 90.0, 2, 2, 0.217, 0.11);
        near(Q.clear_rigid, P.dmin_rigid, 1e-4, "row and pair commands: same rigid tip distance");
    }

    std::printf("\n%d checks, %d failed\n", checks, failures);
    return failures ? 1 : 0;
}
