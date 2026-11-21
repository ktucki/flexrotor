// FlexRotor command-line interface.
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace fr;

namespace {
struct Entry { std::string name, file, setting, variant; double rpm = 0, u = 0; };

std::vector<Entry> read_table(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open run table " + path);
    auto p = path.find_last_of("/\\");
    std::string dir = p == std::string::npos ? "" : path.substr(0, p + 1);
    std::vector<Entry> out;
    std::string line;
    bool header = true;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (header) { header = false; continue; }
        std::stringstream ss(line);
        std::vector<std::string> c;
        std::string cell;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        if (c.size() < 6) throw std::runtime_error("run table: expected name,file,setting,variant,rpm,u");
        Entry e{c[0], dir + c[1], c[2], c[3], std::stod(c[4]), std::stod(c[5])};
        out.push_back(e);
    }
    return out;
}

struct Args {
    std::string cmd, table, out = ".", setting, params, config, tag;
    std::string phase = "";
    int order = 2, z = -1, reps = 40;
    double duration = 30, soil_scale = 1.0;
    double radius = -1, h = -1, b = -1, sigma = 1.5, zq = 1.645, shift = -1, phi = -1, spacing = -1;
    bool help = false;
};
Args parse(int argc, char** argv) {
    Args a;
    if (argc < 2) { a.help = true; return a; }
    a.cmd = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string k = argv[i];
        auto val = [&]() { if (i + 1 >= argc) throw std::runtime_error("missing value for " + k); return std::string(argv[++i]); };
        if (k == "--out") a.out = val();
        else if (k == "--setting") a.setting = val();
        else if (k == "--params") a.params = val();
        else if (k == "--config") a.config = val();
        else if (k == "--phase") a.phase = val();
        else if (k == "--order") a.order = std::stoi(val());
        else if (k == "--z") a.z = std::stoi(val());
        else if (k == "--radius") a.radius = std::stod(val());
        else if (k == "--depth") a.h = std::stod(val());
        else if (k == "--width") a.b = std::stod(val());
        else if (k == "--sigma") a.sigma = std::stod(val());
        else if (k == "--shift") a.shift = std::stod(val());
        else if (k == "--phi") a.phi = std::stod(val());
        else if (k == "--spacing") a.spacing = std::stod(val());
        else if (k == "--reps") a.reps = std::stoi(val());
        else if (k == "--duration") a.duration = std::stod(val());
        else if (k == "--soil-scale") a.soil_scale = std::stod(val());
        else if (k == "--tag") a.tag = val();
        else if (k == "-h" || k == "--help") a.help = true;
        else if (a.table.empty()) a.table = k;
        else throw std::runtime_error("unknown argument " + k);
    }
    return a;
}

// Two processing conventions: 'angle' (identification of the dynamic model on the nominal
// shaft speed, soil phase = disc lifted from the initial stop) and 'travel' (mean values over
// the middle 60% of the carriage travel, measured shaft speed).
RigConfig rig;  // set from --config in main()

ProcessOptions options(const std::string& phase, double rpm) {
    ProcessOptions o;
    o.ch = rig.ch; o.loss = rig.loss; o.lowpass_hz = rig.lowpass_hz;
    o.phase_angle_deg = rig.phase_angle_deg; o.phase_speed = rig.phase_speed; o.phase_trim = rig.phase_trim;
    if (phase == "travel") {
        o.phase = SoilPhase::CarriageSpeed; o.zero_alpha_before_travel = true; o.use_timestamps = true;
    } else {
        o.phase = SoilPhase::DiscAngle; o.nominal_rpm = rpm;
    }
    return o;
}

std::vector<Run> load(const std::vector<Entry>& T, const std::string& setting, const std::string& phase) {
    std::vector<Run> R;
    for (auto& e : T) {
        if (!setting.empty() && e.setting != setting) continue;
        Run r = process(read_csv(e.file), options(phase, e.rpm));
        r.name = e.name;
        R.push_back(std::move(r));
    }
    if (R.empty()) throw std::runtime_error("no runs selected");
    return R;
}
std::vector<const Run*> ptrs(const std::vector<Run>& R) { std::vector<const Run*> p; for (auto& r : R) p.push_back(&r); return p; }

DiscParams read_params(const std::string& path) {
    DiscParams p = rig.disc;
    if (path.empty()) return p;
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::string line;
    while (std::getline(f, line)) {
        auto q = line.find('=');
        if (q == std::string::npos || line[0] == '#') continue;
        std::string k = line.substr(0, q);
        k.erase(std::remove(k.begin(), k.end(), ' '), k.end());
        double v = std::stod(line.substr(q + 1));
        if (k == "k") p.k = v; else if (k == "preload") p.preload = v; else if (k == "mu") p.mu = v;
        else if (k == "Mc") p.Mc = v; else if (k == "J") p.J = v; else if (k == "amax") p.amax = v;
    }
    return p;
}
void write_params(const std::string& path, const DiscParams& p, const Metrics& m) {
    std::ofstream f(path);
    f << "# FlexRotor " << version() << " identified disc model\n";
    f << "k = " << p.k << "        # N*mm/deg\npreload = " << p.preload << "  # deg\nmu = " << p.mu << "      # N*m*s/rad\nMc = " << p.Mc
      << "      # N*m\nJ = " << p.J << "\namax = " << p.amax << "\n";
    f << "# threshold torque M* = " << p.threshold() << " N*m; RMSE " << m.rmse << " deg; R2 " << m.r2 << "\n";
}

// Calibration of the design path from the records ('travel' convention): coupling characteristics in
// soil, operating variants, mean-load model, its leave-one-out chain and the range of bite lengths.
struct Calib {
    std::map<std::string, Characteristic> C;
    std::vector<Variant> V;
    LoadModel L;
    LooResult Lo;
    double smin = 1e9, smax = 0, offm = 0;
    SelectionRule rule(const Args& a, double scale, double amax = 0) const {
        return selection_rule(L.Ma + offm, L.K, smin, smax, a.sigma, a.zq, amax > 0 ? amax : rig.disc.amax, scale);
    }
    int variant(const std::string& setting, const std::string& name) const {
        for (size_t i = 0; i < V.size(); ++i) if (V[i].setting == setting && V[i].name == name) return (int)i;
        throw std::runtime_error("variant " + setting + "/" + name + " not found");
    }
};
Calib calibrate(const std::vector<Entry>& T, const std::vector<Run>& R, const Args& a, const std::string& o) {
    Calib K;
    std::map<std::string, std::vector<const Run*>> S;
    for (auto& r : R) for (auto& e : T) if (e.name == r.name) S[e.setting].push_back(&r);
    for (auto& s : S) {
        std::vector<Bin> bins;
        K.C[s.first] = characteristic_binned(s.second, 5, 10, 85, &bins);
        if (o.empty()) continue;
        std::ofstream fb(o + "characteristic_" + s.first + ".csv");
        fb << "alpha_bin_centre_deg,n,Mg_q25_Nm,Mg_median_Nm,Mg_q75_Nm,used_in_fit,Mg_fit_Nm\n";
        for (auto& b : bins) {
            bool used = b.center - 2.5 >= 10 - 1e-9 && b.center + 2.5 <= 85 + 1e-9 && b.n >= 500;
            fb << b.center << "," << b.n << "," << b.q25 << "," << b.median << "," << b.q75 << "," << int(used) << "," << K.C[s.first].M0 + K.C[s.first].k * b.center << "\n";
        }
    }
    // variants: mean of all runs with the same setting and group
    std::map<std::pair<std::string, std::string>, std::vector<const Run*>> VG;
    for (auto& r : R) for (auto& e : T) if (e.name == r.name) VG[{e.setting, e.variant}].push_back(&r);
    for (auto& g : VG) {
        Variant v; v.setting = g.first.first; v.name = g.first.second;
        for (auto* r : g.second) { v.M += r->Mg_mean; v.alpha += r->alpha_mean; v.omega += r->omega; v.u += r->v; v.Msd += r->Mg_sd; }
        double n = g.second.size(); v.M /= n; v.alpha /= n; v.omega /= n; v.u /= n; v.Msd /= n;
        K.V.push_back(v);
    }
    K.L = fit_load_model(K.V, a.h, a.b, a.z);
    K.Lo = loo_chain(K.V, K.C, a.h, a.b, a.z);
    for (auto& v : K.V) { double s = K.L.bite(v.u, v.omega); K.smin = std::min(K.smin, s); K.smax = std::max(K.smax, s); }
    for (auto& kv : K.L.offset) K.offm += kv.second;
    K.offm /= K.L.offset.size();
    return K;
}

void usage() {
    std::cout <<
R"(FlexRotor )" << version() << R"( - compliant tine disc coupling of a power harrow rotor

Usage: flexrotor <command> runs.csv [options]

Commands
  summary      soil phase statistics of every run (Mg, alpha, end-stop time, speeds)
  identify     identify k, preload, damping and dry friction of the disc model
  cv           leave-one-run-out and leave-one-speed-group-out validation
  landscape    RMSE of the disc angle over stiffness and preload (identifiability)
  orders       angle-synchronous order analysis of torque and disc angle
  kinematics   bite length factor and coupling dissipation
  map          spring selection map over stiffness and threshold torque
  loadmodel    coupling characteristic in soil, mean-load model, LOO chain
               prediction of the disc angle and the stiffness selection rule
  virtual      virtual soil bin test: disc angle of every run predicted from the
               mean-load model and the load fluctuation of the other runs
  spring       torsion spring after EN 13906-3: stress check of the coupling spring
               and the most compact springs that satisfy the selection rule
  row          two adjacent rotors designed from the geometry: tine clearance and
               largest uncut soil element over the spacing and the phase offset
  pair         prediction for two adjacent counter-rotating rotors on a common gear
               train: relative disc angle, tine-order torque of the pair and tip
               distance (each rotor carries the load measured on one rotor)

Run table (CSV): name,file,setting,variant,rpm,u   (file relative to the table;
  runs with the same setting and variant are averaged by 'loadmodel';
  'cv' groups the runs by nominal rpm)

Options
  --setting S   use only runs of setting S
  --config F    rig and rotor description (key = value; see examples/rig_sggw.cfg);
                without it the defaults of the SGGW soil bin rotor are used
  --params F    disc parameters (written by 'identify')
  --phase P     angle | travel  (soil phase convention; default per command)
  --out DIR     output directory (default .)
  --order n     order for 'orders' (default 2 = tine entries for two tines)
  --z n         tines per disc (2)           --radius r  tine radius, m (0.11)
  --depth h     working depth, m (0.14)      --width b   working width, m (0.25)
                (these override the values in --config)
  --sigma s     SD of slow torque fluctuations for the selection rule (1.5 N*m)
  --shift s     'pair': time shift of the record used for rotor 2, s
                (default: 0.5, 1, 2 and 3 s)
  --phi deg     'pair': tine phase offset of adjacent rotors (default 180/z)
  --reps n      'virtual': realisations of the synthetic load per run (40)
  --duration s  'row': simulated time, s (30)
  --soil-scale f 'virtual': mean soil torque multiplied by f (harder soil; 1)
  --tag name    'virtual': output file virtual_name.csv
  --spacing a   'pair', 'row': axis distance of adjacent rotors, m (0.217; tip distance only
                if given and smaller than twice the tip radius)
)";
}
}  // namespace

int main(int argc, char** argv) {
    try {
        Args a = parse(argc, argv);
        if (a.help || a.table.empty()) { usage(); return a.help ? 0 : 1; }
        if (!a.config.empty()) rig = read_config(a.config);
        if (a.z < 0) a.z = rig.z;
        if (a.radius < 0) a.radius = rig.radius;
        if (a.h < 0) a.h = rig.depth;
        if (a.b < 0) a.b = rig.width;
        if (a.spacing < 0) a.spacing = rig.spacing;
        if (a.phi < 0) a.phi = rig.phase_deg > 0 ? rig.phase_deg : 180.0 / a.z;
        auto T = read_table(a.table);
        const std::string o = a.out + "/";

        if (a.cmd == "summary") {
            auto R = load(T, a.setting, a.phase.empty() ? "travel" : a.phase);
            std::ofstream f(o + "summary.csv");
            f << "run,omega_rad_s,v_m_s,Mg_mean_Nm,Mg_sd_Nm,alpha_mean_deg,alpha_sd_deg,end_stop_pct\n";
            std::printf("%-10s %8s %7s %8s %7s %8s %7s %7s\n", "run", "omega", "v", "Mg", "Mg_sd", "alpha", "a_sd", "stop%");
            for (auto& r : R) {
                f << r.name << "," << r.omega << "," << r.v << "," << r.Mg_mean << "," << r.Mg_sd << "," << r.alpha_mean << "," << r.alpha_sd << "," << r.end_stop_pct << "\n";
                std::printf("%-10s %8.2f %7.3f %8.2f %7.2f %8.1f %7.1f %7.1f\n", r.name.c_str(), r.omega, r.v, r.Mg_mean, r.Mg_sd, r.alpha_mean, r.alpha_sd, r.end_stop_pct);
            }
        } else if (a.cmd == "identify") {
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            auto P = ptrs(R);
            FitOptions fo; fo.base = rig.disc;
            auto F = fit(P, fo);
            std::printf("dry friction model: k %.2f N*mm/deg  preload %.2f deg  mu %.3f N*m*s/rad  Mc %.3f N*m\n", F.p.k, F.p.preload, F.p.mu, F.p.Mc);
            std::printf("threshold torque M* = %.2f N*m   RMSE %.2f deg  MAE %.2f  bias %.2f  R2 %.3f\n", F.p.threshold(), F.m.rmse, F.m.mae, F.m.bias, F.m.r2);
            std::printf("parameter correlation (k, preload, mu, Mc):\n");
            for (auto& r : F.corr) { for (double c : r) std::printf(" %7.3f", c); std::printf("\n"); }
            for (auto& kv : F.m.per) std::printf("  %-8s RMSE %5.2f  bias %6.2f  R2 %.3f\n", kv.first.c_str(), kv.second[0], kv.second[2], kv.second[3]);
            FitOptions fv = fo; fv.friction = false; fv.x0 = {135, 86, 1.8, 0};
            auto V = fit(P, fv);
            std::printf("viscous-only model: k %.2f preload %.2f mu %.3f  RMSE %.2f deg  R2 %.3f\n", V.p.k, V.p.preload, V.p.mu, V.m.rmse, V.m.r2);
            write_params(o + "params.txt", F.p, F.m);
            for (auto* r : P) {
                auto m = simulate(r->Mf, r->dt, F.p, 20);
                PlotSpec s; s.xlabel = "Time (s)"; s.ylabel = "Disc lag angle α (°)";
                double xs = r->t.back() > 20 ? 4 : 2;
                s.xmin = 0; s.xstep = xs; s.xmax = std::ceil(r->t.back() / xs) * xs; s.ymin = -10; s.ymax = 100; s.ystep = 10;
                Series me{r->t, r->alpha_deg, "measured", "#4472C4", "", false, false}; Series mo{r->t, m, "model (dry friction)", "#000000", "3,2"};
                s.series = {me, mo};
                write_svg(o + "fit_" + r->name + ".svg", s);
                auto mv = simulate(r->Mf, r->dt, V.p, 20);
                std::ofstream fc(o + "fit_" + r->name + ".csv");
                fc << "t_s,Mg_filtered_Nm,alpha_measured_deg,alpha_model_deg,alpha_viscous_deg,in_soil\n";
                for (size_t i = 0; i < m.size(); i += 5)
                    fc << r->t[i] << "," << r->Mf[i] << "," << r->alpha_deg[i] << "," << m[i] << "," << mv[i] << "," << int(r->in_soil[i]) << "\n";
            }
        } else if (a.cmd == "cv") {
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            auto P = ptrs(R);
            FitOptions fo; fo.base = rig.disc;
            auto sq = [](const Run* r, const DiscParams& p, double& se, double& n) {
                auto m = simulate(r->Mf, r->dt, p, 20);
                for (size_t i = 0; i < m.size(); ++i) if (r->in_soil[i]) { se += (m[i] - r->alpha_deg[i]) * (m[i] - r->alpha_deg[i]); n += 1; }
            };
            double se1 = 0, n1 = 0;
            std::printf("leave one run out\n");
            for (size_t i = 0; i < P.size(); ++i) {
                std::vector<const Run*> tr; for (size_t j = 0; j < P.size(); ++j) if (j != i) tr.push_back(P[j]);
                auto F = fit(tr, fo); auto m = metrics({P[i]}, F.p);
                sq(P[i], F.p, se1, n1);
                std::printf("  %-8s k %.1f preload %.1f mu %.3f Mc %.3f  test RMSE %.2f deg\n", P[i]->name.c_str(), F.p.k, F.p.preload, F.p.mu, F.p.Mc, m.rmse);
            }
            std::printf("  pooled test RMSE %.2f deg\n", std::sqrt(se1 / n1));
            std::map<std::string, std::vector<const Run*>> G;
            for (size_t i = 0; i < P.size(); ++i) for (auto& e : T) if (e.name == P[i]->name) G[std::to_string((int)e.rpm) + " rpm"].push_back(P[i]);
            double se2 = 0, n2 = 0;
            std::printf("leave one group out\n");
            for (auto& g : G) {
                std::vector<const Run*> tr; for (auto* r : P) if (std::find(g.second.begin(), g.second.end(), r) == g.second.end()) tr.push_back(r);
                auto F = fit(tr, fo); auto m = metrics(g.second, F.p);
                for (auto* r : g.second) sq(r, F.p, se2, n2);
                std::printf("  group %-6s test RMSE %.2f deg\n", g.first.c_str(), m.rmse);
            }
            std::printf("  pooled test RMSE %.2f deg\n", std::sqrt(se2 / n2));
        } else if (a.cmd == "landscape") {
            // RMSE of the disc angle over a grid of stiffness and preload (damping and friction
            // kept at the identified values): shows which combination the records determine.
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            auto P = ptrs(R);
            auto p = read_params(a.params);
            std::ofstream f(o + "landscape.csv");
            f << "k_Nmm_deg,preload_deg,Mstar_Nm,rmse_deg\n";
            for (double k = 90; k <= 190 + 1e-9; k += 5)
                for (double b = 60; b <= 130 + 1e-9; b += 3.5) {
                    DiscParams q = p; q.k = k; q.preload = b;
                    double se = 0; size_t n = 0;
                    for (auto* r : P) {
                        auto m = simulate(r->Mf, r->dt, q, 5);
                        for (size_t i = 0; i < m.size(); i += 5) if (r->in_soil[i]) { se += (m[i] - r->alpha_deg[i]) * (m[i] - r->alpha_deg[i]); ++n; }
                    }
                    f << k << "," << b << "," << q.threshold() << "," << std::sqrt(se / n) << "\n";
                }
            std::printf("landscape written\n");
        } else if (a.cmd == "orders") {
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            auto p = read_params(a.params);
            std::printf("%-8s %6s %5s %8s %8s %9s %9s %6s %6s\n", "run", "rpm", "nrev", "|Mt_n|", "|Mg_n|", "|a_n|", "a_n lin", "T_exp", "T_lin");
            for (auto& r : R) {
                auto x = order_analysis(r, p, a.order);
                std::printf("%-8s %6.0f %5d %8.3f %8.3f %9.4f %9.4f %6.3f %6.3f\n", x.name.c_str(), x.rpm, x.nrev, x.Mt_n, x.Mg_n, x.alpha_n_deg, x.alpha_n_lin_deg, x.T_exp, x.T_lin);
            }
        } else if (a.cmd == "kinematics") {
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            auto p = read_params(a.params);
            for (auto& r : R) {
                double u = 0; for (auto& e : T) if (e.name == r.name) u = e.u;
                auto k = kinematics(r, p, a.z, a.radius, u);
                std::printf("%-8s lambda %.1f  dalpha max %.1f p95 %.1f deg  f %.2f-%.2f (2.5-97.5%%: %.2f-%.2f)  Pd/Pg %.3f%%\n", k.name.c_str(), k.lambda,
                            k.dalpha_max, k.dalpha_p95, k.f_min, k.f_max, k.f_p025, k.f_p975, 100 * k.Pd / k.Pg);
            }
        } else if (a.cmd == "map") {
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            auto p = read_params(a.params);
            std::vector<double> ks, Ms;
            for (double k = 50; k <= 300 + 1e-9; k += 10) ks.push_back(k);
            for (double M = 4; M <= 24 + 1e-9; M += 1) Ms.push_back(M);
            auto m = selection_map(ptrs(R), p, ks, Ms);
            std::ofstream f(o + "selection_map.csv");
            f << "k_Nmm_deg,Mstar_Nm,initial_stop,compliant,end_stop\n";
            for (size_t i = 0; i < Ms.size(); ++i) for (size_t j = 0; j < ks.size(); ++j)
                f << ks[j] << "," << Ms[i] << "," << m.low[i][j] << "," << m.mid[i][j] << "," << m.up[i][j] << "\n";
            write_map_svg(o + "selection_map.svg", m, {{p.k, p.threshold()}});
            std::printf("selection map written (%zu x %zu)\n", Ms.size(), ks.size());
        } else if (a.cmd == "pair") {
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            std::vector<double> shifts = a.shift >= 0 ? std::vector<double>{a.shift} : std::vector<double>{0.5, 1, 2, 3};
            const bool sim = !a.params.empty();
            const DiscParams p = read_params(a.params);
            std::printf("rotor pair (PREDICTION): two counter-rotating rotors on a common gear train, tine phase offset %.1f deg,\n"
                        "order %d; rotor 2 = same record shifted in time; each rotor carries the load measured on one rotor\n", a.phi, a.order);
            std::printf("disc angle: %s\n", sim ? "simulated with --params" : "measured");
            std::printf("rho: order-%d torque of the pair relative to the sum of the two rotors (rigid tines: %.3f)\n", a.order,
                        std::fabs(std::cos(0.5 * a.order * a.phi * PI / 180)));
            std::printf("run      shift_s  corr   dalpha_sd  dalpha 2.5..97.5 %%  |dalpha|max  |Mt_n|   rho    rho*2|Mt_n|");
            if (a.spacing > 0) std::printf("  dmin_mm  dmin_rigid_mm");
            std::printf("\n");
            std::ofstream f(o + "pair.csv"), g(o + "pair_phase.csv");
            f << "run,shift_s,corr,dalpha_mean,dalpha_sd,dalpha_p2_5,dalpha_p97_5,dalpha_absmax,Mt_n,rho,rho_rigid,dmin_m,dmin_rigid_m\n";
            g << "run,shift_s,phi_deg,rho_compliant,rho_rigid\n";
            std::ofstream dd(o + "pair_distance.csv");
            dd << "run,shift_s,spacing_m,t_s,dmin_rev_m,dmin_rigid_m\n";
            for (auto& r : R) {
                std::vector<double> asim;
                if (sim) asim = simulate(r.Mf, r.dt, p, 20);
                const double Mt = order_analysis(r, p, a.order).Mt_n;
                for (double sh : shifts) {
                    auto q = pair_analysis(r, sh, a.phi, a.order, a.z, a.spacing, a.radius, sim ? &asim : nullptr);
                    std::printf("%-8s %5.1f   %5.2f  %7.2f    %7.2f..%-7.2f   %7.2f     %6.3f  %6.3f  %6.3f",
                                r.name.c_str(), sh, q.corr, q.da_sd, q.da_lo, q.da_hi, q.da_max, Mt, q.rho, q.rho * 2 * Mt);
                    if (q.dmin >= 0) std::printf("     %6.1f    %6.1f", q.dmin * 1e3, q.dmin_rigid * 1e3);
                    std::printf("\n");
                    f << r.name << "," << sh << "," << q.corr << "," << q.da_mean << "," << q.da_sd << "," << q.da_lo << "," << q.da_hi << ","
                      << q.da_max << "," << Mt << "," << q.rho << "," << q.rho_rigid << "," << q.dmin << "," << q.dmin_rigid << "\n";
                    for (size_t k = 0; k < q.phi_deg.size(); ++k)
                        g << r.name << "," << sh << "," << q.phi_deg[k] << "," << q.ratio[k] << "," << q.ratio_rigid[k] << "\n";
                    for (size_t k = 0; k < q.rev_t.size(); ++k)
                        dd << r.name << "," << sh << "," << a.spacing << "," << q.rev_t[k] << "," << q.rev_dmin[k] << "," << q.dmin_rigid << "\n";
                }
            }
            if (a.spacing <= 0) std::printf("tip distance not computed: rotor spacing not given (--spacing or 'spacing' in --config)\n");
        } else if (a.cmd == "loadmodel") {
            auto R = load(T, "", a.phase.empty() ? "travel" : a.phase);
            Calib K = calibrate(T, R, a, o);
            auto& C = K.C; auto& V = K.V; auto& L = K.L; auto& Lo = K.Lo;
            for (auto& c : C) std::printf("setting %-4s characteristic in soil: M0 %.2f N*m  k %.1f N*mm/deg  preload %.1f deg\n", c.first.c_str(), c.second.M0,
                                          c.second.k * 1e3, c.second.M0 / c.second.k);
            std::printf("mean load: M = %.2f + d_setting + %.4f s_k   (s_k in mm); specific work of cutting e = %.1f kJ/m3\n", L.Ma, L.K, L.e / 1e3);
            for (auto& kv : L.offset) std::printf("  offset %-4s %+.2f N*m\n", kv.first.c_str(), kv.second);
            std::printf("LOO chain prediction of the disc angle: RMSE %.2f deg  R2 %.3f  (torque RMSE %.2f N*m)\n", Lo.rmse_alpha, Lo.r2_alpha, Lo.rmse_M);
            for (double sc : {1.0, 1.5}) {
                auto sr = K.rule(a, sc);
                std::printf("selection rule (soil scale %.1f): torque %.1f-%.1f N*m  ->  k >= %.0f N*mm/deg, threshold torque ~ %.1f N*m\n", sc, sr.M_lo, sr.M_hi, sr.k_min, sr.M_lo);
            }
            std::ofstream f(o + "variants.csv");
            f << "setting,variant,omega_rad_s,u_m_s,bite_mm,M_Nm,alpha_deg,M_pred_loo,alpha_pred_loo\n";
            PlotSpec s; s.xlabel = "Measured disc lag angle α (°)"; s.ylabel = "Predicted disc lag angle α (°)";
            s.xmin = 30; s.xmax = 90; s.xstep = 10; s.ymin = 30; s.ymax = 90; s.ystep = 10;
            s.series.push_back({{30, 90}, {30, 90}, "1:1 line", "#000000", "", false, false});
            std::map<std::string, Series> pts;
            const char* cols[2] = {"#4472C4", "#ED7D31"};
            int ci = 0;
            for (size_t i = 0; i < V.size(); ++i) {
                f << V[i].setting << "," << V[i].name << "," << V[i].omega << "," << V[i].u << "," << L.bite(V[i].u, V[i].omega) << "," << V[i].M << ","
                  << V[i].alpha << "," << Lo.M_pred[i] << "," << Lo.alpha_pred[i] << "\n";
                if (!pts.count(V[i].setting)) { Series se; se.label = "setting " + V[i].setting; se.color = cols[ci++ % 2]; se.markers = true; se.open = ci == 1; pts[V[i].setting] = se; }
                pts[V[i].setting].x.push_back(V[i].alpha); pts[V[i].setting].y.push_back(Lo.alpha_pred[i]);
            }
            for (auto& kv : pts) s.series.push_back(kv.second);
            write_svg(o + "loo_disc_angle.svg", s);
        } else if (a.cmd == "virtual") {
            // Virtual soil bin test: every run of the setting is predicted without its own record. The mean
            // soil torque comes from the leave-one-out mean-load model, the slow fluctuation (Ornstein-Uhlenbeck
            // process) from the other runs; the disc model with the identified parameters gives the disc angle.
            auto Rt = load(T, "", "travel");
            Calib K = calibrate(T, Rt, a, "");
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            const DiscParams p = read_params(a.params);
            std::vector<Fluct> F; for (auto& r : R) F.push_back(fluctuation(r));
            std::ofstream f(o + "virtual" + (a.tag.empty() ? "" : "_" + a.tag) + ".csv");
            f << "run,M_meas_Nm,M_pred_Nm,sd_used_Nm,tau_used_s,a_mean_meas,a_mean_pred,a_sd_meas,a_sd_pred,a_p05_meas,a_p05_pred,a_p95_meas,a_p95_pred,end_meas,end_pred,init_meas,init_pred\n";
            std::printf("virtual soil bin test (%d realisations per run; each run predicted without its own record)\n", a.reps);
            std::printf("%-8s %7s %7s | %15s | %13s | %13s\n", "run", "M meas", "M pred", "alpha mean m/p", "alpha SD m/p", "end stop m/p");
            double se = 0, ss = 0, mm = 0; std::vector<double> am;
            for (size_t j = 0; j < R.size(); ++j) {
                const Run& r = R[j];
                const Entry* e = nullptr; for (auto& x : T) if (x.name == r.name) e = &x;
                const double Mp = a.soil_scale * K.Lo.M_pred[K.variant(e->setting, e->variant)];
                double s2 = 0; std::vector<double> taus;
                for (size_t i = 0; i < R.size(); ++i) if (i != j) { s2 += F[i].sd * F[i].sd; taus.push_back(F[i].tau); }
                Fluct fj; fj.sd = std::sqrt(s2 / taus.size()); fj.tau = percentile(taus, 50);
                size_t s0 = r.in_soil.size(), cnt = 0;
                for (size_t i = 0; i < r.in_soil.size(); ++i) if (r.in_soil[i]) { s0 = std::min(s0, i); ++cnt; }
                const size_t i0 = s0 + (size_t)(0.2 * cnt), i1 = s0 + (size_t)(0.8 * cnt), n = i1 - i0, warm = (size_t)(3.0 / r.dt);
                auto m = disc_stats(r.alpha_deg, p.amax, i0, i1);
                DiscStats q; Rng g(1000 + (unsigned)j);
                for (int k = 0; k < a.reps; ++k) {
                    auto x = ou_process(fj.sd, fj.tau, r.dt, warm + n, g);
                    for (double& v : x) v += Mp;
                    auto st = disc_stats(simulate(x, r.dt, p, 10), p.amax, warm, warm + n);
                    q.mean += st.mean; q.sd += st.sd; q.p05 += st.p05; q.p95 += st.p95; q.low += st.low; q.up += st.up;
                }
                q.mean /= a.reps; q.sd /= a.reps; q.p05 /= a.reps; q.p95 /= a.reps; q.low /= a.reps; q.up /= a.reps;
                std::printf("%-8s %7.2f %7.2f | %6.1f  %6.1f | %5.1f  %5.1f | %5.1f%% %5.1f%%\n", r.name.c_str(), F[j].mean, Mp, m.mean, q.mean, m.sd, q.sd, 100 * m.up, 100 * q.up);
                f << r.name << "," << F[j].mean << "," << Mp << "," << fj.sd << "," << fj.tau << "," << m.mean << "," << q.mean << "," << m.sd << "," << q.sd << ","
                  << m.p05 << "," << q.p05 << "," << m.p95 << "," << q.p95 << "," << m.up << "," << q.up << "," << m.low << "," << q.low << "\n";
                se += (q.mean - m.mean) * (q.mean - m.mean); am.push_back(m.mean); mm += m.mean;
            }
            mm /= am.size(); for (double v : am) ss += (v - mm) * (v - mm);
            std::printf("mean disc angle against the records: RMSE %.2f deg  R2 %.3f\n", std::sqrt(se / am.size()), 1 - se / ss);
        } else if (a.cmd == "spring") {
            // Spring dimensioning (EN 13906-3): check of the prototype spring and the most compact springs
            // of standard wire diameters that satisfy the selection rule in the tested and in a harder soil.
            auto Rt = load(T, "", "travel");
            Calib K = calibrate(T, Rt, a, "");
            const DiscParams p = read_params(a.params);
            double Mlo = 1e9, Mhi = 0;
            for (auto& v : K.V) if (a.setting.empty() || v.setting == a.setting) { Mlo = std::min(Mlo, v.M); Mhi = std::max(Mhi, v.M); }
            const TorsionSpring& sp = rig.spring;
            const double Mfl = a.zq * a.sigma;   // amplitude of the slow fluctuation (z_q sigma of the selection rule)
            auto c = check_spring(sp, p.threshold(), Mlo - Mfl, Mhi + Mfl, rig.disc.amax, rig.stress_factor, Mfl, rig.fatigue_r6, rig.fatigue_r5);
            const double Rm = sp.Rm > 0 ? sp.Rm : wire_rm_dh(sp.d);
            std::printf("prototype spring d %.2f mm  D %.1f mm  n %.2f: rate %.1f N*mm/deg (identified in soil %.1f)\n", sp.d, sp.D, sp.n, c.k, p.k);
            std::printf("  measured mean soil torque %.1f-%.1f N*m: static bending stress %.0f-%.0f N/mm2\n", Mlo, Mhi,
                        Mlo * 1e3 / (PI * std::pow(sp.d, 3) / 32), Mhi * 1e3 / (PI * std::pow(sp.d, 3) / 32));
            std::printf("  end stop (M* %.1f N*m): torque %.1f N*m, stress %.0f N/mm2; permissible %.2f Rm = %.0f N/mm2 (Rm %.0f)  -> %s\n",
                        p.threshold(), c.M_end, c.sigma_end, rig.stress_factor, c.sigma_allow, Rm, c.stress_ok ? "OK" : "EXCEEDED");
            std::printf("  inner diameter at the end stop %.1f mm, outer diameter %.1f mm, body length %.1f mm\n", c.Di_end, c.Do, c.Lk);
            std::printf("  fatigue (Gerber): fluctuation %.0f +- %.0f N/mm2, safety factor %.2f (1e6 cycles); operating range %.0f-%.0f N/mm2, %.2f (1e5 cycles)\n",
                        c.sigma_fm, c.sigma_fa, c.nf_fl, c.sigma_lo, c.sigma_hi, c.nf_range);
            std::ofstream f(o + "spring_design.csv"), g(o + "spring_chart.csv");
            f << "soil_scale,k_min,M_lo,M_hi,d_mm,D_mm,n,k,beta_deg,M_end,sigma_end,sigma_allow,sigma_q_lo,sigma_q_hi,Di_end,Do,Lk,mass_kg,amax_deg,sigma_fa,sigma_fm,nf_fluct_1e6,nf_range_1e5\n";
            g << "soil_scale,D_mm,k_min,M_end,d_mm,n_for_kmin,sigma_end,sigma_allow\n";
            f << "prototype,," << Mlo << "," << Mhi << "," << sp.d << "," << sp.D << "," << sp.n << "," << c.k << "," << c.beta << "," << c.M_end << "," << c.sigma_end << ","
              << c.sigma_allow << "," << c.sigma_lo << "," << c.sigma_hi << "," << c.Di_end << "," << c.Do << "," << c.Lk << "," << c.mass << "," << rig.disc.amax << "," << c.sigma_fa << "," << c.sigma_fm << "," << c.nf_fl << "," << c.nf_range << "\n";
            // shortest spring that meets the rule: largest mean coil diameter within Do_max for every standard wire
            auto shortest = [&](const SelectionRule& sr, double am, std::pair<TorsionSpring, SpringCheck>& best) {
                bool found = false;
                for (double d : {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0}) {
                    TorsionSpring t; t.d = d; t.D = std::floor((rig.limits.Do_max - d) / rig.limits.D_step) * rig.limits.D_step; t.E = sp.E; t.Rm = sp.Rm;
                    t.n = std::floor(coils_for_rate(sr.k_min, d, t.D, sp.E) * 4.0) / 4.0;
                    if (t.n < 2) continue;
                    auto q = check_spring(t, sr.M_lo, sr.M_lo, sr.M_hi, am, rig.stress_factor);
                    if (!q.stress_ok || q.Di_end < rig.limits.Di_min) continue;
                    if (!found || q.Lk < best.second.Lk) { best = {t, q}; found = true; }
                }
                return found;
            };
            // rotation range of the disc against the length of the shortest spring
            std::ofstream fr(o + "spring_range.csv");
            fr << "soil_scale,amax_deg,k_min,M_lo,M_hi,d_mm,D_mm,n,k,body_length_mm,sigma_end,sigma_allow\n";
            for (double sc : {1.0, 1.5})
                for (double am = 20; am <= rig.disc.amax + 1e-9; am += 5) {
                    auto sr = K.rule(a, sc, am);
                    std::pair<TorsionSpring, SpringCheck> best;
                    if (!shortest(sr, am, best)) continue;
                    fr << sc << "," << am << "," << sr.k_min << "," << sr.M_lo << "," << sr.M_hi << "," << best.first.d << "," << best.first.D << "," << best.first.n
                       << "," << best.second.k << "," << best.second.Lk << "," << best.second.sigma_end << "," << best.second.sigma_allow << "\n";
                }
            for (double sc : {1.0, 1.5}) {   // stress at the end stop against the wire diameter (full rotation range)
                auto sr = K.rule(a, sc);
                const double Mend = sr.M_lo + sr.k_min * rig.disc.amax * 1e-3;
                for (double d = 3.0; d <= 8.0 + 1e-9; d += 0.05)
                    g << sc << "," << sp.D << "," << sr.k_min << "," << Mend << "," << d << "," << coils_for_rate(sr.k_min, d, sp.D, sp.E) << ","
                      << Mend * 1e3 / (PI * d * d * d / 32) << "," << rig.stress_factor * (sp.Rm > 0 ? sp.Rm : wire_rm_dh(d)) << "\n";
            }
            auto put = [&](const std::string& tag, double sc, double am, const SelectionRule& sr, std::pair<TorsionSpring, SpringCheck> x) {
                x.second = check_spring(x.first, sr.M_lo, sr.M_lo, sr.M_hi, am, rig.stress_factor, a.zq * a.sigma, rig.fatigue_r6, rig.fatigue_r5);
                std::printf("  %-9s d %.1f  D %.1f  n %.2f  k %.1f  beta %.0f deg  range %.0f deg  M_end %.1f  stress %.0f/%.0f  Di(end) %.1f  Do %.1f  L %.0f mm  %.2f kg\n",
                            tag.c_str(), x.first.d, x.first.D, x.first.n, x.second.k, x.second.beta, am, x.second.M_end, x.second.sigma_end, x.second.sigma_allow,
                            x.second.Di_end, x.second.Do, x.second.Lk, x.second.mass);
                std::printf("            fatigue (Gerber): fluctuation %.0f +- %.0f N/mm2, safety factor %.2f (1e6); operating range %.0f-%.0f N/mm2, %.2f (1e5)\n",
                            x.second.sigma_fm, x.second.sigma_fa, x.second.nf_fl, x.second.sigma_lo, x.second.sigma_hi, x.second.nf_range);
                f << sc << "_" << tag << "," << sr.k_min << "," << sr.M_lo << "," << sr.M_hi << "," << x.first.d << "," << x.first.D << "," << x.first.n << "," << x.second.k << ","
                  << x.second.beta << "," << x.second.M_end << "," << x.second.sigma_end << "," << x.second.sigma_allow << "," << x.second.sigma_lo << ","
                  << x.second.sigma_hi << "," << x.second.Di_end << "," << x.second.Do << "," << x.second.Lk << "," << x.second.mass << "," << am << "," << x.second.sigma_fa << "," << x.second.sigma_fm << "," << x.second.nf_fl << "," << x.second.nf_range << "\n";
            };
            for (double sc : {1.0, 1.5}) {
                // (1) full rotation range, diameters of the envelope, any length; (2) largest range within the whole envelope
                auto sr = K.rule(a, sc);
                std::printf("soil scale %.1f: rotation range %.0f deg -> k >= %.0f N*mm/deg, M* = %.1f N*m, torque up to %.1f N*m\n", sc, rig.disc.amax, sr.k_min, sr.M_lo, sr.M_hi);
                std::pair<TorsionSpring, SpringCheck> full;
                if (shortest(sr, rig.disc.amax, full)) put("full", sc, rig.disc.amax, sr, full);   // full range, any length
                double am = rig.disc.amax; std::vector<std::pair<TorsionSpring, SpringCheck>> D;
                for (; am >= 10 - 1e-9; am -= 1) {
                    sr = K.rule(a, sc, am);
                    D = design_springs(sr.k_min, sr.M_lo, sr.M_lo, sr.M_hi, am, rig.limits, sp.E, rig.stress_factor);
                    if (!D.empty()) break;
                }
                std::printf("  envelope Do <= %.0f, Di(end) >= %.1f, L <= %.0f mm: ", rig.limits.Do_max, rig.limits.Di_min, rig.limits.L_max);
                if (D.empty()) { std::printf("no spring for any rotation range >= 10 deg\n"); continue; }
                std::printf("largest rotation range %.0f deg (k >= %.0f N*mm/deg)\n", am, sr.k_min);
                auto best = *std::min_element(D.begin(), D.end(), [](const auto& x, const auto& y) { return x.second.Do < y.second.Do; });
                put("envelope", sc, am, sr, best);
                DiscParams q = p; q.k = best.second.k; q.preload = best.second.beta; q.amax = am;
                char nm[64]; std::snprintf(nm, sizeof nm, "spring_params_%.1f.txt", sc);
                std::ofstream fq(o + nm);
                fq << "# FlexRotor " << version() << " disc model with the designed spring d " << best.first.d << " mm, D " << best.first.D
                   << " mm, n " << best.first.n << ", rotation range " << am << " deg (soil scale " << sc << ")\n";
                fq << "k = " << q.k << "\npreload = " << q.preload << "\nmu = " << q.mu << "\nMc = " << q.Mc << "\nJ = " << q.J << "\namax = " << q.amax << "\n";
            }
        } else if (a.cmd == "row") {
            // Two adjacent rotors of a row designed from the geometry: synthetic soil loads (mean-load model,
            // fluctuation of the records, correlation load_corr), disc model, tine clearance and uncut element.
            auto Rt = load(T, "", "travel");
            Calib K = calibrate(T, Rt, a, "");
            auto R = load(T, a.setting, a.phase.empty() ? "angle" : a.phase);
            const DiscParams p = read_params(a.params);
            double s2 = 0; std::vector<double> taus;
            for (auto& r : R) { auto x = fluctuation(r); s2 += x.sd * x.sd; taus.push_back(x.tau); }
            Fluct fl; fl.sd = std::sqrt(s2 / R.size()); fl.tau = percentile(taus, 50);
            const double rho = rig.load_corr, dt = R[0].dt;
            RowGeometry G; G.r = a.radius; G.z = a.z; G.phi_deg = a.phi; G.a = a.spacing > 0 ? a.spacing : 2 * a.radius;
            std::printf("row of rotors: r %.3f m, z %d, spacing %.3f m (tine path overlap %.1f mm), phase offset %.1f deg,\n"
                        "load fluctuation SD %.2f N*m, tau %.2f s, correlation %.2f\n", G.r, G.z, G.a, (2 * G.r - G.a) * 1e3, G.phi_deg, fl.sd, fl.tau, rho);
            if (G.a > 2 * G.r) std::printf("WARNING: spacing larger than the tip diameter: a strip of %.1f mm between the rotors is not tilled\n", (G.a - 2 * G.r) * 1e3);
            std::ofstream fs(o + "row_spacing.csv"), fp(o + "row_phase.csv"), fr_(o + "row_runs.csv");
            fs << "run,spacing_m,clear_rigid_m,clear_min_m,clear_p01_m,uncut_m,uncut_rigid_m\n";
            fp << "run,phi_deg,clear_rigid_m,clear_min_m\n";
            fr_ << "run,omega,u,bite_mm,lambda,M_mean,dalpha_sd,dalpha_p2_5,dalpha_p97_5,alpha_mean,end_stop\n";
            std::ofstream fk(o + "row_stiffness.csv");
            fk << "run,k_Nmm_deg,preload_deg,clear_min_m,end_stop\n";
            std::vector<std::string> seen;
            for (size_t j = 0; j < R.size(); ++j) {
                const Entry* e = nullptr; for (auto& x : T) if (x.name == R[j].name) e = &x;
                if (std::find(seen.begin(), seen.end(), e->variant) != seen.end()) continue;
                seen.push_back(e->variant);
                const Variant& v = K.V[K.variant(e->setting, e->variant)];
                const double Mm = K.L.torque(e->setting, v.u, v.omega);
                const size_t warm = (size_t)(3.0 / dt), n = (size_t)(a.duration / dt);
                Rng g(2000 + (unsigned)j);
                auto x1 = ou_process(fl.sd, fl.tau, dt, warm + n, g), x2 = ou_process(fl.sd, fl.tau, dt, warm + n, g);
                for (size_t i = 0; i < x1.size(); ++i) { x2[i] = rho * x1[i] + std::sqrt(1 - rho * rho) * x2[i]; x1[i] += Mm; x2[i] += Mm; }
                auto a1 = simulate(x1, dt, p, 10), a2 = simulate(x2, dt, p, 10);
                a1.erase(a1.begin(), a1.begin() + warm); a2.erase(a2.begin(), a2.begin() + warm);
                std::vector<double> da(a1.size()); for (size_t i = 0; i < da.size(); ++i) da[i] = a1[i] - a2[i];
                auto st = disc_stats(a1, p.amax);
                fr_ << R[j].name << "," << v.omega << "," << v.u << "," << K.L.bite(v.u, v.omega) << "," << G.r * v.omega / v.u << "," << Mm << "," << stdev(da) << ","
                    << percentile(da, 2.5) << "," << percentile(da, 97.5) << "," << st.mean << "," << st.up << "\n";
                for (double sa = 0.180; sa <= 0.240 + 1e-9; sa += 0.005) {
                    RowGeometry Gs = G; Gs.a = sa;
                    auto rr = row_analysis(Gs, v.omega, v.u, a1, a2, dt, 0.001, 0, 2.0);
                    fs << R[j].name << "," << sa << "," << rr.clear_rigid << "," << rr.clear_min << "," << percentile(rr.rev_clear, 1) << "," << rr.uncut << "," << rr.uncut_rigid << "\n";
                }
                {
                    auto rr = row_analysis(G, v.omega, v.u, a1, a2, dt, 0.001, 0.30, 2.0);
                    std::printf("%-8s %5.0f rpm %.3f m/s  bite %4.1f mm  dalpha SD %4.1f deg (%.1f..%.1f)  clearance %5.1f mm (rigid %5.1f)  uncut %4.1f mm (rigid %4.1f)\n",
                                R[j].name.c_str(), v.omega * 30 / PI, v.u, K.L.bite(v.u, v.omega), stdev(da), percentile(da, 2.5), percentile(da, 97.5),
                                rr.clear_min * 1e3, rr.clear_rigid * 1e3, rr.uncut * 1e3, rr.uncut_rigid * 1e3);
                    std::ofstream fpth(o + "row_paths_" + R[j].name + ".csv");
                    fpth << "x_m,y_m,rotor_tine\n";
                    for (auto& q : rr.path) fpth << q[0] << "," << q[1] << "," << q[2] << "\n";
                    fpth << "# largest uncut element: x " << rr.uncut_circle[0] << " y " << rr.uncut_circle[1] << " d " << rr.uncut_circle[2] << "\n";
                }
                // clearance for softer and stiffer springs of the same threshold torque
                for (double kk : {40.0, 60.0, 80.0, 100.0, 134.0, 200.0, 300.0}) {
                    DiscParams q = p; q.k = kk; q.preload = rad2deg(p.threshold() / k_si(kk));
                    auto b1 = simulate(x1, dt, q, 10), b2 = simulate(x2, dt, q, 10);
                    double cm = 1e9, up = 0;
                    for (size_t i = warm; i < b1.size(); ++i) {
                        const double th = v.omega * (i - warm) * dt;
                        cm = std::min(cm, tip_distance(G, th - b1[i] * PI / 180, th + (G.phi_deg - b2[i]) * PI / 180));
                        if (b1[i] > q.amax - 0.5) up += 1;
                    }
                    fk << R[j].name << "," << kk << "," << q.preload << "," << cm << "," << up / (b1.size() - warm) << "\n";
                }
                for (double ph = 0; ph <= 180 + 1e-9; ph += 2) {
                    RowGeometry Gp = G; Gp.phi_deg = ph;
                    double cr = 1e9, cm = 1e9;
                    for (int k = 0; k < 3600; ++k) { const double t1 = 2 * PI * k / 3600.0; cr = std::min(cr, tip_distance(Gp, t1, t1 + ph * PI / 180)); }
                    for (size_t i = 0; i < a1.size(); ++i) {
                        const double th = v.omega * i * dt;
                        cm = std::min(cm, tip_distance(Gp, th - a1[i] * PI / 180, th + (ph - a2[i]) * PI / 180));
                    }
                    fp << R[j].name << "," << ph << "," << cr << "," << cm << "\n";
                }
            }
        } else {
            usage();
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "flexrotor: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
