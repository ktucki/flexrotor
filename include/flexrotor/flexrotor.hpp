// FlexRotor - design of the compliant tine disc coupling of a power harrow rotor from the
// rotor geometry, calibrated on soil bin records.
// Copyright (c) 2026 Karol Tucki. MIT License.
#pragma once
#include <array>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace fr {

constexpr double PI = 3.14159265358979323846;
inline double deg2rad(double d) { return d * PI / 180.0; }
inline double rad2deg(double r) { return r * 180.0 / PI; }
// spring stiffness: N*mm/deg -> N*m/rad
inline double k_si(double k_nmm_deg) { return k_nmm_deg * 1e-3 * 180.0 / PI; }

// ---------------------------------------------------------------- input data
// Column layout of a record. Defaults follow the reduced CSV written by
// tools/myrio_xlsx2csv.py: t_us, shaft_cnt, disc_cnt, cart_cnt, torque_x10, speed_x10.
struct Channels {
    int t = 0, shaft = 1, disc = 2, cart = 3, torque = 4, speed = 5;
    double counts_per_rev = 4096.0;      // quadrature counts per revolution (1024 lines x 4)
    double torque_scale = -0.1;          // N*m per raw unit (logger stores -10 x torque)
    double speed_scale = 0.1;            // rpm per raw unit
    double cart_counts_per_m = 51300.0;  // carriage encoder
    double dt_nominal = 0.002;           // s, logger loop period
};

struct RawRun {
    std::string name;
    std::vector<double> t_us, shaft, disc, cart, torque, speed;
    size_t size() const { return t_us.size(); }
};

RawRun read_csv(const std::string& path, const Channels& ch = Channels());

// Torque lost in the drive train when the rotor turns in air, as a function of the
// shaft angular velocity (rad/s). Default: quadratic fit for the test rig.
struct DriveLoss {
    double c0 = 10.72, c1 = 0.6325, c2 = -1.060e-2;
    double operator()(double om) const { return c0 + c1 * om + c2 * om * om; }
};

enum class SoilPhase { DiscAngle, CarriageSpeed };

struct ProcessOptions {
    Channels ch;
    DriveLoss loss;
    double nominal_rpm = 0;       // >0: use nominal shaft speed instead of the measured mean
    bool use_timestamps = false;  // false: uniform grid with ch.dt_nominal
    SoilPhase phase = SoilPhase::DiscAngle;
    double phase_angle_deg = 1.0;     // DiscAngle: in soil while alpha > threshold
    double phase_speed = 0.058;       // CarriageSpeed: m/s threshold
    double phase_trim = 0.2;          // CarriageSpeed: fraction trimmed at each end
    bool zero_alpha_before_travel = false;
    double lowpass_hz = 1.5;          // zero-phase 4th-order Butterworth; <=0 disables
};

struct Run {
    std::string name;
    double dt = 0.002;
    std::vector<double> t, alpha_deg, theta_rad, Md, Mg, Mf, x_m;
    std::vector<char> in_soil;
    double omega = 0;   // rad/s
    double v = 0;       // m/s (carriage, median in soil)
    double Mg_mean = 0, Mg_sd = 0, alpha_mean = 0, alpha_sd = 0, end_stop_pct = 0;
    size_t i0 = 0, i1 = 0;  // soil phase bounds used for the statistics
};

Run process(const RawRun& raw, const ProcessOptions& opt);

// ---------------------------------------------------------------- signal tools
std::vector<double> butter_lowpass_filtfilt(const std::vector<double>& x, double fc, double fs, int order = 4);
std::vector<double> moving_average(const std::vector<double>& x, int w);
std::vector<double> savgol(const std::vector<double>& x, int window, int poly, int deriv, double delta);
double percentile(std::vector<double> v, double p);
double mean(const std::vector<double>& v);
double stdev(const std::vector<double>& v, int ddof = 0);

// ---------------------------------------------------------------- disc model
struct DiscParams {
    double k = 131.8;        // spring stiffness, N*mm/deg
    double preload = 89.5;   // spring preload angle, deg
    double mu = 0.46;        // viscous damping, N*m*s/rad
    double Mc = 1.42;        // dry friction torque, N*m
    double J = 0.023;        // disc + tines inertia, kg*m^2
    double ko = 1e4;         // end-stop contact stiffness, N*m/rad
    double muo = 20.0;       // end-stop contact damping, N*m*s/rad
    double amax = 90.0;      // rotation range, deg
    double v0 = 0.05;        // friction regularisation velocity, rad/s
    bool friction = true;
    double threshold() const { return k_si(k) * deg2rad(preload); }  // M*, N*m
};

// Disc angle (deg) driven by the soil resistance torque Mg sampled at dt.
// RK4 with nsub substeps per sample, linear interpolation of Mg. Optionally returns
// the torque transmitted to the shaft.
std::vector<double> simulate(const std::vector<double>& Mg, double dt, const DiscParams& p, int nsub = 20,
                             std::vector<double>* Mt = nullptr);

// Quasi-static play model: the disc moves only when |Mg - k(alpha+beta)| > Mc.
std::vector<double> simulate_play(const std::vector<double>& Mg, const DiscParams& p);

// ---------------------------------------------------------------- identification
struct Metrics {
    double rmse = 0, mae = 0, bias = 0, r2 = 0;
    std::map<std::string, std::array<double, 4>> per;  // rmse, mae, bias, r2 per run
};
Metrics metrics(const std::vector<const Run*>& runs, const DiscParams& p, int nsub = 20);

struct FitOptions {
    std::array<double, 4> x0{135, 88, 0.45, 1.4};
    std::array<double, 4> lb{60, 40, 0.001, 0};
    std::array<double, 4> ub{300, 140, 3, 10};
    std::map<int, double> fixed;  // 0 k, 1 preload, 2 mu, 3 Mc
    bool friction = true;
    int step = 5;     // residual decimation
    int nsub = 10;    // substeps during fitting
    int max_iter = 200;
    DiscParams base;
};
struct FitResult {
    DiscParams p;
    Metrics m;
    std::vector<std::vector<double>> corr;  // correlation of the free parameters
    std::vector<int> free_idx;
    double cost = 0;
    int iterations = 0;
};
FitResult fit(const std::vector<const Run*>& runs, const FitOptions& opt);

// ---------------------------------------------------------------- analyses
struct OrderResult {
    std::string name;
    double rpm = 0, Omega = 0;
    int nrev = 0;
    double Mt_n = 0, Mg_n = 0, alpha_n_deg = 0, alpha_n_lin_deg = 0, T_exp = 0, T_lin = 0;
};
// Angle-synchronous Fourier coefficients of order n of shaft torque and disc angle.
OrderResult order_analysis(const Run& r, const DiscParams& p, int n);

struct KinResult {
    std::string name;
    double lambda = 0, k_crit = 0, dalpha_max = 0, dalpha_p95 = 0, f_min = 0, f_max = 0, f_p025 = 0, f_p975 = 0;
    double Pd = 0, Pg = 0;
};
// Bite length factor of the compliant rotor: f = (2pi/z) / (2pi/z - delta alpha over one tine period).
KinResult kinematics(const Run& r, const DiscParams& p, int z, double radius, double u);

struct SelectionMap {
    std::vector<double> ks, Ms;               // N*mm/deg, N*m (threshold torque)
    std::vector<std::vector<double>> low, mid, up;  // time fractions at initial stop / free / end stop
};
SelectionMap selection_map(const std::vector<const Run*>& runs, const DiscParams& p,
                           const std::vector<double>& ks, const std::vector<double>& Ms, int nsub = 5);

// ---------------------------------------------------------------- load model (mean values)
struct Variant {
    std::string name, setting;
    double omega = 0, u = 0, M = 0, alpha = 0, Msd = 0;
};
struct Characteristic { double M0 = 0, k = 0; };  // M = M0 + k*alpha, k in N*m/deg
Characteristic characteristic(const std::vector<Variant>& v, double amin = 1.0, double amax = 88.0);
// Characteristic identified in soil: Mg and alpha smoothed with a 0.2 s moving average,
// grouped into angle bins; a line is fitted to the bin medians of the torque within [lo, hi].
struct Bin { double center, q25, median, q75; size_t n; };
Characteristic characteristic_binned(const std::vector<const Run*>& runs, double width = 5, double lo = 10,
                                     double hi = 85, std::vector<Bin>* bins = nullptr, size_t min_n = 500);

struct LoadModel {
    double Ma = 0, K = 0, e = 0;                 // M = Ma + d_s + K*s_k ; e = pi*K/(h*b)
    std::map<std::string, double> offset;        // per setting
    std::vector<std::string> settings;
    double h = 0.14, b = 0.25;                   // working depth and tine width, m
    int z = 2;                                   // tines per disc
    double bite(double u, double omega) const { return 2 * PI * u / (z * omega) * 1e3; }  // mm
    double torque(const std::string& s, double u, double omega) const;
};
LoadModel fit_load_model(const std::vector<Variant>& v, double h, double b, int z);
struct LooResult { std::vector<double> M_pred, alpha_pred; double rmse_alpha = 0, r2_alpha = 0, rmse_M = 0; };
LooResult loo_chain(const std::vector<Variant>& v, const std::map<std::string, Characteristic>& ch, double h, double b, int z);

struct SelectionRule { double M_lo = 0, M_hi = 0, k_min = 0; };  // k_min in N*mm/deg
SelectionRule selection_rule(double Ma, double K, double s_min, double s_max, double sigma, double z_q,
                             double amax_deg, double soil_scale = 1.0);

// ---------------------------------------------------------------- output
struct Series { std::vector<double> x, y; std::string label, color = "#000000"; std::string dash; bool markers = false; bool open = false; };
struct PlotSpec {
    std::string xlabel, ylabel;
    double xmin = 0, xmax = 1, xstep = 0.2, ymin = 0, ymax = 1, ystep = 0.2;
    double width = 340, height = 250;
    std::vector<Series> series;
    std::vector<std::pair<std::pair<double, double>, std::string>> notes;
};
void write_svg(const std::string& path, const PlotSpec& s);
void write_map_svg(const std::string& path, const SelectionMap& m, const std::vector<std::pair<double, double>>& points);

// ---------------------------------------------------------------- rotor pair (prediction)
// Two adjacent, counter-rotating rotors driven by a common gear train: the shaft of rotor 2 is
// locked to rotor 1 with the tine phase offset phi0. Rotor 1 is the measured record; rotor 2 is
// the same record shifted in time by 'shift' seconds, which sets the correlation of the two soil
// loads (0 s: identical loads). The soil interaction between the rotors is NOT represented:
// every rotor carries the load measured on a single rotor. Results are predictions.
struct PairResult {
    std::string name;
    double shift = 0, corr = 0;                       // s; correlation of the filtered soil torques
    double da_mean = 0, da_sd = 0, da_lo = 0, da_hi = 0, da_max = 0;  // alpha1 - alpha2, deg (2.5/97.5 %)
    int nrev = 0;
    double rho = 0, rho_rigid = 0;                    // order-n torque of the pair / sum of both rotors
    double dmin = -1, dmin_rigid = -1;                // minimum tip-to-tip distance, m (-1: not computed)
    std::vector<double> rev_t, rev_dmin;              // minimum tip distance in each shaft revolution (s from start, m)
    std::vector<double> phi_deg, ratio, ratio_rigid;  // rho over the phase offset 0..180 deg
};
// alpha2 (optional): disc angle of rotor 2 from a simulation; otherwise the measured angle is used
// for both rotors. 'spacing' and 'radius' (tip radius) enable the tip distance (spacing < 2 radius).
PairResult pair_analysis(const Run& r, double shift, double phi0_deg, int n, int z, double spacing, double radius,
                         const std::vector<double>* alpha_sim = nullptr);

// ---------------------------------------------------------------- design path (from geometry)
// Slow fluctuation of the filtered soil torque over the middle 60 % of the soil phase:
// standard deviation and correlation time (lag at which the autocorrelation falls to 1/e).
struct Fluct { double sd = 0, tau = 0, mean = 0; };
Fluct fluctuation(const Run& r);

// Normal random numbers that are identical on every compiler (mt19937 + Box-Muller).
class Rng {
public:
    explicit Rng(unsigned seed = 1);
    double normal();
private:
    std::mt19937 m_;
    double uniform();
    bool has_ = false; double next_ = 0;
};
// Zero-mean Ornstein-Uhlenbeck process with standard deviation sd and correlation time tau.
std::vector<double> ou_process(double sd, double tau, double dt, size_t n, Rng& g);

// Statistics of a disc angle series (deg) and the time fractions at the stops.
struct DiscStats { double mean = 0, sd = 0, p05 = 0, p95 = 0, low = 0, mid = 0, up = 0; };
DiscStats disc_stats(const std::vector<double>& a, double amax, size_t i0 = 0, size_t i1 = 0);

// Torsion spring with straight legs after EN 13906-3 (body only, legs neglected).
struct TorsionSpring {
    double d = 5.0, D = 46.5, n = 5.5;   // wire diameter, mean coil diameter (mm), active coils
    double E = 206000;                   // modulus of elasticity, N/mm^2
    double Rm = 0;                       // tensile strength, N/mm^2 (0: EN 10270-1 SH/DH minimum)
};
double wire_rm_dh(double d);                       // EN 10270-1 grade SH/DH minimum tensile strength, 3-8 mm
double spring_rate(const TorsionSpring& s);        // N*mm/deg, k = E d^4 / (3667.3 D n)
double coils_for_rate(double k, double d, double D, double E);
struct SpringCheck {
    double k = 0, beta = 0;                  // N*mm/deg; preload angle for the threshold torque, deg
    double M_end = 0, sigma_end = 0;         // spring torque at the end stop (N*m) and bending stress (N/mm^2)
    double sigma_lo = 0, sigma_hi = 0;       // bending stress at the ends of the operating torque range
    double q = 1, sigma_allow = 0;           // stress correction factor; permissible static stress 0.7 Rm
    double Di_end = 0, Do = 0, Lk = 0;       // inner diameter at the end stop, outer diameter, body length (mm)
    double mass = 0;                         // kg, active coils
    bool stress_ok = false;
    // Fatigue (Gerber criterion, corrected stresses): (1) slow fluctuation of amplitude M_fl about the highest mean
    // torque M_hi - M_fl, 1e6 cycles; (2) whole operating range M_lo-M_hi, 1e5 cycles. Fatigue strength of a repeated
    // stress S_r = r6 Rm (1e6) or r5 Rm (1e5) of torsion springs from patented wire, not shot peened (SMI handbook).
    double sigma_fa = 0, sigma_fm = 0;       // amplitude and mean stress of the fluctuation
    double nf_fl = 0, nf_range = 0;          // fatigue safety factors (0: not evaluated)
};
SpringCheck check_spring(const TorsionSpring& s, double Mstar, double M_lo, double M_hi, double amax_deg,
                         double stress_factor = 0.7, double M_fl = 0, double r6 = 0.50, double r5 = 0.53);
// Gerber criterion: endurance amplitude from the strength of a repeated stress S_r, and the safety factor
double gerber_se(double Sr, double Sut);
double gerber_nf(double sa, double sm, double Se, double Sut);
// For every standard wire diameter: the spring of the smallest outer diameter within the limits that has a
// rate of at least k_min (active coils in quarter turns), a static bending stress at the end stop within
// stress_factor * Rm, an inner diameter at the end stop of at least Di_min and a body length up to L_max.
struct SpringLimits { double Do_max = 70, Di_min = 37.7, L_max = 40, D_step = 0.5; std::vector<double> d_list; };
std::vector<std::pair<TorsionSpring, SpringCheck>> design_springs(double k_min, double Mstar, double M_lo, double M_hi,
                                                                  double amax_deg, const SpringLimits& lim,
                                                                  double E = 206000, double stress_factor = 0.7);

// Two adjacent counter-rotating rotors of a row (rotor 2 is the mirror image of rotor 1 about the
// line midway between the axes) moving forward at speed u. Tine tip paths in the soil, minimum
// tip-to-tip distance in every shaft revolution and the largest soil element not crossed by any
// tip path (diameter of the largest empty circle; the row is laterally periodic with period 2a).
struct RowGeometry { double r = 0.11, a = 0.217, phi_deg = 90; int z = 2; };
struct RowResult {
    double clear_min = 0, clear_rigid = 0;   // m
    double uncut = 0, uncut_rigid = 0;       // m
    std::vector<double> rev_clear;           // m, per shaft revolution
    std::vector<std::array<double, 3>> path; // x, y (m), rotor + 0.1 tine: tip paths around the largest uncut element
    std::array<double, 3> uncut_circle{0, 0, 0}; // centre x, y and diameter (m) of the largest empty circle
};
RowResult row_analysis(const RowGeometry& g, double omega, double u, const std::vector<double>& alpha1_deg,
                       const std::vector<double>& alpha2_deg, double dt, double res = 0.001, double plot_len = 0,
                       double uncut_len = 2.0);
double tip_distance(const RowGeometry& g, double psi1, double psi2);   // minimum over the tine pairs, m

// Test rig and rotor description read from a plain 'key = value' file, so that records of
// another logger, rotor or coupling can be processed without recompiling. Unknown keys are
// rejected; keys that are not given keep their defaults (the SGGW soil bin rotor).
struct RigConfig {
    Channels ch;
    DriveLoss loss;
    double lowpass_hz = 1.5, phase_angle_deg = 1.0, phase_speed = 0.058, phase_trim = 0.2;
    DiscParams disc;            // J, ko, muo, amax, v0 and the start values of k, preload, mu, Mc
    int z = 2;                  // tines per disc
    double radius = 0.11, depth = 0.14, width = 0.25;  // m
    double spacing = 0.217;     // m, axis distance of adjacent rotors (tested harrow: 1.3 m / 6 rotors; 0: not given)
    double phase_deg = 0;       // deg, tine phase offset of adjacent rotors (0: 180/z)
    TorsionSpring spring;       // spring of the coupling (default: prototype spring of setting N2)
    SpringLimits limits;        // envelope for the spring design
    double stress_factor = 0.7; // permissible static bending stress / Rm
    double fatigue_r6 = 0.50, fatigue_r5 = 0.53;   // strength of a repeated stress / Rm at 1e6 and 1e5 cycles
    double load_corr = 0;       // correlation of the slow soil loads of adjacent rotors (design path)
};
RigConfig read_config(const std::string& path);

std::string version();
}  // namespace fr
