// Reading soil bin records and separating the soil resistance torque.
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fr {

std::string version() { return "1.0.0"; }

static bool parse_row(const std::string& line, std::vector<double>& v) {
    v.clear();
    std::string cell;
    std::string s = line;
    for (char& c : s) if (c == ';' || c == '\t') c = ',';
    std::stringstream ss(s);
    while (std::getline(ss, cell, ',')) {
        // accept decimal comma inside quotes is not supported; logger exports use dots
        if (cell.empty()) { v.push_back(NAN); continue; }
        char* end = nullptr;
        double d = std::strtod(cell.c_str(), &end);
        if (end == cell.c_str()) return false;
        v.push_back(d);
    }
    return !v.empty();
}

RawRun read_csv(const std::string& path, const Channels& ch) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    RawRun r;
    auto p = path.find_last_of("/\\");
    r.name = path.substr(p == std::string::npos ? 0 : p + 1);
    if (r.name.size() > 4 && r.name.substr(r.name.size() - 4) == ".csv") r.name.resize(r.name.size() - 4);
    std::string line;
    std::vector<double> v;
    const int need = std::max({ch.t, ch.shaft, ch.disc, ch.cart, ch.torque, ch.speed});
    while (std::getline(f, line)) {
        if (!parse_row(line, v)) continue;  // header or text
        if ((int)v.size() <= need) continue;
        bool ok = true;
        for (int c : {ch.t, ch.shaft, ch.disc, ch.cart, ch.torque, ch.speed}) if (std::isnan(v[c])) ok = false;
        if (!ok) continue;
        r.t_us.push_back(v[ch.t]); r.shaft.push_back(v[ch.shaft]); r.disc.push_back(v[ch.disc]);
        r.cart.push_back(v[ch.cart]); r.torque.push_back(v[ch.torque]); r.speed.push_back(v[ch.speed]);
    }
    if (r.size() < 1000) throw std::runtime_error("too few samples in " + path);
    return r;
}

Run process(const RawRun& raw, const ProcessOptions& o) {
    const auto& ch = o.ch;
    const size_t n = raw.size();
    Run r;
    r.name = raw.name;
    r.dt = ch.dt_nominal;
    r.t.resize(n); r.alpha_deg.resize(n); r.theta_rad.resize(n); r.Md.resize(n); r.x_m.resize(n);
    // time base: logger timestamps in microseconds with 32-bit wrap-around
    r.t[0] = 0;
    for (size_t i = 1; i < n; ++i) {
        if (!o.use_timestamps) { r.t[i] = i * ch.dt_nominal; continue; }
        double d = raw.t_us[i] - raw.t_us[i - 1];
        if (d < -2147483648.0) d += 4294967296.0;
        if (d <= 0 || d > 1e5) d = ch.dt_nominal * 1e6;
        r.t[i] = r.t[i - 1] + d * 1e-6;
    }
    std::vector<double> rpm(n);
    for (size_t i = 0; i < n; ++i) {
        // the disc encoder counts in the opposite sense; alpha = disc angle relative to the shaft
        r.alpha_deg[i] = -(raw.shaft[i] + raw.disc[i]) / ch.counts_per_rev * 360.0;
        r.theta_rad[i] = std::fabs(raw.shaft[i]) / ch.counts_per_rev * 2 * PI;
        r.Md[i] = raw.torque[i] * ch.torque_scale;
        rpm[i] = raw.speed[i] * ch.speed_scale;
        r.x_m[i] = std::fabs((raw.cart[i] - raw.cart[0]) / ch.cart_counts_per_m);
    }
    // carriage speed and travel phase
    auto xs = moving_average(r.x_m, 250);
    std::vector<double> v(n);
    for (size_t i = 0; i < n; ++i) {
        size_t a = i ? i - 1 : 0, b = i + 1 < n ? i + 1 : n - 1;
        v[i] = std::fabs((xs[b] - xs[a]) / (r.t[b] - r.t[a]));
    }
    size_t t0 = n, t1 = 0, cnt = 0;
    for (size_t i = 0; i < n; ++i) if (v[i] > o.phase_speed) { t0 = std::min(t0, i); t1 = i; ++cnt; }
    const bool travel = cnt >= 500;
    if (o.zero_alpha_before_travel && travel) {
        std::vector<double> pre;
        for (size_t i = 0; i < n; ++i) if (r.t[i] < r.t[t0] - 0.2 && r.t[i] > r.t[t0] - 3.2) pre.push_back(r.alpha_deg[i]);
        if (pre.size() <= 100) pre.assign(r.alpha_deg.begin(), r.alpha_deg.begin() + std::max<size_t>(t0, 1));
        double a0 = percentile(pre, 50);
        for (auto& a : r.alpha_deg) a -= a0;
    }
    r.in_soil.assign(n, 0);
    if (o.phase == SoilPhase::DiscAngle) {
        for (size_t i = 0; i < n; ++i) r.in_soil[i] = r.alpha_deg[i] > o.phase_angle_deg;
        size_t a = n, b = 0;
        for (size_t i = 0; i < n; ++i) if (r.in_soil[i]) { a = std::min(a, i); b = i; }
        r.i0 = a; r.i1 = b;
    } else {
        if (!travel) throw std::runtime_error(r.name + ": no carriage travel found");
        size_t L = t1 - t0;
        r.i0 = t0 + (size_t)(o.phase_trim * L);
        r.i1 = t0 + (size_t)((1 - o.phase_trim) * L);
        for (size_t i = r.i0; i < r.i1; ++i) r.in_soil[i] = 1;
    }
    // shaft speed and drive-train loss
    if (o.nominal_rpm > 0) r.omega = o.nominal_rpm * 2 * PI / 60;
    else {
        double s = 0; size_t c = 0;
        for (size_t i = r.i0; i < r.i1; ++i) { s += rpm[i]; ++c; }
        r.omega = (c ? s / c : 0) * 2 * PI / 60;
    }
    const double Ms = o.loss(r.omega);
    r.Mg.resize(n);
    for (size_t i = 0; i < n; ++i) r.Mg[i] = r.Md[i] - Ms;
    r.Mf = o.lowpass_hz > 0 ? butter_lowpass_filtfilt(r.Mg, o.lowpass_hz, 1.0 / ch.dt_nominal, 4) : r.Mg;
    // statistics in soil
    std::vector<double> mg, al, vv;
    for (size_t i = 0; i < n; ++i) if (r.in_soil[i]) { mg.push_back(r.Mg[i]); al.push_back(r.alpha_deg[i]); vv.push_back(v[i]); }
    if (mg.empty()) throw std::runtime_error(r.name + ": empty soil phase");
    r.Mg_mean = mean(mg); r.Mg_sd = stdev(mg); r.alpha_mean = mean(al); r.alpha_sd = stdev(al);
    r.v = percentile(vv, 50);
    size_t es = 0; for (double a : al) if (a > 88.0) ++es;
    r.end_stop_pct = 100.0 * es / al.size();
    return r;
}
}  // namespace fr
