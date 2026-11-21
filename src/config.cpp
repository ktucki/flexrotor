// Rig and rotor configuration file ('key = value', '#' starts a comment).
#include "flexrotor/flexrotor.hpp"
#include <fstream>
#include <map>
#include <stdexcept>

namespace fr {

static std::string trim(const std::string& s) {
    const char* ws = " \t\r";
    auto a = s.find_first_not_of(ws);
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(ws) - a + 1);
}

RigConfig read_config(const std::string& path) {
    RigConfig c;
    std::map<std::string, double*> d = {
        {"counts_per_rev", &c.ch.counts_per_rev}, {"torque_scale", &c.ch.torque_scale},
        {"speed_scale", &c.ch.speed_scale},       {"cart_counts_per_m", &c.ch.cart_counts_per_m},
        {"dt_nominal", &c.ch.dt_nominal},
        {"loss_c0", &c.loss.c0}, {"loss_c1", &c.loss.c1}, {"loss_c2", &c.loss.c2},
        {"lowpass_hz", &c.lowpass_hz}, {"phase_angle_deg", &c.phase_angle_deg},
        {"phase_speed", &c.phase_speed}, {"phase_trim", &c.phase_trim},
        {"J", &c.disc.J}, {"ko", &c.disc.ko}, {"muo", &c.disc.muo}, {"amax", &c.disc.amax}, {"v0", &c.disc.v0},
        {"k", &c.disc.k}, {"preload", &c.disc.preload}, {"mu", &c.disc.mu}, {"Mc", &c.disc.Mc},
        {"radius", &c.radius}, {"depth", &c.depth}, {"width", &c.width}, {"spacing", &c.spacing},
        {"phase_deg", &c.phase_deg}, {"spring_d", &c.spring.d}, {"spring_D", &c.spring.D}, {"spring_n", &c.spring.n},
        {"spring_E", &c.spring.E}, {"spring_Rm", &c.spring.Rm}, {"Do_max", &c.limits.Do_max}, {"Di_min", &c.limits.Di_min},
        {"L_max", &c.limits.L_max}, {"stress_factor", &c.stress_factor}, {"load_corr", &c.load_corr},
        {"fatigue_r6", &c.fatigue_r6}, {"fatigue_r5", &c.fatigue_r5}};
    std::map<std::string, int*> n = {
        {"col_t", &c.ch.t}, {"col_shaft", &c.ch.shaft}, {"col_disc", &c.ch.disc}, {"col_cart", &c.ch.cart},
        {"col_torque", &c.ch.torque}, {"col_speed", &c.ch.speed}, {"z", &c.z}};
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open config " + path);
    std::string line;
    int ln = 0;
    while (std::getline(f, line)) {
        ++ln;
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        line = trim(line);
        if (line.empty()) continue;
        auto q = line.find('=');
        if (q == std::string::npos) throw std::runtime_error(path + ":" + std::to_string(ln) + ": expected key = value");
        std::string k = trim(line.substr(0, q)), v = trim(line.substr(q + 1));
        try {
            if (d.count(k)) *d[k] = std::stod(v);
            else if (n.count(k)) *n[k] = std::stoi(v);
            else throw std::runtime_error("unknown key '" + k + "'");
        } catch (const std::invalid_argument&) {
            throw std::runtime_error(path + ":" + std::to_string(ln) + ": bad value for '" + k + "'");
        } catch (const std::runtime_error& e) {
            throw std::runtime_error(path + ":" + std::to_string(ln) + ": " + e.what());
        }
    }
    if (c.ch.counts_per_rev <= 0 || c.ch.dt_nominal <= 0 || c.disc.J <= 0 || c.disc.amax <= 0 || c.z < 1)
        throw std::runtime_error(path + ": counts_per_rev, dt_nominal, J, amax and z must be positive");
    if (c.spring.d <= 0 || c.spring.D <= c.spring.d || c.spring.n <= 0 || c.spring.E <= 0 || c.stress_factor <= 0 || c.fatigue_r6 <= 0 || c.fatigue_r5 <= 0 || c.fatigue_r6 >= 1 || c.fatigue_r5 >= 1 ||
        c.load_corr < -1 || c.load_corr > 1)
        throw std::runtime_error(path + ": spring_d, spring_D > spring_d, spring_n, spring_E, stress_factor must be positive, 0 < fatigue_r6, fatigue_r5 < 1, |load_corr| <= 1");
    return c;
}
}  // namespace fr
