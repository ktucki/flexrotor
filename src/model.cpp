// Dynamic model of the tine disc mounted on the shaft through a preloaded torsion
// spring with end stops, viscous damping and dry friction.
//   J a'' = Mg - k (a + beta) - mu a' - Mc tanh(a'/v0) - Mo(a, a')
// a: disc lag angle relative to the shaft, beta: spring preload angle.
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>

namespace fr {

namespace {
struct Coeffs {
    double k, b, mu, Mc, v0, J, ko, muo, amax;
    bool fric;
};
inline double stop_torque(double a, double v, const Coeffs& c) {
    if (a < 0) return c.ko * a + c.muo * v;
    if (a > c.amax) return c.ko * (a - c.amax) + c.muo * v;
    return 0.0;
}
inline double accel(double Mg, double a, double v, const Coeffs& c) {
    double fr = c.fric ? c.Mc * std::tanh(v / c.v0) : 0.0;
    return (Mg - c.k * (a + c.b) - c.mu * v - fr - stop_torque(a, v, c)) / c.J;
}
}  // namespace

std::vector<double> simulate(const std::vector<double>& Mg, double dt, const DiscParams& p, int nsub,
                             std::vector<double>* Mt) {
    Coeffs c{k_si(p.k), deg2rad(p.preload), p.mu, p.Mc, p.v0, p.J, p.ko, p.muo, deg2rad(p.amax), p.friction};
    const size_t n = Mg.size();
    std::vector<double> out(n, 0.0);
    if (Mt) Mt->assign(n, 0.0);
    double a = 0, v = 0;
    const double h = dt / nsub;
    if (Mt) (*Mt)[0] = c.k * (a + c.b) + c.mu * v + stop_torque(a, v, c);
    for (size_t i = 0; i + 1 < n; ++i) {
        const double g = Mg[i], dg = Mg[i + 1] - Mg[i];
        for (int s = 0; s < nsub; ++s) {
            const double g0 = g + dg * s / nsub, gm = g + dg * (s + 0.5) / nsub, g1 = g + dg * (s + 1.0) / nsub;
            const double k1a = v, k1v = accel(g0, a, v, c);
            const double k2a = v + 0.5 * h * k1v, k2v = accel(gm, a + 0.5 * h * k1a, v + 0.5 * h * k1v, c);
            const double k3a = v + 0.5 * h * k2v, k3v = accel(gm, a + 0.5 * h * k2a, v + 0.5 * h * k2v, c);
            const double k4a = v + h * k3v, k4v = accel(g1, a + h * k3a, v + h * k3v, c);
            a += h / 6 * (k1a + 2 * k2a + 2 * k3a + k4a);
            v += h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);
        }
        out[i + 1] = rad2deg(a);
        if (Mt) (*Mt)[i + 1] = c.k * (a + c.b) + c.mu * v + stop_torque(a, v, c);
    }
    return out;
}

std::vector<double> simulate_play(const std::vector<double>& Mg, const DiscParams& p) {
    const double kr = k_si(p.k), br = deg2rad(p.preload), amax = deg2rad(p.amax);
    std::vector<double> q(Mg.size());
    double a = 0;
    for (size_t i = 0; i < Mg.size(); ++i) {
        const double lo = (Mg[i] - p.Mc) / kr - br, hi = (Mg[i] + p.Mc) / kr - br;
        a = std::min(std::max(a, lo), hi);
        a = std::min(std::max(a, 0.0), amax);
        q[i] = rad2deg(a);
    }
    return q;
}
}  // namespace fr
