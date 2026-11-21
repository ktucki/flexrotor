// Minimal SVG plots: thin lines, outward ticks, axes ending on a labelled tick,
// legend above the plot area without a frame.
#include "flexrotor/flexrotor.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace fr {

static std::string fmt(double v, double step) {
    char b[32];
    int dec = 0;
    if (step < 1) dec = (int)std::ceil(-std::log10(step) - 1e-9);
    std::snprintf(b, sizeof b, "%.*f", dec, std::fabs(v) < 1e-12 ? 0.0 : std::fabs(v));
    return (v < -1e-12 ? std::string("\u2212") : std::string()) + b;
}
static std::string esc(const std::string& s) {
    std::string o;
    for (char c : s) { if (c == '<') o += "&lt;"; else if (c == '>') o += "&gt;"; else if (c == '&') o += "&amp;"; else o += c; }
    return o;
}

void write_svg(const std::string& path, const PlotSpec& s) {
    const double L = 52, R = 14, T = 18 + 14.0 * ((s.series.size() + 1) / 2), B = 40;
    const double W = s.width, H = s.height + T;
    const double pw = W - L - R, ph = H - T - B;
    auto X = [&](double x) { return L + (x - s.xmin) / (s.xmax - s.xmin) * pw; };
    auto Y = [&](double y) { return T + ph - (y - s.ymin) / (s.ymax - s.ymin) * ph; };
    std::ostringstream o;
    o << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W << "' height='" << H << "' viewBox='0 0 " << W << " " << H
      << "' font-family='Times New Roman, Liberation Serif, serif' font-size='9'>\n<rect width='100%' height='100%' fill='white'/>\n";
    // light horizontal grid
    for (double y = s.ymin; y <= s.ymax + 1e-9; y += s.ystep)
        o << "<line x1='" << L << "' x2='" << L + pw << "' y1='" << Y(y) << "' y2='" << Y(y) << "' stroke='#d9d9d9' stroke-width='0.4'/>\n";
    // axes
    o << "<line x1='" << L << "' x2='" << L + pw << "' y1='" << T + ph << "' y2='" << T + ph << "' stroke='black' stroke-width='0.6'/>\n";
    o << "<line x1='" << L << "' x2='" << L << "' y1='" << T << "' y2='" << T + ph << "' stroke='black' stroke-width='0.6'/>\n";
    for (double x = s.xmin; x <= s.xmax + 1e-9; x += s.xstep) {
        o << "<line x1='" << X(x) << "' x2='" << X(x) << "' y1='" << T + ph << "' y2='" << T + ph + 3 << "' stroke='black' stroke-width='0.6'/>";
        o << "<text x='" << X(x) << "' y='" << T + ph + 13 << "' text-anchor='middle'>" << fmt(x, s.xstep) << "</text>\n";
    }
    for (double y = s.ymin; y <= s.ymax + 1e-9; y += s.ystep) {
        o << "<line x1='" << L - 3 << "' x2='" << L << "' y1='" << Y(y) << "' y2='" << Y(y) << "' stroke='black' stroke-width='0.6'/>";
        o << "<text x='" << L - 5 << "' y='" << Y(y) + 3 << "' text-anchor='end'>" << fmt(y, s.ystep) << "</text>\n";
    }
    o << "<text x='" << L + pw / 2 << "' y='" << H - 8 << "' text-anchor='middle'>" << esc(s.xlabel) << "</text>\n";
    o << "<text transform='translate(12," << T + ph / 2 << ") rotate(-90)' text-anchor='middle'>" << esc(s.ylabel) << "</text>\n";
    o << "<clipPath id='c'><rect x='" << L << "' y='" << T << "' width='" << pw << "' height='" << ph << "'/></clipPath>\n";
    int li = 0;
    for (auto& se : s.series) {
        if (!se.markers) {
            o << "<polyline clip-path='url(#c)' fill='none' stroke='" << se.color << "' stroke-width='0.6'";
            if (!se.dash.empty()) o << " stroke-dasharray='" << se.dash << "'";
            o << " points='";
            for (size_t i = 0; i < se.x.size(); ++i) o << X(se.x[i]) << "," << Y(se.y[i]) << " ";
            o << "'/>\n";
        } else {
            for (size_t i = 0; i < se.x.size(); ++i)
                o << "<circle cx='" << X(se.x[i]) << "' cy='" << Y(se.y[i]) << "' r='2.2' stroke='" << se.color
                  << "' stroke-width='0.6' fill='" << (se.open ? "white" : se.color) << "'/>\n";
        }
        // legend entry
        const double lx = L - 30 + (li % 2) * (pw / 2 + 20), ly = 10 + 14.0 * (li / 2);
        if (se.markers)
            o << "<circle cx='" << lx + 9 << "' cy='" << ly << "' r='2.2' stroke='" << se.color << "' stroke-width='0.6' fill='"
              << (se.open ? "white" : se.color) << "'/>";
        else {
            o << "<line x1='" << lx << "' x2='" << lx + 18 << "' y1='" << ly << "' y2='" << ly << "' stroke='" << se.color << "' stroke-width='0.6'";
            if (!se.dash.empty()) o << " stroke-dasharray='" << se.dash << "'";
            o << "/>";
        }
        o << "<text x='" << lx + 23 << "' y='" << ly + 3 << "'>" << esc(se.label) << "</text>\n";
        ++li;
    }
    for (auto& n : s.notes) o << "<text x='" << X(n.first.first) << "' y='" << Y(n.first.second) << "' font-size='8'>" << esc(n.second) << "</text>\n";
    o << "</svg>\n";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << o.str();
}

// Spring selection map: each cell is coloured by the state of the disc over the soil phase
// (green: compliant operation, blue: on the initial stop > 5% of the time,
// orange: on the end stop > 5% of the time).
void write_map_svg(const std::string& path, const SelectionMap& m, const std::vector<std::pair<double, double>>& points) {
    PlotSpec s;
    s.xlabel = "Spring stiffness k (N·mm/°)"; s.ylabel = "Threshold torque M* (N·m)";
    s.xstep = 50; s.xmin = std::floor(m.ks.front() / s.xstep) * s.xstep; s.xmax = std::ceil(m.ks.back() / s.xstep) * s.xstep;
    s.ystep = 4; s.ymin = std::floor(m.Ms.front() / s.ystep) * s.ystep; s.ymax = std::ceil(m.Ms.back() / s.ystep) * s.ystep;
    const double L = 52, R = 14, T = 18 + 28, B = 40, W = s.width, H = s.height + T;
    const double pw = W - L - R, ph = H - T - B;
    auto X = [&](double x) { return L + (x - s.xmin) / (s.xmax - s.xmin) * pw; };
    auto Y = [&](double y) { return T + ph - (y - s.ymin) / (s.ymax - s.ymin) * ph; };
    std::ostringstream o;
    o << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W << "' height='" << H << "' font-family='Times New Roman, Liberation Serif, serif' font-size='9'>\n"
      << "<rect width='100%' height='100%' fill='white'/>\n";
    const double dx = (m.ks.size() > 1 ? m.ks[1] - m.ks[0] : 1) / 2, dy = (m.Ms.size() > 1 ? m.Ms[1] - m.Ms[0] : 1) / 2;
    for (size_t i = 0; i < m.Ms.size(); ++i)
        for (size_t j = 0; j < m.ks.size(); ++j) {
            std::string c = "#FFFFFF";
            if (m.up[i][j] > 0.05) c = "#F8CBAD";        // end stop
            else if (m.low[i][j] > 0.05) c = "#BDD7EE";  // initial stop
            else c = "#E2EFDA";                           // compliant operation
            double x0 = std::max(X(m.ks[j] - dx), L), x1 = std::min(X(m.ks[j] + dx), L + pw);
            double y0 = std::max(Y(m.Ms[i] + dy), T), y1 = std::min(Y(m.Ms[i] - dy), T + ph);
            o << "<rect x='" << x0 << "' y='" << y0 << "' width='" << x1 - x0 << "' height='" << y1 - y0 << "' fill='" << c << "' stroke='none'/>\n";
        }
    o << "<line x1='" << L << "' x2='" << L + pw << "' y1='" << T + ph << "' y2='" << T + ph << "' stroke='black' stroke-width='0.6'/>\n";
    o << "<line x1='" << L << "' x2='" << L << "' y1='" << T << "' y2='" << T + ph << "' stroke='black' stroke-width='0.6'/>\n";
    for (double x = s.xmin; x <= s.xmax + 1e-9; x += s.xstep)
        o << "<line x1='" << X(x) << "' x2='" << X(x) << "' y1='" << T + ph << "' y2='" << T + ph + 3 << "' stroke='black' stroke-width='0.6'/><text x='" << X(x)
          << "' y='" << T + ph + 13 << "' text-anchor='middle'>" << fmt(x, s.xstep) << "</text>\n";
    for (double y = s.ymin; y <= s.ymax + 1e-9; y += s.ystep)
        o << "<line x1='" << L - 3 << "' x2='" << L << "' y1='" << Y(y) << "' y2='" << Y(y) << "' stroke='black' stroke-width='0.6'/><text x='" << L - 5
          << "' y='" << Y(y) + 3 << "' text-anchor='end'>" << fmt(y, s.ystep) << "</text>\n";
    o << "<text x='" << L + pw / 2 << "' y='" << H - 8 << "' text-anchor='middle'>" << esc(s.xlabel) << "</text>\n";
    o << "<text transform='translate(12," << T + ph / 2 << ") rotate(-90)' text-anchor='middle'>" << esc(s.ylabel) << "</text>\n";
    const char* lab[3] = {"compliant operation", "initial stop &gt; 5% of time", "end stop &gt; 5% of time"};
    const char* col[3] = {"#E2EFDA", "#BDD7EE", "#F8CBAD"};
    for (int i = 0; i < 3; ++i)
        o << "<rect x='" << L - 30 + i * 105 << "' y='" << 8 << "' width='10' height='8' fill='" << col[i] << "'/><text x='" << L - 16 + i * 105
          << "' y='" << 15 << "'>" << lab[i] << "</text>\n";
    for (auto& p : points)
        o << "<circle cx='" << X(p.first) << "' cy='" << Y(p.second) << "' r='2.5' fill='black'/>\n";
    o << "</svg>\n";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << o.str();
}
}  // namespace fr
