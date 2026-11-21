// Signal processing helpers (no external dependencies).
#include "flexrotor/flexrotor.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

namespace fr {
using cd = std::complex<double>;

static std::vector<double> poly_from_roots(const std::vector<cd>& r) {
    std::vector<cd> c{1.0};
    for (auto z : r) {
        std::vector<cd> n(c.size() + 1, 0.0);
        for (size_t i = 0; i < c.size(); ++i) { n[i] += c[i]; n[i + 1] -= z * c[i]; }
        c = n;
    }
    std::vector<double> out(c.size());
    for (size_t i = 0; i < c.size(); ++i) out[i] = c[i].real();
    return out;
}

// Digital Butterworth low-pass, transfer function form (as scipy.signal.butter).
static void butter_ba(int N, double wn, std::vector<double>& b, std::vector<double>& a) {
    const double fs = 2.0;
    const double warped = 2 * fs * std::tan(PI * wn / fs);
    std::vector<cd> p;
    for (int m = -N + 1; m < N; m += 2) p.push_back(-std::exp(cd(0, PI * m / (2.0 * N))));
    double kgain = 1.0;
    for (auto& z : p) z *= warped;
    kgain = std::pow(warped, N);
    // bilinear
    const double fs2 = 2 * fs;
    std::vector<cd> pz, zz(N, -1.0);
    cd den = 1.0;
    for (auto z : p) { pz.push_back((fs2 + z) / (fs2 - z)); den *= (fs2 - z); }
    double kz = kgain * (1.0 / den).real();
    b = poly_from_roots(zz);
    for (auto& v : b) v *= kz;
    a = poly_from_roots(pz);
}

static std::vector<double> lfilter(const std::vector<double>& b, const std::vector<double>& a,
                                   const std::vector<double>& x, std::vector<double> z) {
    const size_t n = a.size();
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        double yi = b[0] * x[i] + z[0];
        for (size_t j = 1; j < n; ++j) z[j - 1] = (j < n - 1 ? z[j] : 0.0) + b[j] * x[i] - a[j] * yi;
        y[i] = yi;
    }
    return y;
}

// Steady-state initial conditions for a unit step (scipy.signal.lfilter_zi).
static std::vector<double> lfilter_zi(const std::vector<double>& b, const std::vector<double>& a) {
    const int n = (int)a.size() - 1;
    std::vector<std::vector<double>> M(n, std::vector<double>(n + 1, 0.0));
    // (I - companion(a)^T) zi = b[1:] - a[1:] b[0]; companion(a): first row -a[1:]/a[0], subdiagonal 1
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
    for (int j = 0; j < n; ++j) C[0][j] = -a[j + 1] / a[0];
    for (int i = 1; i < n; ++i) C[i][i - 1] = 1.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) M[i][j] = (i == j ? 1.0 : 0.0) - C[j][i];
    for (int i = 0; i < n; ++i) M[i][n] = b[i + 1] - a[i + 1] * b[0];
    // Gaussian elimination with partial pivoting
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
        std::swap(M[c], M[piv]);
        for (int r = 0; r < n; ++r) if (r != c) {
            double f = M[r][c] / M[c][c];
            for (int k = c; k <= n; ++k) M[r][k] -= f * M[c][k];
        }
    }
    std::vector<double> zi(n);
    for (int i = 0; i < n; ++i) zi[i] = M[i][n] / M[i][i];
    return zi;
}

std::vector<double> butter_lowpass_filtfilt(const std::vector<double>& x, double fc, double fs, int order) {
    std::vector<double> b, a;
    butter_ba(order, fc / (fs / 2.0), b, a);
    const size_t padlen = 3 * std::max(a.size(), b.size());
    if (x.size() <= padlen) throw std::runtime_error("filtfilt: signal too short");
    // odd extension
    std::vector<double> ext;
    ext.reserve(x.size() + 2 * padlen);
    for (size_t i = padlen; i >= 1; --i) ext.push_back(2 * x[0] - x[i]);
    ext.insert(ext.end(), x.begin(), x.end());
    const size_t n = x.size();
    for (size_t i = 1; i <= padlen; ++i) ext.push_back(2 * x[n - 1] - x[n - 1 - i]);
    auto zi = lfilter_zi(b, a);
    std::vector<double> z(zi.size());
    for (size_t i = 0; i < zi.size(); ++i) z[i] = zi[i] * ext.front();
    auto y = lfilter(b, a, ext, z);
    std::reverse(y.begin(), y.end());
    for (size_t i = 0; i < zi.size(); ++i) z[i] = zi[i] * y.front();
    y = lfilter(b, a, y, z);
    std::reverse(y.begin(), y.end());
    return std::vector<double>(y.begin() + padlen, y.begin() + padlen + n);
}

// Centred moving average with reflected edges (scipy.ndimage.uniform_filter1d, mode='reflect').
std::vector<double> moving_average(const std::vector<double>& x, int w) {
    const int n = (int)x.size();
    std::vector<double> y(n);
    const int left = w / 2;
    auto at = [&](int i) {
        while (i < 0 || i >= n) { if (i < 0) i = -i - 1; if (i >= n) i = 2 * n - i - 1; }
        return x[i];
    };
    for (int i = 0; i < n; ++i) {
        double s = 0;
        for (int k = i - left; k < i - left + w; ++k) s += at(k);
        y[i] = s / w;
    }
    return y;
}

// Savitzky-Golay filter with polynomial fit at the edges (scipy mode='interp').
std::vector<double> savgol(const std::vector<double>& x, int window, int poly, int deriv, double delta) {
    const int n = (int)x.size(), h = window / 2;
    auto fitcoef = [&](int npts, int center, int p, int d) {
        // least-squares weights for derivative d at position 'center' of npts samples
        std::vector<std::vector<double>> A(npts, std::vector<double>(p + 1));
        for (int i = 0; i < npts; ++i) for (int j = 0; j <= p; ++j) A[i][j] = std::pow(double(i - center), j);
        // normal matrix
        std::vector<std::vector<double>> N(p + 1, std::vector<double>(p + 1 + npts, 0.0));
        for (int r = 0; r <= p; ++r) {
            for (int c = 0; c <= p; ++c) for (int i = 0; i < npts; ++i) N[r][c] += A[i][r] * A[i][c];
            for (int i = 0; i < npts; ++i) N[r][p + 1 + i] = A[i][r];
        }
        for (int c = 0; c <= p; ++c) {
            int piv = c;
            for (int r = c + 1; r <= p; ++r) if (std::fabs(N[r][c]) > std::fabs(N[piv][c])) piv = r;
            std::swap(N[c], N[piv]);
            for (int r = 0; r <= p; ++r) if (r != c) {
                double f = N[r][c] / N[c][c];
                for (size_t k = c; k < N[r].size(); ++k) N[r][k] -= f * N[c][k];
            }
        }
        std::vector<double> w(npts);
        double fact = 1; for (int i = 2; i <= d; ++i) fact *= i;
        for (int i = 0; i < npts; ++i) w[i] = N[d][p + 1 + i] / N[d][d] * fact;
        return w;
    };
    std::vector<double> y(n, 0.0);
    const double sc = std::pow(delta, -deriv);
    auto wc = fitcoef(window, h, poly, deriv);
    for (int i = h; i < n - h; ++i) {
        double s = 0;
        for (int k = 0; k < window; ++k) s += wc[k] * x[i - h + k];
        y[i] = s * sc;
    }
    for (int i = 0; i < h && i < n; ++i) {
        auto w0 = fitcoef(window, i, poly, deriv);
        double s = 0; for (int k = 0; k < window; ++k) s += w0[k] * x[k];
        y[i] = s * sc;
        auto w1 = fitcoef(window, window - 1 - i, poly, deriv);
        s = 0; for (int k = 0; k < window; ++k) s += w1[k] * x[n - window + k];
        y[n - 1 - i] = s * sc;
    }
    return y;
}

double percentile(std::vector<double> v, double p) {
    if (v.empty()) return NAN;
    std::sort(v.begin(), v.end());
    double pos = p / 100.0 * (v.size() - 1);
    size_t i = (size_t)std::floor(pos);
    double f = pos - i;
    return i + 1 < v.size() ? v[i] * (1 - f) + v[i + 1] * f : v[i];
}
double mean(const std::vector<double>& v) { return v.empty() ? NAN : std::accumulate(v.begin(), v.end(), 0.0) / v.size(); }
double stdev(const std::vector<double>& v, int ddof) {
    double m = mean(v), s = 0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (v.size() - ddof));
}
}  // namespace fr
