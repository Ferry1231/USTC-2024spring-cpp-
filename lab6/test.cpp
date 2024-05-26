#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;
const double PI = std::acos(-1);

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;

// Recursive FFT function
ComplexVector FFT(const ComplexVector &f) {
    int n = f.size();
    if (n == 1) {
        return f;
    }
    
    Complex w_n(std::cos(2 * PI / n), -std::sin(2 * PI / n));
    Complex w(1, 0);
    
    ComplexVector f0(n / 2), f1(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        f0[i] = f[2 * i];
        f1[i] = f[2 * i + 1];
    }
    
    ComplexVector g0 = FFT(f0);
    ComplexVector g1 = FFT(f1);
    
    ComplexVector g(n);
    for (int k = 0; k < n / 2; ++k) {
        g[k] = (g0[k] + w * g1[k]) / 2.0;
        g[k + n / 2] = (g0[k] - w * g1[k]) / 2.0;
        w *= w_n;
    }
    
    return g;
}

// Recursive IFFT function
ComplexVector IFFT(const ComplexVector &f) {
    int n = f.size();
    if (n == 1) {
        return f;
    }
    
    Complex w_n(std::cos(2 * PI / n), std::sin(2 * PI / n));
    Complex w(1, 0);
    
    ComplexVector f0(n / 2), f1(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        f0[i] = f[2 * i];
        f1[i] = f[2 * i + 1];
    }
    
    ComplexVector g0 = IFFT(f0);
    ComplexVector g1 = IFFT(f1);
    
    ComplexVector g(n);
    for (int k = 0; k < n / 2; ++k) {
        g[k] = g0[k] + w * g1[k];
        g[k + n / 2] = g0[k] - w * g1[k];
        w *= w_n;
    }
    
    return g;
}

int main() {
    // Example usage
    ComplexVector f = {1, 1, 1, 1, 0, 0, 0, 0};
    
    ComplexVector result = FFT(f);
    
    std::cout << "FFT result:" << std::endl;
    for (const auto &val : result) {
        std::cout << val << std::endl;
    }
    
    ComplexVector inv_result = IFFT(result);
    
    std::cout << "IFFT result:" << std::endl;
    for (const auto &val : inv_result) {
        std::cout << val << std::endl;
    }
    
    return 0;
}