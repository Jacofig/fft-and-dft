#include <iostream>
#include <complex>
#include <cmath>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>

typedef std::complex<double> Complex;

template<typename T>
std::complex<T>* dft(const T* arr, int size) {
    std::complex<T>* c = new std::complex<T>[size];
    for (int i = 0; i < size; i++) {
        c[i] = 0;
        for (int n = 0; n < size; n++) {
            T power = -2.0 * M_PI * i * n / size;
            c[i] += arr[n] * std::exp(std::complex<T>(0, power));
        }
    }
    return c;
}

template<typename T>
std::complex<T>* fft(const T* arr, int size) {
    if (size == 1) {
        return new std::complex<T>[1] { arr[0] };
    }

    int sizeHalf = size / 2;

    T* evenHalf = new T[sizeHalf];
    T* oddHalf = new T[sizeHalf];

    for (int i = 0; i < sizeHalf; i++) {
        evenHalf[i] = arr[2 * i];
        oddHalf[i] = arr[2 * i + 1];
    }

    std::complex<T>* even = fft(evenHalf, sizeHalf);
    std::complex<T>* odd = fft(oddHalf, sizeHalf);

    std::complex<T>* c = new std::complex<T>[size];
    for (int i = 0; i < sizeHalf; i++) {
        std::complex<T> t = std::polar<T>(1.0, -2.0 * M_PI * i / size) * odd[i];
        c[i] = even[i] + t;
        c[i + sizeHalf] = even[i] - t;
    } 

    delete[] evenHalf;
    delete[] oddHalf;
    delete[] even;
    delete[] odd;

    return c;
}

template<typename T>
T err(const std::complex<T>* cDFT, const std::complex<T>* cFFT, int size) {
    T sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += std::abs(cDFT[i] - cFFT[i]);
    }
    return sum / size;
}

int main() {
    const int MAX_ORDER = 13;
    const bool PRINT_COEFS = true;

    for (int o = 4; o <= MAX_ORDER; o++) {
        const int N = 1 << o;
        std::cout << "N: " << N << "    order: " << o << std::endl;

        double* f = new double[N];
        for (int n = 0; n < N; n++) {
            f[n] = n / (double)N;
        }

        clock_t t1 = clock();
        std::complex<double>* cDFT = dft(f, N);
        clock_t t2 = clock();
        double dft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
        std::cout << "DFT time [ms]: " << dft_time << std::endl;

        t1 = clock();
        std::complex<double>* cFFT = fft(f, N);
        t2 = clock();
        double fft_time = (t2 - t1) / (double)CLOCKS_PER_SEC * 1000.0;
        std::cout << "FFT time [ms]: " << fft_time << std::endl;

        std::cout << "mean error: " << err(cDFT, cFFT, N) << std::endl;

        if (PRINT_COEFS) {
            for (int k = 0; k < 16; k++) {
                std::cout << "DFT: " << cDFT[k] << ", FFT: " << cFFT[k] << std::endl;
            }
        }
        std::cout << "--------" << std::endl;

        delete[] f;
        delete[] cDFT;
        delete[] cFFT;
    }
    return 0;
}


