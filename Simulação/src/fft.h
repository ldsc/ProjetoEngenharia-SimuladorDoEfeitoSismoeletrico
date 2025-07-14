// FFTWrapper.h
#ifndef FFTWRAPPER_H
#define FFTWRAPPER_H

#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>

#include <fftw3/fftw3.h>

class FFTWrapper {
public:
    FFTWrapper();
    FFTWrapper(int size);
    ~FFTWrapper();

    // Configura o tamanho da FFT
    void setSize(int n);

    // Executa a transformada de Fourier direta
    std::vector<std::complex<double>> forward(const std::vector<std::complex<double>>& input);

    // Executa a transformada de Fourier inversa
    std::vector<std::complex<double>> inverse(const std::vector<std::complex<double>>& input);

    // Verifica se a FFT foi inicializada corretamente
    bool isInitialized() const;

private:
    int size; // Tamanho da transformada
    fftw_complex *in, *out;
    fftw_plan plan_forward, plan_backward;

    // Funções internas para inicializar e liberar memória
    void initialize();
    void cleanup();

    // Funções auxiliares para simular o comportamento do MATLAB
    std::vector<std::complex<double>> fftshift(const std::vector<std::complex<double>>& input);
    std::vector<std::complex<double>> ifftshift(const std::vector<std::complex<double>>& input);
};

#endif // FFTWRAPPER_H
