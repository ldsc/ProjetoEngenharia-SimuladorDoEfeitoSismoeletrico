#ifndef SEISMIC_SIMULATOR_H
#define SEISMIC_SIMULATOR_H

#include "systemParameters.h"
#include "rockProperties.h"
#include "fluidProperties.h"
#include "eigenValueSolver.h"
#include "gnuPlotter.h"
#include "fft.h"
#include <vector>
#include <complex>
#include <iostream>

class SeismicSimulator {
public:
    SeismicSimulator(const SystemParameters& params);

    void initialize();
    void calculateFrequencies();
    void calculateDynamicParameters();
    void performFFT();
    void constructFrequencySolution();
    void performInverseFFT();
    void normalizeResults();
    void plotResults(const std::string& gnuplotPath);

private:
    SystemParameters systemParams;
    int numPoints;
    RockProperties rock; 
    EigenvalueSolver eigenSolver;
    FFTWrapper fftProcessor;

    std::vector<FluidProperties> fluids;
    std::vector<double> timeVector;
    std::vector<double> frequencyVector;
    std::vector<std::vector<double>> v, x; // Resultados no domínio do tempo
    std::vector<std::vector<std::complex<double>>> vfreq, wfreq; // Resultados no domínio da frequência
    std::vector<std::complex<double>> h; // Vetor h
    std::vector<std::complex<double>> b11, b12, b21, b22; // Parâmetros calculados
    double dominantFrequency;

private:
    std::vector<double> generateRickerWavelet(int numPoints, double dt, double dominantFrequency);
};

#endif // SEISMIC_SIMULATOR_H