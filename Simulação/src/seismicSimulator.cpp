
#include "seismicSimulator.h"
#include <functional>
#include <iostream>
#include <numeric>
#include <cmath>

using namespace std::literals::complex_literals;

SeismicSimulator::SeismicSimulator(const SystemParameters& params)
    : systemParams(params), numPoints(params.getNumPoints()),
      rock(2619.0, 1590000000.0, 0.4, 262000000.0, 2.4, 0.88, 8150000000.0, 1.97e-3),
      fftProcessor(numPoints) { // Inicializa o fftProcessor com numPoints
    
    std::cout << "Construtor do SeismicSimulator chamado com sucesso." << std::endl;

    dominantFrequency = systemParams.getDominantFrequency();
    std::cout << "Frequência dominante configurada: " << dominantFrequency << " Hz." << std::endl;

    // Inicializar os fluidos conforme os parâmetros fornecidos
    fluids.push_back(FluidProperties({1040.0}, {0.0015}, 80.0, -0.044, 7.54e-4));
    fluids.push_back(FluidProperties({844.8}, {0.008989}, 80.0, -0.044, 7.54e-4));
    fluids.push_back(FluidProperties({905.1}, {0.007379}, 80.0, -0.044, 7.54e-4));
    fluids.push_back(FluidProperties({967.2}, {0.0085378}, 80.0, -0.044, 7.54e-4));
    std::cout << "Propriedades dos fluidos inicializadas com sucesso." << std::endl;

    // Inicialização adicional
    initialize();

    std::cout << "Construtor do SeismicSimulator concluído com sucesso." << std::endl;
}

void SeismicSimulator::initialize() {
    std::cout << "Inicializando vetores de tempo e frequência." << std::endl;

    //double T = systemParams.getTotalTime();
    double T = 20;
    int N = numPoints;
    double dt = T / (N - 1); // Passo de tempo
    double df = (N - 1) / (N * T); // Passo de discretização da frequência
    double fNyq = 1.0 / (2.0 * dt); // Frequência de Nyquist

    // Inicializar vetores de tempo e frequência
    timeVector.resize(N);
    frequencyVector.resize(N);
    h.resize(N, 0.0); // Inicializa h com zeros

    // Preencher o vetor de tempo
    for (int j = 0; j < N; ++j) {
        timeVector[j] = j * dt;
    }

    std::cout << "Vetores de tempo e frequência inicializados." << std::endl;

    // Inicializar vfreq e wfreq para armazenar os resultados da FFT
    vfreq.resize(fluids.size(), std::vector<std::complex<double>>(N, 0.0));
    wfreq.resize(fluids.size(), std::vector<std::complex<double>>(N, 0.0));

    std::cout << "Vetores vfreq e wfreq inicializados." << std::endl;

    // Configurar o FFTWrapper
    fftProcessor.setSize(N);
    if (!fftProcessor.isInitialized()) {
        throw std::runtime_error("Falha ao inicializar o FFTProcessor.");
    }

    // Inicializar o vetor v com zeros para todos os fluidos
    v.resize(fluids.size(), std::vector<double>(N, 1.0));
    x.resize(fluids.size(), std::vector<double>(N, 1.0));

    std::cout << "FFTProcessor inicializado com sucesso." << std::endl;
    std::cout << "Inicialização do simulador concluída com " << N << " pontos." << std::endl;
}



void SeismicSimulator::calculateDynamicParameters() {
    std::cout << "Calculando parâmetros dinâmicos." << std::endl;

    size_t numFluids = fluids.size();
    b11.resize(numFluids);
    b12.resize(numFluids);
    b21.resize(numFluids);
    b22.resize(numFluids);

    double wd = 4 * M_PI * dominantFrequency; // Frequência angular dominante
    double C = rock.getBiotCoefficient();
    double M = rock.getBiotModulus();

    for (size_t i = 0; i < numFluids; ++i) {
        const auto& fluid = fluids[i];

        // Propriedades físicas necessárias
        double phi = rock.getPorosity();
        double eta = fluid.getViscosities()[0];
        double rhof = fluid.getDensities()[0];
        double F = rock.getTortuosity() / phi;
        double kk = rock.getAbsolutePermeability();
        double rho = (1.0 - phi) * rock.getDensity() + phi * rhof;
        double lamb = rock.getBulkModulus();
        double G = rock.getShearModulus();
        double beta = 1.0 / (C * C - M * (lamb + 2 * G));

        // Cálculo de kappa conforme o MATLAB
        double kappa = kk * (1.0 - std::sqrt(eta / (F * kk * rhof) * 4.0 / wd));

        // Construir a matriz para cálculo dos autovalores
        Eigen::MatrixXd matrix(2, 2);
        matrix(0, 0) = beta * (2.0 * C * rhof - kappa * (lamb + 2 * G));
        matrix(1, 1) = beta * (2.0 * C * rhof + kappa * (lamb + 2 * G));
        matrix(0, 1) = kappa * M;
        matrix(1, 0) = kappa * M;

        // Calcular autovalores e autovetores
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(matrix);
        if (eigenSolver.info() != Eigen::Success) {
            std::cerr << "Erro ao calcular autovalores para o fluido " << i + 1 << "." << std::endl;
            continue;
        }

        // Obter os autovalores e autovetores
        Eigen::VectorXd eigenvalues = eigenSolver.eigenvalues().real();
        Eigen::MatrixXd eigenvectors = eigenSolver.eigenvectors().real();

        // Armazenar os valores nos vetores bij
        b11[i] = eigenvalues[0];
        b12[i] = eigenvalues[1];
        b21[i] = eigenvectors(0, 0);
        b22[i] = eigenvectors(1, 0);
    }

    std::cout << "Parâmetros dinâmicos calculados com sucesso." << std::endl;
}

void SeismicSimulator::calculateFrequencies() {
    // Configuração das frequências conforme o MATLAB
    std::cout << "Calculando as frequências." << std::endl;

    // Número de pontos
    size_t N = numPoints;

    // Tempo total e intervalo de tempo
    double T = systemParams.getTotalTime();
    double dt = T / N;

    // Frequência de Nyquist
    double fNyq = 1.0 / (2.0 * dt);

    // Passo de discretização da frequência
    double df = 1.0 / T;

    // Definindo o vetor de frequências de acordo com o MATLAB
    frequencyVector.resize(N);
    // Calcular o vetor de frequência conforme o MATLAB
    for (int j = 0.0001; j < N; ++j) {
        frequencyVector[j] = j * df;
    }


    std::cout << "Frequências calculadas com sucesso." << std::endl;
}

void SeismicSimulator::performFFT() {
    // Atualiza o vetor de frequências
    calculateFrequencies();

    // Realiza a FFT para cada fluido
    std::cout << "Realizando a FFT para cada fluido." << std::endl;
    for (size_t i = 0; i < fluids.size(); ++i) {
        std::cout << "Processando FFT para o fluido " << i + 1 << "." << std::endl;

        // Geração do sinal de entrada (Ricker Wavelet)
        double dt = systemParams.getTotalTime() / (numPoints - 1);
        double dominantFrequency = systemParams.getDominantFrequency();
        std::vector<double> inputSignal = generateRickerWavelet(numPoints, dt, dominantFrequency);

        // Garantir que o tamanho do sinal de entrada seja exatamente igual a numPoints
        if (inputSignal.size() != numPoints) {
            std::cerr << "Erro: tamanho do sinal de entrada não corresponde a numPoints." << std::endl;
            inputSignal.resize(numPoints, 0.0); // Redimensiona preenchendo com zeros se necessário
        }

        // Verificar se o sinal gerado não está zerado
        double maxInputSignal = *std::max_element(inputSignal.begin(), inputSignal.end());
        std::cout << "Máxima amplitude do sinal de entrada: " << maxInputSignal << std::endl;
        if (maxInputSignal < 1e-6) {
            std::cerr << "Aviso: o sinal de entrada para o fluido " << i + 1 << " tem amplitude muito baixa." << std::endl;
            continue; // Pular este fluido, pois o sinal está muito baixo
        }

        // Preencher o vetor complexInput com os dados do sinal de entrada
        std::vector<std::complex<double>> complexInput(numPoints);
        for (int j = 0; j < numPoints; ++j) {
            complexInput[j] = std::complex<double>(inputSignal[j], 0.0);
        }

        try {
            // Verificar se o FFTWrapper foi inicializado corretamente
            if (!fftProcessor.isInitialized()) {
                throw std::runtime_error("FFTWrapper não foi inicializado corretamente.");
            }

            // Realizar a FFT para todos os pontos e armazenar o resultado em fftResult
            std::cout << "Executando a FFT para o fluido " << i + 1 << "." << std::endl;
            std::vector<std::complex<double>> fftResult = fftProcessor.forward(complexInput);
            std::cout << "FFT para o fluido " << i + 1 << " concluída." << std::endl;

            // Normalizar o resultado da FFT dividindo pelo número de pontos
            for (auto& val : fftResult) {
                val /= static_cast<double>(numPoints);
            }

            // Atualizar o vetor vfreq com os resultados para cada ponto
            if (vfreq[i].size() != numPoints) {
                vfreq[i].resize(numPoints);
            }
            for (int j = 0; j < numPoints; ++j) {
                vfreq[i][j] = fftResult[j]; // Armazena o resultado da FFT no vetor de frequência
            }

            // Verificar se o resultado da FFT não está zerado
            double maxFFT = std::abs(*std::max_element(vfreq[i].begin(), vfreq[i].end(),
                [](const std::complex<double>& a, const std::complex<double>& b) {
                    return std::abs(a) < std::abs(b);
                }));
            std::cout << "Máxima amplitude da FFT para o fluido " << i + 1 << ": " << maxFFT << std::endl;
            if (maxFFT < 1e-6) {
                std::cerr << "Aviso: o resultado da FFT para o fluido " << i + 1 << " tem amplitude muito baixa." << std::endl;
            }

            std::cout << "FFT para o fluido " << i + 1 << " concluída com sucesso." << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "Erro ao realizar a FFT para o fluido " << i + 1 << ": " << e.what() << std::endl;
        }

        // Mensagem de depuração para verificar se o loop continua
        std::cout << "Continuando para o próximo fluido." << std::endl;
    }
    std::cout << "FFT para todos os fluidos concluída." << std::endl;
}


void SeismicSimulator::performInverseFFT() {
    std::cout << "Realizando a IFFT para cada fluido." << std::endl;

    for (size_t i = 0; i < fluids.size(); ++i) {
        if (vfreq[i].empty() || wfreq[i].empty()) {
            std::cerr << "Aviso: vfreq ou wfreq vazio para o fluido " << i + 1 << ", pulando IFFT." << std::endl;
            continue;
        }

        // Realiza a IFFT usando a classe FFTWrapper para vfreq e wfreq
        std::vector<std::complex<double>> inverseVResult = fftProcessor.inverse(vfreq[i]);
        std::vector<std::complex<double>> inverseWResult = fftProcessor.inverse(wfreq[i]);

        // Multiplica por (4 * pi) / dt para normalizar
        double dt = systemParams.getTotalTime() / (numPoints - 1);
        double normalizationFactor = (4.0 * M_PI) / dt;
        for (size_t j = 0; j < inverseVResult.size(); ++j) {
            inverseVResult[j] *= normalizationFactor / static_cast<double>(numPoints);
            inverseWResult[j] *= normalizationFactor / static_cast<double>(numPoints);
        }

        // Centraliza os dados em torno de zero (subtraindo a média)
        double meanV = std::accumulate(inverseVResult.begin(), inverseVResult.end(), 0.0,
                                       [](double sum, const std::complex<double>& val) {
                                           return sum + val.real();
                                       }) / numPoints;
        double meanW = std::accumulate(inverseWResult.begin(), inverseWResult.end(), 0.0,
                                       [](double sum, const std::complex<double>& val) {
                                           return sum + val.real();
                                       }) / numPoints;
        for (auto& val : inverseVResult) {
            val -= meanV;
        }
        for (auto& val : inverseWResult) {
            val -= meanW;
        }

        // Atualiza os vetores v e x com os valores reais da transformada inversa e inverte a ordem
        for (int j = 0; j < numPoints; ++j) {
            v[i][j] = inverseVResult[numPoints - 1 - j].real();
            x[i][j] = inverseWResult[numPoints - 1 - j].real();
        }

        std::cout << "IFFT para o fluido " << i + 1 << " concluída com sucesso." << std::endl;
    }
    std::cout << "IFFT para todos os fluidos concluída." << std::endl;
}



void SeismicSimulator::normalizeResults() {
    // Calcular o valor máximo de vmax2 e xmax2
    double vmax2 = 0.0;
    double xmax2 = 0.0;

    std::cout << "fluidos: " << fluids.size() <<" e " << vmax2 << std::endl;

    // Encontrar o máximo valor de v e x para normalização
    for (size_t i = 0; i < fluids.size(); ++i) {
        double currentVMax = *std::max_element(v[i].begin(), v[i].end(), 
            [](double a, double b) { return std::abs(a) < std::abs(b); });
        double currentXMax = *std::max_element(x[i].begin(), x[i].end(), 
            [](double a, double b) { return std::abs(a) < std::abs(b); });

        vmax2 = std::max(vmax2, currentVMax);
        xmax2 = std::max(xmax2, currentXMax);
        std::cout << "maxX e maxY: " << xmax2 <<" e " << vmax2 << std::endl;
    }

    // Normaliza os resultados para cada fluido
    for (size_t i = 0; i < fluids.size(); ++i) {
        if (vmax2 > 0) { // Evitar divisão por zero
            for (size_t j = 0; j < v[i].size(); ++j) {
                v[i][j] /= vmax2; // Normalizar v
            }
        }

        if (xmax2 > 0) { // Evitar divisão por zero
            for (size_t j = 0; j < x[i].size(); ++j) {
                x[i][j] /= xmax2; // Normalizar x
            }
        }
    }
}

void SeismicSimulator::constructFrequencySolution() {
    std::cout << "Construindo a solução no domínio da frequência." << std::endl;

    double ezero = 8.85e-12;
    double F = rock.getTortuosity() / rock.getPorosity();
    double wd = 2 * M_PI * dominantFrequency;
    double z = systemParams.getDepth();
    double ni = sqrt(8 * rock.getAbsolutePermeability() * F);
    double sigmab = fluids[0].getConductivity();
    double m = 8.0;

    // Preparar variáveis para cálculo
    std::vector<double> rho(fluids.size());
    std::vector<double> omegat(fluids.size());
    std::vector<double> Lzero(fluids.size());

    for (size_t k = 0; k < fluids.size(); ++k) {
        rho[k] = (1.0 - rock.getPorosity()) * rock.getDensity() + rock.getPorosity() * fluids[k].getDensities()[0];
        omegat[k] = fluids[k].getViscosities()[0] / (F * rock.getAbsolutePermeability() * fluids[k].getDensities()[0]);
        Lzero[k] = -(1.0 / F) * (ezero * fluids[k].getDielectricConstant() * fluids[k].getDielectricConstant() / fluids[k].getViscosities()[0]);
    }

    // Calcular o vetor de frequências conforme o MATLAB
    double T = systemParams.getTotalTime();
    double df = 1.0 / T;
    std::vector<double> f(numPoints);
    for (int j = 0; j < numPoints; ++j) {
        f[j] = (0.0001 + j * 0.5) * df;
    }

    // Vetor de frequências angulares
    std::vector<std::complex<double>> w(numPoints);
    std::complex<double> delta = std::complex<double>(0, M_PI / T);
    for (int j = 0; j < numPoints; ++j) {
        w[j] = 2.0 * M_PI * f[j] + delta;
    }

    for (size_t i = 0; i < fluids.size(); ++i) {
        std::cout << "Construindo solução para o fluido " << i + 1 << "..." << std::endl;

        for (int j = 0; j < numPoints; ++j) {
            std::complex<double> wj = w[j];

            // Calcular o valor de h[j] com base nas condições de wj
            if (wj == wd) {
                h[j] = -1.0i * M_PI / wd
                    - (21.0 / 16.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 4.0 * wd * wd))
                    + (21.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 16.0 * wd * wd))
                    - (1.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 64.0 * wd * wd));
            } else if (wj == 2.0 * wd) {
                h[j] = wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - wd * wd))
                    + 21.0i * M_PI / (32.0 * wd)
                    + (21.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 16.0 * wd * wd))
                    - (1.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 64.0 * wd * wd));
            } else if (wj == 4.0 * wd) {
                h[j] = wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - wd * wd))
                    - (21.0 / 16.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 4.0 * wd * wd))
                    - 21.0i * M_PI / (256.0 * wd)
                    - (1.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 64.0 * wd * wd));
            } else if (wj == 8.0 * wd) {
                h[j] = wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - wd * wd))
                    - (21.0 / 16.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 4.0 * wd * wd))
                    + (21.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 16.0 * wd * wd))
                    + 1.0i * M_PI / (512.0 * wd);
            } else {
                h[j] = wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - wd * wd))
                    - (21.0 / 16.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 4.0 * wd * wd))
                    + (21.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 16.0 * wd * wd))
                    - (1.0 / 64.0) * wd * ((std::exp(-2.0i * M_PI * wj / wd) - 1.0) / (wj * wj - 64.0 * wd * wd));
            }

            // Cálculo de kappa (permeabilidade dinâmica)
            std::complex<double> kappa = rock.getAbsolutePermeability() * std::pow(1.0 - (1.0i * wj / omegat[i]) * 4.0 / m, 0.5) - (1.0i * wj / omegat[i]);

            // Cálculo de L (coeficiente de acoplamento)
            std::complex<double> L = Lzero[i] * std::pow(1.0 - (1.0i * (wj / omegat[i]) * (m / 4.0)) * std::pow(1.0 - 2.0 * z / ni, 2.0) * std::pow(1.0 - 1.0i * std::sqrt(3.0 / 2.0) * z * std::sqrt(wj * fluids[i].getDensities()[0] / fluids[i].getViscosities()[0]), 2.0), -0.5);

            // Cálculo de Q (termo simplificador)
            std::complex<double> Q = -fluids[i].getViscosities()[0] / (1.0i * wj * kappa * (1.0 - (std::pow(L, 2.0) * fluids[i].getViscosities()[0]) / (kappa * sigmab)));

            // Parâmetros para o cálculo dos autovalores
            double C = rock.getBiotCoefficient();
            double M = rock.getBiotModulus();
            double lamb = rock.getBulkModulus();
            double G = rock.getShearModulus();
            double beta = 1.0 / (C * C - M * (lamb + 2.0 * G));

             // Cálculo dos autovalores q1 e q2
            std::complex<double> q1quad = (beta / 2.0) * ((2.0 * C * fluids[i].getDensities()[0] - M * rho[i] - Q * (lamb + 2.0 * G)) +
                                           std::sqrt(std::pow(Q * (lamb + 2.0 * G) - M * rho[i], 2.0) -
                                                     4.0 * ((M * fluids[i].getDensities()[0] - C * Q) * (C * rho[i] - (lamb + 2.0 * G) * fluids[i].getDensities()[0]))));
            std::complex<double> q2quad = (beta / 2.0) * ((2.0 * C * fluids[i].getDensities()[0] - M * rho[i] - Q * (lamb + 2.0 * G)) -
                                           std::sqrt(std::pow(Q * (lamb + 2.0 * G) - M * rho[i], 2.0) -
                                                     4.0 * ((M * fluids[i].getDensities()[0] - C * Q) * (C * rho[i] - (lamb + 2.0 * G) * fluids[i].getDensities()[0]))));
            std::complex<double> q1 = std::sqrt(q1quad);
            std::complex<double> q2 = std::sqrt(q2quad);

            // Cálculo dos autovalores e autovetores
            std::complex<double> a1 = sqrt((M * fluids[i].getDensities()[0] * M * fluids[i].getDensities()[0] - 2.0 * M * fluids[i].getDensities()[0] * C * Q + C * C * Q * Q) /
                (M * M * (fluids[i].getDensities()[0] * fluids[i].getDensities()[0] + rho[i] * rho[i]) + C * C * (fluids[i].getDensities()[0] * fluids[i].getDensities()[0] + Q * Q) - 
                2.0 * M * fluids[i].getDensities()[0] * C * (rho[i] + Q) - 
                2.0 * (fluids[i].getDensities()[0] * C - M * rho[i]) * ((q1 * q1) / beta) + 
                pow((q1 * q1) / beta, 2.0)));

            std::complex<double> a2 = sqrt((M * fluids[i].getDensities()[0] * M * fluids[i].getDensities()[0] - 2.0 * M * fluids[i].getDensities()[0] * C * Q + C * C * Q * Q) /
                (M * M * (fluids[i].getDensities()[0] * fluids[i].getDensities()[0] + rho[i] * rho[i]) + C * C * (fluids[i].getDensities()[0] * fluids[i].getDensities()[0] + Q * Q) - 
                2.0 * M * fluids[i].getDensities()[0] * C * (rho[i] + Q) - 
                2.0 * (fluids[i].getDensities()[0] * C - M * rho[i]) * ((q2 * q2) / beta) + 
                pow((q2 * q2) / beta, 2.0)));

            // Cálculo dos parâmetros xi1 e xi2
            std::complex<double> xi1 = (C * fluids[i].getDensities()[0] - M * rho[i] - (q1 * q1) / beta) / (M * fluids[i].getDensities()[0] - C * Q);
            std::complex<double> xi2 = (C * fluids[i].getDensities()[0] - M * rho[i] - (q2 * q2) / beta) / (M * fluids[i].getDensities()[0] - C * Q);

            // Cálculo dos parâmetros b11, b12, b21, b22
            b11[i] = -a1;
            b21[i] = a1 * xi1;
            b12[i] = -a2;
            b22[i] = a2 * xi2;

            // Calcular os autovetores y11, y12, y21, y22
            std::complex<double> y11 = (a1 / q1) * (-rho[i] - fluids[i].getDensities()[0] * xi1);
            std::complex<double> y21 = (a1 / q1) * (fluids[i].getDensities()[0] + Q * xi1);
            std::complex<double> y12 = (a2 / q2) * (-rho[i] - fluids[i].getDensities()[0] * xi2);
            std::complex<double> y22 = (a2 / q2) * (fluids[i].getDensities()[0] + Q * xi2);

            // Cálculo de vfreq e wfreq utilizando os autovalores e autovetores
            vfreq[i][j] = h[j] * (y11 * (b11[i] - b21[i]) * exp(1i * wj * q1 * z) + y12 * (b12[i] - b22[i]) * exp(1i * wj * q2 * z));
            wfreq[i][j] = -h[j] * (y21 * (b11[i] - b21[i]) * exp(1i * wj * q1 * z) + y22 * (b12[i] - b22[i]) * exp(1i * wj * q2 * z));
        }

        std::cout << "Solução para o fluido " << i + 1 << " construída com sucesso." << std::endl;
    }
}


void SeismicSimulator::plotResults(const std::string& gnuplotPath) {
    GnuplotPlotter plotter(gnuplotPath);
    plotter.plotResults(v, x, timeVector);
    std::cout << "Resultados plotados com sucesso." << std::endl;
}

std::vector<double> SeismicSimulator::generateRickerWavelet(int numPoints, double dt, double dominantFrequency) {
    // Parâmetros para a onda Ricker
    double piSquared = M_PI * M_PI;
    double fSquared = dominantFrequency * dominantFrequency;
    double t0 = (numPoints - 1) * dt / 2.0; // Centralizar a onda Ricker no meio do intervalo de tempo

    std::vector<double> wavelet(numPoints);

    for (int i = 0; i < numPoints; ++i) {
        double t = i * dt - t0; // Centraliza o tempo em torno de zero
        double term = piSquared * fSquared * t * t;
        wavelet[i] = (1.0 - 2.0 * term) * std::exp(-term);
    }

    return wavelet;
}

