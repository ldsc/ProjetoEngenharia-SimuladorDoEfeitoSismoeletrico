// main.cpp
#include "seismicSimulator.h"
#include <iostream>

void printSystemParameters(const SystemParameters& params) {
    std::cout << "Parametros configurados:" << std::endl;
    std::cout << "Frequencia dominante: " << params.getDominantFrequency() << " Hz" << std::endl;
    std::cout << "Profundidade: " << params.getDepth() << " m" << std::endl;
    std::cout << "Tempo total: " << params.getTotalTime() << " s" << std::endl;
    std::cout << "Numero de pontos: " << params.getNumPoints() << std::endl;
    std::cout << "Comprimento de Debye: " << params.getDebyeLength() << " m" << std::endl;
    std::cout << "Permissividade elétrica: " << params.getElectricalPermittivity() << " F/m" << std::endl;
}

int main() {
    // 1. Selecao do modelo de frequencia
    int modelo;
    std::cout << "Informe o MODELO que deseja simular:" << std::endl;
    std::cout << "\t1 -> Baixas frequencias" << std::endl;
    std::cout << "\t2 -> Altas frequencias" << std::endl;
    std::cin >> modelo;

    double dominantFrequency, totalTime;
    if (modelo == 1) {
        dominantFrequency = 20.0; // Hz para baixas frequencias
        totalTime = 2.0;          // Tempo total [s]
        std::cout << "Modelo de baixas frequencias selecionado." << std::endl;
    } else if (modelo == 2) {
        dominantFrequency = 200000; // Hz para altas frequencias
        totalTime = 1;         // Tempo total [s]
        std::cout << "Modelo de altas frequencias selecionado." << std::endl;
    } else {
        std::cerr << "Opcao invalida. Saindo do programa." << std::endl;
        return 1;
    }

    // 2. Configuracao dos parametros do sistema
    double depth = 1000.0;           // Profundidade
    int numPoints = 320;           // Numero de pontos de discretizacao
    double debyeLength = 1e-2;       // Comprimento de Debye
    double electricalPermittivity = 8.85e-12; // Permissividade elétrica

    std::cout << "Configurando os parametros do sistema." << std::endl;
    SystemParameters systemParams(dominantFrequency, depth, totalTime, numPoints, debyeLength, electricalPermittivity);
    printSystemParameters(systemParams);
    std::cout << "Parametros do sistema configurados com sucesso." << std::endl;

    // 3. Inicializacao do simulador
    std::cout << "Inicializando o simulador." << std::endl;
    try {
        SeismicSimulator simulator(systemParams);

        // 4. Calculo das frequencias e propriedades dinamicas
        std::cout << "Calculando as frequencias." << std::endl;
        simulator.calculateFrequencies();
        std::cout << "Frequencias calculadas com sucesso." << std::endl;

        std::cout << "Calculando os parametros dinamicos." << std::endl;
        simulator.calculateDynamicParameters();
        std::cout << "Parametros dinamicos calculados com sucesso." << std::endl;

        // 5. Realizacao da FFT para transformar os dados para o dominio da frequencia
        std::cout << "Realizando a FFT para cada fluido." << std::endl;
        simulator.performFFT(); // A geracao do sinal de entrada esta dentro deste método
        std::cout << "FFT realizada com sucesso." << std::endl;

        // 6. Construcao da solucao no dominio da frequencia (vfreq e wfreq)
        simulator.constructFrequencySolution();
        std::cout << "Solucao no dominio da frequencia construida com sucesso." << std::endl;

        // 7. Realizacao da transformada inversa para obter os sinais no dominio do tempo
        simulator.performInverseFFT();
        std::cout << "Transformada inversa realizada com sucesso." << std::endl;

        // 8. Normalizacao dos resultados para manter a escala de amplitude identica ao MATLAB
        simulator.normalizeResults();
        std::cout << "Resultados normalizados com sucesso." << std::endl;

        // 9. Plotagem dos resultados
        simulator.plotResults("c:/gnuplot");
        std::cout << "Resultados plotados com sucesso." << std::endl;

        std::cout << "Simulacao concluida com sucesso." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Erro durante a execucao do simulador: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Pressione Enter para sair...";
    std::cin.ignore();
    std::cin.get();

    return 0;
}
