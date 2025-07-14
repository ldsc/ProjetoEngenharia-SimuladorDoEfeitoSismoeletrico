// main.cpp
#include "seismicSimulator.h"
#include <iostream>

void printSystemParameters(const SystemParameters& params) {
    std::cout << "Parâmetros configurados:" << std::endl;
    std::cout << "Frequência dominante: " << params.getDominantFrequency() << " Hz" << std::endl;
    std::cout << "Profundidade: " << params.getDepth() << " m" << std::endl;
    std::cout << "Tempo total: " << params.getTotalTime() << " s" << std::endl;
    std::cout << "Número de pontos: " << params.getNumPoints() << std::endl;
    std::cout << "Comprimento de Debye: " << params.getDebyeLength() << " m" << std::endl;
    std::cout << "Permissividade elétrica: " << params.getElectricalPermittivity() << " F/m" << std::endl;
}

int main() {
    // 1. Seleção do modelo de frequência
    int modelo;
    std::cout << "Informe o MODELO que deseja simular:" << std::endl;
    std::cout << "\t1 -> Baixas frequências" << std::endl;
    std::cout << "\t2 -> Altas frequências" << std::endl;
    std::cin >> modelo;

    double dominantFrequency, totalTime;
    if (modelo == 1) {
        dominantFrequency = 20.0; // Hz para baixas frequências
        totalTime = 2.0;          // Tempo total [s]
        std::cout << "Modelo de baixas frequências selecionado." << std::endl;
    } else if (modelo == 2) {
        dominantFrequency = 200000; // Hz para altas frequências
        totalTime = 1;         // Tempo total [s]
        std::cout << "Modelo de altas frequências selecionado." << std::endl;
    } else {
        std::cerr << "Opção inválida. Saindo do programa." << std::endl;
        return 1;
    }

    // 2. Configuração dos parâmetros do sistema
    double depth = 1000.0;           // Profundidade
    int numPoints = 320;           // Número de pontos de discretização
    double debyeLength = 1e-2;       // Comprimento de Debye
    double electricalPermittivity = 8.85e-12; // Permissividade elétrica

    std::cout << "Configurando os parâmetros do sistema." << std::endl;
    SystemParameters systemParams(dominantFrequency, depth, totalTime, numPoints, debyeLength, electricalPermittivity);
    printSystemParameters(systemParams);
    std::cout << "Parâmetros do sistema configurados com sucesso." << std::endl;

    // 3. Inicialização do simulador
    std::cout << "Inicializando o simulador." << std::endl;
    try {
        SeismicSimulator simulator(systemParams);

        // 4. Cálculo das frequências e propriedades dinâmicas
        std::cout << "Calculando as frequências." << std::endl;
        simulator.calculateFrequencies();
        std::cout << "Frequências calculadas com sucesso." << std::endl;

        std::cout << "Calculando os parâmetros dinâmicos." << std::endl;
        simulator.calculateDynamicParameters();
        std::cout << "Parâmetros dinâmicos calculados com sucesso." << std::endl;

        // 5. Realização da FFT para transformar os dados para o domínio da frequência
        std::cout << "Realizando a FFT para cada fluido." << std::endl;
        simulator.performFFT(); // A geração do sinal de entrada está dentro deste método
        std::cout << "FFT realizada com sucesso." << std::endl;

        // 6. Construção da solução no domínio da frequência (vfreq e wfreq)
        simulator.constructFrequencySolution();
        std::cout << "Solução no domínio da frequência construída com sucesso." << std::endl;

        // 7. Realização da transformada inversa para obter os sinais no domínio do tempo
        simulator.performInverseFFT();
        std::cout << "Transformada inversa realizada com sucesso." << std::endl;

        // 8. Normalização dos resultados para manter a escala de amplitude idêntica ao MATLAB
        simulator.normalizeResults();
        std::cout << "Resultados normalizados com sucesso." << std::endl;

        // 9. Plotagem dos resultados
        simulator.plotResults("c:/gnuplot");
        std::cout << "Resultados plotados com sucesso." << std::endl;

        std::cout << "Simulação concluída com sucesso." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Erro durante a execução do simulador: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Pressione Enter para sair...";
    std::cin.ignore();
    std::cin.get();

    return 0;
}
