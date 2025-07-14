// GnuplotPlotter.cpp
#include "gnuPlotter.h"

// Para portabilidade entre sistemas (Windows vs. Linux/macOS)
#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#endif

GnuplotPlotter::GnuplotPlotter(const std::string& gnuplotPath)
    : gnuplotPath(gnuplotPath) {
    // A variável tempDataFileX membro não é mais usada,
    // pois os nomes dos arquivos temporários são definidos localmente
    // dentro de plotResults para 'v' e 'x'.
    // Portanto, não há necessidade de inicializá-la aqui.
}

void GnuplotPlotter::plotResults(const std::vector<std::vector<double>>& resultsV, 
                                 const std::vector<std::vector<double>>& resultsX, 
                                 const std::vector<double>& timeVector) {
    // Escreve os dados em arquivos temporários separados para 'v' e 'x'
    std::string tempDataFileV = "temp_data_v.txt";
    std::string tempDataFileX = "temp_data_x.txt";
    writeDataToFile(resultsV, timeVector, tempDataFileV);
    writeDataToFile(resultsX, timeVector, tempDataFileX);

    // Configura o comando para o Gnuplot
    // Usando popen/pclose para compatibilidade
    std::string command = "\"" + gnuplotPath + "/bin/gnuplot.exe\" -persist";
    FILE* gnuplotPipe = popen(command.c_str(), "w");
    if (!gnuplotPipe) {
        std::cerr << "Erro ao abrir o pipe para o Gnuplot." << std::endl;
        return;
    }

    // Configuraçoes comuns para Gnuplot
    const char* gnuplotSettings = R"(
        set key outside
        unset ytics
        set grid
        set encoding utf8
        set xlabel 'Tempo (s)'
        set ylabel 'Amplitude Normalizada'
        set tics font ',16'
        set title font ',18'
        set term wxt size 800,600
    )";
    fprintf(gnuplotPipe, "%s\n", gnuplotSettings);

    // Plot para 'x' (Velocidade do fluido) na primeira janela
    fprintf(gnuplotPipe, "set term wxt 0\n"); // Abre a primeira janela
    fprintf(gnuplotPipe, "set title 'Velocidade da Fase Fluida'\n");
    fprintf(gnuplotPipe, "plot '%s' using 1:2 with lines title 'Agua' lc rgb 'blue', \\\n", tempDataFileX.c_str());
    fprintf(gnuplotPipe, "     '%s' using 1:($3+2) with lines title 'Oleo Leve' lc rgb 'green', \\\n", tempDataFileX.c_str());
    fprintf(gnuplotPipe, "     '%s' using 1:($4+4) with lines title 'Oleo Medio' lc rgb 'red', \\\n", tempDataFileX.c_str());
    fprintf(gnuplotPipe, "     '%s' using 1:($5+6) with lines title 'Oleo Pesado' lc rgb 'black'\n", tempDataFileX.c_str());
    fflush(gnuplotPipe); // Garante que os comandos para a primeira janela são enviados

 

    // Fecha o pipe do Gnuplot
    fflush(gnuplotPipe); // Garante que os comandos para a segunda janela são enviados
    pclose(gnuplotPipe); // Usa pclose para compatibilidade

    // Remover os arquivos temporários
    //std::remove(tempDataFileV.c_str());
    //std::remove(tempDataFileX.c_str());
}

void GnuplotPlotter::writeDataToFile(const std::vector<std::vector<double>>& results, 
                                     const std::vector<double>& timeVector, 
                                     const std::string& fileName) {
    std::ofstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Erro: Não foi possível abrir o arquivo " << fileName << " para escrita." << std::endl;
        return;
    }

    // Assumimos que results[0] corresponde à coluna de dados para timeVector
    // e as colunas subsequentes são os dados para cada fluido.
    // O formato esperado pelo gnuplot é: Tempo Fluido1 Fluido2 ...
    // timeVector.size() deve ser igual a results[i].size() para todos i.
    if (timeVector.size() != results[0].size()) {
        std::cerr << "Erro: Tamanho do vetor de tempo não corresponde ao tamanho dos resultados." << std::endl;
        return;
    }

    for (size_t i = 0; i < timeVector.size(); ++i) {
        file << timeVector[i]; // Coluna de tempo
        for (size_t j = 0; j < results.size(); ++j) {
            file << " " << results[j][i]; // Colunas de dados para cada fluido
        }
        file << std::endl;
    }

    file.close();
}