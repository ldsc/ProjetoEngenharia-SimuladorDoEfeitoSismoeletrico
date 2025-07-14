// FFTWrapper.cpp
#include "fft.h"
#include <stdexcept>

FFTWrapper::FFTWrapper()
{
    initialize();
}

FFTWrapper::FFTWrapper(int size) : size(size), in(nullptr), out(nullptr)
{
    initialize();
}

FFTWrapper::~FFTWrapper() {
    cleanup();
}

void FFTWrapper::setSize(int n) {
    size = n;
}

void FFTWrapper::initialize() {
    if (size <= 0) {
        throw std::runtime_error("Tamanho inválido para a FFT.");
    }
    std::cout << "Inicializando FFTWrapper com tamanho: " << size << std::endl;
    
    // Não chame cleanup no início. Em vez disso, verifique se os recursos estão alocados.
    // Alocação de memória para FFT
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);

    if (!in || !out) {
        std::cerr << "Falha ao alocar memória para FFTW." << std::endl;
        throw std::runtime_error("Falha ao alocar memória para FFTW.");
    }

    // Criação dos planos FFT
    std::cout << "Criando planos FFT..." << std::endl;
    plan_forward = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    if (!plan_forward) {
        std::cerr << "Falha ao criar plano FFTW_FORWARD." << std::endl;
        cleanup(); // Libere os recursos já alocados
        throw std::runtime_error("Falha ao criar plano FFTW_FORWARD.");
    }

    if (!plan_backward) {
        std::cerr << "Falha ao criar plano FFTW_BACKWARD." << std::endl;
        cleanup(); // Libere os recursos já alocados
        throw std::runtime_error("Falha ao criar plano FFTW_BACKWARD.");
    }

    std::cout << "FFTWrapper inicializado com sucesso." << std::endl;
}
void FFTWrapper::cleanup() {
    std::cout << "Iniciando cleanup." << std::endl;
    if (plan_forward) {
        fftw_destroy_plan(plan_forward);
        plan_forward = nullptr; // Evitar uso posterior
        std::cout << "Plano FFTW forward destruído." << std::endl;
    }
    if (plan_backward) {
        fftw_destroy_plan(plan_backward);
        plan_backward = nullptr; // Evitar uso posterior
        std::cout << "Plano FFTW backward destruído." << std::endl;
    }
    if (in) {
        fftw_free(in);
        in = nullptr; // Evitar uso posterior
        std::cout << "Memória de entrada FFTW liberada." << std::endl;
    }
    if (out) {
        fftw_free(out);
        out = nullptr; // Evitar uso posterior
        std::cout << "Memória de saída FFTW liberada." << std::endl;
    }
    std::cout << "Cleanup concluído." << std::endl;
}

std::vector<std::complex<double>> FFTWrapper::forward(const std::vector<std::complex<double>>& input) {
    if (input.size() != size) {
        throw std::invalid_argument("O tamanho do vetor de entrada não corresponde ao tamanho da FFT.");
    }

    // Copia os dados de entrada para o array 'in'
    for (int i = 0; i < size; ++i) {
        in[i][0] = input[i].real();
        in[i][1] = input[i].imag();
    }

    // Executa a transformada direta
    fftw_execute(plan_forward);

    // Copia os resultados para um vetor de saída
    std::vector<std::complex<double>> output(size);
    for (int i = 0; i < size; ++i) {
        output[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    // Aplica fftshift para reorganizar as frequências como no MATLAB
    output = fftshift(output);

    // Normaliza o resultado como é feito no MATLAB, dividindo pelo tamanho
    for (auto& val : output) {
        val /= static_cast<double>(size);
    }

    return output;
}


std::vector<std::complex<double>> FFTWrapper::inverse(const std::vector<std::complex<double>>& input) {
    if (!isInitialized()) {
        throw std::runtime_error("FFTWrapper não foi inicializado corretamente.");
    }

    if (input.size() != size) {
        throw std::invalid_argument("O tamanho do vetor de entrada não corresponde ao tamanho da FFT.");
    }

    // Aplica ifftshift ao vetor de entrada antes da transformada inversa
    std::vector<std::complex<double>> shiftedInput = ifftshift(input);

    // Copia os dados de entrada (shiftedInput) para o array 'in'
    for (int i = 0; i < size; ++i) {
        in[i][0] = shiftedInput[i].real();
        in[i][1] = shiftedInput[i].imag();
    }

    // Executa a transformada inversa
    fftw_execute(plan_backward);

    // Copia os resultados para um vetor de saída
    std::vector<std::complex<double>> output(size);
    for (int i = 0; i < size; ++i) {
        output[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    // Aplica ifftshift ao vetor de saída para reorganizar as frequências
    output = ifftshift(output);

    return output;
}


// Funções auxiliares fftshift e ifftshift
std::vector<std::complex<double>> FFTWrapper::fftshift(const std::vector<std::complex<double>>& input) {
    std::vector<std::complex<double>> output(size);
    int mid = size / 2;
    for (int i = 0; i < size; ++i) {
        int shiftedIndex = (i + mid) % size;
        output[i] = input[shiftedIndex];
    }
    return output;
}

std::vector<std::complex<double>> FFTWrapper::ifftshift(const std::vector<std::complex<double>>& input) {
    std::vector<std::complex<double>> output(size);
    int mid = size / 2;
    // Desloca os elementos para a esquerda, para que o centro seja movido para a primeira posição
    for (int i = 0; i < size; ++i) {
        int shiftedIndex = (i + mid) % size;
        output[i] = input[shiftedIndex];
    }
    return output;
}

bool FFTWrapper::isInitialized() const {
    // Verifica se os planos foram criados com sucesso
    return plan_forward != nullptr && plan_backward != nullptr && in != nullptr && out != nullptr;
}