// EigenvalueSolver.cpp
#include "eigenValueSolver.h"
#include <Eigen/Eigenvalues>
#include <iostream>


// Construtor
EigenvalueSolver::EigenvalueSolver() : m_isComputed(false) {}

// Define a matriz
void EigenvalueSolver::setMatrix(const Eigen::MatrixXd& matrix) {
    m_matrix = matrix;
    m_isComputed = false; // Marca como não computado
}

// Calcula os autovalores e autovetores
bool EigenvalueSolver::compute() {
    if (m_matrix.size() == 0) {
        std::cerr << "A matriz não foi definida." << std::endl;
        return false;
    }

    // Calcula os autovalores e autovetores
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_matrix);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Falha ao calcular os autovalores." << std::endl;
        return false;
    }

    // Armazena os resultados
    m_eigenvalues = solver.eigenvalues();
    m_eigenvectors = solver.eigenvectors();
    m_isComputed = true;

    return true;
}

// Retorna os autovalores
Eigen::VectorXd EigenvalueSolver::getEigenvalues() const {
    if (!m_isComputed) {
        std::cerr << "Os autovalores não foram calculados." << std::endl;
    }
    return m_eigenvalues;
}

// Retorna os autovetores
Eigen::MatrixXd EigenvalueSolver::getEigenvectors() const {
    if (!m_isComputed) {
        std::cerr << "Os autovetores não foram calculados." << std::endl;
    }
    return m_eigenvectors;
}
