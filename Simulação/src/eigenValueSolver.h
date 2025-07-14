// EigenvalueSolver.h
#ifndef EIGENVALUE_SOLVER_H
#define EIGENVALUE_SOLVER_H

#include <Eigen/Dense>
#include <vector>

class EigenvalueSolver {
public:
    // Construtor
    EigenvalueSolver();

    // Define a matriz para calcular os autovalores e autovetores
    void setMatrix(const Eigen::MatrixXd& matrix);

    // Calcula os autovalores e autovetores
    bool compute();

    // Obtém os autovalores
    Eigen::VectorXd getEigenvalues() const;

    // Obtém os autovetores
    Eigen::MatrixXd getEigenvectors() const;

private:
    Eigen::MatrixXd m_matrix; // Matriz de entrada
    Eigen::VectorXd m_eigenvalues; // Autovalores calculados
    Eigen::MatrixXd m_eigenvectors; // Autovetores calculados
    bool m_isComputed; // Indica se os autovalores e autovetores foram calculados
};

#endif // EIGENVALUE_SOLVER_H
