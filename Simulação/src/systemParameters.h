// SystemParameters.h
#ifndef SYSTEM_PARAMETERS_H
#define SYSTEM_PARAMETERS_H

class SystemParameters {
public:

    SystemParameters(double dominantFrequency, double depth, double totalTime, int numPoints,
                     double debyeLength, double electricalPermittivity)
        : dominantFrequency(dominantFrequency), depth(depth), totalTime(totalTime), numPoints(numPoints),
          debyeLength(debyeLength), electricalPermittivity(electricalPermittivity) {}

    // Getters
    double getDominantFrequency() const { return dominantFrequency; }
    double getDepth() const { return depth; }
    double getTotalTime() const { return totalTime; }
    int getNumPoints() const { return numPoints; }
    double getDebyeLength() const { return debyeLength; }
    double getElectricalPermittivity() const { return electricalPermittivity; }

    // Setters
    void setDominantFrequency(double value) { dominantFrequency = value; }
    void setDepth(double value) { depth = value; }
    void setTotalTime(double value) { totalTime = value; }
    void setNumPoints(int value) { numPoints = value; }
    void setDebyeLength(double value) { debyeLength = value; }
    void setElectricalPermittivity(double value) { electricalPermittivity = value; }

private:
    double dominantFrequency; // Frequência dominante [Hz]
    double depth; // Profundidade da fonte [m]
    double totalTime; // Tempo total da simulação [s]
    int numPoints; // Número de pontos de discretização
    double debyeLength; // Comprimento de Debye
    double electricalPermittivity; // Permissividade elétrica do meio
};

#endif // SYSTEM_PARAMETERS_H
