#ifndef FLUID_PROPERTIES_H
#define FLUID_PROPERTIES_H

#include <vector>

class FluidProperties {
public:
    FluidProperties(const std::vector<double>& densities, const std::vector<double>& viscosities,
                    double dielectricConstant, double zetaPotential, double conductivity)
        : densities(densities), viscosities(viscosities), dielectricConstant(dielectricConstant),
          zetaPotential(zetaPotential), conductivity(conductivity) {}

    std::vector<double> getDensities() const { return densities; }
    std::vector<double> getViscosities() const { return viscosities; }
    double getDielectricConstant() const { return dielectricConstant; }
    double getZetaPotential() const { return zetaPotential; }
    double getConductivity() const { return conductivity; }

private:
    std::vector<double> densities;
    std::vector<double> viscosities;
    double dielectricConstant;
    double zetaPotential;
    double conductivity;
};

#endif // FLUID_PROPERTIES_H