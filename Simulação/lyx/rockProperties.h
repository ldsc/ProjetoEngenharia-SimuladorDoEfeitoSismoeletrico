#ifndef ROCK_PROPERTIES_H
#define ROCK_PROPERTIES_H

class RockProperties {
public:
    RockProperties(double density, double shearModulus, double porosity, double bulkModulus,
                   double tortuosity, double biotCoefficient, double biotModulus, double permeability)
        : density(density), shearModulus(shearModulus), porosity(porosity), bulkModulus(bulkModulus),
          tortuosity(tortuosity), biotCoefficient(biotCoefficient), biotModulus(biotModulus),
          permeability(permeability) {}

    double getDensity() const { return density; }
    double getShearModulus() const { return shearModulus; }
    double getPorosity() const { return porosity; }
    double getBulkModulus() const { return bulkModulus; }
    double getTortuosity() const { return tortuosity; }
    double getBiotCoefficient() const { return biotCoefficient; }
    double getBiotModulus() const { return biotModulus; }
    double getAbsolutePermeability() const { return permeability; }

private:
    double density;
    double shearModulus;
    double porosity;
    double bulkModulus;
    double tortuosity;
    double biotCoefficient;
    double biotModulus;
    double permeability;
};

#endif // ROCK_PROPERTIES_H