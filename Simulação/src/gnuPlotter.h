// GnuplotPlotter.h
#ifndef GNUPLOT_PLOTTER_H
#define GNUPLOT_PLOTTER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

class GnuplotPlotter {
public:
    GnuplotPlotter(const std::string& gnuplotPath = "c:/gnuplot");
    void plotResults(const std::vector<std::vector<double>>& resultsV, 
                                 const std::vector<std::vector<double>>& resultsX, 
                                 const std::vector<double>& timeVector);

private:
    std::string gnuplotPath;
    std::string tempDataFileX;
    void writeDataToFile(const std::vector<std::vector<double>>& results, 
                                     const std::vector<double>& timeVector, 
                                     const std::string& fileName) ;
};

#endif // GNUPLOT_PLOTTER_H