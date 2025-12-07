#ifndef SOLVEZMTX_HPP
#define SOLVEZMTX_HPP


#include <limits>
#include <vector>
#include <complex>
#include "mininec.hpp"


std::vector<std::complex<double>> solveLinearSystem(
    std::vector<std::vector<std::complex<double>>> A,
    const std::vector<std::complex<double>>& b);

std::complex<double> computeInputImpedance(
    const ImpedanceMatrix& Z,
    int feedSegment);  // 1-baserat index

#endif // SOLVEZMTX_HPP
