#ifndef MININEC_HPP
#define MININEC_HPP

#include <vector>
//#include <array>
#include <complex>
//#include <stdexcept>
#include <algorithm>
#include <cmath>
#include "geometry.hpp"

// -------------------------------------------------------------
// Kernel-setup: motsvarar I6!-logiken i BASIC
// -------------------------------------------------------------
struct KernelSetup {
    double F2;  // F2 i BASIC
    int    L;   // Gauss-ordning (1, 3 eller 7)
    double I6;  // I6! i BASIC (0 => ingen exact kernel)
};

// -------------------------------------------------------------
// Enkel impedansmatris med 1-baserad indexering
// -------------------------------------------------------------
class ImpedanceMatrix {
public:
    using Complex = std::complex<double>;

    ImpedanceMatrix() : n_(0), stride_(0) {}
    explicit ImpedanceMatrix(int n) { resize(n); }

    void resize(int n) {
        if (n <= 0) {
            n_ = 0;
            stride_ = 0;
            data_.clear();
            return;
        }
        n_ = n;
        stride_ = n_ + 1;      // 1-baserad
        data_.assign(stride_ * stride_, Complex(0.0, 0.0));
    }

    int size() const { return n_; }

    void zero() {
        std::fill(data_.begin(), data_.end(), Complex(0.0, 0.0));
    }

    Complex& operator()(int i, int j) {
        return data_.at(i * stride_ + j);
    }

    const Complex& operator()(int i, int j) const {
        return data_.at(i * stride_ + j);
    }

private:
    int n_;
    int stride_;
    std::vector<Complex> data_;
};

// -------------------------------------------------------------
// MININEC-liknande impedans-matrislösare
// -------------------------------------------------------------
class MininecImpedanceSolver {
public:
    using Complex = std::complex<double>;

    // k = vågtal = 2*pi / lambda
    // srm = SRM i BASIC ("small radius"-gräns)
    MininecImpedanceSolver(const Geometry& geom, double k, double srm);

    // Bygger impedansmatrisen Z (N x N)
    void build(ImpedanceMatrix& Z) const;

private:
    const Geometry& g_;
    double k_;    // vågtal
    double w_;    // samma som k_ här (BASIC W)
    double srm_;  // SRM
    double pi_;
    double W2_;   // vektorpotentialens skalningsfaktor

    // Q-vektor & elliptiska integral-koefficienter
    double Q_[15];  // 1..14 används
    double C_[10];  // C0..C9

    // Hjälpare
    int sgn(int x) const { return (x > 0) - (x < 0); }

    // Kernel (GOSUB 28)
    void kernel28(
        double& T3, double& T4,
        double X2, double Y2, double Z2,
        double V1, double V2, double V3,
        double T, int P4, double A2,
        double I6   // samma som I6! i BASIC
        ) const;

    KernelSetup setupKernelI6(
        int I, int J,
        double P2, double P3, int P4,
        double D0, double D3, double S4,
        double FVS,    // FVS från GOSUB 87 (oftast 1 för impedans)
        char   Cmode   // C$ ("N" = no exact kernel, annars tillåten)
        ) const;

    // Gemensam Gauss-integrationskärna för GOSUB 87/102
    void computePsiGauss(
        bool   scalar,
        int    I, int J,   // pulsindex
        double P1,         // "M"-plats (mittpunkt eller nod)
        double P2, double P3,  // U- och V-pos (kan vara N+0.5)
        int    P4,         // segment-index för radie A[P4], S[P4]
        double FVS,        // "FVS" (ofta 1.0 för impedansfallet)
        char   Cmode,      // C$ ("N" = no exact kernel)
        double& outT1,     // Re(PSI)
        double& outT2      // Im(PSI)
        ) const;

    // GOSUB 87 – skalärpotential-psi
    void psiScalarKernel(
        int I, int J,
        double P1, double P2, double P3, int P4,
        double& T1, double& T2) const;

    // GOSUB 102 – vektorpotential-psi
    void psiVectorKernel(
        int I, int J,
        double P1, double P2, double P3, int P4,
        double& T1, double& T2) const;

    // Impedansbidrag för ett par (I,J)
    void psiImpedanceKernel(
        int I, int J,
        double& outRe, double& outIm) const;

    // ⬇⬇⬇ Lägg till denna rad:
    void gradPhiContribution(
        int I, int J,
        double& gRe, double& gIm) const;

    // Initierar Q[] och C[] enligt DATA-raderna
    void initConstants();
};

#endif // MININEC_HPP
