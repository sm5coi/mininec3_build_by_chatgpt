#ifndef MININEC_HPP
#define MININEC_HPP

#include <vector>
#include <array>
#include <complex>
#include <stdexcept>

// 3. Exakt C++-implementation av I6-logiken
// Här är en C++-funktion som du kan lägga i din MininecImpedanceSolver
//och kalla från psiImpedanceKernel innan Gauss-loopen:

struct KernelSetup {
    double F2;          // F2 i BASIC
    int    L;           // Gauss-ordning
    double I6;          // I6! i BASIC (0 => ingen exact kernel)
};


// -------------------------------------------------------------
// 1-baserad geometri, så att BASIC-logiken kan följas rakt av.
// Alla vektorer ska ha storlek >= maxIndex+1, index 0 ignoreras.
// -------------------------------------------------------------
struct Geometry {
    int N = 0;   // antal pulser (okända)
    int G = 1;   // antal "images" (1 eller -1, i praktiken 1 eller 2)

    // Nodpunkter (1..MS)
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;

    // Trådegenskaper (1..MW)
    std::vector<double> A;   // radie
    std::vector<double> CA;  // cos(alpha)
    std::vector<double> CB;  // cos(beta)
    std::vector<double> CG;  // cos(gamma)
    std::vector<double> S;   // segmentlängd per tråd

    // Kopplingar (pulser -> trådar / segment)
    // C[puls][0/1] motsvarar C%(I,1/2)
    std::vector<std::array<int,2>> C;
    // W[puls] motsvarar W%(I): vilken tråd pulsen I ligger på
    std::vector<int> W;

    // J2[wire][0/1] motsvarar J2(W,1/2): första/sista nod på tråd
    std::vector<std::array<int,2>> J2;
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
// MININEC impedans-matrislösare (hybrid: modern klass + BASIC-logik)
// -------------------------------------------------------------
class MininecImpedanceSolver {
public:
    double W2;
    // k = vågtal = 2*pi / lambda
    // srm = SRM i BASIC ("small radius"-gräns)
    MininecImpedanceSolver(const Geometry& geom, double k, double srm);

    // Bygger impedansmatrisen Z (N x N) enligt MININEC3
    void build(ImpedanceMatrix& Z) const;

private:
    using Complex = std::complex<double>;

    const Geometry& g_;
    double k_;           // vågtal
    double w_;           // samma som k_ här (BASIC W efter 1146)
    double w2_;          // W2 = W*W/2
    double srm_;         // SRM
    double pi_;

    // Q-vektor & elliptiska integral-koefficienter
    double Q_[15];       // 1..14 används
    double C_[10];       // C0..C9

    // Hjälpare
    int sgn(int x) const { return (x > 0) - (x < 0); }

    // GOSUB 28: kernel evaluation I2 & I3
    void kernel28(
        double& T3, double& T4,
        double X2, double Y2, double Z2,
        double V1, double V2, double V3,
        double T, int P4, double A2,
        double I6   // samma som I6! i BASIC
        ) const;

    // GOSUB 87/102: PSI(P1,P2,P3) = T1 + j*T2
    // scalar = true => GOSUB 87 (skalärpotential)
    // scalar = false => GOSUB 102 (vektorpotential)

    inline int baseP1(int I) const
    {
        return 2 * g_.W[I] + I - 1;
    };

    inline int baseP2(int J) const
    {
        return 2 * g_.W[J] + J - 1;
    };


    void psiImpedanceKernel(
        int I, int J,
        double& outRe, double& outIm) const;


/*
        void psiImpedanceKernel(bool scalar,
                            int I, int J,
                            int P1, int P2, int P3, int P4,
                            double& T1, double& T2) const;
*/

    // Grad(Φ)-delen (rader 274–312) för ett givet I,J
    void scalarGradientContribution(int I, int J,
                                    int J1, int J2,
                                    double F4, double F5,
                                    int P1_base, int P2_base,
                                    double& U1_out, double& U2_out) const;

    // Initierar Q[] och C[] enligt DATA-raderna
    void initConstants();

    KernelSetup setupKernelI6(int I, int J,
                              int P2, int P3, int P4,
                              double D0, double D3, double S4,
                              double FVS,        // FVS från GOSUB 87 (oftast 1 för impedans)
                              char   Cmode       // C$ ("N" = no exact kernel, annars tillåten)
                              ) const;


    void computePsiGauss(
        bool   scalar,
        int    I, int J,    // pulsindex (behövs för setupKernelI6)
        int    P1,          // "M"-plats (mittpunkt eller nod beroende på scalar/vector)
        int    P2, int P3,  // U- och V-pos i originalet
        int    P4,          // trådens index (seg-index för radie A[P4], S[P4])
        double FVS,         // "FVS" från BASIC (ofta 1.0 för impedansfallet)
        char   Cmode,       // C$ ("N" = no exact kernel)
        double& outT1,      // Re(PSI)
        double& outT2       // Im(PSI)
        ) const;

    void psiScalarKernel(
        int I, int J,
        int P1, int P2, int P3, int P4,
        double& T1, double& T2) const;

    void psiVectorKernel(
        int I, int J,
        int P1, int P2, int P3, int P4,
        double& T1, double& T2) const;

};
#endif // MININEC_HPP
