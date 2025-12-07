#include "mininec.hpp"
#include <iostream>

MininecImpedanceSolver::MininecImpedanceSolver(
    const Geometry& geom,
    double k, double srm)
    : g_(geom),
    k_(k),
    w_(k),
    srm_(srm),
    pi_(std::acos(-1.0)),
    W2_(0.0)
{
    initConstants();

    // W2 motsvarar ofta W0/(8*pi) i MININEC. Här tar vi en enkel variant
    // som du kan justera när du jämför mot originalet:
    // t.ex. W0 ≈ 120*pi (frirymdsimpedans), W2 = W0/(8*pi) = 15.
    W2_ = 15.0;
}

void MininecImpedanceSolver::initConstants() {
    // Q(1..14) från MININEC3.BAS
    Q_[1]  = 0.288675135;
    Q_[2]  = 0.5;
    Q_[3]  = 0.430568156;
    Q_[4]  = 0.173927423;
    Q_[5]  = 0.169990522;
    Q_[6]  = 0.326072577;
    Q_[7]  = 0.480144928;
    Q_[8]  = 0.050614268;
    Q_[9]  = 0.398333239;
    Q_[10] = 0.111190517;
    Q_[11] = 0.262766205;
    Q_[12] = 0.156853323;
    Q_[13] = 0.091717321;
    Q_[14] = 0.181341892;

    // C0..C9 från MININEC3.BAS (elliptisk integral approx)
    C_[0] = 1.38629436112;
    C_[1] = 0.09666344259;
    C_[2] = 0.03590092383;
    C_[3] = 0.03742563713;
    C_[4] = 0.01451196212;
    C_[5] = 0.5;
    C_[6] = 0.12498593397;
    C_[7] = 0.06880248576;
    C_[8] = 0.0332835346;
    C_[9] = 0.00441787012;
}

// -------------------------------------------------------------
// kernel28 – exact + exp(-jkR)/R
// -------------------------------------------------------------
void MininecImpedanceSolver::kernel28(
    double& T3, double& T4,
    double X2, double Y2, double Z2,
    double V1, double V2, double V3,
    double T, int P4, double A2,
    double I6   // samma som I6! i BASIC
    ) const
{
    double X3, Y3, Z3;

    // Linjär interpolation beroende på "T"
    X3 = X2 + T * (V1 - X2);
    Y3 = Y2 + T * (V2 - Y2);
    Z3 = Z2 + T * (V3 - Z2);

    double D3 = X3*X3 + Y3*Y3 + Z3*Z3;
    double D;

    if (g_.A[P4] <= srm_) {
        D = std::sqrt(D3);
    } else {
        D = D3 + A2;
        if (D > 0.0) D = std::sqrt(D);
    }

    // ----- exakt kernel (elliptisk integral) om I6 != 0 -----
    if (I6 != 0.0) {
        double B = D3 / (D3 + 4.0 * A2);

        double W0 = C_[0] + B * (C_[1] + B * (C_[2] + B * (C_[3] + B * C_[4])));
        double W1 = C_[5] + B * (C_[6] + B * (C_[7] + B * (C_[8] + B * C_[9])));
        double V0 = (W0 - W1 * std::log(B)) * std::sqrt(1.0 - B);

        T3 += (V0 + std::log(D3 / (64.0 * A2)) / 2.0) / (pi_ * g_.A[P4]) - 1.0 / D;
    }

    // ----- exp(-j*k*R)/R-delen (alltid) -----
    double B1 = D * w_;
    T3 += std::cos(B1) / D;
    T4 -= std::sin(B1) / D;
}

// -------------------------------------------------------------
// setupKernelI6 – bestämmer F2, L och I6 enligt MININEC-logiken
// (här något förenklad: inga J2-beroende specialfall)
// -------------------------------------------------------------
KernelSetup MininecImpedanceSolver::setupKernelI6(
    int /*I*/, int /*J*/,
    double /*P2*/, double /*P3*/, int P4,
    double D0, double D3, double S4,
    double FVS,
    char   Cmode
    ) const
{
    KernelSetup ks;
    ks.F2 = 1.0;
    ks.L  = 7;
    ks.I6 = 0.0;

    double T = (D0 + D3) / g_.S[P4];

    if (T <= 1.1 && Cmode != 'N') {
        if (g_.A[P4] > srm_) {
            ks.F2 = 2.0; // förenklat: P3-P2 ≈ 1
            ks.I6 = (1.0 - std::log(S4 / ks.F2 / 8.0 / g_.A[P4]))
                    / (pi_ * g_.A[P4]);
            return ks;
        }
    }

    if (T > 6.0)  ks.L = 3;
    if (T > 10.0) ks.L = 1;

    return ks;
}

// -------------------------------------------------------------
// computePsiGauss – gemensam Gauss-integral för GOSUB 87/102
// -------------------------------------------------------------
void MininecImpedanceSolver::computePsiGauss(
    bool   scalar,
    int    /*I*/, int /*J*/,
    double P1,
    double P2, double P3,
    int    P4,
    double FVS,
    char   Cmode,
    double& outT1,
    double& outT2
    ) const
{
    (void)FVS; // används ej i denna förenklade version

    // 1) S(M)
    double X1, Y1, Z1;
    if (scalar) {
        // Mittpunkt av segment P1 => nod P1 och P1+1
        int i4 = static_cast<int>(P1);
        int i5 = i4 + 1;
        X1 = 0.5 * (g_.X[i4] + g_.X[i5]);
        Y1 = 0.5 * (g_.Y[i4] + g_.Y[i5]);
        Z1 = 0.5 * (g_.Z[i4] + g_.Z[i5]);
    } else {
        // nod P1
        int i1 = static_cast<int>(P1);
        X1 = g_.X[i1];
        Y1 = g_.Y[i1];
        Z1 = g_.Z[i1];
    }

    // 2) S(U)-S(M)
    double X2, Y2, Z2;
    int i4 = static_cast<int>(P2);
    if (i4 != static_cast<int>(P2)) {
        int i5 = i4 + 1;
        X2 = 0.5 * (g_.X[i4] + g_.X[i5]) - X1;
        Y2 = 0.5 * (g_.Y[i4] + g_.Y[i5]) - Y1;
        Z2 = 0.5 * (g_.Z[i4] + g_.Z[i5]) - Z1;
    } else {
        X2 = g_.X[i4] - X1;
        Y2 = g_.Y[i4] - Y1;
        Z2 = g_.Z[i4] - Z1;
    }

    // 3) S(V)-S(M)
    double V1, V2, V3;
    i4 = static_cast<int>(P3);
    if (i4 != static_cast<int>(P3)) {
        int i5 = i4 + 1;
        V1 = 0.5 * (g_.X[i4] + g_.X[i5]) - X1;
        V2 = 0.5 * (g_.Y[i4] + g_.Y[i5]) - Y1;
        V3 = 0.5 * (g_.Z[i4] + g_.Z[i5]) - Z1;
    } else {
        V1 = g_.X[i4] - X1;
        V2 = g_.Y[i4] - Y1;
        V3 = g_.Z[i4] - Z1;
    }

    // 4) D0, D3, S4
    double D0 = std::sqrt(X2*X2 + Y2*Y2 + Z2*Z2);
    double D3 = std::sqrt(V1*V1 + V2*V2 + V3*V3);
    double S4 = std::fabs(P3 - P2) * g_.S[P4];
    double A2 = g_.A[P4] * g_.A[P4];

    // 5) Kernel-setup
    KernelSetup ks = setupKernelI6(
        0, 0, P2, P3, P4,
        D0, D3, S4,
        FVS,
        Cmode
        );

    if (ks.L == 0) {
        outT1 = 0.0;
        outT2 = 0.0;
        return;
    }

    // 6) Gauss-kvadratur
    double T1 = 0.0, T2 = 0.0;
    double F2 = ks.F2;
    int    L  = ks.L;
    int    I5 = 2 * L;

    while (true) {
        double T3 = 0.0;
        double T4 = 0.0;

        double T = (Q_[L] + 0.5) / F2;
        kernel28(T3, T4, X2, Y2, Z2, V1, V2, V3, T, P4, A2, ks.I6);

        T = (0.5 - Q_[L]) / F2;
        kernel28(T3, T4, X2, Y2, Z2, V1, V2, V3, T, P4, A2, ks.I6);

        L += 1;
        T1 += Q_[L] * T3;
        T2 += Q_[L] * T4;

        L += 1;
        if (L >= I5) break;
    }

    T1 = S4 * (T1 + ks.I6);
    T2 = S4 * T2;

    outT1 = T1;
    outT2 = T2;
}

// -------------------------------------------------------------
// GOSUB 87 – skalärpotential-psi
// -------------------------------------------------------------
void MininecImpedanceSolver::psiScalarKernel(
    int I, int J,
    double P1, double P2, double P3, int P4,
    double& T1, double& T2) const
{
    (void)I; (void)J;
    double FVS  = 1.0;
    char   Cmod = 'Y';
    computePsiGauss(true, I, J, P1, P2, P3, P4, FVS, Cmod, T1, T2);
}

// -------------------------------------------------------------
// GOSUB 102 – vektorpotential-psi
// -------------------------------------------------------------
void MininecImpedanceSolver::psiVectorKernel(
    int I, int J,
    double P1, double P2, double P3, int P4,
    double& T1, double& T2) const
{
    (void)I; (void)J;
    double FVS  = 1.0;
    char   Cmod = 'Y';
    computePsiGauss(false, I, J, P1, P2, P3, P4, FVS, Cmod, T1, T2);
}

// -------------------------------------------------------------
// psiImpedanceKernel – enkel kombination av scalar + vector
// (OBS: förenklad: ingen full grad(Φ)-term ännu)
// -------------------------------------------------------------
void MininecImpedanceSolver::psiImpedanceKernel(
    int I, int J,
    double& outRe, double& outIm) const
{
    // Vi använder segmentindex direkt:
    // P1: observation på segment I (mittpunkt)
    // P2,P3: täcker segment J (nod J..J+1)
    double P1 = static_cast<double>(I);
    double P2 = static_cast<double>(J);
    double P3 = static_cast<double>(J + 1);
    int    P4 = J;

    double T1s = 0.0, T2s = 0.0;
    double T1v = 0.0, T2v = 0.0;

    psiScalarKernel(I, J, P1, P2, P3, P4, T1s, T2s);
    psiVectorKernel(I, J, P1, P2, P3, P4, T1v, T2v);

    // Enkel kombination: Z ∝ k * (scalar + vector)
    outRe = k_ * (T1s + T1v);
    outIm = k_ * (T2s + T2v);
}

// -------------------------------------------------------------
// Bygg impedansmatrisen genom att anropa psiImpedanceKernel
// -------------------------------------------------------------
void MininecImpedanceSolver::build(ImpedanceMatrix& Z) const {
    const int N = g_.N;
    if (N <= 0) {
        throw std::runtime_error("Geometry.N must be > 0");
    }
    if (Z.size() != N) {
        Z.resize(N);
    }
    Z.zero();

    for (int I = 1; I <= N; ++I) {
        for (int J = 1; J <= N; ++J) {
            double zr = 0.0, zi = 0.0;
            psiImpedanceKernel(I, J, zr, zi);
            Z(I, J) = Complex(zr, zi);
        }
    }
}
