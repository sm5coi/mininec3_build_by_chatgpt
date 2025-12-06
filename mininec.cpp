#include "mininec.hpp"
#include <cmath>
#include <iostream>

// GOSUB 87 – skalärpotential-psi
void MininecImpedanceSolver::psiScalarKernel(
    int I, int J,
    int P1, int P2, int P3, int P4,
    double& T1, double& T2) const
{
    // Här kan du lägga in specialfallet 87–93 (super-nära med slutna log-uttryck)
    // och annars bara kalla computePsiGauss(true, ...)

    double FVS  = 1.0;   // i impedansfallet
    char   Cmod = 'Y';   // eller 'N' om du vill stänga exact kernel

    computePsiGauss(true, I, J, P1, P2, P3, P4, FVS, Cmod, T1, T2);
}

// GOSUB 102 – vektorpotential-psi
void MininecImpedanceSolver::psiVectorKernel(
    int I, int J,
    int P1, int P2, int P3, int P4,
    double& T1, double& T2) const
{
    double FVS  = 1.0;   // i impedansfallet normalt
    char   Cmod = 'Y';

    computePsiGauss(false, I, J, P1, P2, P3, P4, FVS, Cmod, T1, T2);
}

MininecImpedanceSolver::MininecImpedanceSolver(
    const Geometry& geom,
    double k, double srm)
    : g_(geom),
    k_(k),
    w_(k),
    w2_(0.5 * k * k),
    srm_(srm),
    pi_(std::acos(-1.0))
{
    initConstants();
}

void MininecImpedanceSolver::initConstants() {
    // Q(1..14) från rader 1506–1508 i MININEC3.BAS
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

    // C0..C9 från rader 1511–1512
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

// Beräkna PSI(P1,P2,P3) = T1 + j*T2 med Gauss + kernel28 + I6-hantering
// scalar = true  => GOSUB 87 (skalärpotential)
// scalar = false => GOSUB 102 (vektorpotential)
void MininecImpedanceSolver::computePsiGauss(
    bool   scalar,
    int    I, int J,    // pulsindex (behövs för setupKernelI6)
    int    P1,          // "M"-plats (mittpunkt eller nod beroende på scalar/vector)
    int    P2, int P3,  // U- och V-pos i originalet
    int    P4,          // trådens index (seg-index för radie A[P4], S[P4])
    double FVS,         // "FVS" från BASIC (ofta 1.0 för impedansfallet)
    char   Cmode,       // C$ ("N" = no exact kernel)
    double& outT1,      // Re(PSI)
    double& outT2       // Im(PSI)
    ) const
{
    // -------------------------
    // 1) S(M) i X1,Y1,Z1
    // -------------------------
    double X1, Y1, Z1;
    if (scalar) {
        // GOSUB 87: mittpunkt av segment P1
        int I4 = static_cast<int>(P1);
        int I5 = I4 + 1;
        X1 = 0.5 * (g_.X[I4] + g_.X[I5]);
        Y1 = 0.5 * (g_.Y[I4] + g_.Y[I5]);
        Z1 = 0.5 * (g_.Z[I4] + g_.Z[I5]);
    } else {
        // GOSUB 102: nod P1
        X1 = g_.X[P1];
        Y1 = g_.Y[P1];
        Z1 = g_.Z[P1];
    }

    // -------------------------
    // 2) S(U)-S(M) → X2,Y2,Z2  (113–123)
    // -------------------------
    double X2, Y2, Z2;
    int I4 = static_cast<int>(P2);
    if (I4 != P2) {
        int I5 = I4 + 1;
        X2 = 0.5 * (g_.X[I4] + g_.X[I5]) - X1;
        Y2 = 0.5 * (g_.Y[I4] + g_.Y[I5]) - Y1;
        Z2 = k_ * 0.5 * (g_.Z[I4] + g_.Z[I5]) - Z1;
    } else {
        X2 = g_.X[P2] - X1;
        Y2 = g_.Y[P2] - Y1;
        Z2 = k_ * g_.Z[P2] - Z1;
    }

    // -------------------------
    // 3) S(V)-S(M) → V1,V2,V3 (124–133)
    // -------------------------
    double V1, V2, V3;
    I4 = static_cast<int>(P3);
    if (I4 != P3) {
        int I5 = I4 + 1;
        V1 = 0.5 * (g_.X[I4] + g_.X[I5]) - X1;
        V2 = 0.5 * (g_.Y[I4] + g_.Y[I5]) - Y1;
        V3 = k_ * 0.5 * (g_.Z[I4] + g_.Z[I5]) - Z1;
    } else {
        V1 = g_.X[P3] - X1;
        V2 = g_.Y[P3] - Y1;
        V3 = k_ * g_.Z[P3] - Z1;
    }

    // -------------------------
    // 4) D0, D3, S4 (135–143)
    // -------------------------
    double D0 = X2*X2 + Y2*Y2 + Z2*Z2;
    if (D0 > 0.0) D0 = std::sqrt(D0);

    double D3 = V1*V1 + V2*V2 + V3*V3;
    if (D3 > 0.0) D3 = std::sqrt(D3);

    double S4 = (P3 - P2) * g_.S[P4];   // (143)

    double A2 = g_.A[P4] * g_.A[P4];

    // -------------------------
    // 5) Kernel-setup (145–167)
    // -------------------------
    KernelSetup ks = setupKernelI6(
        I, J,
        P2, P3, P4,
        D0, D3, S4,
        FVS,
        Cmode
        );

    // Specialfallet A<=SRM & FVS=1 hoppar till slutna uttryck (rad 91/106).
    // Om du redan har portat de slutna uttrycken kan du testa på ks.L == 0 här:
    if (ks.L == 0) {
        // TODO: anropa din "short formula"-funktion här och sätt outT1,outT2.
        outT1 = 0.0;
        outT2 = 0.0;
        return;
    }

    // -------------------------
    // 6) Gauss-kvadratur (168–178)
    // -------------------------
    double T1 = 0.0;
    double T2 = 0.0;

    double F2 = ks.F2;
    int    L  = ks.L;
    int    I5 = 2 * L;

    while (true) {
        double T3 = 0.0;
        double T4 = 0.0;

        // Första Gauss-punkten: T = (Q(L)+0.5)/F2
        double T = (Q_[L] + 0.5) / F2;
        kernel28(T3, T4, X2, Y2, Z2, V1, V2, V3, T, P4, A2, ks.I6);

        // Andra Gauss-punkten: T = (0.5-Q(L))/F2
        T = (0.5 - Q_[L]) / F2;
        kernel28(T3, T4, X2, Y2, Z2, V1, V2, V3, T, P4, A2, ks.I6);

        // Nästa Gaussvikt: Q(L+1)
        L += 1;
        T1 += Q_[L] * T3;
        T2 += Q_[L] * T4;

        // Öka L igen enligt BASIC
        L += 1;
        if (L >= I5) break;
    }

    // -------------------------
    // 7) Slutlig skalning (179–180)
    // -------------------------
    T1 = S4 * (T1 + ks.I6);
    T2 = S4 * T2;

    outT1 = T1;
    outT2 = T2;
    int onodig = 0;
    int onodig1 = 0;
}



// -------------------------------------------------------------
// GOSUB 28 – kernel evaluation I2 & I3 (rad 27–53)
// -------------------------------------------------------------
/*
// Old version
// void MininecImpedanceSolver::kernel28(
//     double& T3, double& T4,
//     double X2, double Y2, double Z2,
//     double V1, double V2, double V3,
//     double T, int P4, double A2) const
// {
//     double X3, Y3, Z3;
//     double D3, D;
//     double B, W0, W1, V0;
//     double B1;

//     // 28–35: linjär interpolation beroende på K (tecken)
//     if (k_ >= 0) {
//         X3 = X2 + T * (V1 - X2);
//         Y3 = Y2 + T * (V2 - Y2);
//         Z3 = Z2 + T * (V3 - Z2);
//     } else {
//         X3 = V1 + T * (X2 - V1);
//         Y3 = V2 + T * (Y2 - V2);
//         Z3 = V3 + T * (Z2 - V3);
//     }

//     D3 = X3*X3 + Y3*Y3 + Z3*Z3;

//     // 38–40: small-radius-mod
//     if (g_.A[P4] <= srm_) {
//         D = std::sqrt(D3);
//     } else {
//         D = D3 + A2;
//         if (D > 0) D = std::sqrt(D);
//     }

//     // Här förenklar vi: vi använder alltid reducerad kernel (I6!=0-del
//     // från BASIC är utelämnad för överskådlighet).
//     // Vill du ha full "exact kernel" med elliptiska integraler:
//     // portera raderna 41–48 rad-för-rad och använd C_[].

//     // Exp(-j*k*r)/r
//     B1 = D * w_;
//     T3 += std::cos(B1) / D;
//     T4 -= std::sin(B1) / D;
// }

// // -------------------------------------------------------------
// // PSI(P1,P2,P3) för impedansmatrisen (GOSUB 87/102)
// // scalar = true  => skalärpotential (rad 84–99 + 113–181)
// // scalar = false => vektorpotential (rad 100–111 + 113–181)
// // OBS: här är "exact kernel"-kriterierna förenklade något.
// // -------------------------------------------------------------
*/
// New version 2025-12-02

void MininecImpedanceSolver::kernel28(
    double& T3, double& T4,
    double X2, double Y2, double Z2,
    double V1, double V2, double V3,
    double T, int P4, double A2,
    double I6   // samma som I6! i BASIC
    ) const
{
    double X3, Y3, Z3;

    if (k_ >= 0) {
        X3 = X2 + T * (V1 - X2);
        Y3 = Y2 + T * (V2 - Y2);
        Z3 = Z2 + T * (V3 - Z2);
    } else {
        X3 = V1 + T * (X2 - V1);
        Y3 = V2 + T * (Y2 - V2);
        Z3 = V3 + T * (Z2 - V3);
    }

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



KernelSetup MininecImpedanceSolver::setupKernelI6(
    int I, int J,
    int P2, int P3, int P4,
    double D0, double D3, double S4,
    double FVS,        // FVS från GOSUB 87 (oftast 1 för impedans)
    char   Cmode       // C$ ("N" = no exact kernel, annars tillåten)
    ) const
{
    KernelSetup ks;
    ks.F2 = 1.0;
    ks.L  = 7;
    ks.I6 = 0.0;

    // 151: T = (D0 + D3) / S(P4)
    double T = (D0 + D3) / g_.S[P4];

    // 152–158: kriterier för exact kernel
    if (T <= 1.1 && Cmode != 'N') {
        int wi = g_.W[I];   // W%(I)
        int wj = g_.W[J];   // W%(J)

        int wi_1 = g_.J2[wi][0];
        int wi_2 = g_.J2[wi][1];
        int wj_1 = g_.J2[wj][0];
        int wj_2 = g_.J2[wj][1];

        bool sameEnd =
            (wi_1 == wj_1) ||
            (wi_1 == wj_2) ||
            (wi_2 == wj_1) ||
            (wi_2 == wj_2);

        if (sameEnd) {
            // Motsvarar rad 160–163
            if (g_.A[P4] <= srm_) {
                // Rad 161: specialfall med slutna formler
                // IF FVS = 1 THEN 91 ELSE GOTO 106
                // => Hanteras *utanför* denna funktion:
                //    här signalerar vi bara att ingen numerisk kernel ska köras.
                // Du kan t.ex. markera detta genom att sätta L=0:
                ks.L  = 0;
                ks.F2 = 1.0;
                ks.I6 = 0.0;
                return ks;
            } else {
                // A(P4) > SRM: exact kernel med elliptisk integral
                ks.F2 = 2.0 * (P3 - P2);

                // 163: I6! = (1 - LOG(S4 / F2 / 8 / A(P4))) / P / A(P4)
                ks.I6 = (1.0 - std::log(S4 / ks.F2 / 8.0 / g_.A[P4]))
                        / (pi_ * g_.A[P4]);

                // L lämnas = 7 (full Gauss)
                return ks;
            }
        }
    }

    // 165–166: ingen exact kernel, justera L beroende på T
    if (T > 6.0)  ks.L = 3;
    if (T > 10.0) ks.L = 1;

    // F2 = 1, I6 = 0 (ingen exact kernel)
    return ks;
}


void MininecImpedanceSolver::psiImpedanceKernel(
    bool scalar,
    int I, int J,
    int P1, int P2, int P3, int P4,
    double& T1, double& T2) const
{
    // Förenklad version: vi behåller specialfallen 91/106 och
    // en "ren" Gauss-integration 145–181, men hoppar över
    // C$- och J2-baserade finfilter (153–160).

    // --- Specialfall för mycket små radier (87–93 och 102–108) ---
    if (k_ >= 1.0 && g_.A[P4] <= srm_) {
        // Impedans-matrisfallet använder P3 = P2+1 resp. .5
        if (scalar) {
            // 90–93: P3=P2+1 och P1=(P2+P3)/2
            if ((P3 == P2 + 1) && (P1 == (P2 + P3) / 2)) {
                T1 = 2.0 * std::log(g_.S[P4] / g_.A[P4]);
                T2 = -w_ * g_.S[P4];
                return;
            }
        } else {
            // 105–108: I=J och P3 = P2 + 0.5
            if (I == J && std::abs(P3 - (P2 + 0.5)) < 1e-9) {
                T1 = std::log(g_.S[P4] / g_.A[P4]);
                T2 = -w_ * g_.S[P4] / 2.0;
                return;
            }
        }
    }

    double FVS = 1.0;
    char Cmode = 'Y';

    computePsiGauss(scalar, I, J, P1, P2, P3,  P4, FVS,  Cmode, T1, T2 );

/*
        // ----- S(M) i (X1,Y1,Z1) -----
    double X1, Y1, Z1;
    if (scalar) {
        // 94–99: mittpunkt av segment P1
        int I4 = static_cast<int>(P1);
        int I5 = I4 + 1;
        X1 = 0.5 * (g_.X[I4] + g_.X[I5]);
        Y1 = 0.5 * (g_.Y[I4] + g_.Y[I5]);
        Z1 = 0.5 * (g_.Z[I4] + g_.Z[I5]);
    } else {
        // 109–111: nod P1
        X1 = g_.X[P1];
        Y1 = g_.Y[P1];
        Z1 = g_.Z[P1];
    }

    // ----- S(U)-S(M) i (X2,Y2,Z2)  (113–123) -----
    double X2, Y2, Z2;
    int I4 = static_cast<int>(P2);
    if (I4 != P2) {
        int I5 = I4 + 1;
        X2 = 0.5 * (g_.X[I4] + g_.X[I5]) - X1;
        Y2 = 0.5 * (g_.Y[I4] + g_.Y[I5]) - Y1;
        Z2 = k_ * 0.5 * (g_.Z[I4] + g_.Z[I5]) - Z1;
    } else {
        X2 = g_.X[P2] - X1;
        Y2 = g_.Y[P2] - Y1;
        Z2 = k_ * g_.Z[P2] - Z1;
    }

    // ----- S(V)-S(M) i (V1,V2,V3) (124–133) -----
    double V1, V2, V3;
    I4 = static_cast<int>(P3);
    if (I4 != P3) {
        int I5 = I4 + 1;
        V1 = 0.5 * (g_.X[I4] + g_.X[I5]) - X1;
        V2 = 0.5 * (g_.Y[I4] + g_.Y[I5]) - Y1;
        V3 = k_ * 0.5 * (g_.Z[I4] + g_.Z[I5]) - Z1;
    } else {
        V1 = g_.X[P3] - X1;
        V2 = g_.Y[P3] - Y1;
        V3 = k_ * g_.Z[P3] - Z1;
    }

    // ----- D0, D3, A2, S4 (135–143) -----
    double D0 = X2*X2 + Y2*Y2 + Z2*Z2;
    if (D0 > 0) D0 = std::sqrt(D0);
    double D3 = V1*V1 + V2*V2 + V3*V3;
    if (D3 > 0) D3 = std::sqrt(D3);

    double A2 = g_.A[P4] * g_.A[P4];
    double S4 = (P3 - P2) * g_.S[P4];

    double FVS = 1.0;
    char Cmode = 'Y';

    KernelSetup ks = setupKernelI6(
        I, J,
        P2, P3, P4,
        D0, D3, S4,
        FVS,    // 1 för impedans (GOSUB 87), annat vid behov
        Cmode   // t.ex. 'Y' för "tillåt exact kernel"
        );


    // ----- Gaussordning L (145–167, förenklad) -----
    T1 = 0.0;
    T2 = 0.0;

    double T = (D0 + D3) / g_.S[P4];

    if (T > 6.0)  ks.L = 3;
    if (T > 10.0) ks.L = 1;
    int I5 = ks.L + ks.L;

    // ----- Gauss-kvadratur (168–178) -----
    while (true) {
        double T3 = 0.0;
        double T4 = 0.0;

        T = (Q_[ks.L] + 0.5) / ks.F2;
        kernel28(T3, T4, X2, Y2, Z2, V1, V2, V3, T, P4, A2, ks.I6);

        T = (0.5 - Q_[ks.L]) / ks.F2;
        kernel28(T3, T4, X2, Y2, Z2, V1, V2, V3, T, P4, A2, ks.I6);

        ks.L += 1;
        T1 += Q_[ks.L] * T3;
        T2 += Q_[ks.L] * T4;
        ks.L += 1;
        if (ks.L >= I5) break;
    }

    T1 = S4 * (T1 + ks.I6);
    T2 = S4 * T2;

*/

}

// -------------------------------------------------------------
// Grad(Φ)-bidrag för ett I,J (rader 273–312), lätt förenklad
// -------------------------------------------------------------
void MininecImpedanceSolver::scalarGradientContribution(
    int I, int J,
    int J1, int J2,
    double F4, double F5,
    int P1_base, int P2_base,
    double& U1_out, double& U2_out) const
{
    // Vi följer strukturen i rader 273–312:
    // beräknar fyra psi-skalärtermer och kombinerar till derivator.

    double T1, T2;
    double U1, U2, U3, U4;
    double U5 = 0.0, U6 = 0.0;

    // ----- psi(M+1/2, N, N+1) -----
    int P1 = P1_base + 0;     // M+1/2 ligger på P1_base+0.5 i BASIC, vi approximerar här
    int P2 = P2_base;
    int P3 = P2_base + 1;
    int P4 = J2;

    // 283 GOSUB 87
    psiImpedanceKernel(true, I, J, P1, P2, P3, P4, T1, T2);
    U5 = T1;
    U6 = T2;

    // psi(M-1/2, N, N+1)
    P1 = P1_base - 1;
    psiImpedanceKernel(true, I, J, P1, P2, P3, P4, T1, T2);
    U1 = (U5 - T1) / g_.S[J2];
    U2 = (U6 - T2) / g_.S[J2];

    // ----- psi(M+1/2, N-1, N) -----
    if (J1 != J2) {
        P1 = P1_base + 0;
        P2 = P2_base - 1;
        P3 = P2_base;
        P4 = J1;

        psiImpedanceKernel(true, I, J, P1, P2, P3, P4, T1, T2);
        U3 = T1;
        U4 = T2;

        // psi(M-1/2, N-1, N)
        P1 = P1_base - 1;
        psiImpedanceKernel(true, I, J, P1, P2, P3, P4, T1, T2);

        U3 = (U3 - T1) / g_.S[J1];
        U4 = (U4 - T2) / g_.S[J1];

        U1 += U3;
        U2 += U4;
    }

    // slutlig gradΦ-projektion längs trådriktningarna (310–311)
    U1_out = F4 * U1;
    U2_out = F4 * U2;
}

// -------------------------------------------------------------
// Själva impedansmatris-beräkningen (rader 195–336)
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

    // Observation-loop (rad 210–221)
    for (int I = 1; I <= N; ++I) {
        int I1 = std::abs(g_.C[I][0]);
        int I2 = std::abs(g_.C[I][1]);

        double F4_obs = sgn(g_.C[I][0]) * g_.S[I1];
        double F5_obs = sgn(g_.C[I][1]) * g_.S[I2];

        // R(M+1/2)-R(M-1/2) = (T5,T6,T7)
        double T5 = F4_obs * g_.CA[I1] + F5_obs * g_.CA[I2];
        double T6 = F4_obs * g_.CB[I1] + F5_obs * g_.CB[I2];
        double T7 = F4_obs * g_.CG[I1] + F5_obs * g_.CG[I2];
        if (g_.C[I][0] == -g_.C[I][1]) {
            T7 = g_.S[I1] * (g_.CG[I1] + g_.CG[I2]);
        }

        // Källa-loop (rad 221–333)
        for (int J = 1; J <= N; ++J) {
            int J1 = std::abs(g_.C[J][0]);
            int J2 = std::abs(g_.C[J][1]);
            double F4 = sgn(g_.C[J][0]);
            double F5 = sgn(g_.C[J][1]);
            double F6 = 1.0;
            double F7 = 1.0;

            // summering över images (rad 230–332)
            double D1_tot = 0.0;  // vektorpotential (Re)
            double D2_tot = 0.0;  // vektorpotential (Im)
            double U1_tot = 0.0;  // gradΦ (Re)
            double U2_tot = 0.0;  // gradΦ (Im)

            for (int Kimg = 1; Kimg >= -1 && Kimg >= 2 - g_.G; Kimg -= 2) {
                int F8 = 0;

                if (g_.C[J][0] == -g_.C[J][1]) {
                    if (Kimg < 0) break;
                    F6 = F4;
                    F7 = F5;
                }

                // Förenklad variant: vi hoppar över F8-optimeringen och
                // symmetrikontrollerna; Z(I,J) beräknas alltid.
                // (Du kan lägga in raderna 235–246 här om du vill ha fullt stöd.)

                // PSI(M,N,N+1/2) (vektorpotential, GOSUB 102)
                int P1 = 2 * g_.W[I] + I - 1;   // 248
                int P2 = 2 * g_.W[J] + J - 1;   // 249
                int P3 = P2 + 1;     // N+1/2 ~ P2+0.5, men vi jobbar med heltal
                int P4 = J2;                    // 251

                // 252 GOSUB 102
                double T1, T2;
                psiImpedanceKernel(false, I, J, P1, P2, P3, P4, T1, T2);
                double U1 = F5 * T1;            // 253
                double U2 = F5 * T2;            // 254

                // PSI(M,N-1/2,N)
                P3 = P2;                        // 256
                P2 = P2 - 1;                    // 257  ?? P2=P2-.5 ??
                P4 = J1;                        // 258

                // 259 IF F8<2 THEN GOSUB 102
                psiImpedanceKernel(false, I, J, P1, P2, P3, P4, T1, T2);
                double V1 = F4 * T1;
                double V2 = F4 * T2;

                // Vektorpotential-bidrag (rad 262–272)
                double X3 = U1 * g_.CA[J2] + V1 * g_.CA[J1];
                double Y3 = U1 * g_.CB[J2] + V1 * g_.CB[J1];
                double Z3 = (F7 * U1 * g_.CG[J2] + F6 * V1 * g_.CG[J1]) * Kimg;
                double D1 = w2_ * (X3 * T5 + Y3 * T6 + Z3 * T7);

                X3 = U2 * g_.CA[J2] + V2 * g_.CA[J1];
                Y3 = U2 * g_.CB[J2] + V2 * g_.CB[J1];
                Z3 = (F7 * U2 * g_.CG[J2] + F6 * V2 * g_.CG[J1]) * Kimg;
                double D2 = w2_ * (X3 * T5 + Y3 * T6 + Z3 * T7);


                // Grad(Φ)-bidrag (274–312, förenklad)
                double gradU1 = 0.0, gradU2 = 0.0;
                scalarGradientContribution(I, J, J1, J2, F4, F5,
                                           P1, P2 + 1, gradU1, gradU2);

                D1_tot += D1;
                D2_tot += D2;
                U1_tot += gradU1;
                U2_tot += gradU2;
            }

            // Summera in i Z (314–315)
            Z(I, J) += Complex(k_ * (D1_tot + U1_tot),
                               k_ * (D2_tot + U2_tot));
        }

        // Ingen tidsmätning/progress här; lägg till vid behov
    }
}
