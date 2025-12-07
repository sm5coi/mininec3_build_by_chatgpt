#include "mininec.hpp"
#include "solveZmtx.hpp"
#include <iostream>
#include <iomanip>

// Bygg en enkel vertikal tråd längs z-axeln
// med Nseg segment, längd L, radie a.
// G = 1: fri rymd, G = 2: perfekt jordplan (ej utnyttjad i denna version).
Geometry buildSimpleVerticalWire(int Nseg, double L, double a, int G = 1)
{
    Geometry g;
    g.N = Nseg;
    g.G = G;

    int Nnodes = Nseg + 1;
    g.X.assign(Nnodes + 1, 0.0);
    g.Y.assign(Nnodes + 1, 0.0);
    g.Z.assign(Nnodes + 1, 0.0);

    double dz = L / Nseg;
    for (int n = 1; n <= Nnodes; ++n) {
        g.X[n] = 0.0;
        g.Y[n] = 0.0;
        g.Z[n] = (n - 1) * dz;
    }

    g.A.assign(Nseg + 1, a);
    g.S.assign(Nseg + 1, dz);
    g.CA.assign(Nseg + 1, 0.0);
    g.CB.assign(Nseg + 1, 0.0);
    g.CG.assign(Nseg + 1, 1.0); // längs z

    g.C.assign(Nseg + 1, {0,0});
    g.W.assign(Nseg + 1, 1);
    g.J2.assign(2, {0,0}); // wire-index 1..1

    // Puls I = segment I (1..Nseg)
    for (int I = 1; I <= Nseg; ++I) {
        g.C[I][0] = I;   // "nedre" segmentindex
        g.C[I][1] = I;   // "övre" segmentindex (förenklat lika)
        g.W[I]     = 1;  // alla pulser på tråd 1
    }

    g.J2[1][0] = 1;         // första nod på tråd 1
    g.J2[1][1] = Nnodes;    // sista nod

    return g;
}

int main()
{
    // Exempel: 11-segments vertikal tråd, 1 m lång, radie 1 mm
    int    Nseg = 11;
    double L    = 15.0;
    double a    = 0.001;
    double freq = 10e6;   // 10 MHz
    double c0   = 299792458.0;
    double lambda = c0 / freq;
    double k     = 2.0 * std::acos(-1.0) / lambda;
    double srm   = a * 2.0;  // enkel SRM-sättning

    Geometry geom = buildSimpleVerticalWire(Nseg, L, a, 1);

    MininecImpedanceSolver solver(geom, k, srm);
    ImpedanceMatrix Z(Nseg);

    solver.build(Z);

    std::cout << "Z-matris (" << Nseg << "x" << Nseg << ")\n";

    int feedSeg = Nseg / 2 + 1;  // mittsegmentet (t.ex. segment 6 av 11)
    auto Zin = computeInputImpedance(Z, feedSeg);

    std::cout << "\nInimpedans i segment " << feedSeg << ":\n";
    std::cout << "Zin = " << Zin.real()
              << (Zin.imag() >= 0 ? " + j" : " - j")
              << std::abs(Zin.imag()) << " ohm\n";


    std::cout << std::fixed << std::setprecision(3);

    for (int i = 1; i <= Nseg; ++i) {
        for (int j = 1; j <= Nseg; ++j) {
            auto Zij = Z(i,j);
            std::cout << std::setw(10) << Zij.real()
                      << (Zij.imag() >= 0 ? "+" : "")
                      << std::setw(8) << Zij.imag() << "j  ";
        }
        std::cout << "\n";
    }

    return 0;
}
