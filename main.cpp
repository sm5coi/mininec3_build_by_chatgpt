#include "geometry.hpp"
#include "mininec.hpp"
#include "solveZmtx.hpp"
#include <iostream>
#include <iomanip>



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

    Geometry geom = Geometry::buildSimpleVerticalWire(Nseg, L, a, 1);


    // From Frequency Input Line 1131 - 1148

        double F = 299.8;
        double W = 299.8/F;
        double S0 = 0.001*W;
        double M = 4.77783352;
        double SRM = 0.0001*W;
        double P = std::acos(-1.0);  // pi
        double P0 = P/180.0;

        W = 2.0*P/W;
        double W2 = W*W/2.0;



 //   Geometry geom = buildMininecInpAntenna();

    // // Geometry geom;

    // // Lägg till noder
    // geom.X = {0, 0, 0,  0, 5};
    // geom.Y = {0, 0, 0,  0, 0};
    // geom.Z = {0, 0,10, 10,10};

    // // Wire1: vertikal (nod1 → nod2)
    // geom.wires.push_back({1, 2, 0.001, 10});

    // // Wire2: horisontell höger (nod2 → nod3)
    // geom.wires.push_back({2, 3, 0.001, 5});

    // // Wire3: horisontell vänster (nod2 → nod4)
    // geom.wires.push_back({2, 4, 0.001, 5});

    // // Generera segment och geometri
    // geom.buildSegments();

    // Nseg = geom.N;

    MininecImpedanceSolver solver(geom, W, SRM);
    ImpedanceMatrix Z(geom.segmentCount());

    solver.build(Z);

    std::cout << "Z-matris (" << geom.segmentCount() << "x" << geom.segmentCount() << ")\n";

    int feedSeg = geom.segmentCount() / 2 + 1;  // mittsegmentet (t.ex. segment 6 av 11)
    auto Zin = computeInputImpedance(Z, feedSeg);

    std::cout << "\nInimpedans i segment " << feedSeg << ":\n";
    std::cout << "Zin = " << Zin.real()
              << (Zin.imag() >= 0 ? " + j" : " - j")
              << std::abs(Zin.imag()) << " ohm\n";


    std::cout << std::fixed << std::setprecision(3);

    for (int i = 1; i <= geom.segmentCount(); ++i) {
        for (int j = 1; j <= geom.segmentCount(); ++j) {
            auto Zij = Z(i,j);
            std::cout << std::setw(10) << Zij.real()
                      << (Zij.imag() >= 0 ? "+" : "")
                      << std::setw(8) << Zij.imag() << "j  ";
        }
        std::cout << "\n";
    }

    return 0;
}
