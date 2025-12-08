#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include <array>
// #include <complex>
// #include <stdexcept>
// #include <algorithm>
#include <cmath>

// -------------------------------------------------------------
// 1-baserad geometri (som i MININEC)
// Alla vektorer ska ha storlek >= maxIndex+1, index 0 ignoreras.
// -------------------------------------------------------------
struct Wire {
    int startNode;    // node index
    int endNode;      // node index
    double radius;    // wire radius
    int segments;     // number of segments
};

struct Geometry {
    int N = 0;  // total number of segments
    int G; // used in buildSimpleVerticalWire()

    std::vector<double> X, Y, Z;         // nodes
    std::vector<Wire> wires;             // wires

    // Per-segment data:
    std::vector<std::array<int,2>> C;    // lower/higher node
    std::vector<int> W;                  // segment→wire index
    std::vector<double> A;               // radius per segment
    std::vector<double> S;               // length per segment
    std::vector<double> CA, CB, CG;      // direction cosines

    // Per-wire:
    std::vector<std::array<int,2>> J2;   // first/last node of wire

    void buildSegments();

    Geometry buildSimpleVerticalWire(int Nseg, double L, double a, int G = 1);
};

#endif // GEOMETRY_HPP
