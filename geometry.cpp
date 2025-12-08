#include "geometry.hpp"


void Geometry::buildSegments()
{
    N = 0;
    C.clear();
    W.clear();
    A.clear();
    S.clear();
    CA.clear();
    CB.clear();
    CG.clear();
    J2.clear();

    J2.resize(wires.size()+1);

    for (int w = 1; w <= (int)wires.size(); ++w) {
        const Wire& Ww = wires[w-1];
        int n1 = Ww.startNode;
        int n2 = Ww.endNode;

        J2[w][0] = n1;
        J2[w][1] = n2;

        double dx = X[n2] - X[n1];
        double dy = Y[n2] - Y[n1];
        double dz = Z[n2] - Z[n1];

        double L = std::sqrt(dx*dx + dy*dy + dz*dz);
        double segLen = L / Ww.segments;

        double ca = dx / L;
        double cb = dy / L;
        double cg = dz / L;

        for (int s = 0; s < Ww.segments; ++s) {
            ++N;
            W.push_back(w);

            int nodeLow  = N * 2 - 1;     // eller gör separat nodtabell
            int nodeHigh = N * 2;

            C.push_back({nodeLow, nodeHigh});

            A.push_back(Ww.radius);
            S.push_back(segLen);
            CA.push_back(ca);
            CB.push_back(cb);
            CG.push_back(cg);
        }
    }
}

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
