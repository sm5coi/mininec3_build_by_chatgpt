#include "mininec.hpp"
#include "GeometryModel.hpp"
#include <iostream>

struct WireSpec {
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    double radius;
};

/*
void buildSimpleGeometry(const std::vector<WireSpec>& wires, Geometry& geom) {
    int NW = static_cast<int>(wires.size());

    // Räkna totala segment
    int totalSegments = 0;
    for (const auto& w : wires) {
        int nodes = static_cast<int>(w.X.size());
        totalSegments += (nodes - 1);
    }

    // Allokera 1-baserat
    geom.N  = totalSegments;
    geom.G  = 1;

    geom.X.assign(1, 0.0); // index 0 oanvänd
    geom.Y.assign(1, 0.0);
    geom.Z.assign(1, 0.0);

    geom.A.assign(totalSegments + 1, 0.0);
    geom.CA.assign(totalSegments + 1, 0.0);
    geom.CB.assign(totalSegments + 1, 0.0);
    geom.CG.assign(totalSegments + 1, 0.0);
    geom.S.assign(totalSegments + 1, 0.0);

    geom.C.assign(geom.N + 1, {0,0});
    geom.W.assign(geom.N + 1, 0);
    geom.J2.assign(NW + 1, {0,0});

    int seg = 0;
    int p   = 0;

    for (int wIndex = 1; wIndex <= NW; ++wIndex) {
        const auto& w = wires[wIndex - 1];
        int nodes = static_cast<int>(w.X.size());

        int firstSeg = seg + 1;

        for (int s = 0; s < nodes - 1; ++s) {
            seg++;
            p++;

            // segmentvektor
            double dx = w.X[s+1] - w.X[s];
            double dy = w.Y[s+1] - w.Y[s];
            double dz = w.Z[s+1] - w.Z[s];
            double L  = std::sqrt(dx*dx + dy*dy + dz*dz);

            geom.S[seg]  = L;
            geom.CA[seg] = (L > 0) ? dx / L : 0.0;
            geom.CB[seg] = (L > 0) ? dy / L : 0.0;
            geom.CG[seg] = (L > 0) ? dz / L : 0.0;
            geom.A[seg]  = w.radius;

            geom.W[p]    = wIndex;
            geom.C[p][0] = seg;
            geom.C[p][1] = seg;
        }

        int lastSeg = seg;
        geom.J2[wIndex][0] = firstSeg;
        geom.J2[wIndex][1] = lastSeg;
    }
}

*/

void buildSimpleGeometry(const std::vector<WireSpec>& wires, Geometry& geom) {
    int NW = static_cast<int>(wires.size());

    // 1) Räkna total antal noder och segment
    int totalNodes = 0;
    int totalSegments = 0;

    for (const auto& w : wires) {
        int nodes = static_cast<int>(w.X.size());
        totalNodes += nodes;

        totalSegments += (nodes - 1);
    }

    // 2) Allokera 1-baserade vektorer
    geom.N  = totalSegments;   // 1 puls per segment
    geom.G  = 1;

    geom.X.assign(totalNodes + 1, 0.0);
    geom.Y.assign(totalNodes + 1, 0.0);
    geom.Z.assign(totalNodes + 1, 0.0);

    geom.A.assign(totalSegments + 1, 0.0);
    geom.CA.assign(totalSegments + 1, 0.0);
    geom.CB.assign(totalSegments + 1, 0.0);
    geom.CG.assign(totalSegments + 1, 0.0);
    geom.S.assign(totalSegments + 1, 0.0);

    geom.C.assign(geom.N + 1, {0,0});
    geom.W.assign(geom.N + 1, 0);
    geom.J2.assign(NW + 1, {0,0});

    // 3) Kopiera alla nodkoordinater till globala listan geom.X/Y/Z
    int nodeIndex = 0;
    for (const auto& w : wires) {
        for (size_t i = 0; i < w.X.size(); ++i) {
            nodeIndex++;
            geom.X[nodeIndex] = w.X[i];
            geom.Y[nodeIndex] = w.Y[i];
            geom.Z[nodeIndex] = w.Z[i];
        }
    }

    // 4) Skapa segment + pulser
    int seg = 0;
    int p   = 0;
    nodeIndex = 0;

    for (int wIndex = 1; wIndex <= NW; ++wIndex) {
        const auto& w = wires[wIndex - 1];
        int nodes = static_cast<int>(w.X.size());

        int firstNode = nodeIndex + 1;
        int firstSeg  = seg + 1;

        // Segment per tråd
        for (int s = 0; s < nodes - 1; ++s) {
            int globalNode1 = firstNode + s;
            int globalNode2 = globalNode1 + 1;

            seg++;
            p++;

            // Beräkna segmentgeometri
            double dx = geom.X[globalNode2] - geom.X[globalNode1];
            double dy = geom.Y[globalNode2] - geom.Y[globalNode1];
            double dz = geom.Z[globalNode2] - geom.Z[globalNode1];
            double L  = std::sqrt(dx*dx + dy*dy + dz*dz);

            geom.S[seg]  = L;
            geom.CA[seg] = (L > 0) ? dx / L : 0.0;
            geom.CB[seg] = (L > 0) ? dy / L : 0.0;
            geom.CG[seg] = (L > 0) ? dz / L : 0.0;
            geom.A[seg]  = w.radius;

            // Puls = segment-index
            geom.W[p]    = wIndex;
            geom.C[p][0] = seg;
            geom.C[p][1] = seg;
        }

        int lastSeg = seg;
        geom.J2[wIndex][0] = firstSeg;
        geom.J2[wIndex][1] = lastSeg;

        nodeIndex += nodes;
    }
}



/*
int main() {
    std::vector<WireSpec> wires;

    //GeometryModel model = readGeometry();

    // Geometri från mininec.doc, L-antenna
    WireSpec w1;
    // Wire 1: lodrät längs z
    w1.X = {0.0, 0.0, 0.0, 0.0, 0.0};
    w1.Y = {0.0, 0.0, 0.0, 0.0, 0.0};
    w1.Z = {0.0, 0.04775, 0.0955, 0.14325, 0.191};
    w1.radius = 0.004;
    wires.push_back(w1);

    WireSpec w2;
    // Wire 2: vågrätt längs y
    w2.X = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    w2.Y = {0.0, 0.0515, 0.103, 0.1545, 0.206, 0.2575, 0.309};
    w2.Z = {0.191, 0.191, 0.191, 0.191, 0.191, 0.191, 0.191};
    w2.radius = 0.004;
    wires.push_back(w2);


    Geometry geom;
    buildSimpleGeometry(wires, geom);

    double lambda = 1.0;
    double k      = 2.0 * 3.141592653589793 / lambda;
    double srm    = 0.01;

    MininecImpedanceSolver solver(geom, k, srm);

    ImpedanceMatrix Z(geom.N);
    solver.build(Z);

    for (int i = 1; i <= geom.N; ++i) {
        for (int j = 1; j <= geom.N; ++j) {
            std::cout << "Z(" << i << "," << j << ") = " << Z(i,j) << "\n";
        }
    }
}
*/

void printGeometry(const Geometry& geom) {
    std::cout << "\n=== Noder (globala) ===\n";
    for (size_t i = 1; i < geom.X.size(); ++i) {
        std::cout << "Node " << i
                  << ": X=" << geom.X[i]
                  << "  Y=" << geom.Y[i]
                  << "  Z=" << geom.Z[i] << "\n";
    }

    std::cout << "\n=== Segment ===\n";
    for (size_t s = 1; s < geom.S.size(); ++s) {
        std::cout << "Seg " << s
                  << ": L=" << geom.S[s]
                  << " CA=" << geom.CA[s]
                  << " CB=" << geom.CB[s]
                  << " CG=" << geom.CG[s]
                  << "  (radie=" << geom.A[s] << ")\n";
    }

    std::cout << "\n=== Pulser (unknowns) ===\n";
    for (int p = 1; p <= geom.N; ++p) {
        std::cout << "Puls " << p
                  << ": C1=" << geom.C[p][0]
                  << " C2=" << geom.C[p][1]
                  << "  wire=" << geom.W[p]
                  << "\n";
    }

    std::cout << "\n=== J2 per tråd (första/sista segment) ===\n";
    for (size_t w = 1; w < geom.J2.size(); ++w) {
        std::cout << "Wire " << w
                  << ": first=" << geom.J2[w][0]
                  << " last="  << geom.J2[w][1]
                  << "\n";
    }
}


int main() {
    std::vector<WireSpec> wires;

/*
    // Wire 1: lodrät tråd (3 segment, 4 noder)
    // z = -0.5, -0.1667, 0.1667, +0.5
    WireSpec w1;
    w1.X = {0.0, 0.0, 0.0, 0.0};
    w1.Y = {0.0, 0.0, 0.0, 0.0};
    w1.Z = {-0.5, -0.1667, 0.1667, 0.5};
    w1.radius = 0.001;
    wires.push_back(w1);

    // Wire 2: lutande tråd (2 segment, 3 noder)
    WireSpec w2;
    w2.X = {1.0, 1.5, 2.0};
    w2.Y = {0.0, 0.2, 0.4};
    w2.Z = {0.0, 0.0, 0.0};
    w2.radius = 0.002;
    wires.push_back(w2);
*/

    // Geometri från mininec.doc, L-antenna
    WireSpec w1;
    // Wire 1: lodrät längs z
    w1.X = {0.0, 0.0, 0.0, 0.0, 0.0};
    w1.Y = {0.0, 0.0, 0.0, 0.0, 0.0};
    w1.Z = {0.0, 0.04775, 0.0955, 0.14325, 0.191};
    w1.radius = 0.004;
    wires.push_back(w1);

    WireSpec w2;
    // Wire 2: vågrätt längs y
    w2.X = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    w2.Y = {0.0, 0.0515, 0.103, 0.1545, 0.206, 0.2575, 0.309};
    w2.Z = {0.191, 0.191, 0.191, 0.191, 0.191, 0.191, 0.191};
    w2.radius = 0.004;
    wires.push_back(w2);

    Geometry geom;
    buildSimpleGeometry(wires, geom);

    printGeometry(geom);
    double lambda = 1.0;
    double k      = 2.0 * 3.141592653589793 / lambda;
    double srm    = 0.01;

    MininecImpedanceSolver solver(geom, k, srm);

    ImpedanceMatrix Z(geom.N);
    solver.build(Z);

    for (int i = 1; i <= geom.N; ++i) {
        for (int j = 1; j <= geom.N; ++j) {
            std::cout << "Z(" << i << "," << j << ") = " << Z(i,j) << "\n";
        }
    }

    return 0;
}
