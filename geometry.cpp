#include "geometry.hpp"
#include <cmath>

Geometry Geometry::buildSimpleVerticalWire(int Nseg, double L, double a, int G)
{
    Geometry geom;
    geom.G_ = G;

    // 1-baserad nodlista
    geom.X_ = {0.0, 0.0};
    geom.Y_ = {0.0, 0.0};
    geom.Z_ = {0.0, L};

    geom.wires_.push_back({
        1, 2, a, Nseg
    });

    geom.buildSegments();
    return geom;
}

void Geometry::buildSegments()
{
    N_ = 0;

    for (size_t w = 0; w < wires_.size(); ++w) {
        const Wire& wire = wires_[w];

        int n1 = wire.startNode;
        int n2 = wire.endNode;
        int ns = wire.segments;

        double dx = (X_[n2] - X_[n1]) / ns;
        double dy = (Y_[n2] - Y_[n1]) / ns;
        double dz = (Z_[n2] - Z_[n1]) / ns;

        J2_.push_back({N_ + 1, N_ + ns});

        for (int i = 0; i < ns; ++i) {
            ++N_;

            C_.push_back({n1, n2});
            W_.push_back(static_cast<int>(w) + 1);
            A_.push_back(wire.radius);

            double len = std::sqrt(dx*dx + dy*dy + dz*dz);
            S_.push_back(len);

            CA_.push_back(dx / len);
            CB_.push_back(dy / len);
            CG_.push_back(dz / len);
        }
    }
}
