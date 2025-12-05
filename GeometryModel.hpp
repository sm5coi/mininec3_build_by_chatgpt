#ifndef GEOMETRY_MODEL_HPP
#define GEOMETRY_MODEL_HPP

#include <vector>

struct WireEnd {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Wire {
    int id = 0;
    int segments = 0;
    double radius = 0.0;

    WireEnd end1;
    WireEnd end2;

    double length = 0.0;
    double segLength = 0.0;

    double dirX = 0.0;
    double dirY = 0.0;
    double dirZ = 0.0;

    int conn1 = 0;
    int conn2 = 0;
};

struct Pulse {
    int id = 0;
    int wireId = 0;
    int connStart = 0;
    int connEnd = 0;
};

struct Breakpoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

class GeometryModel {
public:
    double ground = 0.0;

    std::vector<Wire> wires;
    std::vector<Pulse> pulses;
    std::vector<Breakpoint> points;

    GeometryModel() = default;

    void computeWireGeometry(Wire& w);
    void computeConnections();
    void buildPulsesAndPoints();
};

GeometryModel readGeometry();
void printGeometry(const GeometryModel& model);

#endif
