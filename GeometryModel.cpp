#include "GeometryModel.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;


// ------------------------------------------------------------
// Beräkna riktning & längder för en wire
// ------------------------------------------------------------
void GeometryModel::computeWireGeometry(Wire& w)
{
    double dx = w.end2.x - w.end1.x;
    double dy = w.end2.y - w.end1.y;
    double dz = w.end2.z - w.end1.z;

    w.length = std::sqrt(dx*dx + dy*dy + dz*dz);
    w.segLength = w.length / w.segments;

    w.dirX = dx / w.length;
    w.dirY = dy / w.length;
    w.dirZ = dz / w.length;
}

// ------------------------------------------------------------
// Beräkna anslutningar mellan wires
// ------------------------------------------------------------
void GeometryModel::computeConnections()
{
    for (auto& w : wires)
    {
        // Ground connection
        if (w.end1.z == 0) w.conn1 = -w.id;
        if (w.end2.z == 0) w.conn2 = -w.id;

        // Testa mot alla tidigare wires
        for (auto& other : wires)
        {
            if (other.id == w.id) continue;

            auto same = [](const WireEnd& a, const WireEnd& b) {
                return a.x == b.x && a.y == b.y && a.z == b.z;
            };

            // END1—END1
            if (same(w.end1, other.end1)) w.conn1 = -other.id;
            // END1—END2
            if (same(w.end1, other.end2)) w.conn1 = other.id;

            // END2—END2
            if (same(w.end2, other.end2)) w.conn2 = -other.id;
            // END2—END1
            if (same(w.end2, other.end1)) w.conn2 = other.id;
        }
    }
}

// ------------------------------------------------------------
// Skapa pulser och brytpunkter längs varje wire
// ------------------------------------------------------------
void GeometryModel::buildPulsesAndPoints()
{
    pulses.clear();
    points.clear();

    for (auto& w : wires)
    {
        for (int s = 0; s < w.segments; ++s)
        {
            Pulse p;
            p.id = pulses.size() + 1;
            p.wireId = w.id;

            // Edge connections
            p.connStart = (s == 0) ? w.conn1 : w.id;
            p.connEnd   = (s == w.segments - 1) ? w.conn2 : w.id;

            pulses.push_back(p);

            // Breakpoint
            double t = double(s) / w.segments;

            Breakpoint bp;
            bp.x = w.end1.x + t * (w.end2.x - w.end1.x);
            bp.y = w.end1.y + t * (w.end2.y - w.end1.y);
            bp.z = w.end1.z + t * (w.end2.z - w.end1.z);

            points.push_back(bp);
        }
    }
}

// ------------------------------------------------------------
// Interaktiv läsning av geometri
// ------------------------------------------------------------
GeometryModel readGeometry()
{
    GeometryModel model;
    std::ifstream fid;
    std::string line;
    int NW;

    fid.open("MININEC.INP");

    if (!fid)
    {
        cout << "Could not open infile" << endl;
        return model;
    }
    else
    {
        cout << "    Open MININEC.INP" << endl << endl;
    }
    getline(fid, line);
    istringstream iss(line);
    iss >> NW;



    std::cout << "Number of wires: " << NW << endl;
    //std::cin >> NW;

    int i = 0;
    while (getline(fid, line))
    {
        i++;
        Wire w;
        w.id = i;

        istringstream iss(line);

        iss >> w.segments;
        iss >> w.end1.x >> w.end1.y >> w.end1.z;
        iss >> w.end2.x >> w.end2.y >> w.end2.z;
        iss >> w.radius;

        model.wires.push_back(w);

        //vector<double> row;
        //double value;

        // while (iss >> value)
        //     row.push_back(value);

        // if (!row.empty())
        //     Minp.push_back(row);
    }

    /*
    for (int i = 1; i <= NW; ++i)
    {
        Wire w;
        w.id = i;

        std::cout << "\nWIRE " << i << "\n";

        std::cout << "  Number of segments: ";
        std::cin >> w.segments;

        std::cout << "  End1 (x y z): ";
        std::cin >> w.end1.x >> w.end1.y >> w.end1.z;

        std::cout << "  End2 (x y z): ";
        std::cin >> w.end2.x >> w.end2.y >> w.end2.z;

        std::cout << "  Radius: ";
        std::cin >> w.radius;

        model.wires.push_back(w);
    }
*/

    // Beräkna geometri
    for (auto& w : model.wires)
        model.computeWireGeometry(w);

    // Anslutningar
    model.computeConnections();

    // Puls- och punkt-listor
    model.buildPulsesAndPoints();

    fid.close();

    return model;
}

// ------------------------------------------------------------
// Skriv ut geometridump
// ------------------------------------------------------------
void printGeometry(const GeometryModel& model)
{
    std::cout << "\n===== GEOMETRY =====\n";

    for (auto& w : model.wires)
    {
        std::cout << "WIRE " << w.id << "\n";
        std::cout << "  End1: (" << w.end1.x << ", " << w.end1.y << ", " << w.end1.z << ")\n";
        std::cout << "  End2: (" << w.end2.x << ", " << w.end2.y << ", " << w.end2.z << ")\n";
        std::cout << "  Segments: " << w.segments << "   Radius: " << w.radius << "\n";
        std::cout << "  Connections: " << w.conn1 << "  " << w.conn2 << "\n\n";
    }

    std::cout << "Total pulses: " << model.pulses.size() << "\n";
    std::cout << "Total points: " << model.points.size() << "\n";
}
