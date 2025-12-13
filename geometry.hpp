#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include <array>

struct Wire {
    int startNode;
    int endNode;
    double radius;
    int segments;
};

class Geometry
{
public:
    // ===== Factory / builders =====
    static Geometry buildSimpleVerticalWire(
        int Nseg,
        double L,
        double a,
        int G = 1
        );

    // ===== Read-only API =====
    int segmentCount() const { return N_; }

    // ===== Data access (const) =====
    const std::vector<double>& X()  const { return X_; }
    const std::vector<double>& Y()  const { return Y_; }
    const std::vector<double>& Z()  const { return Z_; }

    const std::vector<Wire>& wires() const { return wires_; }

    const std::vector<std::array<int,2>>& C()  const { return C_; }
    const std::vector<int>& W()   const { return W_; }
    const std::vector<double>& A() const { return A_; }
    const std::vector<double>& S() const { return S_; }
    const std::vector<double>& CA() const { return CA_; }
    const std::vector<double>& CB() const { return CB_; }
    const std::vector<double>& CG() const { return CG_; }

    const std::vector<std::array<int,2>>& J2() const { return J2_; }

private:
    Geometry() = default;

    void buildSegments();

    // ===== Internal data =====
    int N_ = 0;
    int G_ = 0;

    std::vector<double> X_, Y_, Z_;
    std::vector<Wire> wires_;

    std::vector<std::array<int,2>> C_;
    std::vector<int> W_;
    std::vector<double> A_;
    std::vector<double> S_;
    std::vector<double> CA_, CB_, CG_;

    std::vector<std::array<int,2>> J2_;
};

#endif // GEOMETRY_HPP
