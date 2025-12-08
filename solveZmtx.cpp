#include "solveZmtx.hpp"

// Lös A x = b med Gauss-elimination (utan superavancerad pivotering)
// A ändras in-place. Returnerar lösning i x.
std::vector<std::complex<double>> solveLinearSystem(
    std::vector<std::vector<std::complex<double>>> A,
    const std::vector<std::complex<double>>& b)
{
    using Complex = std::complex<double>;
    const std::size_t N = A.size();
    std::vector<Complex> x(b);

    // Framåt-elimination
    for (std::size_t k = 0; k < N; ++k) {
        // Enkel pivot: hitta största element i kolonn k från rad k och nedåt
        std::size_t pivot = k;
        double maxAbs = std::abs(A[k][k]);
        for (std::size_t i = k + 1; i < N; ++i) {
            double val = std::abs(A[i][k]);
            if (val > maxAbs) {
                maxAbs = val;
                pivot = i;
            }
        }
        if (maxAbs < std::numeric_limits<double>::epsilon()) {
            throw std::runtime_error("Singulär matris i solveLinearSystem");
        }

        if (pivot != k) {
            std::swap(A[k], A[pivot]);
            std::swap(x[k], x[pivot]);
        }

        Complex Akk = A[k][k];
        for (std::size_t i = k + 1; i < N; ++i) {
            Complex factor = A[i][k] / Akk;
            A[i][k] = 0.0;
            for (std::size_t j = k + 1; j < N; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            x[i] -= factor * x[k];
        }
    }

    // Bakåtsubstitution
    for (int i = int(N) - 1; i >= 0; --i) {
        Complex sum = x[std::size_t(i)];
        for (std::size_t j = std::size_t(i) + 1; j < N; ++j) {
            sum -= A[std::size_t(i)][j] * x[j];
        }
        x[std::size_t(i)] = sum / A[std::size_t(i)][std::size_t(i)];
    }

    return x;
}

std::complex<double> computeInputImpedance(
    const ImpedanceMatrix& Z,
    int feedSegment)  // 1-baserat index
{
    using Complex = std::complex<double>;
    int N = Z.size();
    if (feedSegment < 1 || feedSegment > N) {
        throw std::runtime_error("feedSegment utanför tillåtet intervall");
    }

    // Bygg A-matris
    std::vector<std::vector<Complex>> A{
        std::size_t(N), std::vector<Complex>(std::size_t(N))};

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) {
            A[std::size_t(i - 1)][std::size_t(j - 1)] = Z(i, j);
        }
    }

    // Exciteringsvektor b: 1 V på feedsegmentet, 0 annars
    std::vector<Complex> b(std::size_t(N), Complex(0.0, 0.0));
    b[std::size_t(feedSegment - 1)] = Complex(1.0, 0.0);

    // Lös A I = b
    std::vector<Complex> I = solveLinearSystem(A, b);

    Complex Ifeed = I[std::size_t(feedSegment - 1)];
    return Complex(1.0, 0.0) / Ifeed;
}
