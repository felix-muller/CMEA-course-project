#include "create_poisson_matrix.hpp"

//! Used for filling the sparse matrix.
using Triplet = Eigen::Triplet<double>;

//! Create the 1D Poisson matrix
//! @param[in] N the number of interior points
//! @param[in] a the coefficient function a
//!
//! @returns the Poisson matrix.
SparseMatrix createPoissonMatrix(int N,
    const std::function<double(double)>& a) {

    SparseMatrix A;
// (write your solution here)
    A.resize(N, N);
    std::vector<Triplet> triplets;
    triplets.reserve(N + 2 * N - 2);

    double h = 1.0 / (N + 1);

    for (int i = 1; i <= N; ++i) {
        double x_i = i * h;

        triplets.push_back(Triplet(i - 1, i - 1, 2.0 * a(x_i) / (h * h)));

        if (i > 1) {
            triplets.push_back(Triplet(i - 1, i - 2, -a(x_i) / (h * h)));
        }

        if (i < N) {
            triplets.push_back(Triplet(i - 1, i, -a(x_i) / (h * h)));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());

    return A;
}

