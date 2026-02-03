#include <Eigen/Sparse>
#include <iostream>
#include "writer.hpp"
#include <cmath>
#include <Eigen/SparseCholesky>
#include <stdexcept>
#include <iostream>

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Used for filling the sparse matrix.
typedef Eigen::Triplet<double> Triplet;


//! Vector type
typedef Eigen::VectorXd Vector;


//! Our function pointer, typedef'd to make it easier to use
typedef double(*FunctionPointer)(double, double);

//----------------poissonBegin----------------
//! Create the Poisson matrix for 2D finite difference.
//! @param[out] A will be the Poisson matrix (as in the exercise)
//! @param[in] N number of elements in the x-direction
void createPorousMediaMatrix2D(SparseMatrix& A, FunctionPointer sigma, int N, double dx) {
    std::vector<Triplet> triplets;
    A.resize(N*N, N*N);
    triplets.reserve(5*N*N-4*N);

    // Fill up triplets

// (write your solution here)
for (int j = 1; j <= N; ++j) {
    for (int i = 1; i <= N; ++i) {
        int index = (j - 1) * N + i - 1;

        triplets.push_back(Triplet(index, index, sigma((i - 0.5) * dx, (j) * dx) + sigma((i + 0.5) * dx, (j) * dx) + sigma((i) * dx, (j - 0.5) * dx) + sigma((i) * dx, (j + 0.5) * dx)));

        if (i > 1) {
            triplets.push_back(Triplet(index, index - 1, -sigma((i - 0.5) * dx, (j) * dx)));
        }

        if (i < N) {
            triplets.push_back(Triplet(index, index + 1, -sigma((i + 0.5) * dx, (j) * dx)));
        }

        if (j > 1) {
            triplets.push_back(Triplet(index, index - N, -sigma((i) * dx, (j - 0.5) * dx)));
        }

        if (j < N) {
            triplets.push_back(Triplet(index, index + N, -sigma((i) * dx, (j + 0.5) * dx)));
        }
    }
}

A.setFromTriplets(triplets.begin(), triplets.end());
}
//----------------poissonEnd----------------

//----------------RHSBegin----------------
//! Create the Right hand side for the 2D finite difference
//! @param[out] rhs will at the end contain the right hand side
//! @param[in] f the right hand side function f
//! @param[in] N the number of points in the x direction
//! @param[in] dx the cell width
void createRHS(Vector& rhs, FunctionPointer f, int N, double dx) {
    rhs.resize(N * N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            const double x = (i + 1) * dx;
            const double y = (j + 1) * dx;
            rhs[j * N + i] = dx * dx * f(x, y);
        }
    }
}
//----------------RHSEnd----------------


//----------------solveBegin----------------
//! Solve the Poisson equation in 2D
//! @param[in] f the function pointer to f
//! @param[in] N the number of points to use (in x direction)
//!
//! @returns the approximate solution u
Vector poissonSolve(FunctionPointer f, FunctionPointer sigma, int N) {
    Vector u;
    double dx = 1.0 / (N + 1);

    // Compute A, rhs and u
// (write your solution here)
    // this is exactly like in series 3, problem 2
    SparseMatrix A;
    createPorousMediaMatrix2D(A, sigma, N, dx);

    Vector rhs;
    createRHS(rhs, f, N, dx);

    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Could not decompose matrix");
    }

    u.resize((N + 2) * (N + 2));
    u.setZero();
    Vector innerU = solver.solve(rhs);

    for (int i = 1; i < N + 1; ++i) {
        for (int j = 1; j < N + 1; ++j) {
            u[i * (N + 2) + j] = innerU[(i - 1) * N + j - 1];
        }
    }

    return u;
}
//----------------solveEnd----------------


//! Gives the Right hand side F at the point x, y
double F(double x, double y) {
  return 4*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*(4*cos(2*M_PI*x)*cos(2*M_PI*y)+M_PI);
}

//----------------convergenceBegin----------------
//! Gives the exact solution at the point x,y
double exactSolution(double x, double y) {
// (write your solution here)
    return sin(2*M_PI*x)*sin(2*M_PI*y);
}



void convergeStudy(FunctionPointer F, FunctionPointer sigma) {
    const int startK = 3;
    const int endK = 9;


    Vector errors(endK - startK);
    Vector resolutions(errors.size());
    for (int k = startK; k < endK; ++k) {
        const int N = 1<<k; // 2^k

// (write your solution here)
        auto u = poissonSolve(F, sigma, N);  // Solve the Poisson equation

        double maxError = 0;  // Initialize maximum error for this resolution

        // Loop over grid points (including boundary points)
        for (int j = 0; j < N + 2; ++j) {
            for (int i = 0; i < N + 2; ++i) {
                // Compute normalized coordinates in [0, 1]
                double x = i / (N + 1.0);
                double y = j / (N + 1.0);

                // Evaluate the exact solution at the current position
                double exact = exactSolution(x, y);

                // Calculate the absolute error and update maxError if needed
                double error = std::abs(u(j * (N + 2) + i) - exact);
                maxError = std::max(maxError, error);
            }
        }
    
        // Store the maximum error and resolution for this iteration
        errors[k - startK] = maxError;
        resolutions[k - startK] = N;
    }

    writeToFile("errors.txt", errors);
    writeToFile("resolutions.txt", resolutions);
}
//----------------convergenceEnd----------------

int main(int, char**) {
    auto sigmaCos = [](double x, double y) {
      return M_PI/2.+cos(2*M_PI*x)*cos(2*M_PI*y);
    };

    auto sigmaConstant = [](double x, double y) {
        return 1.0;
    };
    auto u = poissonSolve(F, sigmaCos, 500);
    writeToFile("u_fd.txt", u);

    convergeStudy(F, sigmaCos);
}




