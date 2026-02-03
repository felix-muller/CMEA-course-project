#include "crank_nicolson.hpp"


//! Uses the Crank-Nicolson method to compute u from time 0 to time T
//!
//! @param[in] u0 the initial data, as column vector (size N+2)
//! @param[in] dt the time step size
//! @param[in] T the final time at which to compute the solution (which we assume to be a multiple of dt)
//! @param[in] N the number of interior grid points
//! @param[in] gL function of time with the Dirichlet condition at left boundary
//! @param[in] gR function of time with the Dirichlet condition at right boundary
//! @param[in] a the coefficient function a
//!
//! @note the vector returned should include the boundary values!
//!
std::pair<Eigen::MatrixXd, Eigen::VectorXd> crankNicolson(
    const Eigen::VectorXd& u0,
    double dt, double T, int N,
    const std::function<double(double)>& gL,
    const std::function<double(double)>& gR,
    const std::function<double(double)>& a) {

    Eigen::VectorXd time;
    Eigen::MatrixXd u;


// (write your solution here)
    // Calculate the number of time steps
    const int nsteps = int(round(T / dt));
    // Calculate the grid spacing
    const double h = 1.0 / (N + 1);

    // Resize the solution matrix and time vector
    u.resize(N + 2, nsteps + 1);
    time.resize(nsteps + 1);

    /* Initialize A (Poisson matrix) */
    auto A = createPoissonMatrix(N, a);

    // Initialize B matrix (Identity + dt/2 * A)
    SparseMatrix B(N, N);
    B.setIdentity();
    B += dt / 2.0 * A;

    /* Initialize u */
    u.col(0) = u0;
    // Initialize time
    time[0] = 0.0;

    /* Initialize solver and compute LU decomposition of B (Note: since dt is constant, the matrix B is the same for all timesteps) */
    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(B);

    // Vectors for storing boundary condition values
    Eigen::VectorXd G1 = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd G2 = Eigen::VectorXd::Zero(N);

    // Time-stepping loop
    for (int k = 0; k < nsteps; k++) {
        // Update time
        time[k + 1] = (k + 1) * dt;

        // Update boundary conditions at current and next time steps
        G1[0] = a(h) * dt * gL(time[k]) / (h * h);
        G1[N - 1] = a(1 - h) * dt * gR(time[k]) / (h * h);
        G2[0] = a(h) * dt * gL(time[k + 1]) / (h * h);
        G2[N - 1] = a(1 - h) * dt * gR(time[k + 1]) / (h * h);

        // Compute right-hand side (rhs) of the linear system
        Eigen::VectorXd rhs =
            u.block(1, k, N, 1) - dt / 2 * A * u.block(1, k, N, 1) + 0.5 * (G1 + G2);

        // Solve the linear system using the SparseLU solver
        u.block(1, k + 1, N, 1) = solver.solve(rhs);

        // Set boundary values
        u(0, k + 1) = gL(time[k + 1]);
        u(N + 1, k + 1) = gR(time[k + 1]);
    }
    
    return std::make_pair(u, time);

}
