#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

#include "stiffness_matrix.hpp"

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Used for filling the sparse matrix.
typedef Eigen::Triplet<double> Triplet;

//----------------AssembleMatrixBegin----------------
//! Assemble the stiffness matrix
//! for the linear system
//!
//! @param[out] A will at the end contain the Galerkin matrix
//! @param[in] vertices a list of triangle vertices
//! @param[in] triangles a list of triangles
//! @param[in] sigma the function sigma as in the exercise
//! @param[in] r the parameter r as in the exercise
template<class Matrix>
void assembleStiffnessMatrix(Matrix& A, const Eigen::MatrixXd& vertices,
			     const Eigen::MatrixXi& triangles,
			     const std::function<double(double, double)>& sigma = constantFunction,
			     double r=0)
{
    
    const int numberOfElements = triangles.rows();
    A.resize(vertices.rows(), vertices.rows());
    
    std::vector<Triplet> triplets;

    triplets.reserve(numberOfElements * 3 * 3);
// (write your solution here)

// very close to Series 5, Exercise 1
for (int i = 0; i < numberOfElements; i++){ // Loop over each element in the mesh

    auto& triangle = triangles.row(i);    // Get the indices of the vertices that form the current triangle

    // Get the coordinates of the vertices that form the current triangle
    const auto& a = vertices.row(triangle(0));
    const auto& b = vertices.row(triangle(1));
    const auto& c = vertices.row(triangle(2));

    Eigen::Matrix3d elementStiffnessMatrix; // Create a 3x3 matrix to store the element stiffness matrix

    computeStiffnessMatrix(elementStiffnessMatrix, a, b, c, sigma, r);  // Compute the element stiffness matrix using the coordinates of the vertices and other parameters

    for (int j = 0; j < 3; j++){    // Loop over the rows of the element stiffness matrix

        for (int k = 0; k < 3; k++){    // Loop over the columns of the element stiffness matrix

            // Add a triplet to the list of non-zero entries in the global stiffness matrix.
            // The triplet consists of the row index, column index, and value of the element stiffness matrix.
            triplets.push_back(Triplet(triangle(j), triangle(k), elementStiffnessMatrix(j, k)));
        }
    }
}

A.setFromTriplets(triplets.begin(), triplets.end());
}
//----------------AssembleMatrixEnd----------------
