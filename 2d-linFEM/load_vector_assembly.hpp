#pragma once
#include <Eigen/Core>
#include "load_vector.hpp"

//----------------AssembleVectorBegin----------------
//! Assemble the load vector into the full right hand side
//! for the linear system
//!
//! @param[out] F will at the end contain the RHS values for each vertex.
//! @param[in] vertices a list of triangle vertices
//! @param[in] triangles a list of triangles
//! @param[in] f the RHS function f.
void assembleLoadVector(Eigen::VectorXd& F,
                           const Eigen::MatrixXd& vertices,
                           const Eigen::MatrixXi& triangles,
                           const std::function<double(double, double)>& f)
{
     const int numberOfElements = triangles.rows();

     F.resize(vertices.rows());
     F.setZero();
// (write your solution here)

// very close to Series 5, Exercise 1
for (int i = 0; i < numberOfElements; i++){ // Loop over each element in the mesh

    auto& triangle = triangles.row(i);    // Get the indices of the vertices that form the current triangle

    // Get the coordinates of the vertices that form the current triangle
    const auto& a = vertices.row(triangle(0));
    const auto& b = vertices.row(triangle(1));
    const auto& c = vertices.row(triangle(2));

    Eigen::Vector3d elementLoadVector; // Create a 3x1 vector to store the element load vector

    computeLoadVector(elementLoadVector, a, b, c, f);  // Compute the element load vector using the coordinates of the vertices and other parameters

    for (int j = 0; j < 3; j++){    // Loop over the rows of the element load vector

        // Add the element load vector to the global load vector.
        F(triangle(j)) += elementLoadVector(j);
    }
}
}
//----------------AssembleVectorEnd----------------
