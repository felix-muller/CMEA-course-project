#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "grad_shape.hpp"
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"
#include "constant_function.hpp"

//----------------compMatrixBegin----------------
//! Evaluate the stiffness matrix on the triangle spanned by
//! the points (a, b, c).
//!
//! Here, the stiffness matrix A is a 3x3 matrix
//! 
//! $$A_{ij} = \int_{K} ( sigma(x,y) \nabla \lambda_i^K(x, y) \cdot  \nabla \lambda_j^K(x, y)\; + r \lambda_i^K(x,y)\lambda_j^K(x,y) )dV$$
//! 
//! where $K$ is the triangle spanned by (a, b, c).
//!
//! @param[out] stiffnessMatrix should be a 3x3 matrix
//!                        At the end, will contain the integrals above.
//!
//! @param[in] a the first corner of the triangle
//! @param[in] b the second corner of the triangle
//! @param[in] c the third corner of the triangle
//! @param[in] sigma the function sigma as in the exercise
//! @param[in] r the parameter r as in the exercise
template<class MatrixType, class Point>
void computeStiffnessMatrix(MatrixType& stiffnessMatrix,
                            const Point& a,
                            const Point& b,
                            const Point& c,
			    const std::function<double(double, double)>& sigma = constantFunction,
                            const double r=0)
{
    Eigen::Matrix2d coordinateTransform = makeCoordinateTransform(b - a, c - a);
    double volumeFactor = std::abs(coordinateTransform.determinant());
    
    Eigen::Matrix2d elementMap = coordinateTransform.inverse().transpose();
// (write your solution here)
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){ // Loop over the rows and columns of the element stiffness matrix
            //integrate.hpp demands a function with two arguments, so we define a lambda function in the loop (Series 5, Execise 2)
            std::function<double(double, double)> integrand = [&](double x, double y) {
                Eigen::Vector2d alpha = coordinateTransform * Eigen::Vector2d(x, y) + Eigen::Vector2d(a(0), a(1)); //mapping local coordinates to global coordinates
                return sigma(alpha(0), alpha(1)) * (elementMap * gradientLambda(i, x, y)).dot(elementMap * gradientLambda(j, x, y)) + r * lambda(i, x, y) * lambda(j, x, y);
            };
            stiffnessMatrix(i, j) = volumeFactor * integrate(integrand); //scaling the integral between local and global coordinates
        }
    }
}
//----------------compMatrixEnd----------------
