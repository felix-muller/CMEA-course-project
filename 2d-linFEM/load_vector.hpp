#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include "coordinate_transform.hpp"
#include "integrate.hpp"
#include "shape.hpp"
#include <functional>

//----------------compVectorBegin----------------
//! Evaluate the load vector on the triangle spanned by
//! the points (a, b, c).
//!
//! Here, the load vector is a vector $(v_i)$ of
//! three components, where 
//! 
//! $$v_i = \int_{K} \lambda_i^K(x, y) f(x, y) \; dV$$
//! 
//! where $K$ is the triangle spanned by (a, b, c).
//!
//! @param[out] loadVector should be a vector of length 3. 
//!                        At the end, will contain the integrals above.
//!
//! @param[in] a the first corner of the triangle
//! @param[in] b the second corner of the triangle
//! @param[in] c the third corner of the triangle
//! @param[in] f the function f (LHS).
template<class Vector, class Point>
void computeLoadVector(Vector& loadVector,
                    const Point& a, const Point& b, const Point& c,
                    const std::function<double(double, double)>& f)
{
    Eigen::Matrix2d coordinateTransform = makeCoordinateTransform(b - a, c - a);
    double volumeFactor = std::abs(coordinateTransform.determinant());
// (write your solution here)
    for (int i = 0; i < 3; i++){
        //integrate.hpp demands a function with two arguments, so we define a lambda function in the loop (Series 5, Exercise 1)
        std::function<double(double, double)> integrand = [&](double x, double y) {
            Eigen::Vector2d alpha = coordinateTransform * Eigen::Vector2d(x, y) + Eigen::Vector2d(a(0), a(1)); //mapping local coordinates to global coordinates (Series 5, Exercise 1)
            return f(alpha(0), alpha(1)) * lambda(i, x, y);
        };
        loadVector(i) = volumeFactor * integrate(integrand); //scaling the integral between local and global coordinates
    }
}
//----------------compVectorEnd----------------
