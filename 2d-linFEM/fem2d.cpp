#include <Eigen/Core>
#include "stiffness_matrix.hpp"
#include "fem_solve.hpp"
#include "writer.hpp"
#include <sstream>
#include "Lshape.hpp"
#include "convergence.hpp"


double sigma(double x, double y){
  return 0.01*(x+2)*(x+2);
}


// Sample solution: square

double uex_square(double x, double y) {
    return sin(2*M_PI *x) * sin(2*M_PI * y);
}

double f_square(double x, double y) {
    return -(x+2)*M_PI/25 * cos(2*M_PI*x)*sin(2*M_PI*y) + 
            (0.5 + 0.08*M_PI*M_PI*(x+2)*(x+2))*uex_square(x,y);
}

double g_square(double x, double y) {
    return 0;
}

Eigen::Vector2d uex_grad_square(double x, double y) {
  Eigen::Vector2d grad;
  grad << 2*M_PI*cos(2*M_PI *x)*sin(2*M_PI*y),
          2*M_PI*sin(2*M_PI *x)*cos(2*M_PI*y);
  return grad;
}

int main(int, char**) {
  try {

    solveL(0.5, sigma);

    /* 
    // Uncomment this block to perform analysis of convergence
    std::cout << "Convergence Analysis" << std::endl;
    convergenceAnalysis("square", 6, f_square, sigma, g_square, 0.5,
			uex_square, uex_grad_square); 
    */

    }
    catch (std::runtime_error& e) {
        std::cerr << "An error occurred. Error message: " << std::endl;
        std::cerr << "    \"" << e.what() << "\"" << std::endl;
        return EXIT_FAILURE;
    } 
    catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
        throw ;
    }

    return EXIT_SUCCESS;

}
