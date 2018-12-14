//header guard
#ifndef METHODS_ADD_H
#define METHODS_ADD_H


#include <cmath>
#include "armadillo"

#define PI 3.14159265

using namespace std;
using namespace arma;
//declarations
double analytical_solution(double x, double y, double t);

void forward_step2d(double alpha, mat& v, mat& vPrev, int n);
//end of header guard
#endif
