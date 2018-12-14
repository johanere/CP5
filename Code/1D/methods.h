//header guard
#ifndef METHODS_ADD_H
#define METHODS_ADD_H


#include <cmath>
#include "armadillo"

#define PI 3.14159265

using namespace std;
using namespace arma;
//declarations
void forward_step(double non_diag,double diag, vec& v, vec& vPrev, int n);
void forward_euler(double alpha, vec& v, int n, double T, double dt);
double analytical_solution(double x, double t);
void backward_euler(double alpha, vec& v, int n, double T, double dt);
void tridiag(double non_diag,double diag, vec& y, vec& v, int n);
void crank_nicolson(double alpha, vec& v, int n, double T, double dt); //
//end of header guard
#endif
