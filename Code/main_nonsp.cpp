#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "armadillo"
using namespace std;
using namespace arma;

void tridiag(double a, double d, vec& v, vec& v_new, int n);

int main()
{

  int n=10; // grid containts n+2 points (999)
  double h= (double) 1/(n+1);   // stepsize
  double h2=h*h;                // stepsize squared
  double delta_t=0.5*h2;        // timestep
  double time=0;

  vec x= zeros<vec>(n+2);       //placeholder grid of size n+2: starts at x_0 and goes to x_(n+1)
  vec v = zeros<vec>(n+2);      //placeholder old solution
  vec v_new = zeros<vec>(n+2);  //placeholder new solution
  vec u = zeros<vec>(n+2);      //placeholder actual solution to PDE

  for (int i = 0; i <=n+1; i++) //setting up grid
  {
    x(i)=(double) i*h;
  }
  u(0)=0; u(n+1)=1.0;// boundry conditions
  v=u-x; //change of variables and setting initial conditions

  double a=h2; double d=(1-2*a); // a=non diagonal elements of matrix,  d=diagonal elements of matrix A


  for (time; time <=1; time+=delta_t) //iterating over timesteps. As time is updated at the end of a loop, time prints in loop must be time+delta_t!
  {
    tridiag(a,d,v,v_new,n); //Solving tridiagonal matrix product: Av=v_new
    v=v_new;  //updating solution
    u=v+x;    //extracting actual updated solution
  }
u.print("U final");
} //end of main

void tridiag(double a, double d, vec& v, vec& v_new, int n)  //Solving tridiagonal matrix product: Av=v_new
{
  v_new(1)=d*v(1)+a*v(2); // find  v_(1,j)
  for (int i = 2; i <=n-1 ; i++) //iterating find  v_(2,j) to v_(n-1,j)
  {
    v_new(i)=a*v(i-1)+d*v(i)+a*v(i+1);
  }
    v_new(n)=a*v(n-1)+d*v(n); // find  v_(n,j)
} //end of func. tridiag
