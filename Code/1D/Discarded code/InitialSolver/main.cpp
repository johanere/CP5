#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>


#include <cmath>
#include "armadillo"

#define PI 3.14159265

using namespace std;
using namespace arma;



// output file
ofstream ofile;



int main()
{
  //read from CL
  int method=1;                 // approximation method, 1=forward euler
  double dx=(double) 1.0/10;               // stepsize of x
  double T= 0.1;                  // total time

  int n=(1/dx)-1.0;             // number of gridpoints between x=0 and x=l
  double dt=0.5*dx*dx;          // stepsize time
  double alpha = dt/(dx*dx);    // coefficient used to determine elements in matrix

  vec x = zeros<vec>(n);
  vec v = zeros<vec>(n);
  for(int i=0;i<n;i++)
  {
    x(i)=(i+1)*dx;
    v(i)=-x(i);
  }

  if (method==1)              // Solving using Forward Euler
  {
    forward_euler(alpha,v,n,T,dt);
    ofile << setw(15) << setprecision(8) <<"Method: Forward Euler" <<endl;
  }
  if (method==2)              // Solving using Forward Euler
  {
    backward_euler(alpha,v,n,T,dt);
    ofile << setw(15) << setprecision(8) <<"Method: Forward Euler" <<endl;
  }
}
