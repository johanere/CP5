#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "methods.h"


// output file
ofstream ofile;



int main(int argc, char* argv[])
{
/*  string filename, argument1,argument2;
  double h,T;
  if (argc <= 2)
  {
    cout << "Bad Usage: " << argv[0] <<
    "Provide: method, dx, T" << endl;
    exit(1);
  }
  if (argc > 1)
  {
    h = atof(argv[2]);
    T = atof(argv[3]);

    argument1=argv[3];
    argument2=argv[2];
  }*/
  string filename, argument1,argument2;
  double h,T;

  h=0.01;
  T=0.05;

  argument1="test1";
  argument1="test2";
  string fileout = "Data/";      // name of output file
  fileout.append("results2D_");
  fileout.append("FW");

  fileout.append("_");
  fileout.append(argument1);
  fileout.append("T_");
  fileout.append(argument2);
  fileout.append("h");


  int n=(1/h)-1.0;             // number of gridpoints between x=0 and x=l
  cout<<n<<endl;
  double dt=0.1*0.5*h*h;          // stepsize time
  double alpha = dt/(h*h);    // coefficient used to determine elements in matrix

  vec x = zeros<vec>(n+2);      // placeholder grid of size n+2: starts at x_0 and goes to x_(n+1)
  vec y = zeros<vec>(n+2);      // placeholder y axis
  mat u = zeros<mat>(n+2,n+2);  // placeholder actual solution to PDE
  for (int i = 0; i <=n+1; i++) // setting up grid
  {
    x(i)=(double) i*h;
    y(i)=(double) i*h;
  }

  mat v = zeros<mat>(n+2,n+2);  // placeholder actual solution to PDE
  mat vprev = zeros<mat>(n+2,n+2);  // placeholder actual solution to PDE


  for (int i = 1; i <=n ; i++) //iterating and calc from v_(2,j) to v_(n-1,j)
  {
    for (int j = 1; j <=n ; j++) //iterating and calc from v_(2,j) to v_(n-1,j)
    {
      v(i,j)=x(i)*y(i);
    }
  }

  ofile.open(fileout);
  ofile << setiosflags(ios::showpoint | ios::uppercase);

  for (double t=0; t <=T; t+=dt)  //iterating over timesteps. Time prints in loop must be time+delta_t!
  {
    vprev=v;
    forward_step2d(alpha,v,vprev,n);
  }
  ofile << setw(15) << setprecision(8) <<"Forward Euler" <<endl;


  v.print("etter");

                           // extracting solution from v
  mat exact = zeros<mat>(n+2,n+2);      // placeholder analytical solution to PDE
  mat relerr = zeros<mat>(n+2,n+2);      // placeholder analytical solution to PDE
  for(int i=1; i<=n; i++)
  {
    for (int j=1;j<=n;j++)
    {
      exact(i,j)=analytical_solution(x(i),y(i),T);   //finding exact solution
      relerr(i,j)=abs((v(i,j)-exact(i,j))/(exact(i,j)));

    }
  }

  double avrerr=accu(relerr)/(n*n);
  cout<<avrerr<<endl;
  ofile << setw(15) << setprecision(8) <<n+2<<" gridpoints" <<endl;


  for (int i = 0; i <= n+1; i++)   // printing solution
  {
    for (int j=0;j<=n+1;j++)
    {
      ofile << setw(15) << setprecision(8) << exact(i,j);
    }
    ofile << setw(15) << setprecision(8) << endl;
  }
  ofile.close();
} //end of main
