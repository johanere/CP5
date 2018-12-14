#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "methods.h"


// output file
ofstream ofile;



int main(int argc, char* argv[])
{
  string filename, argument1,argument2;
  int method;
  double dx,T;
  if (argc <= 3)
  {
    cout << "Bad Usage: " << argv[0] <<
    "Provide: method, dx, T" << endl;
    exit(1);
  }
  if (argc > 1)
  {
    method = atoi(argv[1]);
    dx = atof(argv[2]);
    T = atof(argv[3]);

    argument1=argv[3];
    argument2=argv[2];
  }

  //int method=1;                   // approximation method, 1=forward euler
  //double dx=(double) 1.0/10;     // stepsize of x
  //double T=0.3;                  // total time

  string fileout = "Data/";      // name of output file
  fileout.append("results_");
  if (method==1) {fileout.append("FW");}
  if (method==2) {fileout.append("BW");}
  if (method==3) {fileout.append("CN");}
  if (method==4) {fileout.append("LU");}
  fileout.append("_");
  fileout.append(argument1);
  fileout.append("T_");
  fileout.append(argument2);
  fileout.append("dx");

  int n=(1/dx)-1.0;             // number of gridpoints between x=0 and x=l
  double dt=0.5*dx*dx;          // stepsize time
  double alpha = dt/(dx*dx);    // coefficient used to determine elements in matrix


  vec x = zeros<vec>(n+2);      // placeholder grid of size n+2: starts at x_0 and goes to x_(n+1)
  vec v = zeros<vec>(n+2);      // placeholder solution with variable change
  vec u = zeros<vec>(n+2);      // placeholder actual solution to PDE

  for (int i = 0; i <=n+1; i++) // setting up grid
  {
    x(i)=(double) i*dx;
  }
  u(n+1)=1.0;                   // I.C: u(i)=0 for all i =! n+1.
  v=u-x;                        // change of variables and setting initial conditions
                                // initiate print
  ofile.open(fileout);
  ofile << setiosflags(ios::showpoint | ios::uppercase);

  if (method==1)              // Solving using Forward Euler
  {
    forward_euler(alpha,v,n,T,dt);
    ofile << setw(15) << setprecision(8) <<"Forward Euler" <<endl;
  }
  if (method==2)              // Solving using Backward Euler
  {
    backward_euler(alpha,v,n,T,dt);
    ofile << setw(15) << setprecision(8) <<"Backward Euler" <<endl;
  }

  if (method==3)              // Solving using Crank-Nicolson
  {
    crank_nicolson(alpha,v,n,T,dt);
    ofile << setw(15) << setprecision(8) <<"Crank-Nicolson" <<endl;
  }

  if (method==4)              // Solving using LU
  {
    //initialzing vectors and matrix
    int N=n;
    mat A = zeros<mat>(N,N);
    vec rhs = zeros<vec>(N);
    vec sol = zeros<vec>(N);
    for (int i = 1; i <N-1; i++)
    {
      A(i,i)=1+2*alpha;
      A(i,i-1)=A(i,i+1)=-alpha;
      rhs(i)=v(i+1);
    }
    A(0,1)=A(N-1,N-2)=-alpha;
    A(0,0)=A(N-1,N-1)=1+2*alpha;
    rhs(0)=v(1);
    for (double t=0; t <=T; t+=dt) //iterating over timesteps.
    {
      // solve Ax = b
      vec sol = solve(A,rhs);
      rhs=sol;
    }
    for (int i = 1; i <=N; i++)
    {
      v(i)=rhs(i-1);
    }
    ofile << setw(15) << setprecision(8) <<"LU" <<endl;
  } //end of LU

  u=v+x;                            // extracting solution from v
  vec exact = zeros<vec>(n+2);      // placeholder analytical solution to PDE
  for(int j=0; j<=n+1; j++)
  {
    exact[j]=analytical_solution(x[j],T);   //finding exact solution
  }

  ofile << setw(15) << setprecision(8) <<n+2<<" gridpoints" <<endl;
  ofile << setw(15) << setprecision(8) <<n<<" number of gridpoints (n) between x=0 and x=1 " <<endl;
  ofile << setw(15) << setprecision(8) <<T<<" final temperature, T" <<endl;
  ofile << setw(15) << setprecision(8) <<dx<<" = dx"<<endl;
  ofile << setw(15) << setprecision(8) <<dt<<" = dt"<<endl;
  ofile << setw(15) << setprecision(8) <<"    numerical        exact           rel.err"<<endl;

  for (int i = 0; i <= n+1; i++)   // printing solution
  {
    double RelativeError=fabs((u(i)-exact(i))/exact(i));
    if(i==0){RelativeError=0;} //avoid including division by zero
    ofile << setw(15) << setprecision(8) << u(i);
    ofile << setw(15) << setprecision(8) << exact(i);
    ofile << setw(15) << setprecision(8) << RelativeError<< endl;
  }
  ofile.close();
} //end of main
