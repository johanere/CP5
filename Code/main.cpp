#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "armadillo"
using namespace std;
using namespace arma;

#define PI 3.14159265

// output file
ofstream ofile;

//header
void forward_step(double a, double d, vec& v, vec& vPrev, int n);
void forward_euler(double alpha, vec& v, int n, double T, double dt);
double analytical_solution(double x, double t);
///

int main()
{
  //read from CL
  int method=1;                 // approximation method, 1=forward euler
  double dx=(double) 1.0/100;               // stepsize of x
  double T= 0.1;                  // total time
  string fileout = "dingdong";


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
    ofile << setw(15) << setprecision(8) <<"Method: Forward Euler" <<endl;
  }


  vec exact = zeros<vec>(n+2);      // placeholder analytical solution to PDE
  for(int j=0; j<=n+1; j++)
  {
    exact[j]=analytical_solution(x[j],T);
  }


  u=v+x;  // extracting solution from v

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



void forward_step(double a, double d, vec& v, vec& vPrev, int n)
// takes a=alpha, d=(1+2 alpha), v to fill with new solutions based on previous solutions vPrev
{
  v(1)= d*vPrev(1) +  a*vPrev(2); // calc  v_(1,j)
  for (int i = 2; i <=n-1 ; i++) //iterating and calc from v_(2,j) to v_(n-1,j)
  {
    v(i)= a*(vPrev(i-1)+vPrev(i+1)) +  d*vPrev(i);
  }
  v(n)= a*vPrev(n-1) +  d*vPrev(n); // calc  v_(n,j)
} //end forward_step function



void forward_euler(double alpha, vec& v, int n, double T, double dt) //
{
  double a = alpha;
  double d = 1-2*alpha;
  vec vPrev = zeros<vec>(n+2);
  for (double t=0; t <=T; t+=dt) //iterating over timesteps. As time is updated at the end of a loop, time prints in loop must be time+delta_t!
  {
    vec vPrev = v; //set previous solution to current, and iterate current
    forward_step(a,d,v,vPrev,n); //Solving tridiagonal matrix product: Av=v_new
  }
} //end of forward_euler function

void tridiag(double alpha, vec& v, int n, double T, double dt) //
{
  double a = alpha;
  double d = 1-2*alpha;
  vec vPrev = zeros<vec>(n+2);
  for (double t=0; t <=T; t+=dt) //iterating over timesteps. As time is updated at the end of a loop, time prints in loop must be time+delta_t!
  {
    vec vPrev = v; //set previous solution to current, and iterate current
    forward_step(a,d,v,vPrev,n); //Solving tridiagonal matrix product: Av=v_new
  }
} //end of forward_euler function

double analytical_solution(double x, double t)
{
  int N=800; //chosen 3.5 chapter p.97
  double sum=0;
  for (int i=1; i<=N; i+=2) //negative terms
  {
    sum-=(double) 1/i*sin(i*PI*x)*exp(-i*i*PI*PI*t);
  }
  for (int i=2; i<=N; i+=2) //positive terms
  {
    sum+=(double) 1/i*sin(i*PI*x)*exp(-i*i*PI*PI*t);
  }
double ans=x+2/PI*sum;
return ans;
} //end of analytical_solution function
