#include "methods.h"
void forward_step(double non_diag,double diag, vec& v, vec& vPrev, int n)
// takes a=alpha, d=(1+2 alpha), v to fill with new solutions based on previous solutions vPrev
{
  v(1)= diag*vPrev(1) +  non_diag*vPrev(2); // calc  v_(1,j)
  for (int i = 2; i <=n-1 ; i++) //iterating and calc from v_(2,j) to v_(n-1,j)
  {
    v(i)= non_diag*(vPrev(i-1)+vPrev(i+1)) +  diag*vPrev(i);
  }
  v(n)= non_diag*vPrev(n-1) +  diag*vPrev(n); // calc  v_(n,j)
} //end forward_step function



void tridiag(double non_diag,double diag, vec& y, vec& v, int n) //Tridiagonal toeplitz solver
{
    //initialize vectors
    vec a = zeros<vec>(n+2);
    vec c = zeros<vec>(n+2);
    vec d = zeros<vec>(n+2);

    for (int i = 1; i <=n; i++)
    {
      a(i)=c(i)=non_diag; //  assign non diagonal elements
      d(i)=diag;          //  assign diagonal elements
    }
    a(1)=c(n)=0;          // correct non diagonal elements

    // Forward subsitution
    for (int i = 2; i <= n; i++)
    {
      d(i) -= c(i-1)*a(i)/d(i-1);
      v(i) -= v(i-1)*a(i)/d(i-1);
    }

    // Backward subsitution
    y(n)=v(n)/d(n);
    for (int i = n; i >=2 ; i--)
    {
      y(i-1) = (v(i-1)-c(i-1)*y(i))/d(i-1);
    }
}


void forward_euler(double alpha, vec& v, int n, double T, double dt) //
{
  vec vPrev = zeros<vec>(n+2);    //placeholder v(t-1)
  double non_diag = alpha;        //specify non diagonal elements
  double diag = 1-2*alpha;        //specify diagonal elements
  for (double t=0; t <=T; t+=dt)  //iterating over timesteps. Time prints in loop must be time+delta_t!
  {
    vec vPrev = v; //set previous solution to current
    forward_step(non_diag,diag,v,vPrev,n); //Solving tridiagonal matrix product: Av(t-1)=vt(1)
  }
} //end of forward_euler function


void backward_euler(double alpha, vec& v, int n, double T, double dt) //
{
  double non_diag=-alpha;         // specify non diagonal elements
  double diag=1+2*alpha;          // specify diagonal elements
  vec y = zeros<vec>(n+2);        // Ay=v
  for (double t=0; t <=T; t+=dt)  //iterating over timesteps. Time prints in loop must be time+delta_t!
  {
    tridiag(non_diag,diag,y,v,n); //Solving tridiagonal matrix product: Av_(t)=v(t-1)
    v=y;
  }

} //end of forward_euler function


void crank_nicolson(double alpha, vec& v, int n, double T, double dt) //
{
  vec y = zeros<vec>(n+2); //Ay=v
  vec vPrev = zeros<vec>(n+2);

  double FWnon_diag=alpha;  //set non diagonal elements in explicit step
  double FWdiag=2-2*alpha;  //set diagonal elements in explicit step

  double BWnon_diag=-alpha; //set non diagonal elements in implicit step
  double BWdiag=2+2*alpha;  //set diagonal elements in implicit step

  for (double t=0; t <=T; t+=dt) //iterating over timesteps. Time prints in loop must be time+delta_t!
  {
    vPrev=v;
    forward_step(FWnon_diag,FWdiag,v,vPrev,n);  // Solving tridiagonal matrix product: Av(t-1)=vt(1)
    tridiag(BWnon_diag,BWdiag,y,v,n);           // Solving tridiagonal matrix product: Av_(t)=v(t-1)
    v=y;                                        // Update solution
  }

} //end of crank_nicolson function

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
