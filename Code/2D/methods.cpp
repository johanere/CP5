#include "methods.h"

double analytical_solution(double x, double y, double t)
{
  int N=400; //chosen 3.5 chapter p.97
  int M=400; //chosen 3.5 chapter p.97
  double sum=0;
  for (int i=1; i<=M; i+=2) //negative terms
  {
    for(int j=1; i<=N; i+=2)
    {
    sum+=(double)1/(i*j)*sin(PI*i*x)*sin(PI*j*y)*exp(-(i*i+j*j)*PI*PI*t);
    }
    for(int j=2; i<=N; i+=2)
    {
    sum-=(double)1/(i*j)*sin(PI*i*x)*sin(PI*j*y)*exp(-(i*i+j*j)*PI*PI*t);
    }
  }

  for (int i=2; i<=M; i+=2) //negative terms
  {
    for(int j=1; i<=N; i+=2)
    {
    sum-=(double)1/(i*j)*sin(PI*i*x)*sin(PI*j*y)*exp(-(i*i+j*j)*PI*PI*t);
    }
    for(int j=2; i<=N; i+=2)
    {
    sum+=(double)1/(i*j)*sin(PI*i*x)*sin(PI*j*y)*exp(-(i*i+j*j)*PI*PI*t);
    }
  }
double ans=4/(PI*PI)*sum;
return ans;
} //end of analytical_solution function

void forward_step2d(double alpha, mat& v, mat& vPrev, int n)
{

  for (int i = 1; i <=n ; i++) //iterating and calc from v_(2,j) to v_(n-1,j)
  {
    for (int j = 1; j <=n ; j++) //iterating and calc from v_(2,j) to v_(n-1,j)
    {
      v(i,j)=vPrev(i,j)+alpha*(vPrev(i+1,j)+vPrev(i-1,j)+vPrev(i,j+1)+vPrev(i,j-1)-4*vPrev(i,j));
    }
  }

} //end of analytical_solution function*/
