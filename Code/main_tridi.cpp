#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "armadillo"
using namespace std;
using namespace arma;

int main()
{
  int n=5;

  int alpha=1;

  //initialzing vectors and assigning values to known entities
  vec v=ones<vec>(n+2);       //placeholder grid of size n+2: starts at x_0 and goes to x_(n+1) s/*
  vec a=zeros<vec>(n+2);
  vec c=zeros<vec>(n+2);//non-diagonal elements
  vec d=zeros<vec>(n+2); //diagonal elements
  for (int i = 0; i <n+1; i++)
  {
    a(i)=-alpha;
    d(i)=(1+2*alpha);
  }
  a(0)=a(1)=c(0)=c(n)=d(0)=0;

  // Forward subsitution
  double *m=new double[1];
  m = new double;
  for (int i = 2; i <= n; i++)
  {
    *m=a[i]/d[i-1];
    cout<<*m<<endl;
    d[i] -= *m*c[i-1];
    v[i] -= *m*v[i-1];
  }
  delete m;


  // Backward subsitution
  v[n-1]=b[n-1]/d[n-1];
  for (int i = n-2; i > 0; i--)
  {
    v[i] = (b[i]-c[i-1]*v[i+1])/d[i];*/
}
