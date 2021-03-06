#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "methods.h"


// output file
ofstream ofile;

int main()
{
  string filename;

  string h_s="0.01";
  string T_s="0.2";
  double h=stof(h_s);
  double T=stof(T_s);

  int n=(1/h)-1.0;
  int graph=0;
  int incgam=1;
  int factor=6;

  string fileout = "Data/";      // name of output file
  fileout.append("Avr_");
  fileout.append(h_s);
  fileout.append("_");
  fileout.append(T_s);
  ofile.open(fileout);
  ofile << setiosflags(ios::showpoint | ios::uppercase);

  vec x = zeros<vec>(n+2);      // placeholder grid of size n+2: starts at x_0 and goes to x_(n+1)
  vec y = zeros<vec>(n+2);      // placeholder y axis
  mat v = zeros<mat>(n+2,n+2);  // placeholder actual solution to PDE
  mat vprev = zeros<mat>(n+2,n+2);
  mat exact = zeros<mat>(n+2,n+2);
  mat relerr = zeros<mat>(n+2,n+2);
  double avrerr =0;
  double dt;
  double alpha;
  double dtscale;
  for (int k=1;k<=factor;k++)
  {
    cout<<k<<endl;
    dtscale=(double)(4.2-(1.0/k)));

    dt=dtscale*h*h;          // stepsize time
    alpha = dt/(h*h);    // coefficient used to determine elements in matrix

    for (int i = 0; i <=n+1; i++) // setting up grid
    {
      x(i)=(double) i*h;
      y(i)=(double) i*h;
    }




    for (int i = 1; i <=n ; i++)
    {
      for (int j = 1; j <=n ; j++)
      {
        v(i,j)=x(i)*y(i);
      }
    }

    for (double t=0; t <=T; t+=dt)  //iterating over timesteps.
    {
      vprev=v;
      forward_step2d(alpha,v,vprev,n);
    }


    for(int i=1; i<=n; i++)
    {
      for (int j=1;j<=n;j++)
      {
        exact(i,j)=analytical_solution(x(i),y(j),T);   //finding exact solution
        relerr(i,j)=abs(  (v(i,j)-exact(i,j))/(exact(i,j))  );

      }
    }

    avrerr=accu(relerr)/(n*n);
  if (incgam==1)
  {
    ofile << setw(15) << setprecision(8) <<dtscale<<" "<<avrerr<<endl;
  }
  else
  {
    ofile << setw(15) << setprecision(8) <<avrerr<<endl;
  }
  }
  ofile.close();

  if (graph==1) {
  cout<<"producing data for plotting"<<endl;
  fileout = "Data/";      // name of output file
  fileout.append("exact1");
  ofile.open(fileout);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) <<exact<<endl;
  ofile.close();

  fileout = "Data/";      // name of output file
  fileout.append("numerical1");

  ofile.open(fileout);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) <<v<<endl;
  ofile.close();

  fileout = "Data/";      // name of output file
  fileout.append("relerr");

  ofile.open(fileout);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) <<relerr<<endl;
  ofile.close();
  }

} //end of main
