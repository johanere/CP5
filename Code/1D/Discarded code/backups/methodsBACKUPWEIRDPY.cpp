void tridiag(double alpha, vec& v, int n) //Tridiagonal toeplitz solver
{
  vec d=zeros<vec>(n+2);  //non-diagonal elements
  vec c=zeros<vec>(n+2);  //non-diagonal elements
  vec vNext=zeros<vec>(n+2);  //non-diagonal elements
  for (int i=1;i<=n;i++)
  {
  d(i)=(1+2*alpha);
  c(i)=-alpha;
  }

  for (int i=1;i<=n-1;i++)
  {
    //normalize
    c(i)/=d(i);
    v(i)/=d(i);
    d(i)=1.0;

    d(i+1)+=c(i-1)*alpha;
    v(i+1)+=v(i-1)*alpha;
  }
  vNext(n)=v(n)/d(n);

  for(int i=n-1;i<=1;i--)
  {
    vNext(i)=v(i)-vNext(i+1)*c(i);
  }
  for(int i=0;i<=n+1;i++)
  {
    v(i)=vNext(i);
  }
}






















void tridiag(double alpha, vec& v, int n) //Tridiagonal toeplitz solver
{
 //initialzing vectors and assigning values to known entities
  vec c=zeros<vec>(n+2);  //non-diagonal elements
  vec d=zeros<vec>(n+2);  //diagonal elements
  vec a=zeros<vec>(n+2);  //diagonal elements
  vec vPrev=zeros<vec>(n+2);  //diagonal elements
  for (int i = 0; i <= n; i++)
  {
    c(i)=-alpha;
    d(i)=(1+2*alpha);
    a(i)=-alpha;
    vPrev(i)=v(i);
  }


  c(0)=c(n)=c(n+1)=0;
  d(0)=d(n+1)=0;
  a(0)=a(1)=a(n+1)=0;
  vPrev(n+1)=v(n+1);
  vPrev(0)=v(0);
  // Forward subsitution
  double *m=new double[1];
  m = new double;

  for (int i = 2; i <= n; i++)
  {
    *m=a(i)/d(i-1);
    d(i) -= *m*c(i-1);
    vPrev(i) -= *m*v(i-1);
  }
  delete m;

  // Backward subsitution
  v(n)=vPrev(n)/d(n);
  for (int i = n-1; i >= 1; i--)
  {
    v(i)=(vPrev(i)-c(n)*v(i-1))/d(i);
  }
  v(0)=v(n+1)=0;
}//end of tridiag function



void tridiag(double alpha, vec& v, int n) //Tridiagonal toeplitz solver
{
    vec d=zeros<vec>(n+2);  //non-diagonal elements
    vec c=zeros<vec>(n+2);  //non-diagonal elements

    d+=(1+2*alpha);
    c+=-alpha;

    c(1)/=d(1);
    v(1)/=d(1);
    d(1)=1.0;

    for (int i=2;i<=n-1;i++)
    {
      v(i)+=v(i-1)*alpha;
      d(i)+=c(i-1)*alpha;

      c(i)/=d(i);
      v(i)/=d(i);
      d(i)=1.0;

    }
    v(n)/=d(n);

    for(int i=n-1;i<=1;i--)
    {
      v(i)-=v(i+1)*c(i);
    }
  }



 //initialzing vectors and assigning values to known entities
  vec a=zeros<vec>(n+2);  //non-diagonal elements
  vec c=zeros<vec>(n+2);  //non-diagonal elements
  vec d=zeros<vec>(n+2);  //diagonal elements
  vec vnew=zeros<vec>(n+2);  //diagonal elements
  for (int i = 0; i <= n; i++)
  {
    a(i)=-alpha;
    c(i)=-alpha;
    d(i)=(1+2*alpha);
  }
  a(0)=a(1)=c(0)=c(n)=d(0)=0;
  // Forward subsitution
  double *m=new double[1];
  m = new double;
  for (int i = 2; i <= n; i++)
  {
    *m=a(i)/d(i-1);
    d(i) -= *m*c(i-1);
    v(i) -= *m*v(i-1);
  }
  delete m;

  // Backward subsitution
  vnew(n)=v(n)/d(n);
  for (int i = n-1; i >= 1; i--)
  {
    vnew(i)=(v(i)-c(n)*vnew(i-1))/d(i);
  }
  v=vnew;
  v.print("in triag");*/
}//end of tridiag function
