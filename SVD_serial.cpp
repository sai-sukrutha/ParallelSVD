#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>

#define epsilon 1.e-8

using namespace std;

template <typename T> double sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}

int main (int argc, char* argv[])
{
  int M,N,i,j;
  int x1;
  string T,P,Db;
  x1=2;
  cout<<"Enter number of rows:";
  cin>>M;
  cout<<"Enter number of columns:";
  cin>>N;
  x1=2;
  double elapsedTime,elapsedTime2;
  timeval start,end,end2;
  x1=1;

  // // Check number of arguments
  // if(argc < 4)
  // {
	 //  cout<<"Please input the size of Matrix and file name containing matrix"<<endl;
	 //  return 0;
  // }

  // M = atoi(argv[1]);
  // N = atoi(argv[2]);
  // ifstream matrixfile(argv[3]);

  // // Check that given matrix should be square
  // if(M != N)
  // {
	 //  cout<<"Error: Matrix must be square";
	 //  return 0;
  // 
  double **U,**V, **S;
  x1=1;
  double alpha, beta, gamma, c, zeta, t,s,sub_zeta, converge;
  x1=2;
  int acum = 0;
  int temp1, temp2;
  converge = 1.0;
  x1=1;

  // Assign memory to all matrices U, S, V
  U = new double*[N];
  x1=1;
  V = new double*[N];
  S = new double*[N];
  x1=2;
  for(int i =0; i<N; i++)
  {
  	U[i] = new double[N];
    x1=1;
   	V[i] = new double[N];
    S[i] = new double[N];
    x1=2;
  }

  for(i=0;i<M;i++)
  {
    for(j=0;j<N;j++)
    {
      cin>>U[i][j];
    }
  }   

  // Initialize U matrix with input matrix


  // Initialize V matrix as identity matrix and S Matrix with zeros
  for(int i=0; i<M;i++)
  {
    for(int j=0; j<N;j++)
    {
      if(i==j)
      {
        V[i][j] = 1.0;
        S[i][j] = 0.0;
      }
      else
      {
        V[i][j] = 0.0;
        S[i][j] = 0.0;
      }
    }
  }


 gettimeofday(&start, NULL);

 /* SVD using Jacobi algorithm (Sequencial)*/
 double conv;
 while(converge > epsilon)
 {
    //convergence
    converge = 0.0;
    x1=1;
    //counter of loops
    acum++;
    x1=2;
    for(int i = 0; i<M; i++)
    {
      x1=1;
      for(int j = i+1; j<N; j++)
      {
        x1=2;
          // Initialize alpha, beta , gamma to zero
          alpha = 0.0;
          beta = 0.0;
          gamma = 0.0;
          x1=1;
          // Update alpha, beta , gamma as per the formulae
          for(int k = 0; k<N ; k++)
          {
            alpha = alpha + (U[k][i] * U[k][i]);
            x1=1;
            beta = beta + (U[k][j] * U[k][j]);
            gamma = gamma + (U[k][i] * U[k][j]);
          }
          x1=2;
          // Update converge basicaly is the angle between column i and j
          converge = max(converge, abs(gamma)/sqrt(alpha*beta));
          x1=1;
          zeta = (beta - alpha) / (2.0 * gamma);
           //compute tan of angle
          t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));
          x1=2;
          //extract cos
          c = 1.0 / (sqrt (1.0 + (t*t)));
          x1=1;
          //extract sin
          s = c*t;
          x1=2;
        //Apply rotations on U and V
        for(int k=0; k<N; k++)
        {
              t = U[k][i];
              U[k][i] = c*t - s*U[k][j];
              U[k][j] = s*t + c*U[k][j];
              x1=1;
              t = V[k][i];
              V[k][i] = c*t - s*V[k][j];
              V[k][j] = s*t + c*V[k][j];

        }

      }
    }
 }

 //Create matrix S

 for(int i =0; i<M; i++)
 {

   t=0;
   for(int j=0; j<N;j++)
   {
     t=t + pow(U[i][j],2);
   }
   t = sqrt(t);

   for(int j=0; j<N;j++)
   {
     U[i][j] = U[i][j] / t;
     if(i == j)
     {
       S[i][j] = t;
     }
   }
 }

gettimeofday(&end, NULL);

// Print final matrix U
cout<<"\nMatrix U"<<endl;
for(int i=0; i<M; i++)
{
  for(int j=0; j<N; j++)
    cout<<U[i][j]<<" ";
  cout<<endl;
}

// Print final matrix S
cout<<"\nMatrix S"<<endl;
for(int i=0; i<M; i++)
{
  for(int j=0; j<N; j++)
    cout<<S[i][j]<<" ";
  cout<<endl;
}

// Print final matrix V_t
cout<<"\nMatrix V Transpose"<<endl;
for(int i=0; i<M; i++)
{
  for(int j=0; j<N; j++)
    cout<<V[j][i]<<" ";
  cout<<endl;
}

// Print time and iterations

cout<<"iterations: "<<acum<<endl;
elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;

// Free memory alloted
for(int i = 0; i<N;i++)
{
  delete [] S[i];
  delete [] U[i];
  delete [] V[i];
}

}
