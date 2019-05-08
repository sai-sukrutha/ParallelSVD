#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <cstdlib>

#define epsilon 1.e-8
#define num 16

using namespace std;

template <typename T> double sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}

__global__ void parallel2(double** dev_V, double* dev_s, double* dev_c, double** dev_U, int* dev_i, int* dev_j, int* dev_N){
  for(int k=0; k < *dev_N; k++)
    {
          double temp = dev_U[k][*dev_i];
          dev_U[k][*dev_i] = *(dev_c)*temp - *(dev_s)*dev_U[k][*dev_j];
          dev_U[k][*dev_j] = *(dev_s)*temp + *(dev_c)*dev_U[k][*dev_j];

          temp = dev_V[k][*dev_i];
          dev_V[k][*dev_i] = *(dev_c)*temp - *(dev_s)*dev_V[k][*dev_j];
          dev_V[k][*dev_j] = *(dev_s)*temp + *(dev_c)*dev_V[k][*dev_j];

    }
}

__global__ void parallel1(double* dev_alpha, double* dev_beta, double* dev_gamma, double** dev_U, int* dev_i, int* dev_j, int* dev_N){
  for(int k = 0; k < *dev_N ; k++)
    {
      *dev_alpha = *dev_alpha + (dev_U[k][*dev_i] * dev_U[k][*dev_i]);
      *dev_beta = *dev_beta + (dev_U[k][*dev_j] * dev_U[k][*dev_j]);
      *dev_gamma = *dev_gamma + (dev_U[k][*dev_i] * dev_U[k][*dev_j]);
    }
}

int main (int argc, char* argv[])
{
  int M,N;
  string T,P,Db;

  //double elapsedTime,elapsedTime2;
  //timeval start,end,end2;

  // Check number of arguments
  if(argc < 3)
  {
	  cout<<"Please input the size of Matrix\n";
	  return 0;
  }

  M = atoi(argv[1]);
  N = atoi(argv[2]);

  //printf("N-%d\n",N);

  // Check that given matrix should be square
  if(M != N)
  {
	  cout<<"Error: Matrix must be square";
	  return 0;
  }

  double **U,**V, **S,**A;
  double **U_t,**V_t;
  double *alpha, *beta, *gamma, *c, *zeta, t,*s, converge;
  double *dev_alpha ,*dev_beta,*dev_gamma,*dev_c,*dev_s;
  int *dev_N;
  double **dev_U,**dev_V;
  int *dev_i,*dev_j;

  alpha = (double *)malloc(sizeof(double));
  beta = (double *)malloc(sizeof(double));
  gamma = (double *)malloc(sizeof(double));
  zeta = (double *)malloc(sizeof(double));

  int acum = 0;
  converge = 1.0;

  // Assign memory to all matrices U, S, V
  U = new double*[N];
  V = new double*[N];
  U_t = new double*[N];
  V_t = new double*[N];
  A = new double*[N];
  S = new double*[N];

  for(int i =0; i<N; i++){
	U[i] = new double[N];
 	V[i] = new double[N];
	U_t[i] = new double[N];
	V_t[i] = new double[N];
	A[i] = new double[N];
 	S[i] = new double[N];
  }

  //printf("Initializations done\n");

  //Generate random matrix
 
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            //cout<<"In init "<<i<<" "<<j<<" ";
            A[i][j] = rand()%10+1;
        }
    }
  

    cout<<"A"<<endl<<endl;
    for(int i =0; i<M; i++){
        for(int j =0; j<N; j++){

             cout<<A[i][j]<<"  ";
        }
        cout<<endl;
    }

    //Copy to U_t
    for (int i=0;i<M;i++){
        for (int j=0;j<N;j++){
            U_t[i][j]=A[j][i];
        }
    } 

   //printf("Copy to U_t\n");
  
  // Initialize V matrix as identity matrix and S Matrix with zeros
  for(int i=0; i<M;i++)
  {
    for(int j=0; j<N;j++)
    {
      //printf("i-%d,j-%d\n",i,j);

      if(i==j)
      {
        V[i][j] = (double)1.0;
        S[i][j] = (double)0.0;
      }
      else
      {
        V[i][j] = (double)0.0;
        S[i][j] = (double)0.0;
      }
    }
  }

  //printf("V build\n");


   int * N_temp = (int *)malloc(sizeof(int));
   *N_temp = N;

  //printf("Conv loop start\n");
 //gettimeofday(&start, NULL);

 /* SVD using Jacobi algorithm (Sequencial)*/
 while(converge > epsilon)
 {
    //convergence
    converge = 0.0;

    //counter of loops
    acum++;

    int * i_temp = (int *)malloc(sizeof(int));
    int * j_temp = (int *)malloc(sizeof(int));

    for(int i = 0; i<M; i++)
    {
      for(int j = i+1; j<N; j++)
      {

          i_temp = &i;
          j_temp = &j;
          // Initialize alpha, beta , gamma to zero
          *alpha = 0.0;
          *beta = 0.0;
          *gamma = 0.0;

          // Update alpha, beta , gamma as per the formulae
          cudaMalloc ((void**)&dev_alpha   , sizeof(double));
          cudaMalloc ((void**)&dev_beta   , sizeof(double));
          cudaMalloc ((void**)&dev_gamma   , sizeof(double));
          cudaMalloc ((void**)&dev_i   , sizeof(int));
          cudaMalloc ((void**)&dev_j, sizeof(int));
          cudaMalloc ((void**)&dev_N, sizeof(int));
          cudaMalloc ((void***)&dev_U, sizeof(double)*N*N);

          cudaMemcpy(dev_alpha ,alpha, sizeof(double), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_beta ,beta, sizeof(double), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_gamma ,gamma, sizeof(double), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_i ,i_temp, sizeof(int), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_j ,j_temp, sizeof(int), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_N ,N_temp, sizeof(int), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_U ,U, sizeof(double)*N*N, cudaMemcpyHostToDevice);
          
          //printf("Call to parallel1\n");
          parallel1<<<1,4>>>(dev_alpha, dev_beta, dev_gamma, dev_U, dev_i, dev_j, dev_N);
          //printf("After paralell1\n");

          cudaMemcpy(U, dev_U, sizeof(double)*N*N, cudaMemcpyDeviceToHost);
          //printf("1\n");
          //free(dev_alpha);
          //free(dev_beta);
          //free(dev_gamma);
          //printf("2\n");
          //free(dev_i);
          //free(dev_j);
          //free(dev_N);
          //free(dev_U);
          //printf("3\n");
          

          //printf("Before c,s cal\n");
      
          // Update converge basicaly is the angle between column i and j
          double gamma_t = *gamma;
          //printf("**1\n");
          double alpha_t = *alpha;
          //printf("**1\n");
          double beta_t = *beta;
          //printf("**1\n");
          double zeta_t = *zeta;
          //printf("**1\n");
          double c_t = *c;
          //printf("**1\n");
          //double s_t = *s;
          //printf("**1\n");
          converge = max(converge, abs(gamma_t)/sqrt(alpha_t*beta_t));
          //printf("**2\n");
          zeta_t = (beta_t - alpha_t) / (2.0 * gamma_t);
           //compute tan of angle
          t = sgn(zeta_t) / (abs(zeta_t) + sqrt(1.0 + (zeta_t*zeta_t)));
          //extract cos
          c_t = 1.0 / (sqrt (1.0 + (t*t)));
          //extract sin
          *s = c_t*t;
         
          //printf("**\n");
          *alpha = alpha_t;
          *beta = beta_t;
          *gamma = gamma_t;
          *zeta = zeta_t;
          *c = c_t;
          //*s = s_t;

         
          //printf("After c,s cal\n");

        //Apply rotations on U and V
        cudaMalloc ((void**)&dev_i   , sizeof(int));
        cudaMalloc ((void**)&dev_j   , sizeof(int));
        cudaMalloc ((void**)&dev_N   , sizeof(int));
        cudaMalloc ((void***)&dev_U   , sizeof(double)*N*N);
        cudaMalloc ((void***)&dev_V   , sizeof(double)*N*N);
        cudaMalloc ((void**)&dev_s   , sizeof(double));
        cudaMalloc ((void**)&dev_c   , sizeof(double));

        cudaMemcpy(dev_V ,V, sizeof(double)*N*N, cudaMemcpyHostToDevice);
        cudaMemcpy(dev_s ,s, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_c ,c, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_i ,i_temp, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_j ,j_temp, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_N ,N_temp, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_U ,U, sizeof(double)*N*N, cudaMemcpyHostToDevice);
        

        //printf("Call to parallel2\n");
        parallel2<<<1,4>>>(dev_V, dev_s, dev_c, dev_U, dev_i, dev_j, dev_N);
        //printf("After paralell2\n");

        cudaMemcpy(U, dev_U, sizeof(double)*N*N, cudaMemcpyDeviceToHost);
        cudaMemcpy(V, dev_V, sizeof(double)*N*N, cudaMemcpyDeviceToHost);
        free(dev_V);
        free(dev_s);
        free(dev_c);
        free(dev_i);
        free(dev_j);
        free(dev_N);
        free(dev_U);

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

   //double t_t = *t;
   //t_t = sqrt(t_t);
   //*t = t_t;

   t=sqrt(t);

   for(int j=0; j<N;j++)
   {
     U[i][j] = U[i][j] / t;
     if(i == j)
     {
       S[i][j] = t;
     }
   }
 }

//gettimeofday(&end, NULL);

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
//elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
//elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
//cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;



}
