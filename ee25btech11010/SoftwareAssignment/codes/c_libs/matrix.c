#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

void transpose(double **A,double **B,int m,int n){
for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
    B[j][i]=A[i][j];
}
}
}

double frobeniuserror(double **A,double **Ak,int m,int n){
double sum=0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
           double difference=(A[i][j]-Ak[i][j]);
         sum+=difference*difference;
        }
    }
    return sqrt(sum);
}

void mult(double **A,double **B,int m,int n,int p,double **C){
for(int i=0;i<m;i++){
for(int j=0;j<p;j++){
    C[i][j]=0.0;
    for(int k=0;k<n;k++){
    C[i][j]+=(A[i][k]*B[k][j]);
}
}
}
}

void display(double **A,int m,int n){
for(int i=0;i<m;i++){
for(int j=0;j<n;j++){
    printf("%lf ",A[i][j]);
}
printf("\n");
}
}

void singularvalues(int n,double *a,double **V,int issymmetric){
for(int i=0;i<n;i++){
double lam=V[i][i];
// if (issymmetric == 0)
// {
//     if (lam < 0.0)
//     {
//         lam = 0.0;
//     }
//     a[i] = sqrt(fabs(lam));
// }
// else
// {
    a[i]=sqrt(fabs(lam));

}
}

void matrixU(double **A,double **V,double **U,double *sigma,int m,int n,double **B){
for(int i=0;i<m;i++){
for(int j=0;j<n;j++){
    U[i][j]=0.0;
}
}
for(int j=0;j<n;j++){
if(sigma[j] > 1e-12){
  for(int i=0;i<m;i++){
      U[i][j]=B[i][j]/sigma[j];
  }
}
else{
   for(int i=0;i<m;i++){
       U[i][j]=0.0;
   }
}
}
}

void swap(double *a,double *b){
double t;
t=*a;
*a=*b;
*b=t;
}

void swapcolumn(double **A,int rows,int a,int b){
if(a==b){
    return;
}
for(int i=0;i<rows;i++){
    double temp=A[i][a];
    A[i][a]=A[i][b];
    A[i][b]=temp;
}
}

void sortsigma_andE(double *s,double **E,int n){
for(int i=0;i<n-1;i++){
for(int j=0;j<n-i-1;j++){
if(s[j]<s[j+1]){
        swap(&s[j],&s[j+1]);
        swapcolumn(E,n,j,j + 1);
    }
}
}
}

double **allocatememory(int rows,int columns){
double **A;
A=(double **)malloc(rows * sizeof(double *));
for(int i=0;i<rows;i++){
    A[i]=(double *)malloc(columns * sizeof(double));
}
return A;
}

void freememory(double **A,int rows){
for(int i=0;i<rows;i++){
free(A[i]);
}
free(A);
}

void calcUkandVk(double **U,double **Uk,int m,int k){
for(int i=0;i<m;i++){
for(int j=0;j<k;j++){
Uk[i][j]=U[i][j];
}
}
}

void calcsigmak(double *a,double **A,int k){
for(int i=0;i<k;i++){
for(int j=0;j<k;j++){
if(i==j){
    A[i][j]=a[i];
}
else{
    A[i][j]=0.0;
}
}
}
}

double svderror(double *s,int k,int n){
    double sum=0;
    for(int i=k;i<n;i++){
        sum+=s[i]*s[i];
    }
    return sqrt(sum);
}