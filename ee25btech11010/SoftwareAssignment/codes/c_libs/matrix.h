#ifndef MATRIX_H
#define MATRIX_H

double** allocatememory(int rows,int columns);
void freememory(double **A,int rows);

void transpose(double **A,double **B,int m,int n);
void mult(double **A,double **B,int m,int n,int p,double **C);
void display(double **A,int m,int n);

void singularvalues(int n,double *a,double **V,int issymetric);
void matrixU(double **A,double **V,double **U,double *sigma,int m,int n,double **B);

void swap(double *a,double *b);
void swapcolumn(double **A,int rows,int a,int b);
void sortsigma_andE(double *s,double **E,int n);

void calcUkandVk(double **U,double **Uk,int m,int k);
void calcsigmak(double *a,double **A,int k);
double frobeniuserror(double **A,double **Ak,int m,int n);
double svderror(double *s,int k,int n);

#endif