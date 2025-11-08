#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jacobi.h"


void jacobi(double **V,double **E,int n){
double tol=1e-12; // tolerance 
int max_sweeps=50; 
for(int sweep=0;sweep<max_sweeps;sweep++){          //one sweep means covering every single off diagonal elem in a systematic manner
int rotations_done=0;           //variable to check whether any rotaion was performed in this particular sweep
for(int im=0;im<n;im++){
for(int jm=im+1;jm<n;jm++){             //parsing thru all upper triangular elements
if(fabs(V[im][jm])<tol){
    continue; 
}
//if elem is big enough then rotate matrix
rotations_done++;

// creating rotation parameters
double xt=(2* V[im][jm]);
double t1=(V[im][im] - V[jm][jm]) / xt; 
double xq=1+t1*t1;
double t;
if(t1>=0){
    t= 1 / (t1 + sqrt(xq));
} 
else{
    t= 1 / (t1 - sqrt(xq));
}
double c= 1 / (sqrt(1 + t * t));
double s=t*c;
//update V:the diag matrix with eigen values
double Vii=V[im][im];
double Vjj=V[jm][jm];
double Vij=V[im][jm];
V[im][im]=Vii*c*c+Vjj*s*s+2*Vij*s*c;
V[jm][jm]=Vii*s*s+Vjj*c*c-2*Vij*s*c;
V[im][jm]=V[jm][im]=0.0; // setting rotated element zero

for(int k=0;k<n;k++){
if(k!=im && k!=jm){
    double x=V[k][im]; 
    double y=V[k][jm]; 
    V[k][im]=x*c+y*s;
    V[im][k]=V[k][im];
    V[k][jm]=-x*s+y*c;
    V[jm][k]=V[k][jm];
}
//update E: the matrix with eigen vectors
double a1=E[k][im];
double a2=E[k][jm];
E[k][im]=a1*c+a2*s;
E[k][jm]=-a1*s+a2*c;
    }
} 
} 
if(rotations_done==0){          //convergence condition if this is satisfied means matrix is diagonalized
    break;
}
} 
}





