#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jacobi.h"
#include "matrix.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"          //for loading JPG/PNG images

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"            //for saving PNG images(output)
int main(){
int m,n,k;
int channels;
char image[50];
fgets(image,50,stdin);
int l=strlen(image);
if(l>0 && image[l-1]=='\n'){
image[l-1]='\0';
}
unsigned char *img=stbi_load(image,&n,&m,&channels,1);          //MAIN IMAGE LOADING
// n = width, m = height
if(img==NULL){
printf(" Error: Could not load image");
    return 0;
}
printf("Enter rank k: ");
scanf("%d",&k);
if(k<0){         //restricting the rank 'k' to be between 0 and min(m,n)
k=0;
}
int mn=(m<n?m:n);
if(k>mn){
k=mn;
}
double **A=allocatememory(m,n);
double **B=allocatememory(n,m);          //for At
double **V=allocatememory(n,n);          //will store the eigen values as a diagonal matrix with diag entries as the eigen val
double **E=allocatememory(n,n);          //E matrix is the V matrix in theory formulas
double **U=allocatememory(m,n);
double *singular=(double *)malloc(n * sizeof(double));
double **sigma=allocatememory(n,n);          //diagonal matrix with singular values as entries 
double **AV=allocatememory(m,n);
double **Ak=allocatememory(m,n);
double **Uk=allocatememory(m,k);
double **sigmak=allocatememory(k,k);
double **Vk=allocatememory(n,k);
double **temp=allocatememory(m,k);
double **temp1=allocatememory(k,n);
for(int i=0;i<m;i++){
for(int j=0;j<n;j++){
    A[i][j]=img[(size_t)i*n+j];               //extracting the pixels from the image and storing them to A
}
}
stbi_image_free(img);               //since pixels r now stored we can free the memory used by img array 
for(int i=0;i<n;i++){                                   
for(int j=0;j<n;j++){
if(i==j){                            //initialising E to identity matrix
    E[i][j] = 1;
}
else{
    E[i][j] = 0;
}
}
}
int issymmetric = (m == n);
if(issymmetric){                                               //matrix symmetry check
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
if(fabs(A[i][j] - A[j][i])>1e-12){
    issymmetric=0;
    break;
}
}
if(issymmetric==0){
    break;
}
}
}
if(issymmetric==1){
    for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
    V[i][j]=A[i][j];      //copy A into V
}
}
jacobi(V,E,n);             //if A is symmetric apply jacobi directly on A
}
else{
transpose(A,B,m,n);
mult(B,A,n,m,n,V);                 //if A not symmetric apply jacobi on AtA
jacobi(V,E,n);
}
singularvalues(n,singular,V,issymmetric);      //calculate singular values from eigen values
sortsigma_andE(singular,E,n);                     //sorting singular values in decreasing order and also swapping columns of E so that singular values match with corresponding column of E
for(int i=0;i<n;i++){
for(int j=0;j<n;j++){
    sigma[i][j]=0.0;
}                                                   
    sigma[i][i]=singular[i];
}
mult(A,E,m,n,n,AV);
matrixU(A,E,U,singular,m,n,AV);           //computing U matrix
//image reconstruction
calcUkandVk(U,Uk,m,k);
calcUkandVk(E,Vk,n,k);
calcsigmak(singular,sigmak,k);
transpose(Vk,temp1,n,k);
mult(Uk,sigmak,m,k,k,temp);
mult(temp,temp1,m,k,n,Ak);

// display(Ak,m,n);
unsigned char *out=(unsigned char *)malloc(m * n * sizeof(unsigned char));        //1D array to store reconstructed Ak 
for(int i=0;i<m;i++){
for(int j=0;j<n;j++){
int px=(int)round(Ak[i][j]);      //rounding off double values to whole numbers
if(px<0){            //restricting pixels between 0 and 255
    px=0;
}
if(px>255){
    px=255;
}
out[i*n+j]=(unsigned char)px;
    }
}

stbi_write_png("output.png",n,m,1,out,n);

free(out);
printf("Compression done\n");
printf("Output saved as output.png\n");
printf("Frobenius error is %lf\n",frobeniuserror(A,Ak,m,n));
printf("Theoretical error is %lf\n",svderror(singular,k,n));
freememory(A, m);
freememory(B, n);
freememory(V, n);
freememory(E, n);
freememory(sigma, n);
freememory(AV, m);
free(singular);
freememory(Ak, m);
freememory(Uk, m);
freememory(sigmak, k);
freememory(Vk, n);
freememory(temp, m);
freememory(temp1, k);
freememory(U,m);
return 0;
}
