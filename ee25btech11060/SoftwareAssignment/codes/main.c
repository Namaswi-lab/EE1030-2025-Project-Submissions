#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const int rows=1024;

const int coloumns=1024;




  

double**pixels(int m,int n);
void free_mem(double**ptr,int m);
void transpose(double**ptr,double**ptrT,int m,int n);
void mult(double **A,double**B,double**C,int m,int n,int p);
int pgm_reading(const char*fname,double ***image,int *m,int*n);
void pgm_writing(const char*fname,double **ptr,int m,int n);
void scale_pixel(double **ptr,int m,int n);
void ATA(double**ptr,double**ptrTptr,int m,int n);
void top_k_eigen(double **ptrTptr,double**V,double*S,int n,int k);
void power_iteration(double **ptr,double*eigenvector,double*eigenvalue,int n,int max_iter);
void reconstructed(double**ptr,int m,int n,double**V,double*S,int k,const char*base);
double frobenius_error(double**ptr,double**Ak,int m,int n);


//allocating memory for each number and initilising it to 0
  double **pixels(int m,int n){
        double **ptr=malloc(m*sizeof(double*));
        for(int i=0;i<m;i++){
            ptr[i]=calloc(n,sizeof(double));
        }
        return ptr;
    }


 
void transpose(double**ptr,double**ptrT,int m,int n){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            ptrT[j][i]=ptr[i][j];
    }
}
}



 

void mult(double**A,double**B,double**C,int m,int n,int p){
    for(int i=0;i<m;i++){
        for(int j=0;j<p;j++){
            C[i][j]=0;
            for(int k=0;k<n;k++){
                C[i][j]+=A[i][k]*B[k][j];
            }
        }
    }
}
//function to free memory
void free_mem(double**ptr,int m){
    for(int i=0;i<m;i++){
        free(ptr[i]);
    }
free(ptr);
}


//function to read a pgm file and take pixel values from it

int pgm_reading(const char*fname,double***image,int*m,int*n){
    FILE *f = fopen(fname,"r");
    if(!f){
return 0;
    } 

    char type[3];
    (void)fscanf(f,"%2s",type);
    if(strcmp(type,"P2")!=0){
         fclose(f); 
         return 0;
         }

    int maxval;
    (void)fscanf(f,"%d %d %d", n, m, &maxval);

    *image = pixels(*m,*n);

    for(int i=0;i<*m;i++){ 
        for(int j=0;j<*n;j++){ 
            (void)fscanf(f,"%lf",&((*image)[i][j]));
        }
    }
    fclose(f);
    return 1;
}




//storing those pixel values in a rray
void pgm_writing(const char*fname,double**ptr,int m,int n){
    FILE *f=fopen(fname,"w");
    fprintf(f,"P2\n%d %d\n255\n",n,m);
    for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
            int v=(int)fmin(fmax(ptr[i][j],0.0),255.0);
            fprintf(f,"%d ",v);
        }
fprintf(f,"\n");
    }
    fclose(f); 
}



void scale_pixel(double**ptr,int m,int n){
 double min = ptr[0][0], max = ptr[0][0];
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++){
            if(ptr[i][j]<min){
min = ptr[i][j];
            } 
            if(ptr[i][j]>max){
max = ptr[i][j];
            } 
        }
    double range = max-min;
    if(range <1e-6) {
        range=1.0;
    }
    for(int i=0;i<m;i++){ 
        for(int j=0;j<n;j++){ 
            ptr[i][j] = 255.0*(ptr[i][j]-min)/range;
        }
    }
}

void ATA(double**ptr,double**ptrTptr,int m,int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            ptrTptr[i][j]=0;
            for(int k=0;k<m;k++){ 
                ptrTptr[i][j]+=ptr[k][i]*ptr[k][j];
            }
        }
    }
}


//finding eigenvalues and vectors
 void power_iteration(double **ptr,double *eigenvector,double *eigenvalue,int n,int max_iter){
    for(int i=0;i<n;i++){
        eigenvector[i]=1.0;
    }
    for(int iter=0;iter<max_iter;iter++){
        double *y=calloc(n,sizeof(double));
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                y[i]+=ptr[i][j]*eigenvector[j];
            }
        }
        double norm=0;
        for(int i=0;i<n;i++){
            norm+=y[i]*y[i];
        }
        norm=sqrt(norm);

        for(int i=0;i<n;i++){
            eigenvector[i]=y[i]/norm;
        }
        free(y);
    }

    double lambda=0;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            lambda+=eigenvector[i]*ptr[i][j]*eigenvector[j];
        }
    }
            *eigenvalue=lambda;
        
    
}


void top_k_eigen(double **ptrTptr,double **V,double *S,int n,int k){
    double **B=pixels(n,n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            B[i][j]=ptrTptr[i][j];
        }
    }
    for(int k_idx=0;k_idx<k;k_idx++){
        double *v= calloc(n,sizeof(double));
        double lambda;
        power_iteration(B,v,&lambda,n,100);
        S[k_idx]=sqrt(lambda);

    for(int i=0;i<n;i++){
        V[i][k_idx]=v[i];
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            B[i][j]-=lambda*v[i]*v[j];
        }
    }
    
    free(v);
    }
    free_mem(B,n);
}

//forming reconstructed matrix by taking top k valuess
void reconstructed(double **ptr,int m,int n,double **V,double *S,int k,const char *base){
    double **U=pixels(m,k);

    for(int i=0;i<k;i++){
        for(int j=0;j<m;j++){
            U[j][i]=0;
            for(int l=0;l<n;l++){
                U[j][i]+=ptr[j][l]*V[l][i];
            }
            U[j][i]/=S[i];
        }
    }

    double **VkT=pixels(k,n);

    for(int i=0;i<k;i++){
        for(int j=0;j<n;j++){
            VkT[i][j]=V[j][i];
        }
    }

    double **sigma_U=pixels(m,k);
    for(int i=0;i<m;i++){
        for(int j=0;j<k;j++){
            sigma_U[i][j]=U[i][j]*S[j];
        }
    }

    double **Ak=pixels(m,n);
    mult(sigma_U,VkT,Ak,m,k,n);

    scale_pixel(Ak,m,n);

    double error_value=frobenius_error(ptr,Ak,m,n);

char out[256];
sprintf(out, "../figs/%s_k%d.pgm", base, k);
pgm_writing(out,Ak,m,n);


    printf("For k=%d,Forbenius norm error=%lf\n",k,error_value);

    free_mem(U,m);
    free_mem(VkT,k);
    free_mem(sigma_U,m);
    free_mem(Ak,m);
}
//calclating frobnius error 
double frobenius_error(double **ptr,double **Ak,int m,int n){
    double sum=0.0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            double diff=ptr[i][j]-Ak[i][j];
            sum+=diff*diff;
        }
    }

    return sqrt(sum);
}

int main(int argc,char *argv[]){
    if(argc<2){
        return 1;
    }

    const char *filename=argv[1];
    int im_rows,im_coloumns;
    double **ptr;

    if(!pgm_reading(filename,&ptr,&im_rows,&im_coloumns)){
        return 1;
    }

    char base[128];
    const char *slash=strrchr(filename,'/');
    if(slash){
        strcpy(base,slash+1);
    }else{
        strcpy(base,filename);
    }

    char *dot=strrchr(base,'.');
    if(dot){
        *dot=0;
    }

    printf("Loaded image: %dx%d - %s\n",im_rows,im_coloumns,base);

    double **ptrTptr = pixels(im_coloumns,im_coloumns);
    ATA(ptr,ptrTptr,im_rows,im_coloumns);

    double **V=pixels(im_coloumns,100);
    double *S=calloc(100,sizeof(double));

    top_k_eigen(ptrTptr,V,S,im_coloumns,100);

    int ks[]={5,20,50,100};
    for(int i=0;i<4;i++){
        reconstructed(ptr,im_rows,im_coloumns,V,S,ks[i],base);
        printf("reconstructed image with k=%d\n",ks[i]);
    }


    free_mem(ptr,im_rows);
    free_mem(ptrTptr,im_coloumns);
    free_mem(V,im_coloumns);
    free(S);

    return 0;
}
