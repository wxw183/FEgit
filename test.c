#include<stdio.h>
#include<stdlib.h>

int np=4,nbw=2,jgf=4,jgsm=8,jend=8+2*4;
void dcmpbd(double *a){
    int jj,ii,j1,k,m,n;
    double aa,mk,nk,kk,bb;
    j1=jgsm+(jj-ii)*np+ii-(jj-ii-1)*(jj-ii)/2;
    for(k=1;k<np;k++){
        j1=jgsm+k;
        kk=a[j1];//提取K_kk
        for(m=k;m<np;m++){
            j1=jgsm+(m-k)*np+k-(m-k-1)*(m-k)/2;
            mk=a[j1];//提取K_km
            
            
            bb=a[jgf+m];//提取b_m
            a[jgf+m]-=bb*mk/kk;
            for(n=m;(n<np)&&(n<m+nbw);n++){
                j1=jgsm+(n-k)*np+k-(n-k-1)*(n-k)/2;
                nk=a[j1];//提取K_kn
                
                j1=jgsm+(n-m)*np+m-(n-m-1)*(n-m)/2;
                a[j1]-=mk*nk/kk;
                
            }
        }
    }
}
/*void slvbd(double *a){
    int jj,ii,j1;
    double aa;
    j1=jgsm+(jj-ii)*np+ii-(jj-ii-1)*(jj-ii)/2;
    for(ii=1;ii<np;ii++){
        ;
    }
}*/
int main(){
    double a[16]={0,0,0,0,1,2,3,4,1,1,1,1,2,2,2};
    dcmpbd(a);
    return(0);
}
