#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
double a[22];
double esm[6][6],x[3],y[3],**d;// esm[6][6]—单元刚度矩阵， x[3],y[3]—单元节点坐标， d[3][3]—材料性质矩阵
double b[3][6],ar2;//b[3][6]—几何矩阵， ar2—三角形面积的二倍
double stra[3],stre[3];// a[8500]—存储节点位移向量、节点力向量和总体刚度矩阵的数组 a，STRA(3),STRE(3)—存储单元的应变、应力*/
int np=5,nbw=3,jgf=5,jgsm=10,jend=24;

double *k_get(int ii,int jj){
    int j1;
    if(ii<=jj)j1=jgsm+abs(jj-ii)*np+ii-1-(abs(jj-ii)-1)*abs(jj-ii)/2;
    else j1=jgsm+abs(jj-ii)*np+jj-1-(abs(jj-ii)-1)*abs(jj-ii)/2;
    return &a[j1];
}
//位移向量元素获取函数，例如：d_get(3)返回位移向量的第三分量（从1开始）
double *d_get(int ii){
    return &a[ii-1];
}
//力向量元素获取函数，例如：r_get(3)返回位移向量的第三分量（从1开始）
double *r_get(int ii){
    return &a[jgf+ii-1];
}

void dcmpbd(){
    int j,i,j1,k,m;
    double aa,mk,nk,kk,bb;
    for(int s=1;s<=np;s++){
        //往下打洞
        for(i=s+1;i<=np&&i<=s+nbw-1;i++){
            //r向量打洞
            *r_get(i)-=*r_get(s)*(*k_get(s,i))/(*k_get(s,s));
            //矩阵打洞
            for(j=i;j<=np&&j<=s+nbw-1;j++){
                *k_get(j,i)-=*k_get(s,j)*(*k_get(s,i))/(*k_get(s,s));
            }
        }
        //第一行第一个元素变成1
        *r_get(s)/=*k_get(s,s);
        for((s+nbw-1)<=np?(j=s+nbw-1):(j=np);j>=s;j--){
            *k_get(s,j)/=*k_get(s,s);
        }
    }
    
    //开始回代
    for(int i=np-1;i>=1;i--){
        for(int j=i+1;j<=i+nbw-1&&j<=np;j++){
            *r_get(i)-=*k_get(j,i)*(*r_get(j));
        }
    }
    FILE *fpp;//探针输出文件
    fpp=fopen("probe1.txt","w");
    for(int i=1;i<=np;i++)fprintf(fpp,"%lg\t",*d_get(i));
    fprintf(fpp,"\n");
    fflush(fpp);
    for(int i=1;i<=np;i++)fprintf(fpp,"%lg\t",*r_get(i));
    fprintf(fpp,"\n");
    fflush(fpp);
    for(int i=1;i<=np;i++){
        for(int j=1;j<=np;j++){
            fprintf(fpp,"%lg\t",*k_get(i,j));
        }
        fprintf(fpp,"\n");
        fflush(fpp);
    }
    fclose(fpp);

}

int main(){
    //a=(double *)malloc(jend*sizeof(double));
    a[0]=0;
    a[1]=0;
    a[2]=0;
    a[3]=0;
    a[4]=0;
    a[5]=38;
    a[6]=53;
    a[7]=67;
    a[8]=47;
    a[9]=32;
    a[10]=5;
    a[11]=4;
    a[12]=3;
    a[13]=2;
    a[14]=1;
    a[15]=6;
    a[16]=5;
    a[17]=4;
    a[18]=3;
    a[19]=7;
    a[20]=6;
    a[21]=5;

    dcmpbd();
    
    return 0;
}



