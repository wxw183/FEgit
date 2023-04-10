#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double esm[6][6],x[3],y[3],d[3][3];
double b[3][6],ar2;
double a[8500],stra[3],stre[3];
int np,nbw,jgf,jgsm,jend;
int ns[6],u[6];
double em,pr,th;


/*! esm[6][6]—单元刚度矩阵， x[3],y[3]—单元节点坐标， d[3][3]—材料性质矩阵
! b[3][6]—几何矩阵， ar2—三角形面积的二倍， np—自由度总数， nbw—最大半带宽
! a[8500]—存储节点位移向量、节点力向量和总体刚度矩阵的数组 a
! jgf、 jgsm、 jend 为计数单元
! jgf=np－节点位移向量在数组 a 中的位置 ， jgsm=jgf+np－节点力向量在数组 a 中的位置，
! jend=jgsm+np*nbw—刚度矩阵在数组 a 中的位置，数组 a 总长度
! ns[6]—一个单元的节点自由度编号数组， u[6]—一个单元的节点自由度
! IN/60/,IO/61/—输入输出文件设备号， STRA(3),STRE(3)—存储单元的应变、应力*/
int main(){
    char title[35];
    int i,n;
    double xc[i],yc[i];
    int nel[n][i];
    int nn,ne;
    
    FILE *fpi,*fpo;//fpi指向输入文件，fpo指向输出文件
    fpi=fopen("INPUT.DAT","r+");
    if(fp=NULL){
        printf("找不到输入文件\n");
        return(0);
    }
    fpo=fopen("OUTPUT.DAT","w+");
    fgets(title,40,fpi);
    fscanf(fpi,"%d%d",&nn,&ne);
    fprintf("输入数据为：\n节点数=%d\t单元数=%d\n");
    np=2*nn;
    

}


