#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double esm[6][6],x[3],y[3],d[3][3],b[3][6],ar2;
int np,nbw,jgf,jgsm,jend;

/*! esm[6][6]—单元刚度矩阵， x[3],y[3]—单元节点坐标， d[3][3]—材料性质矩阵
! b[3][6]—几何矩阵， ar2—三角形面积的二倍， np—自由度总数， Nnbw—最大半带宽
! A(8500)—存储节点位移向量、节点力向量和总体刚度矩阵的数组 A
! jgf、 jgsm、 jend 为计数单元
! jgf=np－节点位移向量在数组 A 中的位置 ， jgsm=jgf+np－节点力向量在数组 A 中的位置，
! jendgit=jgsm+np*nbw—刚度矩阵在数组 A 中的位置，数组 A 总长度
! NS(6)—一个单元的节点自由度编号数组， U(6)—一个单元的节点自由度
! IN/60/,IO/61/—输入输出文件设备号， STRA(3),STRE(3)—存储单元的应变、应力*/
int main(){



}


