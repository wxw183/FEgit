#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double esm[6][6],x[3],y[3],d[3][3];// esm[6][6]—单元刚度矩阵， x[3],y[3]—单元节点坐标， d[3][3]—材料性质矩阵
double b[3][6],ar2;//b[3][6]—几何矩阵， ar2—三角形面积的二倍
double a[8500],stra[3],stre[3];// a[8500]—存储节点位移向量、节点力向量和总体刚度矩阵的数组 a，STRA(3),STRE(3)—存储单元的应变、应力*/
int np,nbw,jgf,jgsm,jend;// np—自由度总数， nbw—最大半带宽， jgf、 jgsm、 jend 为计数单元，jgf=np－节点位移向量在数组 a 中的位置 ， jgsm=jgf+np－节点力向量在数组 a 中的位置，jend=jgsm+np*nbw—刚度矩阵在数组 a 中的位置，数组 a 总长度
int ns[6],u[6];//ns[6]—一个单元的节点自由度编号数组， u[6]—一个单元的节点自由度
double em,pr,th;//EM－杨氏模量， PR－泊松比， TH－板的厚度
 
int main(){
    char title[100];//存储计算内容标题的字符数组
    int i,n;//N－单元编号
    
    
    int nn,ne;//NN－节点总数， NE－单元总数 ，

    FILE *fpi,*fpo;//fpi指向输入文件，fpo指向输出文件
    fpi=fopen("INPUT.DAT","r+");
    if(fpi==NULL){
        printf("无输入文件\n");
        exit(0);
    }
    fpo=fopen("OUTPUT.DAT","w+");
    fgets(title,100,fpi);
    fputs(title,fpo);
    fscanf(fpi,"%d\n%d\n",&nn,&ne);
    fprintf(fpo,"输入数据为：\n节点数=%d\t单元数=%d\n",nn,ne);
    np=2*nn;
    i=nn;

    double xc[i],yc[i];// XC(I)－节点的 X 轴的坐标， YC(I)－节点的 Y 轴的坐标
    int nel[ne][3];//NEL(N,I) －组成第 N 个三角形单元的第 I 节点的编号（ I=1,2,3）
    fscanf(fpi,"%lg\n%lg\n%lg\n",&em,&pr,&th);
    for(int k=0;k<nn;k++){
        fscanf(fpi,"%lg",&xc[k]);
    }
    for(int k=0;k<nn;k++){
        fscanf(fpi,"%lg",&yc[k]);
    }

    fprintf(fpo,"\n材料常数为：%g\t泊松比为：%g\t厚度为：%g\n",em,pr,th);
    for(int k=0;k<nn;k++){
        fprintf(fpo,"\n节点%d坐标为:X=%fY=%f\n",k+1,xc[k],yc[k]);
    }
    for(int k=0;k<ne;k++){
        fscanf(fpi,"%d%d%d%d",&n,&nel[k][0],&nel[k][1],&nel[k][2]);
        fprintf(fpo,"\n单元号码为：%d\t组成单元的节点号码为:%d\t%d\t%d\n",n,nel[k][0],nel[k][1],nel[k][2]);
    }
    fclose(fpi);
    fclose(fpo);
   return(0);
}


