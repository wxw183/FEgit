/*coding=utf-8*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

double esm[6][6],x[3],y[3],**d;// esm[6][6]—单元刚度矩阵， x[3],y[3]—单元节点坐标， d[3][3]—材料性质矩阵
double b[3][6],ar2;//b[3][6]—几何矩阵， ar2—三角形面积的二倍
double stra[3],stre[3];// a[8500]—存储节点位移向量、节点力向量和总体刚度矩阵的数组 a，STRA(3),STRE(3)—存储单元的应变、应力*/
int np,nbw,jgf,jgsm,jend;// np—自由度总数， nbw—最大半带宽， jgf、 jgsm、 jend 为计数单元，jgf=np－节点位移向量在数组 a 中的位置 ， jgsm=jgf+np－节点力向量在数组 a 中的位置，jend=jgsm+np*nbw—刚度矩阵在数组 a 中的位置，数组 a 总长度
int ns[6],u[6];//ns[6]—一个单元的节点自由度编号数组，从1开始,例如第一个节点自由度编号为1，2，第4个节点自由度编号为7，8， u[6]—一个单元的节点自由度
//double em=200000000000.0,pr=0.3,th=0.001;//EM－杨氏模量， PR－泊松比， TH－板的厚度
double *a;//存储单元刚度矩阵，由于矩阵是对称带状矩阵，因此存储方式如下：比如最大半带宽是3，5阶矩阵
/*
0   5   9
5   1   6   10
9   6   2   7   11
    10  7   3   8
        11  8   4
数字为数组中括号里面的下标，从0开始*/

void modify();
void elstmx(int kk);
void dcmpbd();
double *k_get(int ii,int jj);
double *d_get(int ii);
double *r_get(int ii);
double **material_set(double em,double pr);

int main(){
    char title[100];//存储计算内容标题的字符数组
    int i,n;//N－单元编号
    int bn=0;//边界节点个数
    
    
    int nn,ne;//NN－节点总数， NE－单元总数 ，

    FILE *fpi,*fpo,*fpd,*fpr;//fpi指向输入文件，fpo指向输出文件,fpd指向输入位移载荷文件，fpr指向力载荷文件
    fpi=fopen("00.txt","r+");//input.txt为comsol导出的网格文件（域内的，非边界）
    fpd=fopen("dc.txt","w");
    fpr=fopen("bl.txt","w");
    if(fpi==NULL){
        printf("无输入文件\n");
        exit(0);
    }
//跳过网格文件前4行，开始读取nodes和elements
    fseek(fpi,206,0);
    fscanf(fpi,"%d",&nn);
    fseek(fpi,252,0);
    fscanf(fpi,"%d",&ne);
    
    fseek(fpi,429,0);
    fpo=fopen("outputnet.txt","w");
    np=2*nn;
    i=nn;
    double *xc,*yc;
    xc=(double *)malloc((i+ne)*sizeof(double));
    yc=(double *)malloc((i+ne)*sizeof(double));
    //double xc[i],yc[i];// XC(I)－节点的 X 轴的坐标，YC(I)－节点的 Y 轴的坐标
    double *xc_nel,*yc_nel;//单元中点的坐标
    double x_nel,y_nel;//本单元的中点坐标
    int n_nel=0;//要加密的单元数
    int *n_nel_list,*pointer;//要加密的单元号码列表 从0开始
    n_nel_list=(int *)malloc(ne*sizeof(int));
    int nel_range=(1);//加密范围
    int **nel;
    nel=(int **)malloc(3*ne*sizeof(int *));
    for(i=0;i<3*ne;i++){
        nel[i]=(int *)malloc(3*sizeof(int));
    }
    //int nel[ne][3];//NEL(N,I) －组成第 N 个三角形单元的第 I 节点的编号，从1开始而不是从0开始（ I=1,2,3）
    /*输入材料的杨氏模量 EM，波松比 PR，平板厚度 TH，节点坐标 XC(I)，YC(I)和组成单元的节点 NEL(N,I)
组成单元的节点的编号都按逆时针顺序输入*/
    for(int k=0;k<nn;k++){
        fscanf(fpi,"%lg%lg",&xc[k],&yc[k]);
    }//输入节点坐标
    
    //fprintf(fpo,"输入数据为：\n节点数=%d\t单元数=%d\n",nn,ne);
//输出材料性质和计算模型拓扑数据便于检查时对照
    //fprintf(fpo,"材料常数为：%g\t泊松比为：%g\t厚度为：%g\n",em,pr,th);
    for(int k=0;k<nn;k++){
        //fprintf(fpo,"节点%d坐标为:X=%f\tY=%f\n",1+k,xc[k],yc[k]);
        fflush(fpo);
    }
    while(fgetc(fpi)==' ')fseek(fpi,1L,1);
    fseek(fpi,24L,1);
    //找到要加密的单元号
    for(int k=0;k<ne;k++){
        fscanf(fpi,"%d%d%d",&nel[k][0],&nel[k][1],&nel[k][2]);
        //fprintf(fpo,"单元号码为：%d\t组成单元的节点号码为:%d\t%d\t%d\n",k+1,nel[k][0],nel[k][1],nel[k][2]);
        fflush(fpo);
        x_nel=0;
        y_nel=0;
        for(i=0;i<3;i++){
            x[i]=xc[nel[k][i]-1];
            x_nel+=x[i];
            y[i]=yc[nel[k][i]-1];
            y_nel+=y[i];
        }
        x_nel/=3.0;
        y_nel/=3.0;
        if(x_nel>0.5){//括号内输入要加密的范围
            n_nel_list[n_nel]=k;
            n_nel++;
            
            
        }
        
    }//节点编号从1开始
    //nel=(int **)realloc(nel,(ne+2*n_nel)*sizeof(int *));
    //yc=(double *)realloc(yc,(nn+n_nel)*sizeof(double));
    //xc=(double *)realloc(xc,(nn+n_nel)*sizeof(double));
    //for(i=0;i<n_nel;i++)nel[ne+i]=(int *)malloc(3*sizeof(int));
    //加密单元
    {
        for(int k=0;k<n_nel;k++){
            x_nel=0;
            y_nel=0;
            for(i=0;i<3;i++){
                x[i]=xc[nel[n_nel_list[k]][i]-1];
                x_nel+=x[i];
                y[i]=yc[nel[n_nel_list[k]][i]-1];
                y_nel+=y[i];
            }
            x_nel/=3.0;
            y_nel/=3.0;
            if(1){//括号内输入要加密的范围
                xc[nn+k]=x_nel;
                yc[nn+k]=y_nel;
                
                nel[ne+k*2][0]=nel[n_nel_list[k]][1];
                nel[ne+k*2][1]=nel[n_nel_list[k]][2];
                nel[ne+k*2][2]=nn+1+k;
                nel[ne+k*2+1][0]=nel[n_nel_list[k]][2];
                nel[ne+k*2+1][1]=nel[n_nel_list[k]][0];
                nel[ne+k*2+1][2]=nn+1+k;
                nel[n_nel_list[k]][2]=nn+1+k;
            }
        
        }
    }
    nn+=n_nel;
    ne+=2*n_nel;
    fseek(fpi,206,0);
    fprintf(fpi,"%d",nn);
    fseek(fpi,252,0);
    fprintf(fpi,"%d",ne);
    fseek(fpi,429,0);
    for(int k=0;k<nn;k++){
        fprintf(fpi,"%-24lg %-25lg\n",xc[k],yc[k]);
    }
    fprintf(fpi,"%% Elements (triangles)\n");
    for(int k=0;k<ne;k++){
        fprintf(fpi,"%-10d%-10d%-10d\n",nel[k][0],nel[k][1],nel[k][2]);
    }

    
    fclose(fpi);

}