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
double em=200000000000.0,pr=0.3,th=1;//EM－杨氏模量， PR－泊松比， TH－板的厚度
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
    
    
    int nn,ne;//NN－节点总数， NE－单元总数 ，

    FILE *fpi,*fpo,*fpd,*fpr;//fpi指向输入文件，fpo指向输出文件,fpd指向输入位移载荷文件，fpr指向力载荷文件
    fpi=fopen("1.txt","r");//input.txt为comsol导出的网格文件（域内的，非边界）
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
    fpo=fopen("output.csv","w");
    np=2*nn;
    i=nn;
    double xc[i],yc[i];// XC(I)－节点的 X 轴的坐标，YC(I)－节点的 Y 轴的坐标
    int nel[ne][3];//NEL(N,I) －组成第 N 个三角形单元的第 I 节点的编号，从1开始而不是从0开始（ I=1,2,3）
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
    for(int k=0;k<ne;k++){
        fscanf(fpi,"%d%d%d",&nel[k][0],&nel[k][1],&nel[k][2]);
        //fprintf(fpo,"单元号码为：%d\t组成单元的节点号码为:%d\t%d\t%d\n",k+1,nel[k][0],nel[k][1],nel[k][2]);
        fflush(fpo);
    }//节点编号从1开始
    fclose(fpi);
    /*计算最大半带宽， B=MAXe(De+1)*F
! De 是一个单元各节编点号之差的最大值， F 是一个节点的自由度数*/
    int inbw=0;
    nbw=0;
    int j,k,m,l,nb;
    for(k=0;k<ne;k++){
        for(i=0;i<3;i++)ns[i]=nel[k][i];
        for(i=0;i<2;i++){
            for(j=i+1;j<3;j++){
                nb=abs(ns[i]-ns[j]);
                if(nb<=nbw)continue;
                inbw=k;
                nbw=nb;
            }
        }

    }
    nbw=(nbw+1)*2;
    
    //生成材料特性矩阵d
    d=material_set(em,pr);
    
    jgf=np;
    jgsm=jgf+np;
    jend=jgsm+np*nbw;
    a=(double *)malloc(jend*sizeof(double));//先存位移再存力最后存刚度阵
    int jl=jend-jgf;
    for(i=0;i<jend;i++)a[i]=0.0;
    
    //单元矩阵循环的开始
    k=0;//k为单元号，从0开始
    while(k<ne){
    for(i=0;i<3;i++){
        j=nel[k][i];//提取k单元的第i个节点的节点号
        ns[2*i]=j*2-1;
        ns[2*i+1]=j*2;
        x[i]=xc[j-1];
        y[i]=yc[j-1];
    }
    //探针
    printf("element%d:node%d:x=%lg,y=%lg\n",k+1,i+1,x[i],y[i]);
    

    elstmx(k);
    
    
    //单元刚度矩阵组装成总体刚度矩阵
    int ii,jj,j1;
    for(i=0;i<6;i++){
        ii=ns[i];
        for(j=0;j<6;j++){
            jj=ns[j];
            if(jj<ii)continue;
            *k_get(ii,jj)+=esm[i][j];
        }
    }
    //探针
    /*for(i=0;i<jend;i++){
        if(i%5==0)printf("\n");
        printf("%lg\t",a[i]);
    }*/
    k++;
    }//单元矩阵循环结束

    //调用子程序 MODIFY 输入载荷节点处的载荷值、位移边界节点处的位移值 ,对总体刚度矩阵、位移数组和节点力数组进行相应的修改
    /*FILE *fpp;//探针输出文件
    fpp=fopen("probe.txt","w");
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
    }*/
    modify();

    //fprintf(fpo,"计算结果为：\n");
    //fprintf(fpo,"最大半带宽为%d\n",nbw);
    //fprintf(fpo,"总数组大小为%d\n",jend);
    dcmpbd(a);
    fprintf(fpo,"节点,X坐标,Y坐标,U位移,V位移\n");
    for(int i=1;i<=np;i+=2){
        fprintf(fpo,"%d,%lg,%lg,%lg,%lg\n",(i+1)/2,xc[(i+1)/2-1],yc[(i+1)/2-1],*r_get(i),*r_get(i+1));
    }
    fflush(fpo);
    int f=0;
    f=fclose(fpo);
    printf("\n\n\n\n%d",f);
    return(0);
}

void elstmx(int kk){
    double c[6][3];
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++)b[i][j]=0.0;
    }
    b[0][0]=y[1]-y[2];
    b[0][2]=y[2]-y[0];
    b[0][4]=y[0]-y[1];
    b[1][1]=x[2]-x[1];
    b[1][3]=x[0]-x[2];
    b[1][5]=x[1]-x[0];
    b[2][0]=b[1][1];
    b[2][1]=b[0][0];
    b[2][2]=b[1][3];
    b[2][3]=b[0][2];
    b[2][4]=b[1][5];
    b[2][5]=b[0][4];
    ar2=x[1]*y[2]+x[2]*y[0]+x[0]*y[1]-x[1]*y[0]-x[2]*y[1]-x[0]*y[2];
    //计算矩阵C
    for(int i=0;i<6;i++){
        for(int j=0;j<3;j++){
            c[i][j]=0.0;
            for(int k=0;k<3;k++)c[i][j]+=b[k][i]*d[k][j];
        }
    }
    // 计算矩阵 ESM=[BT][D][B]=[C][B]
    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            double sum=0.0;
            for(int k=0;k<3;k++){
                sum+=c[i][k]*b[k][j];
                esm[i][j]=sum*th*(0.5*ar2);
            }
        }
    }

}

void modify(){
    /*
    输入文件fpd结构：
    自由度序号，载荷\n

    例如：2，3
    表示1号节点Y方向位移约束为3
    输入文件fpr结构类似同上
    */
   FILE *fpd,*fpr;
    double bv;//节点载荷
    int ib;//自由度序号,从1开始
    int i;
    fpd=fopen("dc.txt","r");//dc.txt为位移约束文件
    fpr=fopen("bl.txt","r");//bl.txt为边界载荷文件
    
    //将读入的位移、力载荷存入数组a中
    
    while(feof(fpr)==0){
        fscanf(fpr,"%d,%lg",&ib,&bv);
        *r_get(ib)=bv;
    }
    while(feof(fpd)==0){
        fscanf(fpd,"%d,%lg",&ib,&bv);
        *d_get(ib)=bv;
        //修改k矩阵和r向量
        //k(ib,ib)所在行列除了k(ib,ib)本身全部置零
        for((ib-nbw+1)>1?(i=ib-nbw+1):(i=1);i<=ib+nbw-1&&i<=np;i++){
            if(i!=ib){
                *k_get(ib,i)=0;
            }
        }
        //修改r向量
        //(s+nbw-1)<=np?(j=s+nbw-1):(j=np);j>=s;j--
        for((ib-nbw+1)>1?(i=ib-nbw+1):(i=1);i<=ib+nbw-1&&i<=np;i++){
            if(i!=ib){
                *r_get(i)-=(*k_get(ib,i))*bv;
            }
            else *r_get(i)=(*k_get(i,i))*bv;
        }
    }
    //探针
    
    
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
    //fpp=fopen("probe1.txt","w");
    for(int i=1;i<=np;i++)fprintf(stdout,"%lg\t",*d_get(i));
    fprintf(stdout,"\n");
    //fflush(fpp);
    for(int i=1;i<=np;i++)fprintf(stdout,"%lg\t",*r_get(i));
    fprintf(stdout,"\n");
    //fflush(fpp);
    for(int i=1;i<=np;i++){
        for(int j=1;j<=np;j++){
            //fprintf(stdout,"%lg\t",*k_get(i,j));
        }
        fprintf(stdout,"\n");
        //fflush(fpp);
    }
    //fclose(fpp);

}

///刚度阵元素获取函数，矩阵行列号ii和jj从1开始编
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

//材料特性矩阵生成函数(平面应力问题)
double **material_set(double em,double pr){//em为杨氏模量，pr为泊松比
    double **pointer;
    pointer=(double **)malloc(3*sizeof(double *));
    for(int i=0;i<3;i++)pointer[i]=(double *)malloc(3*sizeof(double));

    double r=em/(1.0-pr*pr);
    pointer[0][0]=r;
    pointer[1][1]=r;
    pointer[2][2]=r*(1.0-pr)/2.0;
    pointer[0][1]=pr*r;
    pointer[1][0]=pointer[0][1];
    pointer[0][2]=0.0;
    pointer[2][0]=0.0;
    pointer[1][2]=0.0;
    pointer[2][1]=0.0;
    return(pointer);
}