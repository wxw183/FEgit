#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double esm[6][6],x[3],y[3],d[3][3];// esm[6][6]—单元刚度矩阵， x[3],y[3]—单元节点坐标， d[3][3]—材料性质矩阵
double b[3][6],ar2;//b[3][6]—几何矩阵， ar2—三角形面积的二倍
double stra[3],stre[3];// a[8500]—存储节点位移向量、节点力向量和总体刚度矩阵的数组 a，STRA(3),STRE(3)—存储单元的应变、应力*/
int np,nbw,jgf,jgsm,jend;// np—自由度总数， nbw—最大半带宽， jgf、 jgsm、 jend 为计数单元，jgf=np－节点位移向量在数组 a 中的位置 ， jgsm=jgf+np－节点力向量在数组 a 中的位置，jend=jgsm+np*nbw—刚度矩阵在数组 a 中的位置，数组 a 总长度
int ns[6],u[6];//ns[6]—一个单元的节点自由度编号数组， u[6]—一个单元的节点自由度
double em,pr,th;//EM－杨氏模量， PR－泊松比， TH－板的厚度

void modify(double *a,FILE *fpi,FILE *fpo);
void elstmx(int kk);
void dcmpbd(double *a);
void slvbd(double *a);

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
    /*输入材料的杨氏模量 EM，波松比 PR，平板厚度 TH，节点坐标 XC(I)，YC(I)和组成单元的节点 NEL(N,I)
组成单元的节点的编号都按逆时针顺序输入*/
    
    fscanf(fpi,"%lg\n%lg\n%lg\n",&em,&pr,&th);
    for(int k=0;k<nn;k++){
        fscanf(fpi,"%lg",&xc[k]);
    }
    for(int k=0;k<nn;k++){
        fscanf(fpi,"%lg",&yc[k]);
    }
//输出材料性质和计算模型拓扑数据便于检查时对照
    fprintf(fpo,"\n材料常数为：%g\t泊松比为：%g\t厚度为：%g\n",em,pr,th);
    for(int k=0;k<nn;k++){
        fprintf(fpo,"\n节点%d坐标为:X=%fY=%f\n",k+1,xc[k],yc[k]);
    }
    for(int k=0;k<ne;k++){
        fscanf(fpi,"%d%d%d%d",&n,&nel[k][0],&nel[k][1],&nel[k][2]);
        fprintf(fpo,"\n单元号码为：%d\t组成单元的节点号码为:%d\t%d\t%d\n",n,nel[k][0],nel[k][1],nel[k][2]);
    }

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
    jgf=np;
    jgsm=jgf+np;
    jend=jgsm+np*nbw;
    double a[jend];
    int jl=jend-jgf;
    for(i=0;i<jend;i++)a[i]=0.0;

    //生成材料特性矩阵d
    double r=em/(1.0-pr*pr);
    d[0][0]=r;
    d[1][1]=r;
    d[2][2]=r*(1.0-pr)/2.0;
    d[1][2]=pr*r;
    d[2][1]=d[1][2];
    d[1][3]=0.0;
    d[3][1]=0.0;
    d[2][3]=0.0;
    d[3][2]=0.0;
    
    //单元矩阵循环的开始
    k=0;//k为单元号
    while(k<ne){
    for(i=0;i<3;i++){
        j=nel[k][i];
        ns[2*i]=(j-1)*2;
        ns[2*i+1]=j*2-1;
        x[i]=xc[j];
        y[i]=yc[j];
    }

    

    elstmx(k);
    
    //单元刚度矩阵组装成总体刚度矩阵
    int ii,jj;
    for(i=0;i<6;i++){
        ii=ns[i];
        for(j=0;j<6;j++){
            jj=ns[j];
            if(jj<ii)continue;
            int j1=jgsm+(jj-ii)*np+ii+1+(jj-ii-1)*(jj-ii)/2;
            a[j1]=a[j1]+esm[i][j];
        }
    }
    k++;
    }

    //调用子程序 MODIFY 输入载荷节点处的载荷值、位移边界节点处的位移值 ,对总体刚度矩阵、位移数组和节点力数组进行相应的修改
    
    modify(a,fpi,fpo);
    fprintf(fpo,"计算结果为：\n");
    fprintf(fpo,"最大半带宽为%d\n",nbw);
    fprintf(fpo,"总数组大小为%d\n",jend);


    fclose(fpi);
    fclose(fpo);
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
                esm[i][j]=sum*th/(2.0*ar2);
            }
        }
    }

}

void modify(double *a,FILE *fpi,FILE *fpo){
    double bv;
    int ib;
    fscanf(fpi,"%d",&ib);
    for(int i=0;ib!=0;i++){
        fscanf(fpi,"%lg",&bv);

        if(ib%2)fprintf(fpo,"节点%d载荷为:PX=%f\n",(ib+1)/2,bv);
        else fprintf(fpo,"节点%d载荷为:PY=%f\n",ib/2,bv);
        a[jgf+ib]+=bv;
        fscanf(fpi,"%d",&ib);
    }
    fscanf(fpi,"%d",&ib);
    for(int i=0;ib!=0;i++){
        fscanf(fpi,"%lg",&bv);

        if(ib%2)fprintf(fpo,"节点%d位移约束为:U=%f\n",(ib+1)/2,bv);
        else fprintf(fpo,"节点%d位移约束为:V=%f\n",ib/2,bv);
        a[jgsm+ib]+=bv;
        int j1,ii,jj;//ii,jj为总体刚度矩阵的行号列号，j1为a矩阵中iijj元素对应的位置
        jj=i;
        ii=i;
        j1=jgsm+(jj-ii)*np+ii-(jj-ii+1)*(jj-ii)/2;
        double aa=a[j1];//暂存k_ii
        for(int j=0;j<np;j++){
            j1=jgsm+(jj-j)*np+j-(jj-j+1)*(jj-j)/2;
            if(jj>=j)a[j1]=0;
        }//k矩阵第i行变为0
        for(int j=0;j<np;j++){
            jj=j;
            j1=jgsm+(j-ii)*np+ii-(j-ii+1)*(j-ii)/2;
            if(j>=ii)a[j1]=0;
        }//k矩阵第i列变为0
        a[j1]=aa;







        fscanf(fpi,"%d",&ib);
    }

    
}

void dcmpbd(double *a){
    ;
    

}

void slvbd(double *a){
    ;
}