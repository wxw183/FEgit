#include<stdio.h>
#include<stdlib.h>

void aa(double *a);

int main(){
    FILE *fp;
    char s1[100],s2[100];
    fp=fopen("t1.txt","r+");
    /*fgets(s1,100,fp);
    fgets(s2,100,fp);
    fclose(fp);
    fp=fopen("t2.txt","w+");
    fputs(s1,fp);
    fputs(s2,fp);
    fclose(fp);
*/
    double n,i;
    fgets(s1,100,fp);
    fputs(s1,stdout);
    fscanf(fp,"%lg%lg",&n,&i);
    fprintf(stdout,"%f\n%f\n",n,i);
    double *nn;
    nn=&n;
    aa(nn);
    printf("%f",n);
    return(0);
}

void aa(double *a){
    *a=159;
}