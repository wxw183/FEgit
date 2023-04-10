#include<stdio.h>
#include<stdlib.h>
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
    int n,i;
    fscanf(fp,"%d%d",&n,&i);
    printf("%d\n%d\n",n,i);
    return(0);
}