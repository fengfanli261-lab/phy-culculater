#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define MAXLINE 100

#ifndef _CULC_H_
#define _CULC_H_

enum MODE{single,multi,link,differ} mode;

//文件处理
//跳转或跳过若干行
enum WHENCE{LS_SET,LS_CUR,LS_END};
//offset:偏移行数 whence:起始位置
int lineskip(FILE* fp,int offset, int whence){
    if(whence == LS_SET) {
        rewind(fp);
        for(int i=1;i<offset;i++) while (fgetc(fp) != '\n');
        return 0;
    }
    
    if(whence == LS_CUR) {
        for(int i=0;i<offset;i++) while (fgetc(fp) != '\n');
        return 0;
    }
    
    if(whence == LS_END) {
        fseek(fp,0,SEEK_END);
        for(int i=0;i<offset;i++) while (fgetc(fp) != '\n') fseek(fp,-2,SEEK_CUR);
        return 0;
    }
    return 1;
}

//读取n个double型数据到arr数组
//head:行首标识符 tail:行尾标识符 mid:数据间隔符
int farrays(FILE*fp,double* arr,int n,char* head,char tail,char mid){
    int count = 0;
    if(!fscanf_s(fp,head)) {
        for(int i=0;i<n;i++){
            fscanf_s(fp,"%lf",&arr[i]);
            count ++;
            char c = fgetc(fp);
            if (c == mid) continue;
            else if (c == tail) {
                if (count != n) {
                printf("数据文本读取行缺失%d个值\n",n-count);
                break;
                }
                else return 0;
            }
            else{
            printf("读取中止,请检查数据文本第%d处字符:\"%c\"的格式\n",ftell(fp),c);
            return 1;
            }
        }
    }
    else {
        printf("读取失败,请检查数据文本第%d处格式\n",ftell(fp));
        return 1;
    }

    return n-count;
}

//读取n个int型数据到arr数组
int farrays_(FILE*fp,int* arr,int n,char* head,char tail,char mid){
    int count = 0;
    if(!fscanf_s(fp,head)) {
        for(int i=0;i<n;i++){
            fscanf_s(fp,"%d",&arr[i]);
            count ++;
            char c = fgetc(fp);
            if (c == mid) continue;
            else if (c == tail) {
                if (count != n) {
                printf("数据文本读取行缺失%d个值\n",n-count);
                break;
                }
                else return 0;
            }
            else{
            printf("读取中止,请检查数据文本第%d处字符\"%c\"的格式\n",ftell(fp),c);
            return 1;
            }
        }
    }
    else {
        printf("读取失败,请检查数据文本第%d处格式\n",ftell(fp));
        return 1;
    }

    return n-count;
}


//均值
double mean(int n,double*x){
    double expc=0.0;
    for (int i=0;i<n;i++) expc += x[i]/n;
    return expc;
}

//方差
double var(int n,double*x){
    if (!n){
        printf("样本数量为0,无法计算方差\n");
        return NAN;
    }
    double expc=mean(n,x);
    double var_=0.0;
    for (int i=0;i<n;i++) var_ += (x[i]-expc)*(x[i]-expc)/n;
    return var_;
}

//标准差
double stddev(int n,double*x){
    if (n<2) {
        printf("样本数量过少，无法计算标准差\n");
        return NAN;
    }
    else return sqrt(var(n,x)*n/(n-1));
}

//不确定度合成
double unctty(int n,double*x,int n_b,double u_b){
    double stddev_ = stddev(n,x);
    double u_b2 = n_b*u_b*u_b;
    if ((stddev_*stddev_ + u_b2)) return sqrt(stddev_*stddev_ + u_b2);
    else {
        printf("精度过低或数据有误\n");
        return NAN;
    }
}

//相对不确定度
double unctty_r(int n,double*x,int n_b,double u_b){ 
    if(mean(n,x)) return unctty(n,x,n_b,u_b)/mean(n,x);
    else {
        printf("均值为0,无法计算相对不确定度\n");
        return NAN;
    }
}

//独立变量积式不确定度传递
//g=Const*x_1^alpha_1*x_2^alpha_2*...*x_n^alpha_n
double u_mtrans(int n_l,int* n_m,int* n_b,double** x,double* u_b,double* alpha){
    double* ur_x = NULL;
    ur_x = malloc(n_l*sizeof(double));
    for (int i=0;i<n_l;i++) ur_x[i] = unctty_r(n_m[i],x[i],n_b[i],u_b[i]);
    double ur_g=0.0;
    for (int i=0;i<n_l;i++) ur_g += (alpha[i]*ur_x[i])*(alpha[i]*ur_x[i]);
    return sqrt(ur_g);
}

double u_mtrans_d(int n_l,int* n_d,int* n_b,double** dx,double* u_b,double* alpha){
    double* ur_dx = NULL;
    ur_dx = malloc(n_l*sizeof(double));
    for (int i=0;i<n_l;i++) ur_dx[i] = unctty_r(n_d[i],dx[i],2*n_b[i],u_b[i]);
    double ur_g=0.0;
    for (int i=0;i<n_l;i++) ur_g += (alpha[i]*ur_dx[i])*(alpha[i]*ur_dx[i]);
    return sqrt(ur_g);
}

//回归计算
//协方差
double cov(int n,double* x,double* y){
    double expx=mean(n,x);
    double expy=mean(n,y);
    double cov_=0.0;
    for (int i=0;i<n;i++) cov_ += (x[i]-expx)*(y[i]-expy)/n;
    return cov_;
}

//相关系数
double r_link(int n,double* x,double* y){
    return cov(n,x,y)/sqrt(var(n,x)*var(n,y));
}

//Y=aX+b
//斜率a
double k_link(int n,double* x,double* y){
    return cov(n,x,y)/var(n,x);
}
//截距b
double b_link(int n,double* x,double* y){
    return mean(n,y)-k_link(n,x,y)*mean(n,x);
}
//残差平方和
double rss_link(int n,double* x,double* y){
    double a=k_link(n,x,y);
    double b=b_link(n,x,y);
    double rss=0.0;
    for (int i=0;i<n;i++) rss += (y[i]-a*x[i]-b)*(y[i]-a*x[i]-b);
    return rss;
}
//RSE
double rse_link(int n,double* x,double* y){
    return sqrt(rss_link(n,x,y)/(n-2));
}
//斜率不确定度
double u_k_link(int n,double* x,double* y){
    return rse_link(n,x,y)/sqrt(n*var(n,x));
} 
//截距不确定度
double u_b_link(int n,double* x,double* y){
    return u_k_link(n,x,y)*sqrt(mean(n,x)*mean(n,x)+var(n,x)*n);
} 


#endif