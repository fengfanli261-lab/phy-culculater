#include "culc.h"

int main() {
    //设置控制台为utf-8编码以支持中文输出
    system("chcp 65001>nul"); 

    //变量定义
    int mode=0;
    int n_l = 1; //数据组数量
    int* n_m = NULL; //各组测量数据数量
    int* n_b= NULL; //各组B类不确定度数量
    double **x = NULL; //测量数据
    double *u_b = NULL; //B类不确定度
    double* alpha = NULL; //积式指数

    //文件读取
    FILE* fpi=NULL;
    fopen_s(&fpi,"data.txt","r+");
    char str[MAXLINE]={0};
    sscanf_s(fgets(str,MAXLINE,fpi),"mode = %d",&mode); //读取mode
    sscanf_s(fgets(str,MAXLINE,fpi),"num_list = %d",&n_l); //读取数据组数量

    n_m = malloc(n_l*sizeof(int));
    n_b = malloc(n_l*sizeof(int));//动态分配数量数组内存
    x = malloc(n_l*sizeof(double*));
    u_b = malloc(n_l*sizeof(double));//动态分配数据数组内存
    lineskip(fpi,3,LS_SET); //转至num_measure行
    farrays_(fpi,n_m,n_l,"num_measure = [",']',','); //读取各组测量数据数量
    lineskip(fpi,4,LS_SET); //转至num_b行
    farrays_(fpi,n_b,n_l,"num_b = [",']',','); //读取各组B类不确定度计数
    
    char* head_x = NULL;
    char* head_b = NULL;
    for (int i=0;i<n_l;i++) {
        x[i] = malloc(n_m[i]*sizeof(double));//动态分配数据内存
        lineskip(fpi,5+i,LS_SET); //转至数据行
        //生成数据行头
        head_x = malloc((7+(int)n_l/10)*sizeof(char));
        sprintf_s(head_x,7+(int)n_l/10,"x%d = [",i+1);
        //读取数据
        farrays(fpi,x[i],n_m[i],head_x,']',','); //读取测量数据
        fscanf_s(fpi,";u_b = %lf",&u_b[i]); //读取B类不确定度
    }
    alpha = malloc(n_l*sizeof(double));//动态分配指数数组内存
    lineskip(fpi,1,LS_END);//转至alpha行
    farrays(fpi,alpha,n_l,"alpha = [",']',','); //读取指数

    //关闭data文件
    fclose(fpi);
    //打开result文件以写入结果
    FILE* fpo=NULL;
    fopen_s(&fpo,"result.txt","a+");
    
    printf("读取完成, ");
    switch (mode)
    {
    case single:
        printf("mode = single\n");
        fprintf_s(fpo,"\nmode = single\n");
        break;
    case multi:
        printf("mode = multi\n");
        fprintf_s(fpo,"\nmode = multi\n");
        break;
    case link:
        printf("mode = link\n");
        fprintf_s(fpo,"\nmode = link\n");
        break;
    default:
        printf("mode参数错误,请检查data.txt文件\n");
        return 1;
        break;
    }
    printf("\n");
    for (int i=0;i<n_l;i++){
        printf("第%d组数据：\n",i+1);
        printf("测量数据x%d：",i+1);
        for (int j=0;j<n_m[i];j++) printf("%.6g ",x[i][j]);
        printf("\nB类不确定度u_b%d：",i+1);
        printf("%.6g ",u_b[i]);
        printf("\n");
        printf("积式指数：%.6g\n",alpha[i]);
        printf("\n");

        fprintf_s(fpo,"第%d组数据：\n",i+1);
        fprintf_s(fpo,"测量数据x%d：",i+1);
        for (int j=0;j<n_m[i];j++) fprintf_s(fpo,"%.6g ",x[i][j]);
        fprintf_s(fpo,"\nB类不确定度u_b%d：",i+1);
        fprintf_s(fpo,"%.6g ",u_b[i]);
        fprintf_s(fpo,"\n");
        fprintf_s(fpo,"积式指数：%.6g\n",alpha[i]);
        fprintf_s(fpo,"\n");
    }
    printf("\n");

    //单变量统计量计算
    if (mode == single) {
        for(int i=0;i<n_l;i++){
            printf("第%d组数据统计量计算:\n",i+1);
            printf("均值x=∑xi/n:%.6g\n",mean(n_m[i],x[i]));
            printf("方差D(x)=∑(xi-x)^2/n:%.6g\n",var(n_m[i],x[i]));
            printf("标准差s(x)=√((xi-x)^2/(n-1)):%.6g\n",stddev(n_m[i],x[i]));
            printf("合成不确定度u(xi):%.6g\n",unctty(n_m[i],x[i],n_b[i],u_b[i]));
            printf("相对不确定度ur(xi):%.6g\n",unctty_r(n_m[i],x[i],n_b[i],u_b[i]));
            printf("\n");

            fprintf_s(fpo,"第%d组数据统计量计算:\n",i+1);
            fprintf_s(fpo,"均值x=∑xi/n:%.6g\n",mean(n_m[i],x[i]));
            fprintf_s(fpo,"方差D(x)=∑(xi-x)^2/n:%.6g\n",var(n_m[i],x[i]));
            fprintf_s(fpo,"标准差s(x)=√((xi-x)^2/(n-1)):%.6g\n",stddev(n_m[i],x[i]));
            fprintf_s(fpo,"合成不确定度u(xi):%.6g\n",unctty(n_m[i],x[i],n_b[i],u_b[i]));
            fprintf_s(fpo,"相对不确定度ur(xi):%.6g\n",unctty_r(n_m[i],x[i],n_b[i],u_b[i]));
            fprintf_s(fpo,"\n");
        }
    }

    //积式g(c*π(xi^ai))相对不确定度计算
    if (mode == multi){
        printf("相对不确定度ur(g)=%.6g\n",u_mtrans(n_l,n_m,n_b,x,u_b,alpha));

        fprintf_s(fpo,"相对不确定度ur(g)=%.6g\n\n",u_mtrans(n_l,n_m,n_b,x,u_b,alpha));
    }

    if (mode == link){
    printf("一元线性回归计算：\n");
    printf("Y=%.6gX+%.6gf\n",k_link(n_m[0],x[0],x[1]),b_link(n_m[0],x[0],x[1]));
    printf("相关系数:1-%.3g\n",1-r_link(n_m[0],x[0],x[1]));
    printf("斜率k=%.6g,不确定度u(k)=%.6g\n",k_link(n_m[0],x[0],x[1]),u_k_link(n_m[0],x[0],x[1]));
    printf("截距b=%.6g,不确定度u(b)=%.6g\n",b_link(n_m[0],x[0],x[1]),u_b_link(n_m[0],x[0],x[1]));

    fprintf_s(fpo,"一元线性回归计算：\n");
    fprintf_s(fpo,"Y=%.6gX+%.6g\n",k_link(n_m[0],x[0],x[1]),b_link(n_m[0],x[0],x[1]));
    fprintf_s(fpo,"相关系数:1-%.3g\n",1-r_link(n_m[0],x[0],x[1]));
    fprintf_s(fpo,"斜率k=%.6g,不确定度u(k)=%.6g\n",k_link(n_m[0],x[0],x[1]),u_k_link(n_m[0],x[0],x[1]));
    fprintf_s(fpo,"截距b=%.6g,不确定度u(b)=%.6g\n\n",b_link(n_m[0],x[0],x[1]),u_b_link(n_m[0],x[0],x[1]));
    }

    //写入时间戳
    time_t timer;
    struct tm *Now;
    time( &timer );
    Now = localtime( &timer );
    fprintf_s(fpo,"本次记录时间:%s\n", asctime(Now));
    fclose(fpo);

    printf_s("\n本次计算已写入历史结果,可在result.txt文件中查询。\n");
    printf_s("%s\n",asctime(Now));
    getchar();

    return 0;
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


