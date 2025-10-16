#include "culc.h"

int main() {
    //设置控制台为utf-8编码以支持中文输出
    system("chcp 65001>nul"); 

    //变量定义
    int mode=0;
    int n_l = 1; //数据组数量
    int* n_m = NULL; //各组测量数据数量
    int* n_b= NULL; //各组B类不确定度计入数
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
    head_x = malloc((7+(int)n_l/10)*sizeof(char));
    for (int i=0;i<n_l;i++) {
        x[i] = malloc(n_m[i]*sizeof(double));//动态分配数据内存
        lineskip(fpi,5+i,LS_SET); //转至数据行
        //生成数据行头
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
    case differ:
        printf("mode = differ\n");
        fprintf_s(fpo,"\nmode = differ\n");
        break;
    default:
        printf("mode参数错误,请检查data.txt文件\n");
        fclose(fpo);
        printf("按任意键退出...");
        getchar();
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



    switch(mode){
    //单变量统计量计算
    case single:
        {
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

            break;
        }

    //积式g(c*π(xi^ai))相对不确定度计算
    case multi:
        {
            printf("相对不确定度ur(g)=%.6g\n",u_mtrans(n_l,n_m,n_b,x,u_b,alpha));

            fprintf_s(fpo,"相对不确定度ur(g)=%.6g\n\n",u_mtrans(n_l,n_m,n_b,x,u_b,alpha));

            break;
        }

    //一元线性回归计算
    case link:
        {
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

        break;
        }
    
    //逐差法计算
    case differ:
    {
        double **dx = NULL;
        int *n_d = NULL;
        n_d = malloc(n_l*sizeof(int));
        for (int i=0;i<n_l;i++) n_d[i] = (n_m[i]-n_m[i]%2)/2;//计算逐差数据数量

        dx = malloc(n_l*sizeof(double*));
        for (int i=0;i<n_l;i++) dx[i] = malloc((n_d[i])*sizeof(double));//动态分配逐差数据内存
        for (int i=0;i<n_l;i++) {
            for (int j=0;j<n_d[i];j++) dx[i][j] = x[i][j+n_d[i]]-x[i][j];
        }//计算逐差数据

        printf("逐差数据计算完成, 结果如下:\n");
        for (int i=0;i<n_l;i++){
            printf("第%d组逐差数据：\n",i+1);
            printf("逐差数据dx%d：",i+1);
            for (int j=0;j<n_d[i];j++) printf("%.6g ",dx[i][j]);
            printf("\n");
            printf("均值dx=∑dxi/n:%.6g\n",mean(n_d[i],dx[i]));
            printf("方差D(dx)=∑(dxi-dx)^2/n:%.6g\n",var(n_d[i],dx[i]));
            printf("标准差s(dx)=√((dxi-dx)^2/(n-1)):%.6g\n",stddev(n_d[i],dx[i]));
            printf("合成不确定度u(dxi):%.6g\n",unctty(n_d[i],dx[i],2*n_b[i],u_b[i]));
            printf("相对不确定度ur(dxi):%.6g\n",unctty_r(n_d[i],dx[i],2*n_b[i],u_b[i]));
            printf("\n");

            fprintf_s(fpo,"第%d组逐差数据：\n",i+1);
            fprintf_s(fpo,"逐差数据dx%d：",i+1);
            for (int j=0;j<n_d[i];j++) fprintf_s(fpo,"%.6g ",dx[i][j]);
            fprintf_s(fpo,"\n均值dx=∑dxi/n:%.6g\n",mean(n_d[i],dx[i]));
            fprintf_s(fpo,"方差D(dx)=∑(dxi-dx)^2/n:%.6g\n",var(n_d[i],dx[i]));
            fprintf_s(fpo,"标准差s(dx)=√((dxi-dx)^2/(n-1)):%.6g\n",stddev(n_d[i],dx[i]));
            fprintf_s(fpo,"合成不确定度u(dxi):%.6g\n",unctty(n_d[i],dx[i],2*n_b[i],u_b[i]));
            fprintf_s(fpo,"相对不确定度ur(dxi):%.6g\n",unctty_r(n_d[i],dx[i],2*n_b[i],u_b[i]));
            fprintf_s(fpo,"\n");
        }

        //释放内存
        for (int i=0;i<n_l;i++) free(dx[i]);
        free(dx);
        free(n_d);

        break;
    }
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

    //释放内存
    free(n_m);
    free(n_b);
    for (int i=0;i<n_l;i++) free(x[i]);
    free(x);
    free(u_b);
    free(alpha);
    free(head_x);



    printf("按任意键退出...");
    getchar();

    return 0;
}

