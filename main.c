/***************************************************************************
@作者：耿翰林
@学号：1160801224
@时间：2018.06.07
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926        /*定义宏变量PI*/

float the_one_R_xb(float,float,float);      /*一级杆组函数申明*/
float the_one_R_dxb(float,float,float,float);
float the_one_R_d2xb(float,float,float,float,float);
float the_one_R_yb(float,float,float);
float the_one_R_dyb(float,float,float,float);
float the_one_R_d2yb(float,float,float,float,float);
float the_two_RRR_xc(float,float,float);        /*二级杆组(RRR)函数申明*/
float the_two_RRR_dxc(float,float,float,float);
float the_two_RRR_d2xc(float,float,float,float,float);
float the_two_RRR_yc(float,float,float);
float the_two_RRR_dyc(float,float,float,float);
float the_two_RRR_d2yc(float,float,float,float,float);
float the_two_RRR_xg(float,float,float);        /*二级杆组(RRR)函数申明*/
float the_two_RRR_dxg(float,float,float,float);
float the_two_RRR_d2xg(float,float,float,float,float);
float the_two_RRR_yg(float,float,float);
float the_two_RRR_dyg(float,float,float,float);
float the_two_RRR_d2yg(float,float,float,float,float);
float the_two_RRR_fi2(float,float,float,float,float,float);
float the_two_RRR_dfi2(float,float,float,float,float,float, float,float);
float the_two_RRR_d2fi2(float,float,float,float,float,float,float,float,float,float);
float the_two_RRR_fj(float,float,float,float);
float the_two_RRR_dfj(float,float,float,float,float,float, float,float);
float the_two_RRR_d2fj(float,float,float,float,float,float,float,float,float,float);

int main()
{
    float xa = 0, dxa = 0,d2xa = 0;     /*定义A点参数*/
    float ya = 0,dya = 0,d2ya = 0;
    float xb = 0,dxb = 0,d2xb = 0;      /*定义B点参数*/
    float yb = 0,dyb = 0,d2yb = 0;
    float xc = 0,dxc = 0,d2xc = 0;      /*定义C点参数*/
    float yc = 0,dyc = 0,d2yc = 0;
    float xd = -239.32,dxd = 0,d2xd = 0;        /*定义D点参数*/
    float yd = 158.4,dyd = 0,d2yd = 0;
    float xg = 0,dxg = 0,d2xg = 0;      /*定义G点参数*/
    float yg = 0,dyg = 0,d2yg = 0;
    float fi = 0,dfi = 10,d2fi = 0;     /*定义AB杆角度参数*/
    float fi2 = 0,dfi2 = 0,d2fi2 = 0;       /*定义BC杆角度参数*/
    float fj = 0,dfj = 0,d2fj = 0;      /*定义CD杆角度参数*/
    float li = 100,li2 = 273,lj = 136;      /*定义AB、BC、CD杆长*/
    float lcg = 232;        /*定义CG杆长*/
    int sw = 0;     /*定义选择变量*/

    printf("请输入需要计算的数据：\n");
    printf("\n\n");
    printf("1:     G点坐标、速度、加速度\n");
    printf("2:     杆6的角度、角速度、角加速度\n");
    printf("3:     杆8的角度、角速度、角加速度\n");
    printf("\n\n");
    scanf("%d",&sw);

    FILE *fp1Write = fopen("point_g.txt","w");       /*建立文本文档*/
    FILE *fp2Write = fopen("bar_6.txt","w");     /*建立文本文档*/
    FILE *fp3Write = fopen("bar_8.txt","w");     /*建立文本文档*/

    for(fi = 0;fi <= 2 * PI;fi = fi + PI / 180)
    {
         xb = the_one_R_xb(xa,li,fi);       /*计算B点*/
         yb = the_one_R_yb(ya,li,fi);
        dxb = the_one_R_dxb(dxa,li,fi,dfi);
        dyb = the_one_R_dyb(dya,li,fi,dfi);
        d2xb = the_one_R_d2xb(d2xa,fi,dfi,d2fi,li);
        d2yb = the_one_R_d2yb(d2ya,fi,dfi,d2fi,li);
        fi2 = the_two_RRR_fi2(xb,xd,yb,yd,li2,lj);      /*计算BC杆角度*/
        xc = the_two_RRR_xc(xb,fi2,li2);        /*计算C点*/
        yc = the_two_RRR_yc(yb,fi2,li2);
        fj = the_two_RRR_fj(xc,xd,yc,yd);       /*计算CD杆角度*/
        dfi2 = the_two_RRR_dfi2(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
        dfj = the_two_RRR_dfj(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
        d2fi2 = the_two_RRR_d2fi2(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
        d2fj = the_two_RRR_d2fj(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
        dxc = the_two_RRR_dxc(dxb,fi2,dfi2,li2);
        dyc = the_two_RRR_dyc(dyb,fi2,dfi2,li2);
        d2xc = the_two_RRR_d2xc(d2xb,fi2,dfi2,d2fi2,li2);
        d2yc = the_two_RRR_d2yc(d2yb,fi2,dfi2,d2fi2,li2);

        if(sw == 1)     /*计算G点*/
        {
            xg = the_two_RRR_xg(xc,fi2,lcg);
            yg = the_two_RRR_yg(yc,fi2,lcg);
            dxg = the_two_RRR_dxg(dxc,fi2,dfi2,lcg);
            dyg = the_two_RRR_dyg(dyc,fi2,dfi2,lcg);
            d2xg = the_two_RRR_d2xg(d2xc,fi2,dfi2,d2fi2,lcg);
            d2yg = the_two_RRR_d2yg(d2yc,fi2,dfi2,d2fi2,lcg);
            printf("%f      xg = %f,yg = %f\n",fi,xg,yg);
            printf("%f      dxg = %f,dyg = %f\n",fi,dxg,dyg);
            printf("%f      d2xg = %f,d2yg = %f\n",fi,d2xg,d2yg);
            fprintf(fp1Write,"************************************************************************");
            fprintf(fp1Write,"************************************************************\n");
            fprintf(fp1Write,"*          *                  *                  *                     *                 *                   *                     *\n");
            fprintf(fp1Write,"* %f * xg = %f *dxg = %f * d2xg = %f *",fi,xg,dxg,d2xg);
            fprintf(fp1Write," yg = %f * dyg = %f * d2yg = %f *\n",yg,dyg,d2yg);
            fprintf(fp1Write,"*          *                  *                  *                     *                 *                   *                     *\n");
        }
        else if(sw == 2)        /*计算GM杆角度、角速度、角加速度*/
        {
            xg = the_two_RRR_xg(xc,fi2,lcg);
            yg = the_two_RRR_yg(yc,fi2,lcg);
            dxg = the_two_RRR_dxg(dxc,fi2,dfi2,lcg);
            dyg = the_two_RRR_dyg(dyc,fi2,dfi2,lcg);
            d2xg = the_two_RRR_d2xg(d2xc,fi2,dfi2,d2fi2,lcg);
            d2yg = the_two_RRR_d2yg(d2yc,fi2,dfi2,d2fi2,lcg);
            xb = xg;
            yb = yg;
            dxb = dxg;
            dyb = dyg;
            d2xb = d2xg;
            d2yb = d2yg;
            li2 = 136,lj = 191;
            xd = -292.44,yd = 168.07;
            fi2 = the_two_RRR_fi2(xb,xd,yb,yd,li2,lj);
            xc = the_two_RRR_xc(xb,fi2,li2);
            yc = the_two_RRR_yc(yb,fi2,li2);
            fj = the_two_RRR_fj(xc,xd,yc,yd);
            dfi2 = the_two_RRR_dfi2(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
            dfj = the_two_RRR_dfj(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
            d2fi2 = the_two_RRR_d2fi2(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
            d2fj = the_two_RRR_d2fj(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
            dxc = the_two_RRR_dxc(dxb,fi2,dfi2,li2);
            dyc = the_two_RRR_dyc(dyb,fi2,dfi2,li2);
            d2xc = the_two_RRR_d2xc(d2xb,fi2,dfi2,d2fi2,li2);
            d2yc = the_two_RRR_d2yc(d2yb,fi2,dfi2,d2fi2,li2);
            printf("%f      fi = %f,dfi = %f,d2fi = %f\n",fi,fi2,dfi2,d2fi2);
            fprintf(fp2Write,"******************************************************************\n");
            fprintf(fp2Write,"*          *               *                 *                   *\n");
            fprintf(fp2Write,"* %f * fi = %f * dfi = %f * d2fi = %f *\n",fi,fi2,dfi2,d2fi2);
            fprintf(fp2Write,"*          *               *                 *                   *\n");
            li2 = 273,lj = 136;
            xd = -239.32,yd = 158.4;
        }
        else if(sw == 3)        /*计算GE杆角度、角速度、角加速度*/
        {
            xg = the_two_RRR_xg(xc,fi2,lcg);
            yg = the_two_RRR_yg(yc,fi2,lcg);
            dxg = the_two_RRR_dxg(dxc,fi2,dfi2,lcg);
            dyg = the_two_RRR_dyg(dyc,fi2,dfi2,lcg);
            d2xg = the_two_RRR_d2xg(d2xc,fi2,dfi2,d2fi2,lcg);
            d2yg = the_two_RRR_d2yg(d2yc,fi2,dfi2,d2fi2,lcg);
            xb = xg;
            yb = yg;
            dxb = dxg;
            dyb = dyg;
            d2xb = d2xg;
            d2yb = d2yg;
            li2 = 145,lj = 282;
            xd = -232.33,yd = -41.48;
            fi2 = the_two_RRR_fi2(xb,xd,yb,yd,li2,lj);
            xc = the_two_RRR_xc(xb,fi2,li2);
            yc = the_two_RRR_yc(yb,fi2,li2);
            fj = the_two_RRR_fj(xc,xd,yc,yd);
            dfi2 = the_two_RRR_dfi2(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
            dfj = the_two_RRR_dfj(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
            d2fi2 = the_two_RRR_d2fi2(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
            d2fj = the_two_RRR_d2fj(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
            dxc = the_two_RRR_dxc(dxb,fi2,dfi2,li2);
            dyc = the_two_RRR_dyc(dyb,fi2,dfi2,li2);
            d2xc = the_two_RRR_d2xc(d2xb,fi2,dfi2,d2fi2,li2);
            d2yc = the_two_RRR_d2yc(d2yb,fi2,dfi2,d2fi2,li2);
            printf("%f      fi = %f,dfi = %f,d2fi = %f\n",fi,fi2,dfi2,d2fi2);
            fprintf(fp3Write,"******************************************************************\n");
            fprintf(fp3Write,"*          *               *                 *                   *\n");
            fprintf(fp3Write,"* %f * fi = %f * dfi = %f * d2fi = %f *\n",fi,fi2,dfi2,d2fi2);
            fprintf(fp3Write,"*          *               *                 *                   *\n");
            li2 = 273,lj = 136;
            xd = -239.32,yd = 158.4;
        }

    }

    fclose(fp1Write);
    fclose(fp2Write);
    fclose(fp3Write);

    return 0;

}

/***************************************************************************
函数名称：the_one_R_xb
函数描述：计算一级杆组中B点的X坐标
输入参数：xa,li,fi
返回值：B点横坐标
***************************************************************************/
float the_one_R_xb(float xa,float li,float fi)
{
    float r = 0;
    r = xa + li * cos(fi);
    return r;
}

/***************************************************************************
函数名称：the_one_R_dxb
函数描述：计算一级杆组中B点速度的X方向分量
输入参数：dxa,li,fi,dfi
返回值：B点速度的X方向分量
***************************************************************************/
float the_one_R_dxb(float dxa,float li,float fi,float dfi)
{
    float r = 0;
    r = dxa - dfi * li * sin(fi);
    return r;
}

/***************************************************************************
函数名称：the_one_R_d2xb
函数描述：计算一级杆组中B点加速度的X方向分量
输入参数：d2xa,fi,dfi,d2fi,li
返回值：B点加速度的X方向分量
***************************************************************************/
float the_one_R_d2xb(float d2xa,float fi,float dfi,float d2fi,float li)
{
    float r = 0;
    r = d2xa - pow(dfi,2) * li * cos(fi) - d2fi * li * sin(fi);
    return r;
}

/***************************************************************************
函数名称：the_one_R_yb
函数描述：计算一级杆组中B点的Y坐标
输入参数：ya,li,fi
返回值：B点的Y坐标
***************************************************************************/
float the_one_R_yb(float ya,float li,float fi)
{
    float r = 0;
    r = ya + li * sin(fi);
    return r;
}

/***************************************************************************
函数名称：the_one_R_dyb
函数描述：计算一级杆组中B点速度的Y方向分量
输入参数：dya,li,fi,dfi
返回值：B点速度的Y方向分量
***************************************************************************/
float the_one_R_dyb(float dya,float li,float fi,float dfi)
{
    float r = 0;
    r = dya + dfi * li * cos(fi);
    return r;
}

/***************************************************************************
函数名称：the_one_R_d2yb
函数描述：计算一级杆组中B点速度的Y方向分量
输入参数：d2ya,fi,dfi,d2fi,li
返回值：B点速度的Y方向分量
***************************************************************************/
float the_one_R_d2yb(float d2ya,float fi,float dfi,float d2fi,float li)
{
    float r = 0;
    r = d2ya - pow(dfi,2) * li * sin(fi) + d2fi * li * cos(fi);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_xc
函数描述：利用B点求解出的参数计算二级杆组(RRR)中C点的X坐标
输入参数：xb,fi2,li2
返回值：C点的X坐标
***************************************************************************/
float the_two_RRR_xc(float xb,float fi2,float li2)
{
    float r = 0;
    r = xb + li2 * cos(fi2);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_dxc
函数描述：利用B点求解出的参数计算二级杆组(RRR)中C点速度的X方向分量
输入参数：dxb,fi2,dfi2,li2
返回值：C点速度的X方向分量
***************************************************************************/
float the_two_RRR_dxc(float dxb,float fi2,float dfi2,float li2)
{
    float r = 0;
    r = dxb - dfi2 * li2 * sin(fi2);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_d2xc
函数描述：利用B点求解出的参数计算二级杆组(RRR)中C点加速度的X方向分量
输入参数：d2xb,fi2,dfi2,d2fi2,li2
返回值：C点加速度的X方向分量
***************************************************************************/
float the_two_RRR_d2xc(float d2xb,float fi2,float dfi2,float d2fi2,float li2)
{
    float r = 0;
    r = d2xb - d2fi2 * li2 * sin(fi2) - pow(dfi2,2) * li2 * cos(fi2);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_yc
函数描述：利用B点求解出的参数计算二级杆组(RRR)中C点的Y坐标
输入参数：yb,fi2,li2
返回值：C点的Y坐标
***************************************************************************/
float the_two_RRR_yc(float yb,float fi2,float li2)
{
    float r = 0;
    r = yb + li2 * sin(fi2);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_dyc
函数描述：利用B点求解出的参数计算二级杆组(RRR)中C点速度的Y方向分量
输入参数：dyb,fi2,dfi2,li2
返回值：C点速度的Y方向分量
***************************************************************************/
float the_two_RRR_dyc(float dyb,float fi2,float dfi2,float li2)
{
    float r = 0;
    r = dyb + dfi2 * li2 * cos(fi2);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_d2yc
函数描述：利用B点求解出的参数计算二级杆组(RRR)中C点加速度的Y方向分量
输入参数：d2yb,fi2,dfi2,d2fi2,li2
返回值：C点加速度的Y方向分量
***************************************************************************/
float the_two_RRR_d2yc(float d2yb,float fi2,float dfi2,float d2fi2,float li2)
{
    float r = 0;
    r = d2yb + d2fi2 * li2 * cos(fi2) - pow(dfi2,2) * li2 * sin(fi2);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_fi2
函数描述：利用B点求解出的参数和D点参数计算二级杆组(RRR)中BC杆的角度
输入参数：xb,xd,yb,yd,li2,lj
返回值：BC杆的角度
***************************************************************************/
float the_two_RRR_fi2(float xb,float xd,float yb,float yd,float li2,float lj)
{
    float a0 = 0,b0 = 0,c0 = 0,lbd = 0;
    float r = 0;
    lbd = sqrt(pow((xd - xb),2) + pow((yd - yb),2));
    a0 = 2 * li2 * (xd - xb);
    b0 = 2 * li2 * (yd - yb);
    c0 = pow(li2,2) + pow(lbd,2) - pow(lj,2);
    r = 2 * atan((b0 + sqrt(pow(a0,2) + pow(b0,2) - pow(c0,2))) / (a0 + c0));
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_dfi2
函数描述：利用B点求解出的参数、D点参数、BC杆角度、CD杆角度计算二级杆组(R             RR)中BC杆的角速度
输入参数：dxb,dxd,dyb,dyd,fi2,fj,li2,lj
返回值：BC杆的角速度
***************************************************************************/
float the_two_RRR_dfi2(float dxb,float dxd,float dyb,float dyd,float fi2,float fj, float li2,float lj)
{
    float ci = 0,si = 0,cj = 0,sj = 0,g1 = 0;
    float r = 0;
    ci = li2 * cos(fi2);
    si = li2 * sin(fi2);
    cj = lj * cos(fj);
    sj = lj * sin(fj);
    g1 = ci * sj - cj * si;
    r = (cj * (dxd - dxb) + sj * (dyd - dyb)) / g1;
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_d2fi2
函数描述：利用B点求解出的参数、D点参数、BC杆角度和角速度、CD杆角度和角速            度计算二级杆组(RRR)中BC杆的角速度
输入参数：d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj
返回值：BC杆的角速度
***************************************************************************/
float the_two_RRR_d2fi2(float d2xb,float d2xd,float d2yb,float d2yd,float fi2,float dfi2,float fj,float dfj,float li2,float lj)
{
    float ci = 0,si = 0,cj = 0,sj = 0,g1 = 0,g2 = 0,g3 = 0;
    float r = 0;
    ci = li2 * cos(fi2);
    si = li2 * sin(fi2);
    cj = lj * cos(fj);
    sj = lj * sin(fj);
    g1 = ci * sj - cj * si;
    g2 = d2xd - d2xb + pow(dfi2,2) * ci - pow(dfj,2) * cj;
    g3 = d2yd - d2yb + pow(dfi2,2) * si - pow(dfj,2) * sj;
    r = (g2 * cj + g3 * sj) / g1;
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_fj
函数描述：利用C点求解出的参数、D点参数计算二级杆组(RRR)中CD杆的角度
输入参数：xc,xd,yc,yd
返回值：CD杆的角度
***************************************************************************/
float the_two_RRR_fj(float xc,float xd,float yc,float yd)
{
    float r = 0;
    r = atan((yc - yd) / (xc - xd));
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_dfj
函数描述：利用B点求解出的参数、D点参数、BC杆角度、CD杆角度计算二级杆组(R             RR)中CD杆的角速度
输入参数：dxb,dxd,dyb,dyd,fi2,fj,li2,lj
返回值：CD杆的角速度
***************************************************************************/
float the_two_RRR_dfj(float dxb,float dxd,float dyb,float dyd,float fi2,float fj, float li2,float lj)
{
    float ci = 0,si = 0,cj = 0,sj = 0,g1 = 0;
    float r = 0;
    ci = li2 * cos(fi2);
    si = li2 * sin(fi2);
    cj = lj * cos(fj);
    sj = lj * sin(fj);
    g1 = ci * sj - cj * si;
    r = (ci * (dxd - dxb) + si * (dyd - dyb)) / g1;
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_d2fj
函数描述：利用B点求解出的参数、D点参数、BC杆角度和角速度、CD杆角度和角速            度计算二级杆组(RRR)中CD杆的角速度
输入参数：d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj
返回值：CD杆的角速度
***************************************************************************/
float the_two_RRR_d2fj(float d2xb,float d2xd,float d2yb,float d2yd,float fi2,float dfi2,float fj,float dfj,float li2,float lj)
{
    float ci = 0,si = 0,cj = 0,sj = 0,g1 = 0,g2 = 0,g3 = 0;
    float r = 0;
    ci = li2 * cos(fi2);
    si = li2 * sin(fi2);
    cj = lj * cos(fj);
    sj = lj * sin(fj);
    g1 = ci * sj - cj * si;
    g2 = d2xd - d2xb + pow(dfi2,2) * ci - pow(dfj,2) * cj;
    g3 = d2yd - d2yb + pow(dfi2,2) * si - pow(dfj,2) * sj;
    r = (g2 * ci + g3 * si) / g1;
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_xg
函数描述：利用C点求解出的参数计算二级杆组(RRR)中G点的X坐标
输入参数：xc,fi2,lcg
返回值：G点的X坐标
***************************************************************************/
float the_two_RRR_xg(float xc,float fi2,float lcg)
{
    float r = 0;
    r = xc + lcg * cos(fi2 - PI / 4);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_dxg
函数描述：利用C点求解出的参数计算二级杆组(RRR)中G点速度的X方向分量
输入参数：dxc,fi2,dfi2,lcg
返回值：G点速度的X方向分量
***************************************************************************/
float the_two_RRR_dxg(float dxc,float fi2,float dfi2,float lcg)
{
    float r = 0;
    r = dxc - dfi2 * lcg * sin(fi2 - PI / 4);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_d2xg
函数描述：利用C点求解出的参数计算二级杆组(RRR)中G点加速度的X方向分量
输入参数：d2xc,fi2,dfi2,d2fi2,lcg
返回值：G点加速度的X方向分量
***************************************************************************/
float the_two_RRR_d2xg(float d2xc,float fi2,float dfi2,float d2fi2,float lcg)
{
    float r = 0;
    r = d2xc - d2fi2 * lcg * sin(fi2 - PI / 4) - pow(dfi2,2) * lcg * cos(fi2 - PI / 4);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_yg
函数描述：利用C点求解出的参数计算二级杆组(RRR)中G点的Y坐标
输入参数：yc,fi2,lcg
返回值：G点的Y坐标
***************************************************************************/
float the_two_RRR_yg(float yc,float fi2,float lcg)
{
    float r = 0;
    r = yc + lcg * sin(fi2 - PI / 4);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_dyg
函数描述：利用C点求解出的参数计算二级杆组(RRR)中G点速度的Y方向分量
输入参数：dyc,fi2,dfi2,lcg
返回值：G点速度的Y方向分量
***************************************************************************/
float the_two_RRR_dyg(float dyc,float fi2,float dfi2,float lcg)
{
    float r = 0;
    r = dyc + dfi2 * lcg * cos(fi2 - PI / 4);
    return r;
}

/***************************************************************************
函数名称：the_two_RRR_d2yg
函数描述：利用C点求解出的参数计算二级杆组(RRR)中G点加速度的Y方向分量
输入参数：d2yc,fi2,dfi2,d2fi2,lcg
返回值：G点加速度的Y方向分量
***************************************************************************/
float the_two_RRR_d2yg(float d2yc,float fi2,float dfi2,float d2fi2,float lcg)
{
    float r = 0;
    r = d2yc + d2fi2 * lcg * cos(fi2 - PI / 4) - pow(dfi2,2) * lcg * sin(fi2 - PI / 4);
    return r;
}
