/***************************************************************************
@���ߣ�������
@ѧ�ţ�1160801224
@ʱ�䣺2018.06.07
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926        /*��������PI*/

float the_one_R_xb(float,float,float);      /*һ�����麯������*/
float the_one_R_dxb(float,float,float,float);
float the_one_R_d2xb(float,float,float,float,float);
float the_one_R_yb(float,float,float);
float the_one_R_dyb(float,float,float,float);
float the_one_R_d2yb(float,float,float,float,float);
float the_two_RRR_xc(float,float,float);        /*��������(RRR)��������*/
float the_two_RRR_dxc(float,float,float,float);
float the_two_RRR_d2xc(float,float,float,float,float);
float the_two_RRR_yc(float,float,float);
float the_two_RRR_dyc(float,float,float,float);
float the_two_RRR_d2yc(float,float,float,float,float);
float the_two_RRR_xg(float,float,float);        /*��������(RRR)��������*/
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
    float xa = 0, dxa = 0,d2xa = 0;     /*����A�����*/
    float ya = 0,dya = 0,d2ya = 0;
    float xb = 0,dxb = 0,d2xb = 0;      /*����B�����*/
    float yb = 0,dyb = 0,d2yb = 0;
    float xc = 0,dxc = 0,d2xc = 0;      /*����C�����*/
    float yc = 0,dyc = 0,d2yc = 0;
    float xd = -239.32,dxd = 0,d2xd = 0;        /*����D�����*/
    float yd = 158.4,dyd = 0,d2yd = 0;
    float xg = 0,dxg = 0,d2xg = 0;      /*����G�����*/
    float yg = 0,dyg = 0,d2yg = 0;
    float fi = 0,dfi = 10,d2fi = 0;     /*����AB�˽ǶȲ���*/
    float fi2 = 0,dfi2 = 0,d2fi2 = 0;       /*����BC�˽ǶȲ���*/
    float fj = 0,dfj = 0,d2fj = 0;      /*����CD�˽ǶȲ���*/
    float li = 100,li2 = 273,lj = 136;      /*����AB��BC��CD�˳�*/
    float lcg = 232;        /*����CG�˳�*/
    int sw = 0;     /*����ѡ�����*/

    printf("��������Ҫ��������ݣ�\n");
    printf("\n\n");
    printf("1:     G�����ꡢ�ٶȡ����ٶ�\n");
    printf("2:     ��6�ĽǶȡ����ٶȡ��Ǽ��ٶ�\n");
    printf("3:     ��8�ĽǶȡ����ٶȡ��Ǽ��ٶ�\n");
    printf("\n\n");
    scanf("%d",&sw);

    FILE *fp1Write = fopen("point_g.txt","w");       /*�����ı��ĵ�*/
    FILE *fp2Write = fopen("bar_6.txt","w");     /*�����ı��ĵ�*/
    FILE *fp3Write = fopen("bar_8.txt","w");     /*�����ı��ĵ�*/

    for(fi = 0;fi <= 2 * PI;fi = fi + PI / 180)
    {
         xb = the_one_R_xb(xa,li,fi);       /*����B��*/
         yb = the_one_R_yb(ya,li,fi);
        dxb = the_one_R_dxb(dxa,li,fi,dfi);
        dyb = the_one_R_dyb(dya,li,fi,dfi);
        d2xb = the_one_R_d2xb(d2xa,fi,dfi,d2fi,li);
        d2yb = the_one_R_d2yb(d2ya,fi,dfi,d2fi,li);
        fi2 = the_two_RRR_fi2(xb,xd,yb,yd,li2,lj);      /*����BC�˽Ƕ�*/
        xc = the_two_RRR_xc(xb,fi2,li2);        /*����C��*/
        yc = the_two_RRR_yc(yb,fi2,li2);
        fj = the_two_RRR_fj(xc,xd,yc,yd);       /*����CD�˽Ƕ�*/
        dfi2 = the_two_RRR_dfi2(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
        dfj = the_two_RRR_dfj(dxb,dxd,dyb,dyd,fi2,fj,li2,lj);
        d2fi2 = the_two_RRR_d2fi2(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
        d2fj = the_two_RRR_d2fj(d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj);
        dxc = the_two_RRR_dxc(dxb,fi2,dfi2,li2);
        dyc = the_two_RRR_dyc(dyb,fi2,dfi2,li2);
        d2xc = the_two_RRR_d2xc(d2xb,fi2,dfi2,d2fi2,li2);
        d2yc = the_two_RRR_d2yc(d2yb,fi2,dfi2,d2fi2,li2);

        if(sw == 1)     /*����G��*/
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
        else if(sw == 2)        /*����GM�˽Ƕȡ����ٶȡ��Ǽ��ٶ�*/
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
        else if(sw == 3)        /*����GE�˽Ƕȡ����ٶȡ��Ǽ��ٶ�*/
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
�������ƣ�the_one_R_xb
��������������һ��������B���X����
���������xa,li,fi
����ֵ��B�������
***************************************************************************/
float the_one_R_xb(float xa,float li,float fi)
{
    float r = 0;
    r = xa + li * cos(fi);
    return r;
}

/***************************************************************************
�������ƣ�the_one_R_dxb
��������������һ��������B���ٶȵ�X�������
���������dxa,li,fi,dfi
����ֵ��B���ٶȵ�X�������
***************************************************************************/
float the_one_R_dxb(float dxa,float li,float fi,float dfi)
{
    float r = 0;
    r = dxa - dfi * li * sin(fi);
    return r;
}

/***************************************************************************
�������ƣ�the_one_R_d2xb
��������������һ��������B����ٶȵ�X�������
���������d2xa,fi,dfi,d2fi,li
����ֵ��B����ٶȵ�X�������
***************************************************************************/
float the_one_R_d2xb(float d2xa,float fi,float dfi,float d2fi,float li)
{
    float r = 0;
    r = d2xa - pow(dfi,2) * li * cos(fi) - d2fi * li * sin(fi);
    return r;
}

/***************************************************************************
�������ƣ�the_one_R_yb
��������������һ��������B���Y����
���������ya,li,fi
����ֵ��B���Y����
***************************************************************************/
float the_one_R_yb(float ya,float li,float fi)
{
    float r = 0;
    r = ya + li * sin(fi);
    return r;
}

/***************************************************************************
�������ƣ�the_one_R_dyb
��������������һ��������B���ٶȵ�Y�������
���������dya,li,fi,dfi
����ֵ��B���ٶȵ�Y�������
***************************************************************************/
float the_one_R_dyb(float dya,float li,float fi,float dfi)
{
    float r = 0;
    r = dya + dfi * li * cos(fi);
    return r;
}

/***************************************************************************
�������ƣ�the_one_R_d2yb
��������������һ��������B���ٶȵ�Y�������
���������d2ya,fi,dfi,d2fi,li
����ֵ��B���ٶȵ�Y�������
***************************************************************************/
float the_one_R_d2yb(float d2ya,float fi,float dfi,float d2fi,float li)
{
    float r = 0;
    r = d2ya - pow(dfi,2) * li * sin(fi) + d2fi * li * cos(fi);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_xc
��������������B�������Ĳ��������������(RRR)��C���X����
���������xb,fi2,li2
����ֵ��C���X����
***************************************************************************/
float the_two_RRR_xc(float xb,float fi2,float li2)
{
    float r = 0;
    r = xb + li2 * cos(fi2);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_dxc
��������������B�������Ĳ��������������(RRR)��C���ٶȵ�X�������
���������dxb,fi2,dfi2,li2
����ֵ��C���ٶȵ�X�������
***************************************************************************/
float the_two_RRR_dxc(float dxb,float fi2,float dfi2,float li2)
{
    float r = 0;
    r = dxb - dfi2 * li2 * sin(fi2);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_d2xc
��������������B�������Ĳ��������������(RRR)��C����ٶȵ�X�������
���������d2xb,fi2,dfi2,d2fi2,li2
����ֵ��C����ٶȵ�X�������
***************************************************************************/
float the_two_RRR_d2xc(float d2xb,float fi2,float dfi2,float d2fi2,float li2)
{
    float r = 0;
    r = d2xb - d2fi2 * li2 * sin(fi2) - pow(dfi2,2) * li2 * cos(fi2);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_yc
��������������B�������Ĳ��������������(RRR)��C���Y����
���������yb,fi2,li2
����ֵ��C���Y����
***************************************************************************/
float the_two_RRR_yc(float yb,float fi2,float li2)
{
    float r = 0;
    r = yb + li2 * sin(fi2);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_dyc
��������������B�������Ĳ��������������(RRR)��C���ٶȵ�Y�������
���������dyb,fi2,dfi2,li2
����ֵ��C���ٶȵ�Y�������
***************************************************************************/
float the_two_RRR_dyc(float dyb,float fi2,float dfi2,float li2)
{
    float r = 0;
    r = dyb + dfi2 * li2 * cos(fi2);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_d2yc
��������������B�������Ĳ��������������(RRR)��C����ٶȵ�Y�������
���������d2yb,fi2,dfi2,d2fi2,li2
����ֵ��C����ٶȵ�Y�������
***************************************************************************/
float the_two_RRR_d2yc(float d2yb,float fi2,float dfi2,float d2fi2,float li2)
{
    float r = 0;
    r = d2yb + d2fi2 * li2 * cos(fi2) - pow(dfi2,2) * li2 * sin(fi2);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_fi2
��������������B�������Ĳ�����D����������������(RRR)��BC�˵ĽǶ�
���������xb,xd,yb,yd,li2,lj
����ֵ��BC�˵ĽǶ�
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
�������ƣ�the_two_RRR_dfi2
��������������B�������Ĳ�����D�������BC�˽Ƕȡ�CD�˽Ƕȼ����������(R             RR)��BC�˵Ľ��ٶ�
���������dxb,dxd,dyb,dyd,fi2,fj,li2,lj
����ֵ��BC�˵Ľ��ٶ�
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
�������ƣ�the_two_RRR_d2fi2
��������������B�������Ĳ�����D�������BC�˽ǶȺͽ��ٶȡ�CD�˽ǶȺͽ���            �ȼ����������(RRR)��BC�˵Ľ��ٶ�
���������d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj
����ֵ��BC�˵Ľ��ٶ�
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
�������ƣ�the_two_RRR_fj
��������������C�������Ĳ�����D����������������(RRR)��CD�˵ĽǶ�
���������xc,xd,yc,yd
����ֵ��CD�˵ĽǶ�
***************************************************************************/
float the_two_RRR_fj(float xc,float xd,float yc,float yd)
{
    float r = 0;
    r = atan((yc - yd) / (xc - xd));
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_dfj
��������������B�������Ĳ�����D�������BC�˽Ƕȡ�CD�˽Ƕȼ����������(R             RR)��CD�˵Ľ��ٶ�
���������dxb,dxd,dyb,dyd,fi2,fj,li2,lj
����ֵ��CD�˵Ľ��ٶ�
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
�������ƣ�the_two_RRR_d2fj
��������������B�������Ĳ�����D�������BC�˽ǶȺͽ��ٶȡ�CD�˽ǶȺͽ���            �ȼ����������(RRR)��CD�˵Ľ��ٶ�
���������d2xb,d2xd,d2yb,d2yd,fi2,dfi2,fj,dfj,li2,lj
����ֵ��CD�˵Ľ��ٶ�
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
�������ƣ�the_two_RRR_xg
��������������C�������Ĳ��������������(RRR)��G���X����
���������xc,fi2,lcg
����ֵ��G���X����
***************************************************************************/
float the_two_RRR_xg(float xc,float fi2,float lcg)
{
    float r = 0;
    r = xc + lcg * cos(fi2 - PI / 4);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_dxg
��������������C�������Ĳ��������������(RRR)��G���ٶȵ�X�������
���������dxc,fi2,dfi2,lcg
����ֵ��G���ٶȵ�X�������
***************************************************************************/
float the_two_RRR_dxg(float dxc,float fi2,float dfi2,float lcg)
{
    float r = 0;
    r = dxc - dfi2 * lcg * sin(fi2 - PI / 4);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_d2xg
��������������C�������Ĳ��������������(RRR)��G����ٶȵ�X�������
���������d2xc,fi2,dfi2,d2fi2,lcg
����ֵ��G����ٶȵ�X�������
***************************************************************************/
float the_two_RRR_d2xg(float d2xc,float fi2,float dfi2,float d2fi2,float lcg)
{
    float r = 0;
    r = d2xc - d2fi2 * lcg * sin(fi2 - PI / 4) - pow(dfi2,2) * lcg * cos(fi2 - PI / 4);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_yg
��������������C�������Ĳ��������������(RRR)��G���Y����
���������yc,fi2,lcg
����ֵ��G���Y����
***************************************************************************/
float the_two_RRR_yg(float yc,float fi2,float lcg)
{
    float r = 0;
    r = yc + lcg * sin(fi2 - PI / 4);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_dyg
��������������C�������Ĳ��������������(RRR)��G���ٶȵ�Y�������
���������dyc,fi2,dfi2,lcg
����ֵ��G���ٶȵ�Y�������
***************************************************************************/
float the_two_RRR_dyg(float dyc,float fi2,float dfi2,float lcg)
{
    float r = 0;
    r = dyc + dfi2 * lcg * cos(fi2 - PI / 4);
    return r;
}

/***************************************************************************
�������ƣ�the_two_RRR_d2yg
��������������C�������Ĳ��������������(RRR)��G����ٶȵ�Y�������
���������d2yc,fi2,dfi2,d2fi2,lcg
����ֵ��G����ٶȵ�Y�������
***************************************************************************/
float the_two_RRR_d2yg(float d2yc,float fi2,float dfi2,float d2fi2,float lcg)
{
    float r = 0;
    r = d2yc + d2fi2 * lcg * cos(fi2 - PI / 4) - pow(dfi2,2) * lcg * sin(fi2 - PI / 4);
    return r;
}
