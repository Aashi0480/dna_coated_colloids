# include <stdio.h>
# include <math.h>
# define SQ(x) ((x)*(x))
int main()
{
    int i,j,k;
    FILE *in;
    FILE *out1;
    FILE *out2;
    FILE *out3;
    FILE *out4;
    char name1[50],name2[50],name3[50],name4[50],name5[50],name6[50],name7[50],name8[50],name9[50],name10[50],name11[50],name12[50],name13[50],name14[50],name15[50],name16[50],name17[50],name18[50],name19[50],name20[50],name21[50];
    int timestep,atoms,a,b,c,d,e,f,g,h,p;
    double type,id,x,y,z,fx,fy,fz,l,na1,na2,nb1,nb2,b2s,b3s,b2ssur,b3ssur,nn,freesur;
    double l_1,l_2,b2s_1,b2s_2,b3s_1,b3s_2,b2ssur_1,b2ssur_2,b3ssur_1,b3ssur_2,na1_1,na1_2,na2_1,na2_2,nb1_1,nb1_2,nb2_1,nb2_2;
    double x1,y1,z1,x2,y2,z2;
    in=fopen("surB_BDon_10k_30.dump","r");
    out3=fopen("distance_sur_B_BDon_10k_1st_30.txt","w");
    out4=fopen("distance_sur_B_BDon_10k_2nd_30.txt","w");
    double suml=0.0;
    double sumna1=0.0;
    double sumna2=0.0;
    double sumnb1=0.0;
    double sumnb2=0.0;
    double sumb2s=0.0;
    double sumb3s=0.0;
    double sumb2ssur=0.0;
    double sumb3ssur=0.0;
    int t = 10000;
    for (i=0;i<t;i++)
    {
        fscanf(in,"%s %s",name1,name2);
        fscanf(in,"%d",&timestep);
        fscanf(in,"%s %s %s %s",name1,name2,name3,name4);
        fscanf(in,"%d",&atoms);
        fscanf(in,"%s %s %s %s %s %s",name1,name2,name3,name4,name5,name6);
        fscanf(in,"%d %d %d",&a,&b,&c);
        fscanf(in,"%d %d %d",&d,&e,&f);
        fscanf(in,"%d %d %d",&g,&h,&p);
        fscanf(in,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14,name15,name16,name17,name18,name19,name20,name21);

        for(k=0;k<1;k++)
        {
            fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&type,&id,&x,&y,&z,&fx,&fy,&fz,&l,&na1,&na2,&nb1,&nb2,&b2s,&b3s,&b2ssur,&b3ssur,&nn,&freesur);
//             printf("id=%lf\nx=%lf\n",id,x);
            if(id==1)
            {
                x1=x;
                y1=y;
                z1=z;
                l_1=l;
                na1_1=na1;
                na2_1=na2;
                nb1_1=nb1;
                nb2_1=nb2;
                b2s_1=b2s;
                b3s_1=b3s;
                b2ssur_1=b2ssur;
                b3ssur_1=b3ssur;
            }
            if(id==2)
            {
                x2=x;
                y2=y;
                z2=z;
                l_2=l;
                na1_2=na1;
                na2_2=na2;
                nb1_2=nb1;
                nb2_2=nb2;
                b2s_2=b2s;
                b3s_2=b3s;
                b2ssur_2=b2ssur;
                b3ssur_2=b3ssur;
            }
        }
        if(i<=t/2) {fprintf(out3,"%lf\n",z1);}
        if(i>t/2) {fprintf(out4,"%lf\n",z1);}
    }

    fclose(in);
    fclose(out3);
    fclose(out4);
    exit(0);
}
