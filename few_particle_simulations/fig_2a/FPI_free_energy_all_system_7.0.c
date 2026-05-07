# include <stdio.h>
# include <math.h>
# define SQ(x) ((x)*(x))

int main()
{
    int Np;
    double timestep;
    int Lx,Ly,Lz;
    int i,j,k,a,b,c,d,e,f;
    int id,type;
    double x,y,z,Fx,Fy,Fz,nbs,nl,a1s2,a2s2,a1s3,a2s3,b1s2,b2s2,b1s3,b2s3,nn,nfa1,nfa2,nfb1,nfb2,lfbs,omega_i,omega_S;
    double x1,y1,z1,x2,y2,z2,displacement;
    double F_i_ss,F_i_aa,F_i_sa,F_i_sb,F_i_bb,F_ij_ab,F_ij_aa,F_ij_bb,F_gas;
    double nbsA,nlA,a1s2A,a2s2A,a1s3A,a2s3A,free_a1,free_a2,A3S_surface,nbsB,nlB,b1s2B,b2s2B,b1s3B,b2s3B,free_b1,free_b2,nb3s,nb3s_A,nb3s_B,N_S_b1,nfs_c;

    double nbs_1,nl_1,a1s2_1,a2s2_1,a1s3_1,a2s3_1,b1s2_1,b2s2_1,b1s3_1,b2s3_1,nb3s_1,free_a1_1,free_a2_1,free_b1_1,free_b2_1,A3S_surface_1,omega_i_1;
    double nbs_2,nl_2,a1s2_2,a2s2_2,a1s3_2,a2s3_2,b1s2_2,b2s2_2,b1s3_2,b2s3_2,nb3s_2,free_a1_2,free_a2_2,free_b1_2,free_b2_2,A3S_surface_2,omega_i_2;
    double F_i_aa_1,F_i_aa_2,F_T_inf_1,F_T_inf_2;

    double NR=5000.0;
//     double N_S_b1=500.0;
    double NL=50;
//     double expo=0.798838;
    double omega_0 = 4 * 3.14 * 5 * 5 * 1;
    double omega_0s = 50 * 50 * 1;
//     double omega_is = 17.802358;
    double omega_i_A,omega_i_B,F_T_inf_sur,dist,F_T_inf_A,F_T_inf_B,e_ij,v_pp,e_is,v_ps;
    double F;
    char name1[50],name2[50],name3[50],name4[50],name5[50],name6[50],name7[50],name8[50],name9[50],name10[50],name11[50],name12[50],name13[50],name14[50],name15[50],name16[50],name17[50],name18[50],name19[50],name20[50],name21[50],name22[50],name23[50],name24[50],name25[50],name26[50],name27[50],name28[50],name29[50],name30[50];
    FILE *in;
    FILE *out;
    in=fopen("k0k1_blue_clean_7.0.dump","r");
    //out=fopen("free_energy_analytical.txt","w");

   // for(i=0;i<1;i++)
    //for(i=0;i<1;i++)
    {
        fscanf(in,"%s %s",name1,name2);
        fscanf(in,"%lf",&timestep);
        fscanf(in,"%s %s %s %s",name1,name2,name3,name4);
        fscanf(in,"%d",&Np);
        fscanf(in,"%s %s %s %s %s %s",name1,name2,name3,name4,name5,name6);
        fscanf(in,"%d %d %d",&a,&Lx,&b);
        fscanf(in,"%d %d %d",&c,&Ly,&d);
        fscanf(in,"%d %d %d",&e,&Lz,&f);
        fscanf(in,"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",name1,name2,name3,name4,name5,name6,name7,name8,name9,name10,name11,name12,name13,name14,name15,name16,name17,name18,name19,name20,name21,name22,name23,name24,name25,name26,name27,name28,name29,name30);

        for(j=0;j<1;j++)
        {
            fscanf(in,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&type,&id,&x,&y,&z,&Fx,&Fy,&Fz,&nbs,&nl,&a1s2,&a2s2,&a1s3,&a2s3,&b1s2,&b2s2,&b1s3,&b2s3,&nb3s,&nn,&nfa1,&nfa2,&nfb1,&nfb2,&lfbs,&nfs_c,&omega_i,&omega_S);
            if(id==1)
            {
                    x1=x;
                    y1=y;
                    z1=z;
                    nbs_1=nbs;
                    nl_1=nl;
                    a1s2_1=a1s2;
                    a2s2_1=a2s2;
                    a1s3_1=a1s3;
                    a2s3_1=a2s3;
                    b1s2_1=b1s2;
                    b2s2_1=b2s2;
                    b1s3_1=b1s3;
                    b2s3_1=b2s3;
                    nb3s_1=nb3s;
                    free_a1_1=nfa1;
                    free_a2_1=nfa2;
                    free_b1_1=nfb1;
                    free_b2_1=nfb2;
                    A3S_surface_1=lfbs;
                    N_S_b1=nfs_c;
                    omega_i_1 = omega_i;
            }
            if(id==2)
            {
                    x2=x;
                    y2=y;
                    z2=z;
                    nbs_2=nbs;
                    nl_2=nl;
                    a1s2_2=a1s2;
                    a2s2_2=a2s2;
                    a1s3_2=a1s3;
                    a2s3_2=a2s3;
                    b1s2_2=b1s2;
                    b2s2_2=b2s2;
                    b1s3_2=b1s3;
                    b2s3_2=b2s3;
                    nb3s_2=nb3s;
                    free_a1_2=nfa1;
                    free_a2_2=nfa2;
                    free_b1_2=nfb1;
                    free_b2_2=nfb2;
                    A3S_surface_2=lfbs;
                    N_S_b1=nfs_c;
                    omega_i_2 = omega_i;
            }
        }

//         dist = sqrt(SQ(x2-x1)+SQ(y2-y1)+SQ(z2-z1));
        F_gas= -458.726474;
//         F_gas= -967.143248;
//         F_i_ss= NR*log10(N_S_b1/NR);
        F_i_ss= NR*log(N_S_b1/NR);
        F_i_aa_1= NL*log(free_a1_1*free_a2_1/SQ(NL))+nl_1;
//         F_i_aa_2=NL*log(free_a1_2*free_a2_2/SQ(NL))+nl_2;
//         F_i_bb=NL*log(free_b1*free_b2/SQ(NL))+nlB;
        F_i_sa=nbs_1+2*A3S_surface_1;
//         F_i_sb=nbsB;
//         F_ij_ab=a1s2A+a2s2A+2*nb3s_A;
//         F_ij_aa=a1s2_1+a2s2_1;
//         F_ij_bb=b1s2B+b2s2B;
//         F_ij_aa=0.0;
//         F_ij_bb=0.0;
        F = F_i_ss + F_i_aa_1 + F_i_sa - F_gas ;

//         printf("F_ij_aa=%lf\n",a1s2_1);
        if (z1<6)
        {
                F_T_inf_sur = NR * log(omega_0s/omega_S) + 2*NL * log(omega_0/omega_i_1);
                F = F + F_T_inf_sur;
        }

//         if(dist<11.0)
//         {
//                 F_T_inf_1 = 2*NL * log(omega_0/omega_i_1);
//                 F_T_inf_2 = 2*NL * log(omega_0/omega_i_2);
//                 F = F + F_T_inf_1 + F_T_inf_2;
// //                 printf("%lf %lf\n",F_T_inf_A,F_T_inf_B);
//         }
//
//         if (dist<10.75)
//         {
//                 e_ij = (3.14/(12*dist)) * (10.75-dist) * (10.75-dist) * (dist*dist+11.5*dist+10*dist-1.6875);
//                 v_pp = 2*500*log(1-(e_ij/235.5));
//                 F = F - v_pp;
//         }

        if (z1 < 5.75)
        {
                e_is=(3.14 * ((5.75-z1)*(5.75-z1)) * (2*5.75+z1))/3;
//                 printf("dist=%lf e_is=%lf",dist,e_is);
                v_ps=2*500*(log(1-(e_is/(4*3.14*5*5*0.75))));
//                 printf("v_pp=%lf",v_pp);
                F = F - v_ps;
        }

        printf("%lf %lf\n",z1,F);
//         printf("");


        //printf("nbsA=%lf nlA=%lf a1s2A=%lf a2s2A=%lf a1s3A=%lf a2s3A=%lf free_a1=%lf free_a2=%lf A3S_surface=%lf\n",nbsA,nlA,a1s2A,a2s2A,a1s3A,a2s3A,free_a1,free_a2,A3S_surface);
        //printf("nbsB=%lf nlB=%lf b1s2B=%lf b2s2B=%lf b1s3B=%lf b2s3B=%lf free_b1=%lf free_b2=%lf\n",nbsB,nlB,b1s2B,b2s2B,b1s3B,b2s3B,free_b1,free_b2);

    }
    fclose(in);
    //fclose(out);
return 0;
}
