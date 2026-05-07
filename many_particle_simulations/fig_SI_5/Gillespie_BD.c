# include <stdio.h>
# include <unistd.h>
# include <sys/types.h>
# include <math.h>
# include <stdlib.h>
# include <time.h>
# include <fstream>
# include <iostream>
# include <vector>
# include <algorithm>
# define maxNp 500
# define maxN 75
# define SQ(x) ((x)*(x))
# define CUB(x) ((x)*(x)*(x))
using namespace std;
void distance_ads(double,double,double,int,double,double,double,double *,double *,double *,double *);
void overlap(double,double,double,double *);
void overlap_prime(double,double,double,double *);
void overlap_surface(double,double,double,double,double *);
void overlap_surface_prime(double,double,double,double,double *);
void overlapvolume(int,double,double,double,double);
void overlapvolume_surface(int ,double ,double ,double ,double , double *);
void normal(double *);
void neighborlist(int Np,int *ns,int *ns_upper,double Lx,double Ly,double Lz,double Sz, double Sz_upper, double Rc,double L,double R2s);
void distance(int,int,double,double,double,double *,double *,double *,double *);
void initialisation(int Np,int ns,int ns_upper);
void initialisation_after_removal(int Np,int ns,int ns_upper);
void Gillespie(int Np,double rho0,double Delta_tB, double kon0, double beta, double DG0, double DG0_L, double DG0_S,int ns,double *nfS_C,double omega_S, double tau0, double Delta_G0,double Delta_GT, double K0, double K1, double n_alpha, double n_T);
void Force(int Np, int ns, int ns_upper, double Rc, double L, double R2s,double Sz, double Sz_upper, double omega_S, double nfS_C_old, double beta,double Lx, double Ly, double Lz, double omega_0_inert);
double Isolated_colloids(int NA, int NB, double DG0_L_star,double beta, double rho0, double omega_0);
void Affinity(int Np,double rho0,double Delta_tB, double kon0, double beta, double DG0, double DG0_L, double DG0_S,int ns,double *nfS_C,double omega_S,double tau0,double Delta_G0,double Delta_GT,double K0,double K1, double n_alpha, double n_T);

typedef struct 
{	
	int nn_index; 			              //index of neighbors
    int nb;
	int nb_AB,nb_BA;				      //number of bridges
	double nb_a1a2,nb_a1b1,nb_a2b2,nb_a2a1,nb_b2b1,nb_a1a2b1,nb_a1a2b2,nb_b1b2,nb_b1b2a1,aon_b1b2_a2,nb_b1b2a2;
    double nl_a1a2,nl_b1b2;
	int n_index;
	double aoff,aon,aon_AB,aon_BA,aoff_AB,aoff_BA,aoff_L,aon_L,aon_a1a2,aoff_a1a2,aon_a1b1,aoff_a1b1,aon_b1b2,aoff_b1b2,aon_L_a1a2,aon_L_b1b2,aoff_L_a1a2,aoff_L_b1b2,aon_a2b2,aoff_a2b2,aon_a1a2_b1,aoff_a1a2_b1,aon_a1_b1_a2,aoff_a1_b1_a2,aon_a1a2_b2,aoff_a1a2_b2,aon_a2_b2_a1,aoff_a2_b2_a1,aon_b1b2_a1,aoff_b1b2_a1,aon_b1_a1_b2,aoff_b1_a1_b2,aoff_b1b2_a2,aon_b2_a2_b1,aoff_b2_a2_b1;
	double poff,pon,poff_AB,pon_AB,poff_BA,pon_BA,poff_L,pon_L;
    double aon_a2a1,aoff_a2a1,aon_b2b1,aoff_b2b1,pon_a2a1,poff_a2a1,pon_b2b1,poff_b2b1;

    double pon_a1a2,poff_a1a2,pon_a1b1,poff_a1b1,pon_a2b2,poff_a2b2,pon_b1b2,poff_b1b2,pon_a1a2_b1,poff_a1a2_b1,pon_a1_b1_a2,poff_a1_b1_a2,pon_a1a2_b2,poff_a1a2_b2,pon_a2_b2_a1,poff_a2_b2_a1,pon_b1b2_a1,poff_b1b2_a1,pon_b1_a1_b2,poff_b1_a1_b2,pon_b1b2_a2,poff_b1b2_a2,pon_b2_a2_b1,poff_b2_a2_b1;

    double tot_2strand_bridges,tot_3strand_bridges;
    double d;
    double e_ij;
    double e_ij_prime;
    double omega_ij;
    double omega_ij_prime;
    double regu_e_ij,regu_e_ij_prime;
    double nl;
}Pairs;
Pairs Pn[maxNp];
typedef struct 
{
	double x,y,z;
    int nf;
	int nf_A,nf_B,maxn,maxnn,count_regularization;			//Maximum number of free linkers, Maximum number of neighbors
    double nf_a1,nf_a2,nf_b1,nf_b2;
    int Sc,Sc_upper;
	double Sd,Sd_upper; 
	Pairs Pn[maxNp];
    double omega_i,omega_i_inert;
    int type;
    double Fx,Fy,Fz;
    double e_is;
    double e_is_prime;
    double e_si;
    double e_si_prime;
    double omega_is;
    double omega_is_prime;
    double aon_AC_S,aoff_AC_S,pon_AC_S,poff_AC_S,aon_a1C_S,aoff_a1C_S,aon_b2C_S,aoff_b2C_S,aon_a1a2_C_S,aoff_a1a2_C_S;
    double nbS_a1C,nbS_a1a2C,nbS_b2C;
    double nbS_AC;
    double aon_a1_C_a2_S,aoff_a1_C_a2_S;
    double total_free_loops;
    double total_bridge_WC;
    double tot_2strand_bridges_surface,tot_3strand_bridges_surface;
    double A_on_a1a2_sum,A_off_a1a2_sum,A_on_a1b1_sum,A_off_a1b1_sum,A_on_a2b2_sum,A_off_a2b2_sum,A_on_b1b2_sum,A_off_b1b2_sum,A_on_a1a2_b1_sum,A_off_a1a2_b1_sum,A_on_a1_b1_a2_sum,A_off_a1_b1_a2_sum,A_on_a1a2_b2_sum,A_off_a1a2_b2_sum,A_on_a2_b2_a1_sum,A_off_a2_b2_a1_sum,A_on_b1b2_a1_sum,A_off_b1b2_a1_sum,A_on_b1_a1_b2_sum,A_off_b1_a1_b2_sum,A_on_b1b2_a2_sum,A_off_b1b2_a2_sum,A_on_b2_a2_b1_sum,A_off_b2_a2_b1_sum,A_on_sum_L_a1a2,A_off_sum_L_a1a2,A_on_sum_L_b1b2,A_off_sum_L_b1b2;
    double A_on_a2a1_sum,A_off_a2a1_sum,A_on_b2b1_sum,A_off_b2b1_sum,P_on_a2a1_sum,P_off_a2a1_sum,P_on_b2b1_sum,P_off_b2b1_sum;
    double A_on_a1C_sum_S,A_off_a1C_sum_S,A_on_b2C_sum_S,A_off_b2C_sum_S,A_on_a1a2C_sum_S,A_off_a1a2C_sum_S,A_on_a1_C_a2_sum_S,A_off_a1_C_a2_sum_S;
    double A_tot_PP_2strand,A_tot_PP_3strand,A_tot_loop,A_tot_PP;
    double A_tot_sur_2strand,A_tot_sur_3strand,A_tot_PS,A_tot_PPS;
    double P_on_a1a2_sum,P_off_a1a2_sum,P_on_a1b1_sum,P_off_a1b1_sum,P_on_a2b2_sum,P_off_a2b2_sum,P_on_b1b2_sum,P_off_b1b2_sum,P_on_a1a2_b1_sum,P_off_a1a2_b1_sum,P_on_a1_b1_a2_sum,P_off_a1_b1_a2_sum,P_on_a1a2_b2_sum,P_off_a1a2_b2_sum,P_on_a2_b2_a1_sum,P_off_a2_b2_a1_sum,P_on_b1b2_a1_sum,P_off_b1b2_a1_sum,P_on_b1_a1_b2_sum,P_off_b1_a1_b2_sum,P_on_b1b2_a2_sum,P_off_b1b2_a2_sum,P_on_b2_a2_b1_sum,P_off_b2_a2_b1_sum,P_on_sum_L_a1a2,P_off_sum_L_a1a2,P_on_sum_L_b1b2,P_off_sum_L_b1b2;
    double P_on_a1C_sum_S,P_off_a1C_sum_S,P_on_b2C_sum_S,P_off_b2C_sum_S,P_on_a1a2C_sum_S,P_off_a1a2C_sum_S,P_on_a1_C_a2_sum_S,P_off_a1_C_a2_sum_S;
    double P_tot_PP_2strand,P_tot_PP_3strand,P_tot_loop;
    double P_tot_sur_2strand,P_tot_sur_3strand;
    double tot_a1a2_brg,tot_a2a1_brg,tot_a1b1_brg,tot_a2b2_brg,tot_b1b2_brg,tot_b2b1_brg,tot_a1a2b1_brg,tot_a1a2b2_brg,tot_b1b2a1_brg,tot_b1b2a2_brg,tot_a1,tot_a2,tot_b1,tot_b2;
}Particles;
Particles P[maxNp];
Particles P_cpy[maxNp];

int main()
{   
    int Np=0; ////currently checking with 2 particles
// 	int N=80;
    int N=150;
    int NA=75;
    int NB=75;
    int Na1=75;
    int Na2=75;
    int Nb1=75;
    int Nb2=75;
//     double NsC=0;
    double NsC=1000.0;
    double Rc=5.0;
    double R2s=1.0;
	double L=1.0;
	double Lx=60.0;
	double Ly=60.0;
	double Lz=60.0;
    double beta=1.0;
	double Delta_tB=0.001;
//     double Delta_tB=1.0;
//     double t_tot=1000000000.0;
    double t_tot=5000001.0;
//     double t_tot=500001.0;
    double kon0=1.0;
    double DG0=-9.0;
    double DG0_L=-9.0;
    double DG0_S=-10.0;
    double DG0_surface=-10.0;
    double rho_ideal=1.0e-05;
    double rho0=1.0;
	double dx,dy,dz,d,X,Y,Z;
	double t,t_bar,tao;

    double tau0=1.0;
    double Delta_G0=-20.0;
    double Delta_GT=-2.0;
    double K0=exp(-(beta*Delta_G0)-(beta*Delta_GT));
    double K1=exp(-(beta*Delta_GT));
    printf("K0=%lf K1=%lf\n",K0,K1);
    double n_alpha = 3.0;
    double n_T = 1.0;

    int stop=1;//////to stop a calculation
	
    int n_deposited_particles, n_checked_particles; 
    double tot_bridges;
    
    int count,count_nn,count_n;
    int seed = time(NULL)*getpid(); //596462648; //
    srand (seed);
    printf("seed = %d\n",seed);
	
	int i,j,k,E,q;
	bool flag;
    ofstream myfile;
    ofstream myfile1;
    ofstream myfile2;
//     myfile.open ("Particle_gilles.txt");
    myfile1.open ("long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_modified.dump");
//     myfile2.open ("Bridge_loop.dump");

    double nfS_C=NsC;
    double e_is,e_is_prime,e_si,e_si_prime,w_is,w_is_prime,e_si_sum,omega_S;
    double omega_S_0=Lx*Ly*R2s;
    double omega_0=4*M_PI*SQ(Rc)*L; 
    double omega_0_inert=4*M_PI*SQ(Rc)*0.75*L;
    
    double V=Lx*Ly*Lz;
    double DG0_surface_star=beta*DG0_surface-log(1/rho0*CUB(L));
    double DG0_star=beta*DG0-log(1/rho0*CUB(L));
    double DG0_L_star=beta*DG0_L-log(1/rho0*CUB(L));
    
    double FreeEnergy_bulk;
    FreeEnergy_bulk = Isolated_colloids(NA, NB, DG0_L_star,beta,rho0,omega_0);
    //exit(0);
    double e=2.71828;
    double nfS_C_cpy, omega_S_cpy, ns_cpy, ns_upper_cpy;
    for(i=0; i<maxNp; i++)
    {
        for(j=0; j<maxNp; j++)
        {
            P[i].Pn[j].nn_index=-1;
        }
        P[i].Sc=-1;
        P[i].Sc_upper=-1;
    }

    for(i=0; i<Np; i++)
    {
        if(i%2==0)
        {
            P[i].type=0;
        }
        else
        {
            P[i].type=1;
        }

    }

/***********Particle position *************************************/

    n_deposited_particles=0;
    flag=false;
    while (flag==false)
    {
       for (i=0; i<1; i++)
       {
            P[i].x =(rand()/(double)RAND_MAX)*Lx;  //Lx
            P[i].y =(rand()/(double)RAND_MAX)*Ly;
            P[i].z =(rand()/(double)RAND_MAX)*Lz;
            if (P[i].z > Rc && P[i].z <(Lz-Rc)) flag=true;
        }
    }


    n_deposited_particles=1;
//     n_deposited_particles=Np;

    while (n_deposited_particles < Np)
    {
        X =(rand()/(double)RAND_MAX)*Lx;
        Y =(rand()/(double)RAND_MAX)*Ly;
        Z =(rand()/(double)RAND_MAX)*Lz;

        n_checked_particles = 0;
        for (i=0; i<n_deposited_particles; i++)
        {
            distance_ads(X,Y,Z,i,Lx,Ly,Lz,&dx,&dy,&dz,&d);

            if (d > 2*Rc && Z>Rc && Z<(Lz-Rc)) {
                n_checked_particles = n_checked_particles + 1;
            }
            else if (d <= 2*Rc || Z <= Rc || Z>=(Lz-Rc) ) break;
        }
        
        if (n_checked_particles == n_deposited_particles){

            P[n_deposited_particles].x=X;
            P[n_deposited_particles].y=Y;
            P[n_deposited_particles].z=Z;

            n_deposited_particles = n_deposited_particles + 1;
        }
    }
    
//     for (i=0; i<Np; i++)
//     {
//         P[0].x = 19;
//         P[0].y = 30;
//         P[0].z = 6;
//         P[0].type = 0;
//
//         P[1].x = 30;
//         P[1].y = 30;
//         P[1].z = 6.0;
//         P[1].type = 0;
//
//         P[2].x = 19;
//         P[2].y = 41;
//         P[2].z = 6;
//         P[2].type = 0;
// 
//         P[3].x = 30;
//         P[3].y = 41;
//         P[3].z = 6;
//         P[3].type = 0;
//
//         P[4].x = 19;
//         P[4].y = 30;
//         P[4].z = 17;
//         P[4].type = 1;
//
//         P[5].x = 30;
//         P[5].y = 30;
//         P[5].z = 17;
//         P[5].type = 1;
//
//         P[6].x = 19;
//         P[6].y = 41;
//         P[6].z = 17;
//         P[6].type = 1;
// 
//         P[7].x = 30;
//         P[7].y = 41;
//         P[7].z = 17;
//         P[7].type = 1;
//
//         P[8].x = 19;
//         P[8].y = 30;
//         P[8].z = 28;
//         P[8].type = 0;
//
//         P[9].x = 30;
//         P[9].y = 30;
//         P[9].z = 28;
//         P[9].type = 0;
//
//         P[10].x = 19;
//         P[10].y = 41;
//         P[10].z = 28;
//         P[10].type = 0;
// 
//         P[11].x = 30;
//         P[11].y = 41;
//         P[11].z = 28;
//         P[11].type = 0;
//     }
    
    
    
//    for (i=0; i<Np; i++)
//     {
//         P[0].x = 19;
//         P[0].y = 30;
//         P[0].z = 6;
//         P[0].type = 0; 
//        
//         P[1].x = 30;
//         P[1].y = 41;
//         P[1].z = 28;
//         P[1].type = 1; //delete
//
//         P[2].x = 30;
//         P[2].y = 30;
//         P[2].z = 6.0;
//         P[2].type = 1;
//
//         
//     }
    
    for (i=0; i<Np; i++) ////commemnted out for timew being
    {
        if(P[i].type==0)
        {
//             P[i].nf_a1=Na1;
//             P[i].nf_a2=Na2;
            P[i].Pn[i].nl_a1a2=75.0;
            P[i].nf_a1=0.0;
            P[i].nf_a2=0.0;
        }
        if(P[i].type==1)
        {
//             P[i].nf_b1=Nb1;
//             P[i].nf_b2=Nb2;
            P[i].Pn[i].nl_b1b2=75.0;
            P[i].nf_b1=0.0;
            P[i].nf_b2=0.0;
//             P[i].nbS_b2C=0.0;
        }
    }

    double Sx=Lx/2.0;
    double Sy=Ly/2.0;
    double Sz=0.0;
    double Sz_upper=Lz;
    double ov; //   , ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
/********neighbor calculations*********************************************************************/
    t=0.0;
    int count_BD_steps=0;
	while (t < t_tot)  //For Brownian dynamics; for the time being it is commented out
	{
        int ns;
        int ns_upper;
        neighborlist(Np,&ns,&ns_upper,Lx,Ly,Lz,Sz, Sz_upper, Rc,L,R2s);
        overlapvolume(Np,Rc,L,omega_0,omega_0_inert);
        overlapvolume_surface(ns,Rc,L,R2s,omega_S_0,&omega_S);
        Gillespie(Np,rho0,Delta_tB, kon0, beta, DG0, DG0_L, DG0_S, ns, &nfS_C, omega_S, tau0, Delta_G0, Delta_GT, K0, K1, n_alpha, n_T);
        //cout << "ns= " << ns << endl;
        double Free_linkers=0.0;
        double total_surfacebridge=0.0;
        
/********Force Calculations************************************************************************/

        int i_new,Np_new;
        double insetion_prob,remove_prob;
        double rand_BD_GC=rand()/(double)RAND_MAX;
//         if (rand_BD_GC>1/1000.0)
        if (rand_BD_GC>1/1000.0 && stop == 1)
        {
            Force(Np,ns, ns_upper, Rc, L, R2s, Sz, Sz_upper, omega_S, nfS_C, beta, Lx, Ly, Lz, omega_0_inert);
            double Rx,Ry,Rz,Rn;
            double X1_new,X2_new,Y1_new,Y2_new,Z1_new,Z2_new;
            double D;
            double X_old,Y_old, Z_old;
            for (i=0; i<Np; i++)
            {
                normal(&Rx);
                normal(&Ry);
                normal(&Rz);
//                P[i].Fx=0.0;
//                P[i].Fy=0.0;
//                P[i].Fz=0.0;
//                Rx=0;Ry=0;Rz=0;
//                X_old=P[i].x;
//                Y_old=P[i].y;
//                Z_old=P[i].z;
                double X_new = P[i].x + beta*P[i].Fx*Delta_tB*L + sqrt(2*Delta_tB)*Rx;
                double Y_new = P[i].y + beta*P[i].Fy*Delta_tB*L + sqrt(2*Delta_tB)*Ry;
                double Z_new = P[i].z + beta*P[i].Fz*Delta_tB*L + sqrt(2*Delta_tB)*Rz;

                if      (X_new <  0)    X_new = X_new + Lx;
                else if (X_new >  Lx)   X_new = X_new - Lx;
            
                if      (Y_new <  0)    Y_new = Y_new + Ly;
                else if (Y_new >  Ly)   Y_new = Y_new - Ly;
                
                P[i].x=X_new;
                P[i].y=Y_new;
                P[i].z=Z_new;
            }
            
            for(i=0; i<Np; i++)
            {
                for (j=0; j<P[i].maxnn; j++)
                {
                    int index=P[i].Pn[j].nn_index;
                    distance(i,index,Lx,Ly,Lz,&dx,&dy,&dz,&d);
                    if (d>12.0)
                    {
                        if (P[i].type == 0 && P[index].type == 0 && i<index)
                        {
                        P[i].nf_a1=P[i].nf_a1+P[i].Pn[index].nb_a1a2;
                        P[index].nf_a2=P[index].nf_a2+P[i].Pn[index].nb_a1a2;

//                         P[i].nf_a2=P[i].nf_a2+P[i].Pn[index].nb_a2a1;
//                         P[index].nf_a1=P[index].nf_a1+P[i].Pn[index].nb_a2a1;
                        }
                        
                        
                        if (P[i].type == 0 && P[index].type == 0 && i>index)
                        {
//                        P[i].nf_a1=P[i].nf_a1+P[i].Pn[index].nb_a1a2;
//                        P[index].nf_a2=P[index].nf_a2+P[i].Pn[index].nb_a1a2;

                         P[i].nf_a1=P[i].nf_a1+P[i].Pn[index].nb_a2a1;
                         P[index].nf_a2=P[index].nf_a2+P[i].Pn[index].nb_a2a1;
                        }
                        
                        if (P[i].type == 0 && P[index].type == 1)
                        {
                        P[i].nf_a1=P[i].nf_a1+P[i].Pn[index].nb_a1b1;
                        P[index].nf_b1=P[index].nf_b1+P[i].Pn[index].nb_a1b1;

                        P[i].nf_a2=P[i].nf_a2+P[i].Pn[index].nb_a2b2;
                        P[index].nf_b2=P[index].nf_b2+P[i].Pn[index].nb_a2b2;

                        P[i].Pn[i].nl_a1a2=P[i].Pn[i].nl_a1a2+P[i].Pn[index].nb_a1a2b1;
                        P[index].nf_b1=P[index].nf_b1+P[i].Pn[index].nb_a1a2b1;

                        P[i].Pn[i].nl_a1a2=P[i].Pn[i].nl_a1a2+P[i].Pn[index].nb_a1a2b2;
                        P[index].nf_b2=P[index].nf_b2+P[i].Pn[index].nb_a1a2b2;
                        }

                        if (P[i].type == 1 && P[index].type == 0)
                        {
                        P[i].Pn[i].nl_b1b2=P[i].Pn[i].nl_b1b2+P[i].Pn[index].nb_b1b2a1;
                        P[index].nf_a1=P[index].nf_a1+P[i].Pn[index].nb_b1b2a1;

                        P[i].Pn[i].nl_b1b2=P[i].Pn[i].nl_b1b2+P[i].Pn[index].nb_b1b2a2;
                        P[index].nf_a2=P[index].nf_a2+P[i].Pn[index].nb_b1b2a2;
                        }

//                        if (P[i].type == 1 && P[index].type == 1)
//                        {
//                        P[i].nf_b1=P[i].nf_b1+P[i].Pn[index].nb_b1b2;
//                        P[index].nf_b2=P[index].nf_b2+P[i].Pn[index].nb_b1b2;
//                        }
                        
                        
                        if (P[i].type == 1 && P[index].type == 1 && i<index)
                        {
                        P[i].nf_b1=P[i].nf_b1+P[i].Pn[index].nb_b1b2;
                        P[index].nf_b2=P[index].nf_b2+P[i].Pn[index].nb_b1b2;

//                         P[i].nf_a2=P[i].nf_a2+P[i].Pn[index].nb_a2a1;
//                         P[index].nf_a1=P[index].nf_a1+P[i].Pn[index].nb_a2a1;
                        }
                        
                        
                        if (P[i].type == 1 && P[index].type == 1 && i>index)
                        {
//                        P[i].nf_a1=P[i].nf_a1+P[i].Pn[index].nb_a1a2;
//                        P[index].nf_a2=P[index].nf_a2+P[i].Pn[index].nb_a1a2;

                         P[i].nf_b1=P[i].nf_b1+P[i].Pn[index].nb_b2b1;
                         P[index].nf_b2=P[index].nf_b2+P[i].Pn[index].nb_b2b1;
                        }
                        
                        

                        if (P[i].type == 0 && P[index].type == 0)
                        {
                        P[i].Pn[index].nb_a1a2=0.0;
                        //P[index].Pn[i].nb_a1a2=0.0;
                         P[i].Pn[index].nb_a2a1=0.0;
//                         P[index].Pn[i].nb_a2a1=0.0;
                        }

                        if (P[i].type == 0 && P[index].type == 1)
                        {
                        P[i].Pn[index].nb_a1b1=0.0;
                        P[index].Pn[i].nb_a1b1=0.0;
                        P[i].Pn[index].nb_a2b2=0.0;
                        P[index].Pn[i].nb_a2b2=0.0;
                        P[i].Pn[index].nb_a1a2b1=0.0;
                        P[index].Pn[i].nb_a1a2b1=0.0;
                        P[i].Pn[index].nb_a1a2b2=0.0;
                        P[index].Pn[i].nb_a1a2b2=0.0;
                        }

                        if (P[i].type == 1 && P[index].type == 0)
                        {
                        P[i].Pn[index].nb_b1b2a1=0.0;
                        P[index].Pn[i].nb_b1b2a1=0.0;
                        P[i].Pn[index].nb_b1b2a2=0.0;
                        P[index].Pn[i].nb_b1b2a2=0.0;
                        }

                        if (P[i].type == 1 && P[index].type == 1)
                        {
                        P[i].Pn[index].nb_b1b2=0.0;
                        //P[index].Pn[i].nb_b1b2=0.0;
                         P[i].Pn[index].nb_b2b1=0.0;
//                         P[index].Pn[i].nb_b2b1=0.0;
                        }


//
// //                         P[i].nf_B=P[i].nf_B+P[i].Pn[index].nb_BA;
// //                         P[index].nf_A=P[index].nf_A+P[i].Pn[index].nb_BA;
//
//
// //                         P[i].Pn[index].nb_BA=0.0;
// //                         P[index].Pn[i].nb_BA=0.0;

                        if (P[i].type == 1 && P[index].type == 0)
//                         if (P[i].type == 1)
                        {
                            if (P[i].nf_b1 < P[i].nf_b2)
                            {
                                int n = P[i].nf_b1;
                                P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 + n;
                                P[i].nf_b1 = P[i].nf_b1 - n;
                                P[i].nf_b2 = P[i].nf_b2 - n;
                            }

                            if (P[i].nf_b1 > P[i].nf_b2)
                            {
                                int n = P[i].nf_b2;
                                P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 + n;
                                P[i].nf_b1 = P[i].nf_b1 - n;
                                P[i].nf_b2 = P[i].nf_b2 - n;
                            }
                        }
                    }
                }
            }
//
            for(q=0; q<ns; q++)
            {
                int index_cs=P[q].Sc;
                if ( P[index_cs].z > 7.0)
                {
//                     P[index_cs].nf_A = P[index_cs].nf_A+P[index_cs].nbS_AC;
//                     nfS_C = nfS_C+P[index_cs].nbS_AC;
//                     P[index_cs].nbS_AC=0.0;

                    if (P[index_cs].type == 0)
                    {
                    P[index_cs].nf_a1 = P[index_cs].nf_a1+P[index_cs].nbS_a1C;
                    nfS_C = nfS_C+P[index_cs].nbS_a1C;

                    P[index_cs].Pn[index_cs].nl_a1a2 = P[index_cs].Pn[index_cs].nl_a1a2+P[index_cs].nbS_a1a2C;
                    nfS_C = nfS_C+P[index_cs].nbS_a1a2C;
                    }

                    if (P[index_cs].type == 1)
                    {
                    P[index_cs].nf_b2 = P[index_cs].nf_b2+P[index_cs].nbS_b2C;
                    nfS_C = nfS_C+P[index_cs].nbS_b2C;
                    }

                    if (P[index_cs].type == 0)
                    {
                    P[index_cs].nbS_a1C=0.0;
                    P[index_cs].nbS_a1a2C=0.0;
                    }
                    if (P[index_cs].type == 1)
                    {
                    P[index_cs].nbS_b2C=0.0;
                    }
                }
            } 
        }
    
//         else
        else if(stop==1 && count_BD_steps % 100==0)
        {
            double Free_energy_old=0.0;
            double Free_energy_regularization=0.0;
            double Free_energy_regularization_surface=0.0;
            double omega_0_inert = 4*M_PI*SQ(Rc)*0.75;

            double rand_insertion_delete=rand()/(double)RAND_MAX;

            /************insertion *****************************************************/
            int m_test;
// 	    if (rand_insertion_delete <0.5) //insertion
            if (rand_insertion_delete < 0.5 && stop == 1)
            //if (stop == 0)
            {
                X =(rand()/(double)RAND_MAX)*Lx;
                Y =(rand()/(double)RAND_MAX)*Ly;
                Z =(rand()/(double)RAND_MAX)*Lz;
                
                
                int count_deposition_ne=0;
                for (i=0; i<Np; i++)
                {
                    distance_ads(X,Y,Z,i,Lx,Ly,Lz,&dx,&dy,&dz,&d);
//                     if ( ( d>2*(Rc+1.0) ) && Z>(Rc+1.0+0.0) && Z<(Lz-Rc))    count_deposition_ne=count_deposition_ne+1;
                    if ( ( d > 2*(Rc+0.5) ) && Z > (Rc+1.0+0.0) && Z < (Lz-Rc-5) )    count_deposition_ne = count_deposition_ne+1;
                    else break;
                }
      		    if (Np==0 && (Z<=(Rc+1) || Z>=(Lz-Rc-5))) m_test=1;
                else m_test=0;

                if (count_deposition_ne==Np && m_test==0)
                {

                    Np_new = Np+1;
                    i_new=Np_new-1;
                    P[i_new].x=X;P[i_new].y=Y;P[i_new].z=Z;

                    P[i_new].type = rand() % 2; /////type of inserted particles
		    //P[i_new].type = 1;

                    double Free_energy_new=0.0;
                    double Free_energy_regularization=0.0;
                    double omega_0_inert = 4*M_PI*SQ(Rc)*0.75;
                    double Free_energy_diff=0;
                    double insertion_prob=(V/(double) Np_new)*rho_ideal*e*exp(-Free_energy_diff);
                    double insertion_rand=rand()/(double)RAND_MAX;

                    if (insertion_prob > insertion_rand)
                    {
                        Np=Np_new;
//                         P[Np-1].nf_A=1;P[Np-1].nf_B=1;P[Np-1].Pn[Np-1].nl=39.0;P[Np-1].nbS_AC=0.0;
//                         if (P[Np-1].type == 0)
                        if (P[i_new].type == 0)
                        {
//                         P[Np-1].nf_a1=0;
//                         P[Np-1].nf_a2=0;
//                         P[Np-1].nf_b1=0;
//                         P[Np-1].nf_b2=0;
//                         P[Np-1].Pn[Np-1].nl_a1a2=50.0;
//                         P[Np-1].nbS_a1C=0.0;
//                         P[Np-1].nbS_a1a2C=0.0;

                            P[i_new].nf_a1=0;
                            P[i_new].nf_a2=0;
                            P[i_new].nf_b1=0;
                            P[i_new].nf_b2=0;
                            P[i_new].Pn[i_new].nl_a1a2=75.0;
                            P[i_new].nbS_a1C=0.0;
                            P[i_new].nbS_a1a2C=0.0;
                        }
//                         if (P[Np-1].type == 1)
                        if (P[i_new].type == 1)
                        {
//                             P[Np-1].nf_a1=0;
//                             P[Np-1].nf_a2=0;
//                             P[Np-1].nf_b1=0;
//                             P[Np-1].nf_b2=0;
//                             P[Np-1].Pn[Np-1].nl_b1b2=50.0;
//                             P[Np-1].nbS_b2C=0.0;

                            P[i_new].nf_a1=0;
                            P[i_new].nf_a2=0;
                            P[i_new].nf_b1=0;
                            P[i_new].nf_b2=0;
                            P[i_new].Pn[i_new].nl_b1b2=75.0;
                            P[i_new].nbS_b2C=0.0;
                        }
                    }
                    else
                    {
                        P[i_new].x=0.0;P[i_new].y=0.0;P[i_new].z=0.0;
                    }
                }
            }
            /************remove ***********************************/
//             else if (rand_insertion_delete >=0.5 && Np>0)
            else if (rand_insertion_delete >=0.5 && Np>0 && stop == 1)
            //if (count_BD_steps>10000 && stop == 1)
            {
                int Rand_particle=rand() % Np;
//                if (P[Np-1].type==0)
//                {
//                    flag=false;
//                    while (flag==false)
//                    {
//                        Rand_particle=rand() % Np;
//                        if (P[Rand_particle].type==1) flag=true;
//                    }
//                }
//                
//                if (P[Np-1].type==1)
//                {
//                    flag=false;
//                    while (flag==false)
//                    {
//                        Rand_particle=rand() % Np;
//                        if (P[Rand_particle].type==0) flag=true;
//                    }
//                }
                //Rand_particle=rand() % Np;
                
                
//                cout << " P[0].Pn[2].nb_a2b2 = " << " " << P[0].Pn[2].nb_a2b2 << " " << endl;
//                cout << " P[2].Pn[0].nb_a2b2 = " << " " << P[2].Pn[0].nb_a2b2 << " " << endl;
//                cout << " P[0].Pn[1].nb_a2b2 = " << " " << P[0].Pn[1].nb_a2b2 << " " << endl;
//                cout << " P[1].Pn[0].nb_a2b2 = " << " " << P[1].Pn[0].nb_a2b2 << " " << endl;
                
                //Rand_particle=1;
                
                X=P[Rand_particle].x; Y=P[Rand_particle].y; Z=P[Rand_particle].z;
                
                
                
                
//                P[0].x = 19;
//         P[0].y = 30;
//         P[0].z = 6;
//         P[0].type = 0; 
//        
//         P[1].x = 30;
//         P[1].y = 41;
//         P[1].z = 28;
//         P[1].type = 1; //delete
//
//         P[2].x = 30;
//         P[2].y = 30;
//         P[2].z = 6.0;
//         P[2].type = 1;
//                
                
                
                
                //cout << "Rand_particle= " << " " << Rand_particle << " " << X << " " << Y << " " << Z << " " << P[Rand_particle].type << endl;
                
                int count_deposition_ne=0;
                for (i=0; i<Np; i++)
                {
                    if ( i!=Rand_particle )
                    {
                        distance_ads(X,Y,Z,i,Lx,Ly,Lz,&dx,&dy,&dz,&d);
//                         if (  d>2*(Rc+0.0) && Z>(Rc+0.0+0.0) && Z<(Lz-Rc) )  count_deposition_ne=count_deposition_ne+1;
                        if (  d>2*(Rc+0.5) && Z>(Rc+1.0+0.0) && Z<(Lz-Rc-5) )  count_deposition_ne=count_deposition_ne+1;
                        else break;
                    }
                }

                int tot_bridges_removal=0;
                for (j=0; j<P[Rand_particle].maxnn; j++)
                {
                    
                    int index=P[Rand_particle].Pn[j].nn_index;
                    //                     tot_bridges_removal=tot_bridges_removal+P[Rand_particle].Pn[index].nb_AB+P[Rand_particle].Pn[index].nb_BA;
                    //tot_bridges_removal=tot_bridges_removal+P[Rand_particle].Pn[index].nb_a1a2+P[Rand_particle].Pn[index].nb_b1b2+P[Rand_particle].Pn[index].nb_a1b1+P[Rand_particle].Pn[index].nb_a2b2+P[Rand_particle].Pn[index].nb_a1a2b1+P[Rand_particle].Pn[index].nb_a1a2b2+P[Rand_particle].Pn[index].nb_b1b2a1+P[Rand_particle].Pn[index].nb_b1b2a2;
                    
                    
                    //P[0].Pn[1].nb_a1a2
                    //int index=P[Rand_particle].Pn[j].nn_index;
                    double aa=0.0;
                    if (index>Rand_particle)
                    {
                        aa=P[Rand_particle].Pn[index].nb_a1a2;
                    }
                    if (Rand_particle>index)
                    {
                        aa=P[Rand_particle].Pn[index].nb_a2a1;
                    }
                    
                    double bb=0.0;
                    if (index>Rand_particle)
                    {
                        bb=P[Rand_particle].Pn[index].nb_b1b2;
                    }
                    if (Rand_particle>index)
                    {
                        bb=P[Rand_particle].Pn[index].nb_b2b1;
                    }
//                    tot_2strand_bridges = tot_2strand_bridges
//                        + aa 
//                        + bb 
//                        + P[Rand_particle].Pn[index].nb_a1b1 
//                        + P[Rand_particle].Pn[index].nb_a2b2;
//                        // + P[i].Pn[index].nb_a2a1 + P[i].Pn[index].nb_b2b1;
//                    tot_3strand_bridges=tot_3strand_bridges+P[Rand_particle].Pn[index].nb_a1a2b1 + P[Rand_particle].Pn[index].nb_a1a2b2 + P[Rand_particle].Pn[index].nb_b1b2a1 + P[Rand_particle].Pn[index].nb_b1b2a2;

                    //tot_bridges = tot_2strand_bridges + tot_3strand_bridges;
                    
                    //tot_bridges_removal=tot_bridges_removal+tot_bridges;
                    tot_bridges_removal=tot_bridges_removal + aa + bb + P[Rand_particle].Pn[index].nb_a1b1 + P[Rand_particle].Pn[index].nb_a2b2+P[Rand_particle].Pn[index].nb_a1a2b1 + P[Rand_particle].Pn[index].nb_a1a2b2 + P[Rand_particle].Pn[index].nb_b1b2a1 + P[Rand_particle].Pn[index].nb_b1b2a2;
                
                }
                
                cout << "Rand_particle = " <<  Rand_particle << " " << "P[Rand_particle].nbS_a1C = " << " " << P[Rand_particle].nbS_a1C << " " <<  "P[Rand_particle].nbS_a1a2C = " << " " << P[Rand_particle].nbS_a1a2C << endl;

//                 if (count_deposition_ne==Np-1 && tot_bridges_removal==0 && P[Rand_particle].nbS_AC==0)
                if (count_deposition_ne==Np-1 && tot_bridges_removal==0 && P[Rand_particle].nbS_a1C==0 && P[Rand_particle].nbS_a1a2C==0 && P[Rand_particle].nbS_b2C==0)

                {
                    //printf("totalbridge=%d\n",tot_bridges_removal);
                    
                    Np_new=Np-1;
                    i_new=Np-1;
                    double Free_energy_new=0.0;
                    double Free_energy_diff=0.0;
                    double remove_prob = (Np/V)*(1/(rho_ideal*e))*exp(-Free_energy_diff);
                    double remove_rand=rand()/(double)RAND_MAX;
                    //remove_prob=1.0;
                    //remove_rand=0.1;
                    
                    
                    cout << "i_new = " <<  i_new << " " << "P[i_new].nbS_a1C = " << " " << P[i_new].nbS_a1C << " " <<  "P[i_new].nbS_a1a2C = " << " " << P[i_new].nbS_a1a2C << endl;
                    
                
                    cout << "remove_prob = " << " " << remove_prob << " " <<  "remove_rand = " << " " << remove_rand << endl;
                    if (remove_prob>remove_rand)
                    {
                        initialisation_after_removal(Np,ns,ns_upper);
                        
                        cout << "Rand_particle= " << " " << Rand_particle << " " << X << " " << Y << " " << Z << " " << P[Rand_particle].type << endl;
                        cout << "i_new= " << " " << i_new << " " << P[i_new].x << " " << P[i_new].y << " " << P[i_new].z << " " << P[i_new].type << endl;
                        P[Rand_particle].x=P[i_new].x;
                        P[Rand_particle].y=P[i_new].y;
                        P[Rand_particle].z=P[i_new].z;
                        P[Rand_particle].type=P[i_new].type;
                        cout << "Rand_particle= " << " " << Rand_particle << " " << P[Rand_particle].x << " " << P[Rand_particle].y << " " << P[Rand_particle].z << " " << P[Rand_particle].type << endl;
                        cout << "i_new= " << " " << i_new << " " << P[i_new].x << " " << P[i_new].y << " " << P[i_new].z << " " << P[i_new].type << endl;
//                         P[Rand_particle].type = rand() % 2;
//                         P[Rand_particle].nbS_AC=P[i_new].nbS_AC;
//                         if (P[Rand_particle].type==0)
                        //exit(0);
                        if (P[i_new].type==0)
                        {
                            P[Rand_particle].nf_a1=P[i_new].nf_a1;
                            P[Rand_particle].nf_a2=P[i_new].nf_a2;
                            P[Rand_particle].nf_b1=0;
                            P[Rand_particle].nf_b2=0;
                            P[Rand_particle].nbS_a1C=P[i_new].nbS_a1C;
                            P[Rand_particle].nbS_a1a2C=P[i_new].nbS_a1a2C;
                            //P[Rand_particle].Pn[Rand_particle].nl_a1a2=P[i_new].Pn[i_new].nl_a1a2
                            
                        }
//                         if (P[Rand_particle].type==1)
                        if (P[i_new].type==1)
                        {
                            P[Rand_particle].nf_a1=0;
                            P[Rand_particle].nf_a2=0;
                            P[Rand_particle].nf_b1=P[i_new].nf_b1;
                            P[Rand_particle].nf_b2=P[i_new].nf_b2;
                            P[Rand_particle].nbS_b2C=P[i_new].nbS_b2C;
                        }
//                         P[Rand_particle].nf_A=P[i_new].nf_A;
//                         P[Rand_particle].nf_B=P[i_new].nf_B;
                         if (P[i_new].type==0)
                         {
                         P[Rand_particle].Pn[Rand_particle].nl_a1a2=P[i_new].Pn[i_new].nl_a1a2;
                         }
                         if (P[i_new].type==1)
                         {
                         P[Rand_particle].Pn[Rand_particle].nl_b1b2=P[i_new].Pn[i_new].nl_b1b2;
                         }
                        P[Rand_particle].maxn=P[i_new].maxn;
                        P[Rand_particle].maxnn=P[i_new].maxnn;
                        for (j=0; j<P[Rand_particle].maxnn; j++)
                        {
                            int index=P[i_new].Pn[j].nn_index;		     
                            //index of neighbors
                            //                             P[Rand_particle].Pn[index].nb_AB = P[i_new].Pn[index].nb_AB;
//                             P[index].Pn[Rand_particle].nb_AB = P[index].Pn[i_new].nb_AB;
//
//                             P[Rand_particle].Pn[index].nb_BA = P[i_new].Pn[index].nb_BA;
//                             P[index].Pn[Rand_particle].nb_BA = P[index].Pn[i_new].nb_BA;

                            if (P[Rand_particle].type == 0 && P[index].type == 0) // && index>Rand_particle)
                            {
                                P[Rand_particle].Pn[index].nb_a1a2 = P[i_new].Pn[index].nb_a1a2;
                                P[index].Pn[Rand_particle].nb_a1a2 = P[index].Pn[i_new].nb_a1a2;
                            }
                            
                             if (P[Rand_particle].type == 0 && P[index].type == 0 ) //&& index<Rand_particle)
                            {
                                P[Rand_particle].Pn[index].nb_a2a1 = P[i_new].Pn[index].nb_a2a1;
                                P[index].Pn[Rand_particle].nb_a2a1 = P[index].Pn[i_new].nb_a2a1;
                            }
                            
                            if (P[Rand_particle].type == 1 && P[index].type == 1 ) //&& index>Rand_particle)
                            {
                                P[Rand_particle].Pn[index].nb_b1b2 = P[i_new].Pn[index].nb_b1b2;
                                P[index].Pn[Rand_particle].nb_b1b2 = P[index].Pn[i_new].nb_b1b2;
                            }
                            
                            if (P[Rand_particle].type == 1 && P[index].type == 1 ) //&& index<Rand_particle)
                            {
                                P[Rand_particle].Pn[index].nb_b2b1 = P[i_new].Pn[index].nb_b2b1;
                                P[index].Pn[Rand_particle].nb_b2b1 = P[index].Pn[i_new].nb_b2b1;
                            }

                            if ( (P[Rand_particle].type == 0 && P[index].type == 1) || (P[Rand_particle].type == 1 && P[index].type == 0))
                            {
                                P[Rand_particle].Pn[index].nb_a1b1 = P[i_new].Pn[index].nb_a1b1;
                                P[index].Pn[Rand_particle].nb_a1b1 = P[index].Pn[i_new].nb_a1b1;

                                P[Rand_particle].Pn[index].nb_a2b2 = P[i_new].Pn[index].nb_a2b2;
                                P[index].Pn[Rand_particle].nb_a2b2 = P[index].Pn[i_new].nb_a2b2;

                                P[Rand_particle].Pn[index].nb_a1a2b1 = P[i_new].Pn[index].nb_a1a2b1;
                                P[index].Pn[Rand_particle].nb_a1a2b1 = P[index].Pn[i_new].nb_a1a2b1;

                                P[Rand_particle].Pn[index].nb_a1a2b2 = P[i_new].Pn[index].nb_a1a2b2;
                                P[index].Pn[Rand_particle].nb_a1a2b2 = P[index].Pn[i_new].nb_a1a2b2;
                            }

                            if ((P[Rand_particle].type == 1 && P[index].type == 0) || (P[Rand_particle].type == 0 && P[index].type == 1))
                            {
                                P[Rand_particle].Pn[index].nb_b1b2a1 = P[i_new].Pn[index].nb_b1b2a1;
                                P[index].Pn[Rand_particle].nb_b1b2a1 = P[index].Pn[i_new].nb_b1b2a1;

                                P[Rand_particle].Pn[index].nb_b1b2a2 = P[i_new].Pn[index].nb_b1b2a2;
                                P[index].Pn[Rand_particle].nb_b1b2a2 = P[index].Pn[i_new].nb_b1b2a2;
                            }
                        }

                        P[i_new].x=1000.0; P[i_new].y=1000.0; P[i_new].z=1000.0;
                        for (j=0; j<P[i_new].maxnn; j++)
                        {
                            int index=P[i_new].Pn[j].nn_index;		     //index of neighbors
//                             P[i_new].Pn[index].nb_AB = 0.0;
//                             P[index].Pn[i_new].nb_AB = 0.0;
//
//                             P[i_new].Pn[index].nb_BA = 0.0;
//                             P[index].Pn[i_new].nb_BA = 0.0;

                            P[i_new].Pn[index].nb_a1a2 = 0.0;
                            P[index].Pn[i_new].nb_a1a2 = 0.0;
                            
                            P[i_new].Pn[index].nb_a2a1 = 0.0;
                            P[index].Pn[i_new].nb_a2a1 = 0.0;

                            P[i_new].Pn[index].nb_b1b2 = 0.0;
                            P[index].Pn[i_new].nb_b1b2 = 0.0;
                            
                            P[i_new].Pn[index].nb_b2b1 = 0.0;
                            P[index].Pn[i_new].nb_b2b1 = 0.0;

                            P[i_new].Pn[index].nb_a1b1 = 0.0;
                            P[index].Pn[i_new].nb_a1b1 = 0.0;

                            P[i_new].Pn[index].nb_a2b2 = 0.0;
                            P[index].Pn[i_new].nb_a2b2 = 0.0;

                            P[i_new].Pn[index].nb_a1a2b1 = 0.0;
                            P[index].Pn[i_new].nb_a1a2b1 = 0.0;

                            P[i_new].Pn[index].nb_a1a2b2 = 0.0;
                            P[index].Pn[i_new].nb_a1a2b2 = 0.0;

                            P[i_new].Pn[index].nb_b1b2a1 = 0.0;
                            P[index].Pn[i_new].nb_b1b2a1 = 0.0;

                            P[i_new].Pn[index].nb_b1b2a2 = 0.0;
                            P[index].Pn[i_new].nb_b1b2a2 = 0.0;
                        }

                        P[i_new].nbS_a1C=0.0;
                        P[i_new].nbS_a1a2C=0.0;
                        P[i_new].nbS_b2C=0.0;

                        P[i_new].maxnn=0;
                        P[i_new].maxn=0;
                        P[i_new].Pn[i_new].nl_a1a2=0.0;
                        P[i_new].Pn[i_new].nl_b1b2=0.0;
                        Np=Np_new;
                        
                        
                        for(i=0; i<Np; i++)
                        {
                        for(j=0; j<Np; j++)
                        {
                        P[i].Pn[j].nn_index=-1;
                        }
                        P[i].Sc=-1;
                        P[i].Sc_upper=-1;
                        }
                        neighborlist(Np,&ns,&ns_upper,Lx,Ly,Lz,Sz, Sz_upper, Rc,L,R2s);
                    }
                    else
                    {

                    }
                    
                    
                    cout << "Rand_particle = " <<  Rand_particle << " " << "P[Rand_particle].nbS_a1C = " << " " << P[Rand_particle].nbS_a1C << " " <<  "P[Rand_particle].nbS_a1a2C = " << " " << P[Rand_particle].nbS_a1a2C << endl;
                    cout << "i_new = " <<  i_new << " " << "P[i_new].nbS_a1C = " << " " << P[i_new].nbS_a1C << " " <<  "P[i_new].nbS_a1a2C = " << " " << P[i_new].nbS_a1a2C << endl;
                    
                    
//                    printf("index=%d\n",P[0].Pn[0].nn_index);
//                    printf("Np=%d %lf %lf %lf\n",Np,P[0].x,P[0].y,P[0].z);
//                    printf("Np=%d %lf %lf %lf\n",Np,P[1].x,P[1].y,P[1].z);
//                    printf("Np=%d %lf %lf %lf\n",Np,P[2].x,P[2].y,P[2].z);
//                    cout << " P[0].Pn[2].nb_a2b2 = " << " " << P[0].Pn[2].nb_a2b2 << " " << endl;
//                cout << " P[2].Pn[0].nb_a2b2 = " << " " << P[2].Pn[0].nb_a2b2 << " " << endl;
//                cout << " P[0].Pn[1].nb_a2b2 = " << " " << P[0].Pn[1].nb_a2b2 << " " << endl;
//                cout << " P[1].Pn[0].nb_a2b2 = " << " " << P[1].Pn[0].nb_a2b2 << " " << endl;
                    
                    
                    //                P[0].x = 19;
//         P[0].y = 30;
//         P[0].z = 6;
//         P[0].type = 0; 
//        
//         P[1].x = 30;
//         P[1].y = 41;
//         P[1].z = 28;
//         P[1].type = 1; //delete
//
//         P[2].x = 30;
//         P[2].y = 30;
//         P[2].z = 6.0;
//         P[2].type = 1;
                    
                    //exit(0);
                    
                
                }
            }
        }
        //printf("count_BD_steps=%d\n",count_BD_steps);
        //printf("hii1\n");
        if (t >0.0 && count_BD_steps % 100000==0)
        //if (t >0.0 && count_BD_steps % 100==0)
        {
            //printf("hii\n");
            myfile1 << "ITEM: TIMESTEP" <<  endl;
            myfile1 << t <<  endl;
            myfile1 << "ITEM: NUMBER OF ATOMS" <<  endl;
            myfile1 << Np <<  endl;
            myfile1 << "ITEM: BOX BOUNDS xy xz yz" <<  endl;
            myfile1 << "0 " << Lx << " 0" << endl;
            myfile1 << "0 " << Ly << " 0" << endl;
            myfile1 << "0 " << Lz << " 0" << endl;
//             myfile1 << "ITEM: ATOMS type id x y z Fx Fy Fz l na1 na2 nb1 nb2 2sbridge 3sbridge sur_2sbrg sur_3sbrg nn freesurface" <<  endl;
            myfile1 << "ITEM: ATOMS type id x y z l na1 na2 nb1 nb2 2sbridge a1a2 a2a1 a1b1 a2b2 b1b2 b2b1 3sbridge a1a2b1 a1a2b2 b1b2a1 b1b2a2 sur_2sbrg sur_3sbrg tot_a1 tot_a2 tot_b1 tot_b2 nn" <<  endl;
            double total_surfacebridge = 0.0;
            double total_surface_2strand_bridge=0.0;
            double total_surface_3strand_bridge=0.0;
            for (i=0; i<Np; i++)
            {
                total_surface_2strand_bridge=total_surface_2strand_bridge+P[i].nbS_a1C + P[i].nbS_b2C;
                total_surface_3strand_bridge=total_surface_3strand_bridge+P[i].nbS_a1a2C;
                total_surfacebridge = total_surface_2strand_bridge + total_surface_3strand_bridge;
//                 P[i].tot_2strand_bridges_surface = P[i].nbS_a1C + P[i].nbS_b2C;
                if (P[i].type==0){P[i].tot_2strand_bridges_surface = P[i].nbS_a1C;}
                if (P[i].type==1){P[i].tot_2strand_bridges_surface = P[i].nbS_b2C;}
                P[i].tot_3strand_bridges_surface = P[i].nbS_a1a2C;
//                 printf("P[i].nbS_a1C=%lf\n",P[i].nbS_a1C);
//                 printf("P[i].nbS_b2C=%lf\n",P[i].nbS_b2C);
//                 printf("P[i].nbS_a1a2C=%lf\n",P[i].nbS_a1a2C);
//                 printf("P[i].tot_2strand_bridges_surface=%lf\n",P[i].tot_2strand_bridges_surface);
//                 printf("P[i].tot_3strand_bridges_surface=%lf\n",P[i].tot_3strand_bridges_surface);
//                 printf("total_surface_2strand_bridge=%lf\n",total_surface_2strand_bridge);
//                 printf("total_surface_3strand_bridge=%lf\n",total_surface_3strand_bridge);
            }
            for(i=0; i<Np; i++)
            {
                tot_bridges = 0.0;
                double tot_2strand_bridges=0.0;
                double tot_3strand_bridges=0.0;

                double tot_a1a2_brg = 0.0;
                double tot_a2a1_brg = 0.0;
                double tot_a1b1_brg = 0.0;
                double tot_a2b2_brg = 0.0;
                double tot_b1b2_brg = 0.0;
                double tot_b2b1_brg = 0.0;
                double tot_a1a2b1_brg = 0.0;
                double tot_a1a2b2_brg = 0.0;
                double tot_b1b2a1_brg = 0.0;
                double tot_b1b2a2_brg = 0.0;

//                 P[i].Pn[i].nl = P[i].Pn[i].nl_a1a2 + P[i].Pn[i].nl_b1b2;
//                 printf("P[i].Pn[i].nl_b1b2=%lf\n",P[i].Pn[i].nl_b1b2);
                if (P[i].type == 0)
                {
                    P[i].Pn[i].nl = P[i].Pn[i].nl_a1a2;
                }
                if (P[i].type == 1)
                {
                    P[i].Pn[i].nl = P[i].Pn[i].nl_b1b2;
                }
                for (j=0; j<P[i].maxnn; j++)
                {
                    //P[0].Pn[1].nb_a1a2
                    int index=P[i].Pn[j].nn_index;
                    double aa=0.0;
                    if (index>i)
                    {
                        aa=P[i].Pn[index].nb_a1a2;
                    }
                    if (i>index)
                    {
                        aa=P[i].Pn[index].nb_a2a1;
                    }
                    
                    double bb=0.0;
                    if (index>i)
                    {
                        bb=P[i].Pn[index].nb_b1b2;
                    }
                    if (i>index)
                    {
                        bb=P[i].Pn[index].nb_b2b1;
                    }
                    //whatever the name a1a2 or a2a1, always i particle is giving the a1 and j particle is giving the a2 when P[i].Pn[j].nb_a1a2 or  P[i].Pn[j].nb_a2a1
                    tot_2strand_bridges=tot_2strand_bridges
                        + aa 
                        + bb 
                        + P[i].Pn[index].nb_a1b1 
                        + P[i].Pn[index].nb_a2b2;
                        // + P[i].Pn[index].nb_a2a1 + P[i].Pn[index].nb_b2b1;
                    tot_3strand_bridges = tot_3strand_bridges + P[i].Pn[index].nb_a1a2b1 + P[i].Pn[index].nb_a1a2b2 + P[i].Pn[index].nb_b1b2a1 + P[i].Pn[index].nb_b1b2a2;

                    tot_bridges = tot_2strand_bridges + tot_3strand_bridges;
                    double a111=0.0;
                    double a222=0.0;
                    
                    if (index>i)
                    {
                        a111 = P[i].Pn[index].nb_a1a2;
                        a222 = P[index].Pn[i].nb_a2a1;
                        //tot_a1a2_brg = tot_a1a2_brg + aa;
                    }
                    if (i>index)
                    {
                        a111=P[i].Pn[index].nb_a2a1;
                        a222=P[index].Pn[i].nb_a1a2;
                        //tot_a2a1_brg = tot_a2a1_brg + aa;
                    }
                    
                    double b111=0.0;
                    double b222=0.0;
                    
                    if (index>i)
                    {
                        b111 = P[i].Pn[index].nb_b1b2;
                        b222 = P[index].Pn[i].nb_b2b1;
                        //tot_a1a2_brg = tot_a1a2_brg + aa;
                    }
                    if (i>index)
                    {
                        b111=P[i].Pn[index].nb_b2b1;
                        b222=P[index].Pn[i].nb_b1b2;
                        //tot_a2a1_brg = tot_a2a1_brg + aa;
                    }
                    
                    
                    
//                    P[i].Pn[i].nl+P[i].nf_a1+P[0].Pn[1].nb_a1a2    a1 0particle
//                    P[i].Pn[i].nl+P[i].nf_a2+P[1].Pn[0].nb_a2a1   a2  0particle
//            
//                    P[i].Pn[i].nl+P[i].nf_a1+ P[1].Pn[0].nb_a2a1    a1 1particle
//                    P[i].Pn[i].nl+P[i].nf_a2+ P[0].Pn[1].nb_a1a2  a2  1particle
                    
                    
                    tot_a1a2_brg = tot_a1a2_brg + a111;
                    tot_a2a1_brg = tot_a2a1_brg + a222;
                    tot_a1b1_brg = tot_a1b1_brg + P[i].Pn[index].nb_a1b1;
                    tot_a2b2_brg = tot_a2b2_brg + P[i].Pn[index].nb_a2b2;
                    tot_b1b2_brg = tot_b1b2_brg + b111;
                    tot_b2b1_brg = tot_b2b1_brg + b222;

                    tot_a1a2b1_brg = tot_a1a2b1_brg + P[i].Pn[index].nb_a1a2b1;
                    tot_a1a2b2_brg = tot_a1a2b2_brg + P[i].Pn[index].nb_a1a2b2;
                    tot_b1b2a1_brg = tot_b1b2a1_brg + P[i].Pn[index].nb_b1b2a1;
                    tot_b1b2a2_brg = tot_b1b2a2_brg + P[i].Pn[index].nb_b1b2a2;

                    //printf("P[i].Pn[index].nb_a1b1=%lf\n",P[i].Pn[index].nb_a1b1);

                    //printf("P=%d nn=%d P[i].Pn[index].nb_a1a2b1=%lf\t P[i].Pn[index].nb_a1a2b2=%lf\t P[i].Pn[index].nb_b1b2a1=%lf\t P[i].Pn[index].nb_b1b2a2=%lf\n",i,index,P[i].Pn[index].nb_a1a2b1,P[i].Pn[index].nb_a1a2b2,P[i].Pn[index].nb_b1b2a1,P[i].Pn[index].nb_b1b2a2);
                    //cout << P[i].Pn[index].nb_AB << " " << P[i].Pn[index].nb_BA << endl;
                }
                
                
                //printf("%lf %lf %lf\n",tot_a1b1_brg,P[0].Pn[1].nb_a1b1,P[1].Pn[0].nb_a1b1);
                //exit(0);

                P[i].tot_a1a2_brg = tot_a1a2_brg;
                P[i].tot_a2a1_brg = tot_a2a1_brg;
                P[i].tot_a1b1_brg = tot_a1b1_brg;
                P[i].tot_a2b2_brg = tot_a2b2_brg;
                P[i].tot_b1b2_brg = tot_b1b2_brg;
                P[i].tot_b2b1_brg = tot_b2b1_brg;
                P[i].tot_a1a2b1_brg = tot_a1a2b1_brg;
                P[i].tot_a1a2b2_brg = tot_a1a2b2_brg;
                P[i].tot_b1b2a1_brg = tot_b1b2a1_brg;
                P[i].tot_b1b2a2_brg = tot_b1b2a2_brg;
                
                if (P[i].type == 0)
                {
                    P[i].tot_a1 = P[i].Pn[i].nl_a1a2 + P[i].tot_a1a2_brg + P[i].tot_a1b1_brg + P[i].tot_a1a2b1_brg + P[i].tot_a1a2b2_brg + P[i].tot_b1b2a1_brg + P[i].nf_a1 + P[i].nbS_a1C + P[i].nbS_a1a2C;
                    
                    P[i].tot_a2 = P[i].Pn[i].nl_a1a2 + P[i].tot_a2a1_brg + P[i].tot_a2b2_brg + P[i].tot_a1a2b1_brg + P[i].tot_a1a2b2_brg + P[i].tot_b1b2a2_brg + P[i].nf_a2 + P[i].nbS_a1a2C;
                    P[i].tot_b1 = 0;
                    P[i].tot_b2 = 0;
                }
                if (P[i].type == 1)
                {
                    P[i].tot_b1 = P[i].Pn[i].nl_b1b2 + P[i].tot_b1b2_brg + P[i].tot_a1b1_brg + P[i].tot_a1a2b1_brg + P[i].tot_b1b2a1_brg + P[i].tot_b1b2a2_brg + P[i].nf_b1;
                    
                    P[i].tot_b2 = P[i].Pn[i].nl_b1b2 + P[i].tot_b2b1_brg + P[i].tot_a2b2_brg + P[i].tot_a1a2b2_brg + P[i].tot_b1b2a1_brg + P[i].tot_b1b2a2_brg + P[i].nf_b2 + P[i].nbS_b2C;
                    P[i].tot_a1 = 0;
                    P[i].tot_a2 = 0;
                }
                
                myfile1 << " " << P[i].type << " " << i+1 << " " << P[i].x << " " << P[i].y << " " << P[i].z  << " " << P[i].Pn[i].nl << " " << P[i].nf_a1 << " " << P[i].nf_a2 << " " << P[i].nf_b1 << " " << P[i].nf_b2 << " " << tot_2strand_bridges << " " << P[i].tot_a1a2_brg << " " << P[i].tot_a2a1_brg << " " << P[i].tot_a1b1_brg << " " << P[i].tot_a2b2_brg << " " << P[i].tot_b1b2_brg << " " << P[i].tot_b2b1_brg << " " << tot_3strand_bridges << " " << P[i].tot_a1a2b1_brg << " " << P[i].tot_a1a2b2_brg << " " << P[i].tot_b1b2a1_brg << " " << P[i].tot_b1b2a2_brg << " " << P[i].tot_2strand_bridges_surface << " " << P[i].tot_3strand_bridges_surface << " " << P[i].tot_a1 << " " << P[i].tot_a2 << " " << P[i].tot_b1 << " " << P[i].tot_b2 << " " << P[i].maxnn << endl;


        if (P[i].type==0)
        {
            if (P[i].tot_a1 != 75 || P[i].tot_a2 != 75)
            {
                printf("total linkers of A colloid not matching at time = %lf\n",t);
                printf("P[i].tot_a1=%lf\n",P[i].tot_a1);
                printf("P[i].tot_a2=%lf\n",P[i].tot_a2);
                cout << " i= " << " " << i << endl;  
                exit(0);
            }
        }
        if (P[i].type==1)
        {
            if (P[i].tot_b1 != 75 || P[i].tot_b2 != 75)
            {
                printf("total linkers of B colloid not matching at time = %lf\n",t);
                printf("P[i].tot_b1=%lf\n",P[i].tot_b1);
                printf("P[i].tot_b2=%lf\n",P[i].tot_b2);
                exit(0);
            }
        }
    } 
//            cout << t << endl;
//            if (t>350)
//            {
//                exit(0);
//            }
        }
        
        t = t + Delta_tB;
        count_BD_steps=count_BD_steps+1;

/********************************************************************************************************************************/   
        initialisation(Np,ns,ns_upper);
        
/*********************************************************************************************************************************************/        
        
    }

//     myfile.close();
    myfile1.close();
//     myfile2.close();
    return 0;
}

void neighborlist(int Np,int *ns, int *ns_upper,double Lx,double Ly,double Lz,double Sz, double Sz_upper,double Rc,double L,double R2s)
{
    *ns=0;
    *ns_upper=0;
    int i,j,count_n,count_nn,count_regularization;
    double dx,dy,dz,d;
    //cout << "NP= " << " " << Np << endl ; 
    for(i=0; i<Np; i++)
    {
        count_regularization=0;
        count_n=0;
        for (j=i+1; j<Np; j++)
        {
            distance(i,j,Lx,Ly,Lz,&dx,&dy,&dz,&d);
            //cout << "i= " << " " << i << " " << "j= " << " " <<j << endl ; 
            
	    if (d < 2*(Rc+L) )                         //Rc=5 && L=1.0
	    {
                count_nn=P[i].maxnn;                     //counter
                P[i].Pn[count_nn].nn_index=j;		     //index of neighbors
                P[i].Pn[j].d=d; 
                P[i].maxnn=P[i].maxnn+1;			     //Maximum number of neighbors
                                                               
                    
                count_nn=P[j].maxnn;                     //counter
                P[j].Pn[count_nn].nn_index=i;		     //index of neighbors
                P[j].Pn[i].d=d; 
                P[j].maxnn=P[j].maxnn+1;			     //Maximum number of neighbors
                    
                    
/********************* To calculate neigbors only one time **********************************************************/              
                P[i].Pn[count_n].n_index=j;				 //index of neighbors without repeat
                count_n=count_n+1;                       //counter
                P[i].maxn=P[i].maxn+1;			         //Maximum number of neighbors without repeat
/*******************************************************************************************************************/
                if (d < 2*Rc+0.75*L) 
                {
                    P[i].count_regularization = P[i].count_regularization+1;
                    P[j].count_regularization = P[j].count_regularization+1;
                }
            }
            //P[j].Pn[count_nn].nn_index
                
            //printf("index=%d\n",P[0].Pn[1].nn_index);
            if (d<9.0) {cout << "Error in Calculation_PP" << endl; exit(0);}
        }

        dz = P[i].z- Sz;
        d = sqrt(SQ(dz) );

        double dz_upper=P[i].z-Sz_upper;
        double d_upper = sqrt( SQ(dz_upper) );
        if (d<(Rc+L+R2s))
        {
            P[*ns].Sc=i;
            P[i].Sd=d;
            *ns=*ns+1;

            if (d<4.0) {cout << "Error in Calculation_PS" << endl;exit(0);}
            
        }
        else if (d_upper<Rc+0.75*L)
        {
            P[*ns_upper].Sc_upper=i;
            P[i].Sd_upper=d_upper;
            *ns_upper=*ns_upper+1;
            
        }
    } 
}

void distance(int i, int j, double Lx, double Ly, double Lz,double *dx,double *dy,double *dz, double *d)
{
    *dx = P[i].x-P[j].x;
    *dx = *dx - Lx * round(*dx/Lx);
            
    *dy = P[i].y-P[j].y;
    *dy = *dy - Ly * round(*dy/Ly);
            
    *dz = P[i].z-P[j].z;
            
    *d=sqrt( SQ(*dx) + SQ(*dy) + SQ(*dz));

}
void distance_ads(double X, double Y, double Z,int i, double Lx, double Ly, double Lz,double *dx,double *dy,double *dz, double *d)
{
    *dx = P[i].x-X;
    *dx = *dx - Lx * round(*dx/Lx);
            
    *dy = P[i].y-Y;
    *dy = *dy - Ly * round(*dy/Ly);
            
    *dz = P[i].z-Z;
            
    *d=sqrt( SQ(*dx) + SQ(*dy) + SQ(*dz));

}
void overlapvolume(int Np,double Rc,double L,double omega_0,double omega_0_inert)
{
    int i,j;
    double R1  =Rc+L;
    double R2  =Rc;
    double D1  =R1+R1;
    double D1_ =R1-R1;
    double D2  =R2+R2;
    double D2_ =R2-R2;
    double D12 =R1+R2;
    double D12_=R1-R2;
    double ov, ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
    double e_ij,w_ij,e_ij_prime,w_ij_prime;
    double regu_e_ij,regu_e_ij_prime;
    double e_ij_sum;
    double regu_e_ij_sum;

    for (i=0; i<Np; i++)
    {
        e_ij_sum=0.0;
        regu_e_ij_sum=0.0;
        for (j=0; j<P[i].maxnn; j++)
        {
            int index = P[i].Pn[j].nn_index;
            double r=P[i].Pn[index].d;
//             printf("r=%lf",r);
            if (r>=(2*(Rc+L))) 
            {
                e_ij=0.0;                                       // corona-sphere overlap
                w_ij=0.0;                                        //corona-corona overlap
                e_ij_prime=0.0;
                w_ij_prime=0.0;
                regu_e_ij=0.0;
                regu_e_ij_prime=0.0;
            }
            else if (r<(2*(Rc+L)) && r>=2*Rc+1.9*L)
            {
                e_ij=0.0;
                e_ij_prime=0.0;
                overlap(2*Rc+1.9*L,D1,D1_,&ov);
                overlap_prime(2*Rc+1.9*L,D1,D1_,&ov_prime);
                w_ij=ov;
                w_ij_prime=ov_prime;
                regu_e_ij=0.0;
                regu_e_ij_prime=0.0;
            }
            else if (r<(2*Rc+1.9*L) && r>=2*Rc+L)
            {
                e_ij=0.0;
                e_ij_prime=0.0;
                    
                overlap(r,D1,D1_,&ov);
                overlap_prime(r,D1,D1_,&ov_prime);
                w_ij=ov;
                w_ij_prime=ov_prime;
                regu_e_ij=0.0;
                regu_e_ij_prime=0.0;
            }
            else if (r<2*Rc+L && r>=2*Rc)
            {
                overlap(r,D12,D12_,&ov);
                overlap_prime(r,D12,D12_,&ov_prime);
                    
                e_ij=ov;
                e_ij_prime=ov_prime;
                    
                overlap(r,D1,D1_,&ov);
                overlap_prime(r,D1,D1_,&ov_prime);
                w_ij=ov-2.0*e_ij;
                w_ij_prime=ov_prime-2.0*e_ij_prime;
                regu_e_ij=0.0;
                regu_e_ij_prime=0.0;
            }
                
            else if (r<2*Rc)
            {
                overlap(r,D12,D12_,&ov1);
                overlap(r,D2,D2_,&ov2);
                overlap_prime(r,D12,D12_,&ov1_prime);
                overlap_prime(r,D2,D2_,&ov2_prime);
                e_ij=ov1-ov2;                                   //overlap(Rc+L,Rc)-overlap(Rc,Rc);
                e_ij_prime=ov1_prime-ov2_prime;
                    
                overlap(r,D1,D1_,&ov1);
                overlap(r,D12,D12_,&ov2);
                overlap(r,D2,D2_,&ov3);
                overlap_prime(r,D1,D1_,&ov1_prime);
                overlap_prime(r,D12,D12_,&ov2_prime);
                overlap_prime(r,D2,D2_,&ov3_prime);
                w_ij=ov1 - 2*ov2 + ov3;                         //overlap(Rc+L,Rc+L)-2*overlap(Rc+L,Rc)+overlap(Rc,Rc);
                w_ij_prime=ov1_prime - 2*ov2_prime + ov3_prime;
                regu_e_ij=0.0;
                regu_e_ij_prime=0.0;
            }
            if (r<2*Rc+0.75*L) {
                overlap(r,Rc+Rc+0.75*L,0.75*L,&ov1);
                overlap_prime(r,Rc+Rc+0.75*L,0.75*L,&ov1_prime);
                
                regu_e_ij=ov1;
                regu_e_ij_prime=ov1_prime;
            }
            e_ij_sum=e_ij_sum+e_ij;
            regu_e_ij_sum=regu_e_ij_sum+regu_e_ij;
            
            P[i].Pn[index].e_ij = e_ij;
            P[i].Pn[index].e_ij_prime = e_ij_prime;
            P[i].Pn[index].omega_ij = w_ij;
            P[i].Pn[index].omega_ij_prime = w_ij_prime;  
            P[i].Pn[index].regu_e_ij=regu_e_ij;
            P[i].Pn[index].regu_e_ij_prime=regu_e_ij_prime;
        }
        P[i].omega_i=omega_0-e_ij_sum;
        P[i].omega_i_inert=omega_0_inert-regu_e_ij_sum;
//         printf(" w_ij=%lf e_ij_sum=%lf\n",w_ij,e_ij_sum);
    }    
}
void overlapvolume_surface(int ns,double Rc,double L,double R2s,double omega_S_0, double *omega_S)
{
    int i;
    double e_is, e_si, w_is, e_is_prime, e_si_prime, w_is_prime;
    double ov, ov1, ov2, ov3,ov4, ov_prime, ov1_prime, ov2_prime, ov3_prime,ov4_prime;
    double R1,R2,Dcs,dcs; 
    double e_si_sum=0.0;
    for (i=0;i<ns; i++)
    {
        int index_cs=P[i].Sc;
        double r=P[index_cs].Sd;
    
        if (r>=Rc+L+R2s) {
            e_is=0.0;                                       //corona-sphere overlap
            w_is=0.0;                                          //corona-corona overlap
            e_si=0.0;
            e_is_prime=0.0;
            w_is_prime=0.0;
            e_si_prime=0.0;
        }
        else if ( r<(Rc+L+R2s) && r >= (Rc+L+0.9*R2s) ){
             R1=Rc+L; R2=R2s; Dcs=R1+R2; dcs=R1-R2;
             r=Rc+L+0.9*R2s;
             e_is=0.0;
             e_is_prime=0.0;
             overlap_surface(r,Dcs,dcs,R2,&ov);
             overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
             w_is=ov;
             w_is_prime=ov_prime;
             e_si=0.0;
             e_si_prime=0.0;   
        }
        else if ( r < (Rc+L+0.9*R2s) && r >= (Rc+L) ){
            R1=Rc+L; R2=R2s; Dcs=R1+R2; dcs=R1-R2;
            e_is=0.0;
            e_is_prime=0.0;
            overlap_surface(r,Dcs,dcs,R2,&ov);
            overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
            w_is=ov;
            w_is_prime=ov_prime;
            e_si=0.0;
            e_si_prime=0.0;
        }
    
        else if (r<Rc+L && r>=Rc){   //previusly I used r<Rc+R2s && r>=Rc
            R1=Rc;R2=R2s;Dcs=R1+R2;dcs=R1-R2;
            overlap_surface(r,Dcs,dcs,R2,&ov);
            overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
            e_si=ov;       //surface linkes are blocked by the sphere
            e_si_prime=ov_prime;   

            R1=Rc+L;R2=0.0;Dcs=R1+R2;dcs=R1-R2;
            overlap_surface(r,Dcs,dcs,R2,&ov);
            overlap_surface_prime(r,Dcs,dcs,R2,&ov_prime);
            e_is=ov;       // colloids linkers are blocked by the surface
            e_is_prime=ov_prime;
                
            R1=Rc+L;R2=R2s;Dcs=R1+R2;dcs=R1-R2; 
            overlap_surface(r,Dcs,dcs,R2s,&ov);
            overlap_surface_prime(r,Dcs,dcs,R2s,&ov_prime);
            w_is=ov-e_is-e_si;
            w_is_prime=ov_prime-e_is_prime-e_si_prime;
        }
    
        else if (r<Rc){
            R1=Rc+L;R2=R2s;Dcs=R1+R2;dcs=R1-R2;
            overlap_surface(r,Dcs,dcs,R2,&ov1); //overlap between (Rc+L) and corona of surface
            overlap_surface_prime(r,Dcs,dcs,R2,&ov1_prime);

            R1=Rc;R2=R2s;Dcs=R1+R2;dcs=R1-R2;
            overlap_surface(r,Dcs,dcs,R2,&ov2); //overlap between (Rc) and corona of surface
            overlap_surface_prime(r,Dcs,dcs,R2,&ov2_prime);

            R1=Rc+L;R2=0.0;Dcs=R1+R2;dcs=R1-R2;
            overlap_surface(r,Dcs,dcs,R2,&ov3); //overlap between (Rc+L) and hard core of the surface
            overlap_surface_prime(r,Dcs,dcs,R2,&ov3_prime);
            
            R1=Rc;R2=0.0;Dcs=R1+R2;dcs=R1-R2;
            overlap_surface(r,Dcs,dcs,R2,&ov4); //overlap between (Rc) and hard core of the surface
            overlap_surface_prime(r,Dcs,dcs,R2,&ov4_prime);
            
            w_is=ov1-ov2-ov3+ov4;
            w_is_prime=ov1_prime-ov2_prime-ov3_prime+ov4_prime;

            e_is=ov3-ov4;       // colloids linkers are blocked by the surface
            e_is_prime=ov3_prime-ov4_prime;
       
       
            e_si=ov2-ov4;       //surface linkes are blocked by the sphere
            e_si_prime=ov2_prime-ov4_prime;
        }   
    
        P[index_cs].e_is = e_is;
        P[index_cs].e_is_prime = e_is_prime;
        P[index_cs].e_si = e_si;
        P[index_cs].e_si_prime = e_si_prime;
        P[index_cs].omega_is = w_is;
        P[index_cs].omega_is_prime = w_is_prime;
        P[index_cs].omega_i = P[index_cs].omega_i-e_is;
        e_si_sum=e_si_sum+e_si;
//         printf("%d P[index_cs].e_is %lf P[index_cs].e_is_prime %lf P[index_cs].e_si %lf P[index_cs].e_si_prime %lf P[index_cs].omega_is %lf  P[index_cs].omega_is_prime %lf P[index_cs].omega_i %lf\n",index_cs,P[index_cs].e_is,P[index_cs].e_is_prime,P[index_cs].e_si,P[index_cs].e_si_prime,P[index_cs].omega_is,P[index_cs].omega_is_prime,P[index_cs].omega_i);
    }
    
    *omega_S=omega_S_0-e_si_sum;

}

void overlap(double r,double D,double D_,double *ov)
{
    *ov=M_PI/(12.0*r) * SQ(D-r) * (SQ(r) + 2*r*D - 3*SQ(D_));
}
void overlap_prime(double r,double D,double D_,double *ov_prime)
{
    *ov_prime  = (M_PI/(4.0*SQ(r))) * (D-r) * (D*SQ(D_) - D*SQ(r) + SQ(D_)*r - CUB(r));
}
void overlap_surface(double r,double D,double D_,double R2s, double *ov)
{
    *ov=M_PI/(3.0) * SQ(D-r) * (R2s + r + 2*D_);
}
void overlap_surface_prime(double r,double D,double D_,double R2s,double *ov_prime)
{
    *ov_prime=M_PI/(3.0) * (D-r) * (D-2*R2s-3*r-4*D_);
}
void normal(double *Rn)
{
    double x1, x2, w, y1, y2;
    do {

        x1 = 2.0 * rand()/(double)RAND_MAX - 1.0;
        x2 = 2.0 * rand()/(double)RAND_MAX - 1.0;
        w = x1 * x1 + x2 * x2;

    } while ( w >= 1.0 );
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w; 
    *Rn=y1;
}


void Gillespie(int Np,double rho0,double Delta_tB, double kon0, double beta, double DG0, double DG0_L, double DG0_S,int ns,double *nfS_C,double omega_S, double tau0, double Delta_G0, double Delta_GT, double K0, double K1, double n_alpha, double n_T)
{

    double a_tot;
    double tao;
    double t_bar=0.0;
    int count_t=0;
    bool flag;
    int i,j,q;

//     double A_tot_Part_Part = 0.0;
//     double A_tot_Part_Sur = 0.0;

    while (t_bar < Delta_tB)
    {
        //////////////Affinity calculations (I will create a subroutine. In the present version, I keep the calculations in the gillespie subroutine.)/////////////

        double A_tot_Part_Part = 0.0;
        double A_tot_Part_Sur = 0.0;

        double a_on_sum,a_off_sum;
        double a_on_AB_sum, a_on_BA_sum, a_off_AB_sum, a_off_BA_sum, a_on_sum_L, a_off_sum_L;
        double a_on_sum_L_a1a2,a_on_sum_L_b1b2,a_off_sum_L_a1a2,a_off_sum_L_b1b2,a_on_a1a2_sum,a_off_a1a2_sum,a_on_a1b1_sum,a_off_a1b1_sum,a_on_b1b2_sum,a_off_b1b2_sum,a_on_a2b2_sum,a_off_a2b2_sum,a_on_a2a1_sum,a_off_a2a1_sum,a_on_b2b1_sum,a_off_b2b1_sum;
        double a_on_a1a2_b1_sum,a_off_a1a2_b1_sum,a_on_a1_b1_a2_sum,a_off_a1_b1_a2_sum,a_on_a1a2_b2_sum,a_off_a1a2_b2_sum,a_on_a2_b2_a1_sum,a_off_a2_b2_a1_sum;
        double a_on_b1b2_a1_sum,a_off_b1b2_a1_sum,a_on_b1_a1_b2_sum,a_off_b1_a1_b2_sum,a_on_b1b2_a2_sum,a_off_b1b2_a2_sum,a_on_b2_a2_b1_sum,a_off_b2_a2_b1_sum;
        double a_on_AC_sum_S, a_off_AC_sum_S;
        double a_on_a1C_sum_S,a_off_a1C_sum_S,a_on_b2C_sum_S,a_off_b2C_sum_S;
        double a_on_a1a2C_sum_S,a_off_a1a2C_sum_S,a_on_a1_C_a2_sum_S,a_off_a1_C_a2_sum_S;
        double kon,koff,kon_L,koff_L,kon_S,koff_S,koff_3s,koff_S_3s;

        ///////////////initialisation of the Affinity terms///////////////////////

        for (i=0; i<Np; i++)
        {
            a_on_a1a2_sum = 0.0;
            a_off_a1a2_sum = 0.0;
            a_on_a1b1_sum = 0.0;
            a_off_a1b1_sum = 0.0;
            a_on_a2b2_sum = 0.0;
            a_off_a2b2_sum = 0.0;
            a_on_b1b2_sum = 0.0;
            a_off_b1b2_sum = 0.0;

            a_on_a2a1_sum = 0.0;
            a_off_a2a1_sum = 0.0;
            a_on_b2b1_sum = 0.0;
            a_off_b2b1_sum = 0.0;

            a_on_a1a2_b1_sum = 0.0;
            a_off_a1a2_b1_sum = 0.0;
            a_on_a1_b1_a2_sum = 0.0;
            a_off_a1_b1_a2_sum = 0.0;
            a_on_a1a2_b2_sum = 0.0;
            a_off_a1a2_b2_sum = 0.0;
            a_on_a2_b2_a1_sum = 0.0;
            a_off_a2_b2_a1_sum = 0.0;

            a_on_b1b2_a1_sum = 0.0;
            a_off_b1b2_a1_sum = 0.0;
            a_on_b1_a1_b2_sum = 0.0;
            a_off_b1_a1_b2_sum = 0.0;
            a_on_b1b2_a2_sum = 0.0;
            a_off_b1b2_a2_sum = 0.0;
            a_on_b2_a2_b1_sum = 0.0;
            a_off_b2_a2_b1_sum = 0.0;

            a_on_sum_L_a1a2 = 0.0;
            a_on_sum_L_b1b2 = 0.0;
            a_off_sum_L_a1a2 = 0.0;
            a_off_sum_L_b1b2 = 0.0;

            kon_L  = tau0 * (1/P[i].omega_i * rho0);
            koff_L = tau0/K0;

//             printf("kon_L=%lf\tkoff_L=%lf\n",kon_L,koff_L);
// printf("type=%d\n",P[i].type);
            for (j=0; j<P[i].maxnn; j++)
            {
//                 printf("maxn=%d\n",P[i].maxn);
// printf("type=%d\n",P[i].type);
                int index = P[i].Pn[j].nn_index;

                kon = tau0 * ( P[i].Pn[index].omega_ij/(P[i].omega_i * P[index].omega_i * rho0) );
//                 printf("P[%d].Pn[%d].omega_ij=%lf\n",i,index,P[i].Pn[index].omega_ij);
                koff = tau0/K0;
                koff_3s = tau0/(K1*(n_alpha+1));

                ////////////***********2strand complex***************///////////////////////////
//                 printf("type=%d Na1=%lf\n",P[i].type,P[i].nf_a1);
                if(P[i].type==0 && P[index].type==0 && i<index)
                {
                    //printf("i,index=%d %d\n",i,index);
                    P[i].Pn[index].aon_a1a2  = P[i].nf_a1 * P[index].nf_a2 * kon;
                    a_on_a1a2_sum = a_on_a1a2_sum + P[i].Pn[index].aon_a1a2;

                    P[i].Pn[index].aoff_a1a2  = P[i].Pn[index].nb_a1a2 * koff;
                    a_off_a1a2_sum = a_off_a1a2_sum + P[i].Pn[index].aoff_a1a2;
                    
                    
//                    P[i].Pn[i].nl+P[i].nf_a1+P[0].Pn[1].nb_a1a2    a1 0particle
//                    P[i].Pn[i].nl+P[i].nf_a2+P[1].Pn[0].nb_a2a1   a2  0particle
//            
//                    P[i].Pn[i].nl+P[i].nf_a1+ P[1].Pn[0].nb_a2a1    a1 1particle
//                    P[i].Pn[i].nl+P[i].nf_a2+ P[0].Pn[1].nb_a1a2  a2  1particle
                    
                    
                      //whatever the name a1a2 or a2a1, always i particle is giving the a1 and j particle is giving the a2 when P[i].Pn[j].nb_a1a2 or  P[i].Pn[j].nb_a2a1
                    
//                    P[i].Pn[index].aon_a2a1  = P[i].nf_a2 * P[index].nf_a1 * kon;
//                    a_on_a2a1_sum = a_on_a2a1_sum + P[i].Pn[index].aon_a2a1;
//
//                    P[i].Pn[index].aoff_a2a1  = P[i].Pn[index].nb_a2a1 * koff;
//                    a_off_a2a1_sum = a_off_a2a1_sum + P[i].Pn[index].aoff_a2a1;
                }
                
                
                
                if(P[i].type==0 && P[index].type==0 && i>index)
                {
                    //printf("i,index=%d %d\n",i,index);
//                    P[i].Pn[index].aon_a1a2  = P[i].nf_a1 * P[index].nf_a2 * kon;
//                    a_on_a1a2_sum = a_on_a1a2_sum + P[i].Pn[index].aon_a1a2;
//
//                    P[i].Pn[index].aoff_a1a2  = P[i].Pn[index].nb_a1a2 * koff;
//                    a_off_a1a2_sum = a_off_a1a2_sum + P[i].Pn[index].aoff_a1a2;
                    
                    P[i].Pn[index].aon_a2a1  = P[i].nf_a1 * P[index].nf_a2 * kon;
                    a_on_a2a1_sum = a_on_a2a1_sum + P[i].Pn[index].aon_a2a1;

                    P[i].Pn[index].aoff_a2a1  = P[i].Pn[index].nb_a2a1 * koff;
                    a_off_a2a1_sum = a_off_a2a1_sum + P[i].Pn[index].aoff_a2a1;
                }
                
                
                if(P[i].type==0 && P[index].type==1)
                {
                    P[i].Pn[index].aon_a1b1  = P[i].nf_a1 * P[index].nf_b1 * kon;
                    a_on_a1b1_sum = a_on_a1b1_sum + P[i].Pn[index].aon_a1b1;

                    P[i].Pn[index].aoff_a1b1  = P[i].Pn[index].nb_a1b1 * koff;
                    a_off_a1b1_sum = a_off_a1b1_sum + P[i].Pn[index].aoff_a1b1;

                    P[i].Pn[index].aon_a2b2  = P[i].nf_a2 * P[index].nf_b2 * kon;
                    a_on_a2b2_sum = a_on_a2b2_sum + P[i].Pn[index].aon_a2b2;

                    P[i].Pn[index].aoff_a2b2  = P[i].Pn[index].nb_a2b2 * koff;
                    a_off_a2b2_sum = a_off_a2b2_sum + P[i].Pn[index].aoff_a2b2;

                ////////////***********3strand complex***************///////////////////////////
                ////////(here bridge is shown using dash e.g. a1a2_b1 = a1a2 are in loop formation first then make 3strand bridge with b1)////////
// P[0].Pn[1].nl_a1a2 = 5;
                    P[i].Pn[index].aon_a1a2_b1 = P[i].Pn[i].nl_a1a2* P[index].nf_b1 * kon;
                    a_on_a1a2_b1_sum = a_on_a1a2_b1_sum + P[i].Pn[index].aon_a1a2_b1;
// printf("P[0].Pn[1].aon_a1a2_b1=%lf\n",P[0].Pn[1].aon_a1a2_b1);
                    P[i].Pn[index].aoff_a1a2_b1 = P[i].Pn[index].nb_a1a2b1 * koff_3s;
                    a_off_a1a2_b1_sum = a_off_a1a2_b1_sum + P[i].Pn[index].aoff_a1a2_b1;

                    P[i].Pn[index].aon_a1_b1_a2 = P[i].Pn[index].nb_a1b1 * P[i].nf_a2 * kon_L;
                    a_on_a1_b1_a2_sum = a_on_a1_b1_a2_sum + P[i].Pn[index].aon_a1_b1_a2;

                    P[i].Pn[index].aoff_a1_b1_a2 = P[i].Pn[index].nb_a1a2b1 * koff_3s;
                    a_off_a1_b1_a2_sum = a_off_a1_b1_a2_sum + P[i].Pn[index].aoff_a1_b1_a2;

                    P[i].Pn[index].aon_a1a2_b2 = P[i].Pn[i].nl_a1a2* P[index].nf_b2 * kon;
                    a_on_a1a2_b2_sum = a_on_a1a2_b2_sum + P[i].Pn[index].aon_a1a2_b2;

                    P[i].Pn[index].aoff_a1a2_b2 = P[i].Pn[index].nb_a1a2b2 * koff_3s;
                    a_off_a1a2_b2_sum = a_off_a1a2_b2_sum + P[i].Pn[index].aoff_a1a2_b2;

                    P[i].Pn[index].aon_a2_b2_a1 = P[i].Pn[index].nb_a2b2 * P[i].nf_a1 * kon_L;
                    a_on_a2_b2_a1_sum = a_on_a2_b2_a1_sum + P[i].Pn[index].aon_a2_b2_a1;

                    P[i].Pn[index].aoff_a2_b2_a1 = P[i].Pn[index].nb_a1a2b2 * koff_3s;
                    a_off_a2_b2_a1_sum = a_off_a2_b2_a1_sum + P[i].Pn[index].aoff_a2_b2_a1;
                }

                if(P[i].type==1 && P[index].type==1 && i<index)
                {
                    P[i].Pn[index].aon_b1b2  = P[i].nf_b1 * P[index].nf_b2 * kon;
                    a_on_b1b2_sum = a_on_b1b2_sum + P[i].Pn[index].aon_b1b2;

// printf("kon=%lf\n",kon);

                    P[i].Pn[index].aoff_b1b2  = P[i].Pn[index].nb_b1b2 * koff;
                    a_off_b1b2_sum = a_off_b1b2_sum + P[i].Pn[index].aoff_b1b2;

                    //P[i].Pn[index].aon_b2b1  = P[i].nf_b2 * P[index].nf_b1 * kon;
                    //a_on_b2b1_sum = a_on_b2b1_sum + P[i].Pn[index].aon_b2b1;

                    //P[i].Pn[index].aoff_b2b1  = P[i].Pn[index].nb_b2b1 * koff;
                    //a_off_b2b1_sum = a_off_b2b1_sum + P[i].Pn[index].aoff_b2b1;

                }
                
                if(P[i].type==1 && P[index].type==1 && i>index)
                {
//                    P[i].Pn[index].aon_b1b2  = P[i].nf_b1 * P[index].nf_b2 * kon;
//                    a_on_b1b2_sum = a_on_b1b2_sum + P[i].Pn[index].aon_b1b2;
//
//                    // printf("kon=%lf\n",kon);
//
//                    P[i].Pn[index].aoff_b1b2  = P[i].Pn[index].nb_b1b2 * koff;
//                    a_off_b1b2_sum = a_off_b1b2_sum + P[i].Pn[index].aoff_b1b2;

                    P[i].Pn[index].aon_b2b1  = P[i].nf_b1 * P[index].nf_b2 * kon;
                    a_on_b2b1_sum = a_on_b2b1_sum + P[i].Pn[index].aon_b2b1;

                     P[i].Pn[index].aoff_b2b1  = P[i].Pn[index].nb_b2b1 * koff;
                     a_off_b2b1_sum = a_off_b2b1_sum + P[i].Pn[index].aoff_b2b1;        
                }
                
                
                if(P[i].type==1 && P[index].type==0)
                {
                    ////////////*******3strand complex***************///////////////////////////

                    P[i].Pn[index].aon_b1b2_a1 = P[i].Pn[i].nl_b1b2* P[index].nf_a1 * kon;
                    a_on_b1b2_a1_sum = a_on_b1b2_a1_sum + P[i].Pn[index].aon_b1b2_a1;
//                     printf("a_on_b1b2_a1_sum=%lf\n",a_on_b1b2_a1_sum);

                    P[i].Pn[index].aoff_b1b2_a1 = P[i].Pn[index].nb_b1b2a1 * koff_3s;
                    a_off_b1b2_a1_sum = a_off_b1b2_a1_sum + P[i].Pn[index].aoff_b1b2_a1;

                    P[i].Pn[index].aon_b1_a1_b2 = P[i].Pn[index].nb_a1b1 * P[i].nf_b2 * kon_L;
                    a_on_b1_a1_b2_sum = a_on_b1_a1_b2_sum + P[i].Pn[index].aon_b1_a1_b2;

                    P[i].Pn[index].aoff_b1_a1_b2 = P[i].Pn[index].nb_b1b2a1 * koff_3s;
                    a_off_b1_a1_b2_sum = a_off_b1_a1_b2_sum + P[i].Pn[index].aoff_b1_a1_b2;
                   ////////////
                    P[i].Pn[index].aon_b1b2_a2 = P[i].Pn[i].nl_b1b2* P[index].nf_a2 * kon;
                    a_on_b1b2_a2_sum = a_on_b1b2_a2_sum + P[i].Pn[index].aon_b1b2_a2;

                    P[i].Pn[index].aoff_b1b2_a2 = P[i].Pn[index].nb_b1b2a2 * koff_3s;
                    a_off_b1b2_a2_sum = a_off_b1b2_a2_sum + P[i].Pn[index].aoff_b1b2_a2;

                    P[i].Pn[index].aon_b2_a2_b1 = P[i].Pn[index].nb_a2b2 * P[i].nf_b1 * kon_L;
                    a_on_b2_a2_b1_sum = a_on_b2_a2_b1_sum + P[i].Pn[index].aon_b2_a2_b1;

                    P[i].Pn[index].aoff_b2_a2_b1 = P[i].Pn[index].nb_b1b2a2 * koff_3s;
                    a_off_b2_a2_b1_sum = a_off_b2_a2_b1_sum + P[i].Pn[index].aoff_b2_a2_b1;
                }
            }
            
            
            P[i].A_on_a1a2_sum = a_on_a1a2_sum;
            P[i].A_off_a1a2_sum = a_off_a1a2_sum;
            P[i].A_on_a1b1_sum = a_on_a1b1_sum;
            P[i].A_off_a1b1_sum = a_off_a1b1_sum;
            P[i].A_on_a2b2_sum = a_on_a2b2_sum;
            P[i].A_off_a2b2_sum = a_off_a2b2_sum;
            P[i].A_on_b1b2_sum = a_on_b1b2_sum;
            P[i].A_off_b1b2_sum = a_off_b1b2_sum;

            P[i].A_on_a2a1_sum = a_on_a2a1_sum;
            P[i].A_off_a2a1_sum = a_off_a2a1_sum;
            P[i].A_on_b2b1_sum = a_on_b2b1_sum;
            P[i].A_off_b2b1_sum = a_off_b2b1_sum;

            P[i].A_tot_PP_2strand = a_on_a1a2_sum + a_off_a1a2_sum + a_on_a1b1_sum + a_off_a1b1_sum + a_on_a2b2_sum + a_off_a2b2_sum + a_on_b1b2_sum + a_off_b1b2_sum + a_on_a2a1_sum + a_off_a2a1_sum + a_on_b2b1_sum + a_off_b2b1_sum; ///this type of terms will save the calculated affinities for each particle, here it is calculating all the particle-particle 2strand complex affinities of each particle over all of its neighbors
            // printf("a_on_b1b2_sum=%lf\n",a_on_b1b2_sum );
            P[i].A_on_a1a2_b1_sum = a_on_a1a2_b1_sum;
            P[i].A_off_a1a2_b1_sum = a_off_a1a2_b1_sum;
            P[i].A_on_a1_b1_a2_sum = a_on_a1_b1_a2_sum;
            P[i].A_off_a1_b1_a2_sum = a_off_a1_b1_a2_sum;
            P[i].A_on_a1a2_b2_sum = a_on_a1a2_b2_sum;
            P[i].A_off_a1a2_b2_sum = a_off_a1a2_b2_sum;
            P[i].A_on_a2_b2_a1_sum = a_on_a2_b2_a1_sum;
            P[i].A_off_a2_b2_a1_sum = a_off_a2_b2_a1_sum;
            P[i].A_on_b1b2_a1_sum = a_on_b1b2_a1_sum;
            P[i].A_off_b1b2_a1_sum = a_off_b1b2_a1_sum;
            P[i].A_on_b1_a1_b2_sum = a_on_b1_a1_b2_sum;
            P[i].A_off_b1_a1_b2_sum = a_off_b1_a1_b2_sum;
            P[i].A_on_b1b2_a2_sum = a_on_b1b2_a2_sum;
            P[i].A_off_b1b2_a2_sum = a_off_b1b2_a2_sum;
            P[i].A_on_b2_a2_b1_sum = a_on_b2_a2_b1_sum;
            P[i].A_off_b2_a2_b1_sum = a_off_b2_a2_b1_sum;

            P[i].A_tot_PP_3strand = a_on_a1a2_b1_sum + a_off_a1a2_b1_sum + a_on_a1_b1_a2_sum + a_off_a1_b1_a2_sum + a_on_a1a2_b2_sum + a_off_a1a2_b2_sum + a_on_a2_b2_a1_sum + a_off_a2_b2_a1_sum + a_on_b1b2_a1_sum + a_off_b1b2_a1_sum + a_on_b1_a1_b2_sum + a_off_b1_a1_b2_sum + a_on_b1b2_a2_sum + a_off_b1b2_a2_sum + a_on_b2_a2_b1_sum + a_off_b2_a2_b1_sum; ///calculating all the particle-particle 3strand complex affinities of each particle over all of its neighbors

            //printf("%d,a_on_a1a2_b1_sum=%lf\n",P[i].type,a_on_a1a2_b1_sum);

            if(P[i].type==0)
            {
                P[i].Pn[i].aon_L_a1a2  = P[i].nf_a1 * P[i].nf_a2 * kon_L;
                a_on_sum_L_a1a2 = a_on_sum_L_a1a2 + P[i].Pn[i].aon_L_a1a2;

                P[i].Pn[i].aoff_L_a1a2 = P[i].Pn[i].nl_a1a2 * koff_L;
                a_off_sum_L_a1a2 = a_off_sum_L_a1a2 + P[i].Pn[i].aoff_L_a1a2;
            }
            
            
            //cout << "i = " <<  i << " " << "P[i].Pn[i].aon_L_a1a2 = " << " " << P[i].Pn[i].aon_L_a1a2 << " " <<  "P[i_new].nbS_a1a2C = " << " " << P[i_new].nbS_a1a2C << endl;

            if(P[i].type==1)
            {
                P[i].Pn[i].aon_L_b1b2  = P[i].nf_b1 * P[i].nf_b2 * kon_L;
                a_on_sum_L_b1b2 = a_on_sum_L_b1b2 + P[i].Pn[i].aon_L_b1b2;

                P[i].Pn[i].aoff_L_b1b2 = P[i].Pn[i].nl_b1b2 * koff_L;
                a_off_sum_L_b1b2 = a_off_sum_L_b1b2 + P[i].Pn[i].aoff_L_b1b2;
//                 printf("P[i].Pn[i].nl_b1b2=%lf\n",P[i].Pn[i].nl_b1b2);
            }

            P[i].A_on_sum_L_a1a2 = a_on_sum_L_a1a2;
            P[i].A_off_sum_L_a1a2 = a_off_sum_L_a1a2;
            P[i].A_on_sum_L_b1b2 = a_on_sum_L_b1b2;
            P[i].A_off_sum_L_b1b2 = a_off_sum_L_b1b2;

            P[i].A_tot_loop = a_on_sum_L_a1a2 + a_off_sum_L_a1a2 + a_on_sum_L_b1b2 + a_off_sum_L_b1b2;////Affinity of making loops for each particle
//              printf("P[i].Pn[i].nl_b1b2=%lf\n",P[i].Pn[i].nl_b1b2);
//             printf("a_off_sum_L_b1b2=%lf\n",a_off_sum_L_b1b2);
            P[i].A_tot_PP = P[i].A_tot_PP_2strand + P[i].A_tot_PP_3strand + P[i].A_tot_loop;////total affinity of each particle of making loops and bridges with all of its neighbours
//          printf("P[i].A_tot_PP_2strand=%lf\n",P[i].A_tot_PP_2strand);
//             A_tot_Part_Part = A_tot_Part_Part + P[i].A_tot_PP;////affinity sum over all the particles
            A_tot_Part_Part = A_tot_Part_Part + P[i].A_tot_PP_2strand + P[i].A_tot_PP_3strand + P[i].A_tot_loop;
//             printf("Ptype=%d\tA_tot_PP_3strand[i]=%lf\n",P[i].type,A_tot_PP_3strand[i]);

        }
//         exit(0);

//         double A_tot_Part_Sur = 0.0;

        for (q=0; q<ns; q++)
        {
//             printf("hiihii\n");
//             printf("ns=%d\n",ns);
            a_on_a1C_sum_S = 0.0;
            a_off_a1C_sum_S = 0.0;
            a_on_b2C_sum_S = 0.0;
            a_off_b2C_sum_S = 0.0;

            a_on_a1a2C_sum_S = 0.0;
            a_off_a1a2C_sum_S = 0.0;
            a_on_a1_C_a2_sum_S = 0.0;
            a_off_a1_C_a2_sum_S = 0.0;

            int index_cs=P[q].Sc;
            //printf("index_cs=%d\n",index_cs);
            //if (index_cs==1) exit(0);
            kon_S = tau0 * ( P[index_cs].omega_is/(P[index_cs].omega_i * omega_S * rho0) );
            koff_S = tau0/K0;
            kon_L  = tau0 * (1/P[index_cs].omega_i * rho0);
            //koff_S_3s = tau0/(K1*(n_alpha+1));
            koff_S_3s = tau0/(K1*(n_alpha+1));
            //printf("koff_S_3s=%lf\n",koff_S_3s);

            if(P[index_cs].type==0)
            {
                //////////////********2 strand complexes with surface*********//////////////////////

                P[index_cs].aon_a1C_S = P[index_cs].nf_a1 * *nfS_C * kon_S;
                a_on_a1C_sum_S = a_on_a1C_sum_S + P[index_cs].aon_a1C_S;

                P[index_cs].aoff_a1C_S=P[index_cs].nbS_a1C * koff_S;
                a_off_a1C_sum_S = a_off_a1C_sum_S + P[index_cs].aoff_a1C_S;
                // if(P[index_cs].nbS_a1C<0){printf("a_off_a1C_sum_S=%lf  P[%d].aoff_a1C_S=%lf\n",a_off_a1C_sum_S,index_cs,P[index_cs].aoff_a1C_S);}
        //////////////////***************3 strand complexes with surface**************//////////////////////

                P[index_cs].aon_a1a2_C_S = P[index_cs].Pn[index_cs].nl_a1a2 * *nfS_C * kon_S;
                a_on_a1a2C_sum_S = a_on_a1a2C_sum_S + P[index_cs].aon_a1a2_C_S;

                P[index_cs].aoff_a1a2_C_S = P[index_cs].nbS_a1a2C * koff_S_3s;
                a_off_a1a2C_sum_S = a_off_a1a2C_sum_S + P[index_cs].aoff_a1a2_C_S;

                P[index_cs].aon_a1_C_a2_S = P[index_cs].nbS_a1C * P[index_cs].nf_a2 * kon_L;
                a_on_a1_C_a2_sum_S = a_on_a1_C_a2_sum_S + P[index_cs].aon_a1_C_a2_S;

                P[index_cs].aoff_a1_C_a2_S = P[index_cs].nbS_a1a2C * koff_S_3s;
                a_off_a1_C_a2_sum_S = a_off_a1_C_a2_sum_S + P[index_cs].aoff_a1_C_a2_S;
                //printf("P[index_cs].nbS_a1a2C=%lf\n",P[index_cs].nbS_a1a2C);
// printf("P[index_cs].aoff_a1_C_a2_S=%lf P[index_cs].nf_a2=%lf\n",P[index_cs].aoff_a1_C_a2_S,P[index_cs].nf_a2);
//                 printf("P[index_cs].nbS_a1a2C=%lf P[index_cs].nbS_a1C=%lf P[index_cs].nf_a2=%lf\n",P[index_cs].nbS_a1a2C,P[index_cs].nbS_a1C,P[index_cs].nf_a2);
            }

            if(P[index_cs].type==1)
            {
                P[index_cs].aon_b2C_S = P[index_cs].nf_b2 * *nfS_C * kon_S;
                a_on_b2C_sum_S = a_on_b2C_sum_S + P[index_cs].aon_b2C_S;

                P[index_cs].aoff_b2C_S=P[index_cs].nbS_b2C * koff_S;
                a_off_b2C_sum_S = a_off_b2C_sum_S + P[index_cs].aoff_b2C_S;
            }

            P[index_cs].A_on_a1C_sum_S = a_on_a1C_sum_S;
            P[index_cs].A_off_a1C_sum_S = a_off_a1C_sum_S;
            P[index_cs].A_on_b2C_sum_S = a_on_b2C_sum_S;
            P[index_cs].A_off_b2C_sum_S = a_off_b2C_sum_S;

            P[index_cs].A_tot_sur_2strand = a_on_a1C_sum_S + a_off_a1C_sum_S + a_on_b2C_sum_S + a_off_b2C_sum_S;////affinity of making 2 strand bridges with surface for each particle

            P[index_cs].A_on_a1a2C_sum_S = a_on_a1a2C_sum_S;
            P[index_cs].A_off_a1a2C_sum_S = a_off_a1a2C_sum_S;
            P[index_cs].A_on_a1_C_a2_sum_S = a_on_a1_C_a2_sum_S;
            P[index_cs].A_off_a1_C_a2_sum_S = a_off_a1_C_a2_sum_S;

            P[index_cs].A_tot_sur_3strand = a_on_a1a2C_sum_S + a_off_a1a2C_sum_S + a_on_a1_C_a2_sum_S + a_off_a1_C_a2_sum_S;////affinity of making 3 strand bridges with surface for each particle
// printf("Ptype=%d,a_on_a1_C_a2_sum_S=%lf\n",P[index_cs].type,a_on_a1_C_a2_sum_S);
            P[index_cs].A_tot_PS = P[index_cs].A_tot_sur_2strand + P[index_cs].A_tot_sur_3strand;////affinity of making all type of bridges with surface for each particle

//             A_tot_Part_Sur = A_tot_Part_Sur + P[index_cs].A_tot_PS;////sum over all Particles close to surface
            A_tot_Part_Sur = A_tot_Part_Sur + P[index_cs].A_tot_sur_2strand + P[index_cs].A_tot_sur_3strand;
// printf("Ptype=%d,P[index_cs].A_tot_sur_3strand=%lf\n",P[index_cs].type,P[index_cs].A_tot_sur_3strand);
            P[index_cs].A_tot_PPS = P[index_cs].A_tot_PP + P[index_cs].A_tot_PS;////total affinity of each particle making brigdes with other Particles as well as surface

        }

        a_tot = A_tot_Part_Part + A_tot_Part_Sur;////total affinity sum over all the particles
//         printf("A_tot_Part_Part=%lf\n",A_tot_Part_Part);
//         printf("A_tot_Part_Sur=%lf\n",A_tot_Part_Sur);
//         printf("a_tot=%lf\n",a_tot);
//      if(a_tot==50.0) {exit(0);}
        double Rtao=rand()/(double)RAND_MAX;
        tao=(1.0/a_tot) * log (1/Rtao);

        t_bar = t_bar+tao;
        count_t=count_t+1;

//         printf("While loop working confirmed\n");



        if (t_bar < Delta_tB)
        {     
//             printf("if t_bar < Delta_tB loop working confirmed\n");

//             int reaction=0;
//             double sum = 0.0;
            double Random = rand()/(double)RAND_MAX;
//             double Random = 0.414019;


            double sum = 0.0;
            flag = true;
            for (i=0; i<Np; i++)
            {
//              printf("entered\n");
                // double sum = 0.0;
                int reaction=0;
                int click = 0;
                // flag = true;

            ///////////////////listing down all the probabilities here////////////////
// printf("sum1=%lf\trandom=%lf\n",sum,Random);
                P[i].P_tot_PP_2strand = P[i].A_tot_PP_2strand/a_tot;
                sum = sum + P[i].P_tot_PP_2strand;
// printf("sum#=%lf\n",sum);
//                 printf("P[i].P_tot_PP_2strand=%lf\n",P[i].P_tot_PP_2strand);
                if(sum > Random && flag==true)
                {
                    sum = sum - P[i].P_tot_PP_2strand;
                    P[i].P_on_a1a2_sum=P[i].A_on_a1a2_sum/a_tot;
                    sum = sum + P[i].P_on_a1a2_sum;
// printf("sum1-8=%lf\n",sum);
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=1 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1a2_sum=P[i].A_off_a1a2_sum/a_tot;
                        sum = sum + P[i].P_off_a1a2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=2;
                        }
                    
                    
                    if(flag==true)
                    {
                        P[i].P_on_a2a1_sum=P[i].A_on_a2a1_sum/a_tot;
                        sum = sum + P[i].P_on_a2a1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=111;
                        }
                    if(flag==true)
                    {
                        P[i].P_off_a2a1_sum=P[i].A_off_a2a1_sum/a_tot;
                        sum = sum + P[i].P_off_a2a1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=222;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_a1b1_sum=P[i].A_on_a1b1_sum/a_tot;
                        sum = sum + P[i].P_on_a1b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=3;
//                             printf("sum3=%lf\n",sum);
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1b1_sum=P[i].A_off_a1b1_sum/a_tot;
                        sum = sum + P[i].P_off_a1b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=4 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_a2b2_sum=P[i].A_on_a2b2_sum/a_tot;
                        sum = sum + P[i].P_on_a2b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=5 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a2b2_sum=P[i].A_off_a2b2_sum/a_tot;
                        sum = sum + P[i].P_off_a2b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=6 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_b1b2_sum=P[i].A_on_b1b2_sum/a_tot;
                        sum = sum + P[i].P_on_b1b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=7 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b1b2_sum=P[i].A_off_b1b2_sum/a_tot;
                        sum = sum + P[i].P_off_b1b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=8 ;
                        }
                    
                    
                    
                    if(flag==true)
                    {
                        P[i].P_on_b2b1_sum=P[i].A_on_b2b1_sum/a_tot;
                        sum = sum + P[i].P_on_b2b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=777 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b2b1_sum=P[i].A_off_b2b1_sum/a_tot;
                        sum = sum + P[i].P_off_b2b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=888 ;
                        }

//                     if(flag==true)
//                     {
//                         P[i].P_on_a2a1_sum=P[i].A_on_a2a1_sum/a_tot;
//                         sum = sum + P[i].P_on_a2a1_sum;
//                     }
//                     if(sum > Random && flag==true)
//                         {
//                             flag=false;
//                             reaction=111;
//                         }
//                     if(flag==true)
//                     {
//                         P[i].P_off_a2a1_sum=P[i].A_off_a2a1_sum/a_tot;
//                         sum = sum + P[i].P_off_a2a1_sum;
//                     }
//                     if(sum > Random && flag==true)
//                         {
//                             flag=false;
//                             reaction=222;
//                         }
//                     if(flag==true)
//                     {
//                         P[i].P_on_b2b1_sum=P[i].A_on_b2b1_sum/a_tot;
//                         sum = sum + P[i].P_on_b2b1_sum;
//                     }
//                     if(sum > Random && flag==true)
//                         {
//                             flag=false;
//                             reaction=1111;
//                         }
//                     if(flag==true)
//                     {
//                         P[i].P_off_b2b1_sum=P[i].A_off_b2b1_sum/a_tot;
//                         sum = sum + P[i].P_off_b2b1_sum;
//                     }
//                     if(sum > Random && flag==true)
//                         {
//                             flag=false;
//                             reaction=2222;
//                         }
                }

                if(flag==true)
                {
                    P[i].P_tot_PP_3strand = P[i].A_tot_PP_3strand/a_tot;
                    sum = sum + P[i].P_tot_PP_3strand;
                }
// printf("sum$$=%lf\n",sum);
// printf("random=%lf sum=%lf\n",Random,sum);
                if(sum > Random && flag==true)
                {
//                     printf("sum9-24=%lf\n",sum);
                    sum = sum - P[i].P_tot_PP_3strand;

                    P[i].P_on_a1a2_b1_sum=P[i].A_on_a1a2_b1_sum/a_tot;
                    sum = sum + P[i].P_on_a1a2_b1_sum;
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=9 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1a2_b1_sum=P[i].A_off_a1a2_b1_sum/a_tot;
                        sum = sum + P[i].P_off_a1a2_b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=10 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_a1_b1_a2_sum=P[i].A_on_a1_b1_a2_sum/a_tot;
                        sum = sum + P[i].P_on_a1_b1_a2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=11 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1_b1_a2_sum=P[i].A_off_a1_b1_a2_sum/a_tot;
                        sum = sum + P[i].P_off_a1_b1_a2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=12 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_a1a2_b2_sum=P[i].A_on_a1a2_b2_sum/a_tot;
                        sum = sum + P[i].P_on_a1a2_b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=13 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1a2_b2_sum=P[i].A_off_a1a2_b2_sum/a_tot;
                        sum = sum + P[i].P_off_a1a2_b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=14 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_a2_b2_a1_sum=P[i].A_on_a2_b2_a1_sum/a_tot;
                        sum = sum + P[i].P_on_a2_b2_a1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=15 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a2_b2_a1_sum=P[i].A_off_a2_b2_a1_sum/a_tot;
                        sum = sum + P[i].P_off_a2_b2_a1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=16 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_b1b2_a1_sum=P[i].A_on_b1b2_a1_sum/a_tot;
//                         printf("P[i].P_on_b1b2_a1_sum=%lf\n",P[i].P_on_b1b2_a1_sum);
                        sum = sum + P[i].P_on_b1b2_a1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
//                             printf("random=%lf sum=%lf\n",Random,sum);
//                             printf("a_on_b1b2_a1_sum=%lf\n",a_on_b1b2_a1_sum);
                            flag=false;
                            reaction=17;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b1b2_a1_sum=P[i].A_off_b1b2_a1_sum/a_tot;
                        sum = sum + P[i].P_off_b1b2_a1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=18 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_b1_a1_b2_sum=P[i].A_on_b1_a1_b2_sum/a_tot;
                        sum = sum + P[i].P_on_b1_a1_b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=19 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b1_a1_b2_sum=P[i].A_off_b1_a1_b2_sum/a_tot;
                        sum = sum + P[i].P_off_b1_a1_b2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=20 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_b1b2_a2_sum=P[i].A_on_b1b2_a2_sum/a_tot;
                        sum = sum + P[i].P_on_b1b2_a2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=21 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b1b2_a2_sum=P[i].A_off_b1b2_a2_sum/a_tot;
                        sum = sum + P[i].P_off_b1b2_a2_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=22 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_b2_a2_b1_sum=P[i].A_on_b2_a2_b1_sum/a_tot;
                        sum = sum + P[i].P_on_b2_a2_b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=23 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b2_a2_b1_sum=P[i].A_off_b2_a2_b1_sum/a_tot;
                        sum = sum + P[i].P_off_b2_a2_b1_sum;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=24 ;
                        }
                }

                if(flag==true)
                {
                    P[i].P_tot_loop = P[i].A_tot_loop/a_tot;
                    sum = sum + P[i].P_tot_loop;
//                     printf("sum25-28=%lf\n\n",sum);
                }

                if(sum > Random && flag==true)
                {
                    sum = sum - P[i].P_tot_loop;
// printf("sum25-28**=%lf\n\n",sum);
                    P[i].P_on_sum_L_a1a2=P[i].A_on_sum_L_a1a2/a_tot;
                    sum = sum + P[i].P_on_sum_L_a1a2;
                    
                    cout << "i = " <<  i << " " << "P[i].P_on_sum_L_a1a2 = " << " " << P[i].P_on_sum_L_a1a2  << endl;
// printf("sum25-28&&=%lf\n\n",sum);
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=25 ;
//                             printf("sum25=%lf\n\n",sum);
                        }

                    if(flag==true)
                    {
                        P[i].P_off_sum_L_a1a2=P[i].A_off_sum_L_a1a2/a_tot;
                        sum = sum + P[i].P_off_sum_L_a1a2;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=26 ;
//                             printf("sum26=%lf\n\n",sum);
                        }

                    if(flag==true)
                    {
                        P[i].P_on_sum_L_b1b2=P[i].A_on_sum_L_b1b2/a_tot;
                        sum = sum + P[i].P_on_sum_L_b1b2;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=27 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_sum_L_b1b2=P[i].A_off_sum_L_b1b2/a_tot;
                        sum = sum + P[i].P_off_sum_L_b1b2;
//                         printf("P_off_sum_L_b1b2=%lf\n",P_off_sum_L_b1b2);
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=28 ;
                        }
                }

                if(flag==true)
                {
                    P[i].P_tot_sur_2strand = P[i].A_tot_sur_2strand/a_tot;
                    sum = sum + P[i].P_tot_sur_2strand;
//                     printf("sumsur2=%lf\n",sum);
                }

                if(sum > Random && flag==true)
                {
                    sum = sum - P[i].P_tot_sur_2strand;

                    P[i].P_on_a1C_sum_S=P[i].A_on_a1C_sum_S/a_tot;
                    sum = sum + P[i].P_on_a1C_sum_S;
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=29 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1C_sum_S=P[i].A_off_a1C_sum_S/a_tot;
                        sum = sum + P[i].P_off_a1C_sum_S;
//                         printf("P[i].A_off_a1C_sum_S=%lf\n",P[i].A_off_a1C_sum_S);
                    }
                    if(sum > Random && flag==true)
                        {
//                             printf("P[%d].A_off_a1C_sum_S=%lf\tP[%d].nbS_a1C=%lf\n",i,P[i].A_off_a1C_sum_S,i,P[i].nbS_a1C);
                            flag=false;
                            reaction=30;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_b2C_sum_S=P[i].A_on_b2C_sum_S/a_tot;
                        sum = sum + P[i].P_on_b2C_sum_S;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=31 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_b2C_sum_S=P[i].A_off_b2C_sum_S/a_tot;
                        sum = sum + P[i].P_off_b2C_sum_S;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=32 ;
                        }
                }

                if(flag==true)
                {
                    P[i].P_tot_sur_3strand = P[i].A_tot_sur_3strand/a_tot;
                    sum = sum + P[i].P_tot_sur_3strand;
                }

                if(sum > Random && flag==true)
                {
//                     printf("3strand_sur\n");
                    sum = sum - P[i].P_tot_sur_3strand;

                    P[i].P_on_a1a2C_sum_S=P[i].A_on_a1a2C_sum_S/a_tot;
                    sum = sum + P[i].P_on_a1a2C_sum_S;
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=33 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1a2C_sum_S=P[i].A_off_a1a2C_sum_S/a_tot;
                        sum = sum + P[i].P_off_a1a2C_sum_S;
                        cout << "i = " <<  i << " " << "P[i].P_off_a1a2C_sum_S = " << " " << P[i].P_off_a1a2C_sum_S  << endl;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=34 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_on_a1_C_a2_sum_S=P[i].A_on_a1_C_a2_sum_S/a_tot;
                        sum = sum + P[i].P_on_a1_C_a2_sum_S;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=35 ;
                        }

                    if(flag==true)
                    {
                        P[i].P_off_a1_C_a2_sum_S=P[i].A_off_a1_C_a2_sum_S/a_tot;
                        sum = sum + P[i].P_off_a1_C_a2_sum_S;
                    }
                    if(sum > Random && flag==true)
                        {
                            flag=false;
                            reaction=36 ;
                        }
                }
//                 printf("sum_final=%lf\n\n",sum);
//                 exit(0);
//              printf("random=%lf\tPtype=%d\treaction=%d\n",Random,P[i].type,reaction);
//              exit(0);
//              if(reaction==17){printf("sum0=%lf\n",sum);}
// printf("sum0=%lf\trandom=%lf\n",sum,Random);
        ///////////////////////*********Reaction Firing************///////////////////////
//              printf("type=%d Nb1=%lf\n",P[i].type,P[i].nf_b1);

                if(flag==false && reaction>0)
                {
//                     double Pij= 0.0;
//                  printf("hii\n");
                    if (reaction==1)
                    {
                        sum = sum - P[i].P_on_a1a2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==0 && i<index) //////Suggested by AKJ
                            {
                                P[i].Pn[index].pon_a1a2=P[i].Pn[index].aon_a1a2/a_tot;
                                sum=sum+P[i].Pn[index].pon_a1a2;
//                                 printf("yes**\n");
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2 = P[i].Pn[index].nb_a1a2 + 1;
//                                     P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 + 1;
                                    P[index].Pn[i].nb_a1a2 = P[index].Pn[i].nb_a1a2 + 1;
                                    P[i].nf_a1 = P[i].nf_a1 - 1;
                                    P[index].nf_a2 = P[index].nf_a2 - 1;
//                                     printf("yes\n");
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==2)
                    {
                        sum = sum - P[i].P_off_a1a2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==0 && i<index) //////Suggested by AKJ
                            {
                                P[i].Pn[index].poff_a1a2=P[i].Pn[index].aoff_a1a2/a_tot;
                                sum=sum+P[i].Pn[index].poff_a1a2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2 = P[i].Pn[index].nb_a1a2 - 1;
//                                     P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 - 1;
                                    P[index].Pn[i].nb_a1a2 = P[index].Pn[i].nb_a1a2 - 1;
                                    P[i].nf_a1 = P[i].nf_a1 + 1;
                                    P[index].nf_a2 = P[index].nf_a2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }
                    
                    
                    if (reaction==111)
                    {
                        sum = sum - P[i].P_on_a2a1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==0 && i>index)
                            {
                                P[i].Pn[index].pon_a2a1=P[i].Pn[index].aon_a2a1/a_tot;
                                sum=sum+P[i].Pn[index].pon_a2a1;
//                                 printf("yes**\n");
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a2a1 = P[i].Pn[index].nb_a2a1 + 1;
//                                     P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 + 1;
                                    P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 + 1;
                                    P[i].nf_a1 = P[i].nf_a1 - 1;
                                    P[index].nf_a2 = P[index].nf_a2 - 1;
//                                     printf("yes\n");
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }
                    
                    if (reaction==222)
                    {
                        sum = sum - P[i].P_off_a2a1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==0 && i>index)
                            {
                                P[i].Pn[index].poff_a2a1=P[i].Pn[index].aoff_a2a1/a_tot;
                                sum=sum+P[i].Pn[index].poff_a2a1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a2a1 = P[i].Pn[index].nb_a2a1 - 1;
//                                  P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 - 1;
                                    P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 - 1;
                                    P[i].nf_a1 = P[i].nf_a1 + 1;
                                    P[index].nf_a2 = P[index].nf_a2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }
                    

                    if (reaction==3)
                    {
                        sum = sum - P[i].P_on_a1b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].pon_a1b1=P[i].Pn[index].aon_a1b1/a_tot;
                                sum=sum+P[i].Pn[index].pon_a1b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1b1 = P[i].Pn[index].nb_a1b1 + 1;
                                    P[index].Pn[i].nb_a1b1 = P[index].Pn[i].nb_a1b1 + 1;
                                    P[i].nf_a1 = P[i].nf_a1 - 1;
                                    P[index].nf_b1 = P[index].nf_b1 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==4)
                    {
                        sum = sum - P[i].P_off_a1b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].poff_a1b1=P[i].Pn[index].aoff_a1b1/a_tot;
                                sum=sum+P[i].Pn[index].poff_a1b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1b1 = P[i].Pn[index].nb_a1b1 - 1;
                                    P[index].Pn[i].nb_a1b1 = P[index].Pn[i].nb_a1b1 - 1;
                                    P[i].nf_a1 = P[i].nf_a1 + 1;
                                    P[index].nf_b1 = P[index].nf_b1 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==5)
                    {
                        sum = sum - P[i].P_on_a2b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].pon_a2b2=P[i].Pn[index].aon_a2b2/a_tot;
                                sum=sum+P[i].Pn[index].pon_a2b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a2b2 = P[i].Pn[index].nb_a2b2 + 1;
                                    P[index].Pn[i].nb_a2b2 = P[index].Pn[i].nb_a2b2 + 1;
                                    P[i].nf_a2 = P[i].nf_a2 - 1;
                                    P[index].nf_b2 = P[index].nf_b2 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==6)
                    {
                        sum = sum - P[i].P_off_a2b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].poff_a2b2=P[i].Pn[index].aoff_a2b2/a_tot;
                                sum=sum+P[i].Pn[index].poff_a2b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a2b2 = P[i].Pn[index].nb_a2b2 - 1;
                                    P[index].Pn[i].nb_a2b2 = P[index].Pn[i].nb_a2b2 - 1;
                                    P[i].nf_a2 = P[i].nf_a2 + 1;
                                    P[index].nf_b2 = P[index].nf_b2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==7)
                    {
                        sum = sum - P[i].P_on_b1b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==1 && i<index) //////////////Suggested By AKJ
                            {
                                P[i].Pn[index].pon_b1b2=P[i].Pn[index].aon_b1b2/a_tot;
                                sum=sum+P[i].Pn[index].pon_b1b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2 = P[i].Pn[index].nb_b1b2 + 1;
//                                     P[index].Pn[i].nb_b2b1 = P[index].Pn[i].nb_b2b1 + 1;
                                    P[index].Pn[i].nb_b1b2 = P[index].Pn[i].nb_b1b2 + 1;
                                    P[i].nf_b1 = P[i].nf_b1 - 1;
                                    P[index].nf_b2 = P[index].nf_b2 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==8)
                    {
                        sum = sum - P[i].P_off_b1b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==1 && i<index) //////////////Suggested By AKJ
                            {
                                P[i].Pn[index].poff_b1b2=P[i].Pn[index].aoff_b1b2/a_tot;
                                sum=sum+P[i].Pn[index].poff_b1b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2 = P[i].Pn[index].nb_b1b2 - 1;
//                                     P[index].Pn[i].nb_b2b1 = P[index].Pn[i].nb_b2b1 - 1;
                                    P[index].Pn[i].nb_b1b2 = P[index].Pn[i].nb_b1b2 - 1;
                                    P[i].nf_b1 = P[i].nf_b1 + 1;
                                    P[index].nf_b2 = P[index].nf_b2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }
                    
                    
                    
                    if (reaction==777)
                    {
                        sum = sum - P[i].P_on_b2b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==1 && i>index)
                            {
                                P[i].Pn[index].pon_b2b1=P[i].Pn[index].aon_b2b1/a_tot;
                                sum=sum+P[i].Pn[index].pon_b2b1;
//                                 printf("yes**\n");
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b2b1 = P[i].Pn[index].nb_b2b1 + 1;
//                                     P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 + 1;
                                    P[index].Pn[i].nb_b2b1 = P[index].Pn[i].nb_b2b1 + 1;
                                    P[i].nf_b1 = P[i].nf_b1 - 1;
                                    P[index].nf_b2 = P[index].nf_b2 - 1;
//                                     printf("yes\n");
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }
                    
                    if (reaction==888)
                    {
                        sum = sum - P[i].P_off_b2b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==1 && i>index)
                            {
                                P[i].Pn[index].poff_b2b1=P[i].Pn[index].aoff_b2b1/a_tot;
                                sum=sum+P[i].Pn[index].poff_b2b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b2b1 = P[i].Pn[index].nb_b2b1 - 1;
//                                  P[index].Pn[i].nb_a2a1 = P[index].Pn[i].nb_a2a1 - 1;
                                    P[index].Pn[i].nb_b2b1 = P[index].Pn[i].nb_b2b1 - 1;
                                    P[i].nf_b1 = P[i].nf_b1 + 1;
                                    P[index].nf_b2 = P[index].nf_b2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }



//                     if (reaction==111)
//                     {
//                         sum = sum - P[i].P_on_a2a1_sum;
//                         for (j=0; j<P[i].maxn; j++)
//                         {
// //                             printf("i=%d P[i].maxnn=%d\n",i,P[i].maxnn);
//                             int index=P[i].Pn[j].n_index;
//                             if(P[i].type==0 && P[index].type==0)
//                             {
//                                 P[i].Pn[index].pon_a2a1=P[i].Pn[index].aon_a2a1/a_tot;
//                                 sum=sum+P[i].Pn[index].pon_a2a1;
// //                                 printf("yes**\n");
//                                 if (sum > Random)
//                                 {
//                                     P[i].Pn[index].nb_a2a1 = P[i].Pn[index].nb_a2a1 + 1;
//                                     P[index].Pn[i].nb_a1a2 = P[index].Pn[i].nb_a1a2 + 1;
//                                     P[i].nf_a2 = P[i].nf_a2 - 1;
//                                     P[index].nf_a1 = P[index].nf_a1 - 1;
// //                                     flag=true;
// //                                     printf("yes\n");
//                                     break;
//                                 }
//                             }
//                         }
// //                         if (flag==true) break;
//                     }
//
//                     if (reaction==222)
//                     {
//                         sum = sum - P[i].P_off_a2a1_sum;
//                         for (j=0; j<P[i].maxn; j++)
//                         {
// //                             printf("i=%d P[i].maxnn=%d\n",i,P[i].maxnn);
//                             int index=P[i].Pn[j].n_index;
//                             if(P[i].type==0 && P[index].type==0)
//                             {
//                                 P[i].Pn[index].poff_a2a1=P[i].Pn[index].aoff_a2a1/a_tot;
//                                 sum=sum+P[i].Pn[index].poff_a2a1;
//                                 if (sum > Random)
//                                 {
//                                     P[i].Pn[index].nb_a2a1 = P[i].Pn[index].nb_a2a1 - 1;
//                                     P[index].Pn[i].nb_a1a2 = P[index].Pn[i].nb_a1a2 - 1;
//                                     P[i].nf_a2 = P[i].nf_a2 + 1;
//                                     P[index].nf_a1 = P[index].nf_a1 + 1;
// //                                     flag=true;
//                                     break;
//                                 }
//                             }
//                         }
// //                         if (flag==true) break;
//                     }
//                     if (reaction==1111)
//                     {
//                         sum = sum - P[i].P_on_b2b1_sum;
//                         for (j=0; j<P[i].maxn; j++)
//                         {
// //                             printf("i=%d P[i].maxnn=%d\n",i,P[i].maxnn);
//                             int index=P[i].Pn[j].n_index;
//                             if(P[i].type==1 && P[index].type==1)
//                             {
//                                 P[i].Pn[index].pon_b2b1=P[i].Pn[index].aon_b2b1/a_tot;
//                                 sum=sum+P[i].Pn[index].pon_b2b1;
// //                                 printf("yes**\n");
//                                 if (sum > Random)
//                                 {
//                                     P[i].Pn[index].nb_b2b1 = P[i].Pn[index].nb_b2b1 + 1;
//                                     P[index].Pn[i].nb_b1b2 = P[index].Pn[i].nb_b1b2 + 1;
//                                     P[i].nf_b2 = P[i].nf_b2 - 1;
//                                     P[index].nf_b1 = P[index].nf_b1 - 1;
// //                                     flag=true;
// //                                     printf("yes\n");
//                                     break;
//                                 }
//                             }
//                         }
// //                         if (flag==true) break;
//                     }
//
//                     if (reaction==2222)
//                     {
//                         sum = sum - P[i].P_off_b2b1_sum;
//                         for (j=0; j<P[i].maxn; j++)
//                         {
// //                             printf("i=%d P[i].maxnn=%d\n",i,P[i].maxnn);
//                             int index=P[i].Pn[j].n_index;
//                             if(P[i].type==1 && P[index].type==1)
//                             {
//                                 P[i].Pn[index].poff_b2b1=P[i].Pn[index].aoff_b2b1/a_tot;
//                                 sum=sum+P[i].Pn[index].poff_b2b1;
//                                 if (sum > Random)
//                                 {
//                                     P[i].Pn[index].nb_b2b1 = P[i].Pn[index].nb_b2b1 - 1;
//                                     P[index].Pn[i].nb_b1b2 = P[index].Pn[i].nb_b1b2 - 1;
//                                     P[i].nf_b2 = P[i].nf_b2 + 1;
//                                     P[index].nf_b1 = P[index].nf_b1 + 1;
// //                                     flag=true;
//                                     break;
//                                 }
//                             }
//                         }
// //                         if (flag==true) break;
//                     }




                    if (reaction==9)
                    {
                        sum = sum - P[i].P_on_a1a2_b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].pon_a1a2_b1=P[i].Pn[index].aon_a1a2_b1/a_tot;
                                sum=sum+P[i].Pn[index].pon_a1a2_b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b1 = P[i].Pn[index].nb_a1a2b1 + 1;
                                    P[index].Pn[i].nb_a1a2b1 = P[index].Pn[i].nb_a1a2b1 + 1;
                                    P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 -1;
                                    P[index].nf_b1 = P[index].nf_b1 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==10)
                    {
                        sum = sum - P[i].P_off_a1a2_b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].poff_a1a2_b1=P[i].Pn[index].aoff_a1a2_b1/a_tot;
                                sum=sum+P[i].Pn[index].poff_a1a2_b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b1 = P[i].Pn[index].nb_a1a2b1 - 1;
                                    P[index].Pn[i].nb_a1a2b1 = P[index].Pn[i].nb_a1a2b1 - 1;
                                    P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 + 1;
                                    P[index].nf_b1 = P[index].nf_b1 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==11)
                    {
                        sum = sum - P[i].P_on_a1_b1_a2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].pon_a1_b1_a2=P[i].Pn[index].aon_a1_b1_a2/a_tot;
                                sum=sum+P[i].Pn[index].pon_a1_b1_a2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b1 = P[i].Pn[index].nb_a1a2b1 + 1;
                                    P[index].Pn[i].nb_a1a2b1 = P[index].Pn[i].nb_a1a2b1 + 1;
                                    P[i].Pn[index].nb_a1b1 = P[i].Pn[index].nb_a1b1 -1;
                                    P[index].Pn[i].nb_a1b1 = P[index].Pn[i].nb_a1b1 -1;
                                    P[i].nf_a2 = P[i].nf_a2-1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==12)
                    {
                        sum = sum - P[i].P_off_a1_b1_a2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].poff_a1_b1_a2=P[i].Pn[index].aoff_a1_b1_a2/a_tot;
                                sum=sum+P[i].Pn[index].poff_a1_b1_a2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b1 = P[i].Pn[index].nb_a1a2b1 - 1;
                                    P[index].Pn[i].nb_a1a2b1 = P[index].Pn[i].nb_a1a2b1 - 1;
                                    P[i].Pn[index].nb_a1b1 = P[i].Pn[index].nb_a1b1 + 1;
                                    P[index].Pn[i].nb_a1b1 = P[index].Pn[i].nb_a1b1 + 1;
                                    P[i].nf_a2 = P[i].nf_a2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==13)
                    {
                        sum = sum - P[i].P_on_a1a2_b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].pon_a1a2_b2=P[i].Pn[index].aon_a1a2_b2/a_tot;
                                sum=sum+P[i].Pn[index].pon_a1a2_b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b2 = P[i].Pn[index].nb_a1a2b2 + 1;
                                    P[index].Pn[i].nb_a1a2b2 = P[index].Pn[i].nb_a1a2b2 + 1;
                                    P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 -1;
                                    P[index].nf_b2 = P[index].nf_b2 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==14)
                    {
                        sum = sum - P[i].P_off_a1a2_b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].poff_a1a2_b2=P[i].Pn[index].aoff_a1a2_b2/a_tot;
                                sum=sum+P[i].Pn[index].poff_a1a2_b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b2 = P[i].Pn[index].nb_a1a2b2 - 1;
                                    P[index].Pn[i].nb_a1a2b2 = P[index].Pn[i].nb_a1a2b2 - 1;
                                    P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 + 1;
                                    P[index].nf_b2 = P[index].nf_b2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==15)
                    {
                        sum = sum - P[i].P_on_a2_b2_a1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].pon_a2_b2_a1=P[i].Pn[index].aon_a2_b2_a1/a_tot;
                                sum=sum+P[i].Pn[index].pon_a2_b2_a1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b2 = P[i].Pn[index].nb_a1a2b2 + 1;
                                    P[index].Pn[i].nb_a1a2b2 = P[index].Pn[i].nb_a1a2b2 + 1;
                                    P[i].Pn[index].nb_a2b2 = P[i].Pn[index].nb_a2b2 - 1;
                                    P[index].Pn[i].nb_a2b2 = P[index].Pn[i].nb_a2b2 - 1;
                                    P[i].nf_a1 = P[i].nf_a1 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==16)
                    {
                        sum = sum - P[i].P_off_a2_b2_a1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==0 && P[index].type==1)
                            {
                                P[i].Pn[index].poff_a2_b2_a1=P[i].Pn[index].aoff_a2_b2_a1/a_tot;
                                sum=sum+P[i].Pn[index].poff_a2_b2_a1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_a1a2b2 = P[i].Pn[index].nb_a1a2b2 - 1;
                                    P[index].Pn[i].nb_a1a2b2 = P[index].Pn[i].nb_a1a2b2 - 1;
                                    P[i].Pn[index].nb_a2b2 = P[i].Pn[index].nb_a2b2 + 1;
                                    P[index].Pn[i].nb_a2b2 = P[index].Pn[i].nb_a2b2 + 1;
                                    P[i].nf_a1 = P[i].nf_a1 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==17)
                    {
                        sum = sum - P[i].P_on_b1b2_a1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].pon_b1b2_a1=P[i].Pn[index].aon_b1b2_a1/a_tot;
                                sum=sum+P[i].Pn[index].pon_b1b2_a1;

                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a1 = P[i].Pn[index].nb_b1b2a1 + 1;
                                    P[index].Pn[i].nb_b1b2a1 = P[index].Pn[i].nb_b1b2a1 + 1;
                                    P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 - 1;
                                    P[index].nf_a1 = P[index].nf_a1 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==18)
                    {
                        sum = sum - P[i].P_off_b1b2_a1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].poff_b1b2_a1=P[i].Pn[index].aoff_b1b2_a1/a_tot;
                                sum=sum+P[i].Pn[index].poff_b1b2_a1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a1 = P[i].Pn[index].nb_b1b2a1 - 1;
                                    P[index].Pn[i].nb_b1b2a1 = P[index].Pn[i].nb_b1b2a1 - 1;
                                    P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 + 1;
                                    P[index].nf_a1 = P[index].nf_a1 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==19)
                    {
                        sum = sum - P[i].P_on_b1_a1_b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].pon_b1_a1_b2=P[i].Pn[index].aon_b1_a1_b2/a_tot;
                                sum=sum+P[i].Pn[index].pon_b1_a1_b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a1 = P[i].Pn[index].nb_b1b2a1 + 1;
                                    P[index].Pn[i].nb_b1b2a1 = P[index].Pn[i].nb_b1b2a1 + 1;
                                    P[i].Pn[index].nb_a1b1 = P[i].Pn[index].nb_a1b1 - 1;
                                    P[index].Pn[i].nb_a1b1 = P[index].Pn[i].nb_a1b1 - 1;
                                    P[i].nf_b2 = P[i].nf_b2 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==20)
                    {
                        sum = sum - P[i].P_off_b1_a1_b2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].poff_b1_a1_b2=P[i].Pn[index].aoff_b1_a1_b2/a_tot;
                                sum=sum+P[i].Pn[index].poff_b1_a1_b2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a1 = P[i].Pn[index].nb_b1b2a1 - 1;
                                    P[index].Pn[i].nb_b1b2a1 = P[index].Pn[i].nb_b1b2a1 - 1;
                                    P[i].Pn[index].nb_a1b1 = P[i].Pn[index].nb_a1b1 + 1;
                                    P[index].Pn[i].nb_a1b1 = P[index].Pn[i].nb_a1b1 + 1;
                                    P[i].nf_b2 = P[i].nf_b2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==21)
                    {
                        sum = sum - P[i].P_on_b1b2_a2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].pon_b1b2_a2=P[i].Pn[index].aon_b1b2_a2/a_tot;
                                sum=sum+P[i].Pn[index].pon_b1b2_a2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a2 = P[i].Pn[index].nb_b1b2a2 + 1;
                                    P[index].Pn[i].nb_b1b2a2 = P[index].Pn[i].nb_b1b2a2 + 1;
                                    P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 - 1;
                                    P[index].nf_a2 = P[index].nf_a2 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==22)
                    {
                        sum = sum - P[i].P_off_b1b2_a2_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].poff_b1b2_a2=P[i].Pn[index].aoff_b1b2_a2/a_tot;
                                sum=sum+P[i].Pn[index].poff_b1b2_a2;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a2 = P[i].Pn[index].nb_b1b2a2 - 1;
                                    P[index].Pn[i].nb_b1b2a2 = P[index].Pn[i].nb_b1b2a2 - 1;
                                    P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 + 1;
                                    P[index].nf_a2 = P[index].nf_a2 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==23)
                    {
                        sum = sum - P[i].P_on_b2_a2_b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].pon_b2_a2_b1=P[i].Pn[index].aon_b2_a2_b1/a_tot;
                                sum=sum+P[i].Pn[index].pon_b2_a2_b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a2 = P[i].Pn[index].nb_b1b2a2 + 1;
                                    P[index].Pn[i].nb_b1b2a2 = P[index].Pn[i].nb_b1b2a2 + 1;
                                    P[i].Pn[index].nb_a2b2 = P[i].Pn[index].nb_a2b2 - 1;
                                    P[index].Pn[i].nb_a2b2 = P[index].Pn[i].nb_a2b2 - 1;
                                    P[i].nf_b1 = P[i].nf_b1 - 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }

                    if (reaction==24)
                    {
                        sum = sum - P[i].P_off_b2_a2_b1_sum;
                        for (j=0; j<P[i].maxnn; j++)
                        {
                            int index=P[i].Pn[j].nn_index;
                            if(P[i].type==1 && P[index].type==0)
                            {
                                P[i].Pn[index].poff_b2_a2_b1=P[i].Pn[index].aoff_b2_a2_b1/a_tot;
                                sum=sum+P[i].Pn[index].poff_b2_a2_b1;
                                if (sum > Random)
                                {
                                    P[i].Pn[index].nb_b1b2a2 = P[i].Pn[index].nb_b1b2a2 - 1;
                                    P[index].Pn[i].nb_b1b2a2 = P[index].Pn[i].nb_b1b2a2 - 1;
                                    P[i].Pn[index].nb_a2b2 = P[i].Pn[index].nb_a2b2 + 1;
                                    P[index].Pn[i].nb_a2b2 = P[index].Pn[i].nb_a2b2 + 1;
                                    P[i].nf_b1 = P[i].nf_b1 + 1;
                                    click = 1;
                                    break;
                                }
                            }
                        }
                    }
// printf("P[i].Pn[i].nl_a1a2=%lf\tP[i].nf_a1=%lf\tP[i].nf_a2=%lf\n",P[i].Pn[i].nl_a1a2,P[i].nf_a1,P[i].nf_a2);
                    if (reaction==25 && P[i].type==0)
                    {
                        cout << "i = " << i << " " << " P[i].Pn[i].nl_a1a2 " << " " <<  P[i].Pn[i].nl_a1a2 << endl;
                        P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 + 1;
                        P[i].nf_a1 = P[i].nf_a1-1;
                        P[i].nf_a2 = P[i].nf_a2-1;
                        click = 1;
//                         printf("hii\n");
                    }

// exit(0);
                    if (reaction==26 && P[i].type==0)
                    {
                        P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 - 1;
                        P[i].nf_a1 = P[i].nf_a1 + 1;
                        P[i].nf_a2 = P[i].nf_a2 + 1;
                        click = 1;
                    }

                    if (reaction==27 && P[i].type==1)
                    {
                        P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 + 1;
                        P[i].nf_b1 = P[i].nf_b1 - 1;
                        P[i].nf_b2 = P[i].nf_b2 - 1;
                        click = 1;
                    }
                    if (reaction==28 && P[i].type==1)
                    {
                        P[i].Pn[i].nl_b1b2 = P[i].Pn[i].nl_b1b2 - 1;
                        P[i].nf_b1 = P[i].nf_b1 + 1;
                        P[i].nf_b2 = P[i].nf_b2 + 1;
                        click = 1;
                    }
//                     printf("P[i].Pn[i].nl_b1b2=%lf\tP[i].nf_b1=%lf\tP[i].nf_b2=%lf\n",P[i].Pn[i].nl_b1b2,P[i].nf_b1,P[i].nf_b2);

                    if (reaction==29)
                    {
                        P[i].nbS_a1C = P[i].nbS_a1C + 1;
                        P[i].nf_a1 = P[i].nf_a1-1;
                        *nfS_C = *nfS_C-1;
                        click = 1;
                    }

                    if (reaction==30)
                    {
//                         printf("P[%d].nbS_a1C=%lf\n",i,P[i].nbS_a1C);
                        P[i].nbS_a1C = P[i].nbS_a1C - 1;
                        P[i].nf_a1 = P[i].nf_a1 + 1;
                        *nfS_C = *nfS_C + 1;
                        click = 1;
                    }

                    if (reaction==31)
                    {
                        P[i].nbS_b2C = P[i].nbS_b2C + 1;
                        P[i].nf_b2 = P[i].nf_b2-1;
                        *nfS_C = *nfS_C-1;
                        click = 1;
                    }

                    if (reaction==32)
                    {
                        P[i].nbS_b2C = P[i].nbS_b2C - 1;
                        P[i].nf_b2 = P[i].nf_b2 + 1;
                        *nfS_C = *nfS_C + 1;
                        click = 1;
                    }

                    if (reaction==33)
                    {
                        P[i].nbS_a1a2C = P[i].nbS_a1a2C + 1;
                        P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 - 1;
                        *nfS_C = *nfS_C-1;
                        click = 1;
                    }

                    if (reaction==34)
                    {
                        P[i].nbS_a1a2C = P[i].nbS_a1a2C - 1;
                        P[i].Pn[i].nl_a1a2 = P[i].Pn[i].nl_a1a2 + 1;
                        
                        cout << "34i = " << i << " " << " P[i].Pn[i].nl_a1a2 " << " " <<  P[i].Pn[i].nl_a1a2 << endl;
                        *nfS_C = *nfS_C + 1;
                        click = 1;
                    }

                    if (reaction==35)
                    {
                        P[i].nbS_a1a2C = P[i].nbS_a1a2C + 1;
                        P[i].nbS_a1C = P[i].nbS_a1C - 1;
                        P[i].nf_a2 = P[i].nf_a2 - 1;
                        click = 1;
                    }

                    if (reaction==36)
                    {
                        P[i].nbS_a1a2C = P[i].nbS_a1a2C - 1;
                        P[i].nbS_a1C = P[i].nbS_a1C + 1;
                        P[i].nf_a2 = P[i].nf_a2 + 1;
                        click = 1;
                    }
                    if (click == 1){break;}
                }
            }
//             printf("P[i].Pn[i].nl_a1a2**=%lf\n",P[i].Pn[i].nl_a1a2);
        }
//         exit(0);
    }

    for(i=0; i<Np; i++)
    {
        //tot_bridges=0.0;
        if(P[i].type==0)
        {
        P[i].total_free_loops=P[i].nf_a1+P[i].nf_a2+P[i].Pn[i].nl_a1a2;
//         printf("P[i].Pn[i].nl_a1a2=%lf\n",P[i].Pn[i].nl_a1a2);
        }
        if(P[i].type==1)
        {
        P[i].total_free_loops=P[i].nf_b1+P[i].nf_b2+P[i].Pn[i].nl_b1b2;
//         printf("P[i].Pn[i].nl_b1b2=%lf\n",P[i].Pn[i].nl_b1b2);
        }
//         P[i].Pn[i].nl = P[i].Pn[i].nl_a1a2 + P[i].Pn[i].nl_b1b2;

        if(P[i].type==0)
        {
        P[i].Pn[i].nl = P[i].Pn[i].nl_a1a2;
        }
        if(P[i].type==1)
        {
        P[i].Pn[i].nl = P[i].Pn[i].nl_b1b2;
        }
//         printf("%d,P[i].nf_a1=%lf\tP[i].nf_a2=%lf\tP[i].nf_b1=%lf\t P[i].nf_b2=%lf\t P[i].Pn[i].nl_a1a2=%lf\tP[i].Pn[i].nl_b1b2=%lf\n\n\n",i,P[i].nf_a1,P[i].nf_a2,P[i].nf_b1,P[i].nf_b2,P[i].Pn[i].nl_a1a2,P[i].Pn[i].nl_b1b2);
        for (j=0; j<P[i].maxnn; j++)
        {
            int index = P[i].Pn[j].nn_index;
            P[i].Pn[index].tot_2strand_bridges = P[i].Pn[index].nb_a1a2 + P[i].Pn[index].nb_b1b2 + P[i].Pn[index].nb_a1b1 + P[i].Pn[index].nb_a2b2 + P[i].Pn[index].nb_a2a1 + P[i].Pn[index].nb_b2b1;
//             printf("type=%d P[i].Pn[index].nb_a2b2=%lf\n",P[i].type,P[i].Pn[index].nb_a2b2);
            P[i].Pn[index].tot_3strand_bridges = P[i].Pn[index].nb_a1a2b1 + P[i].Pn[index].nb_a1a2b2 + P[i].Pn[index].nb_b1b2a1 + P[i].Pn[index].nb_b1b2a2;
//             printf("type=%d P[i].Pn[index].tot_3strand_bridges=%lf\n",P[i].type,P[i].Pn[index].tot_3strand_bridges);
//             printf("P=%d\tn=%d\n",P[i].type,P[index].type);
//             printf("P[i].Pn[index].nb_a1a2=%lf\tP[i].Pn[index].nb_b1b2=%lf\t P[i].Pn[index].nb_a1b1=%lf\t P[i].Pn[index].nb_a2b2=%lf\t P[i].Pn[index].tot_2strand_bridges=%lf\n\n\n",P[i].Pn[index].nb_a1a2,P[i].Pn[index].nb_b1b2,P[i].Pn[index].nb_a1b1,P[i].Pn[index].nb_a2b2,P[i].Pn[index].tot_2strand_bridges);

//             printf("type=%d ntype=%d P[i].Pn[index].nb_a1a2b1=%lf\tP[i].Pn[index].nb_a1a2b2=%lf\t P[i].Pn[index].nb_b1b2a1=%lf\t P[i].Pn[index].nb_b1b2a2=%lf\t P[i].Pn[index].tot_3strand_bridges=%lf\n\n\n",P[i].type,P[index].type,P[i].Pn[index].nb_a1a2b1,P[i].Pn[index].nb_a1a2b2,P[i].Pn[index].nb_b1b2a1,P[i].Pn[index].nb_b1b2a2,P[i].Pn[index].tot_3strand_bridges);
//             P[i].total_bridge_WC= P[i].total_bridge_WC + P[i].Pn[index].nb_AB + P[i].Pn[index].nb_BA;
            //cout << P[i].Pn[index].nb_AB << " " << P[i].Pn[index].nb_BA << endl;
        }
    }
    for (q=0; q<ns; q++)
    {
//         printf("ns=%d\n",ns);
        int index_cs=P[q].Sc;
        P[index_cs].tot_2strand_bridges_surface = P[index_cs].nbS_a1C + P[index_cs].nbS_b2C;
        P[index_cs].tot_3strand_bridges_surface = P[index_cs].nbS_a1a2C;
//         total_surfacebridge = total_surface_2strand_bridge + total_surface_3strand_bridge;
    }
    
}


void Force(int Np, int ns, int ns_upper, double Rc, double L, double R2s,double Sz, double Sz_upper, double omega_S, double nfS_C, double beta,double Lx, double Ly, double Lz, double omega_0_inert)
{
    double dx,dy,dz,d,ov,ov_prime;
    double term1x, term1y, term1z; 
    double term2x, term2y, term2z;
    double term3x, term3y, term3z;
    double sum_term1x, sum_term1y, sum_term1z; 
    double sum_term2x, sum_term2y, sum_term2z;
    double sum_term3x, sum_term3y, sum_term3z;
    int i,j,q;
    for (i=0; i<Np; i++)
    {
        sum_term1x=0.0;
        sum_term1y=0.0;
        sum_term1z=0.0;
            
        sum_term2x=0.0;
        sum_term2y=0.0;
        sum_term2z=0.0;
            
        sum_term3x=0.0;
        sum_term3y=0.0;
        sum_term3z=0.0;
            
        //sum_term3=0.0;
        for (j=0; j<P[i].maxnn; j++)
        {
            int index=P[i].Pn[j].nn_index;
            distance(i,index,Lx,Ly,Lz,&dx,&dy,&dz,&d);
                
            double ux=dx/P[i].Pn[index].d;
            double uy=dy/P[i].Pn[index].d;
            double uz=dz/P[i].Pn[index].d;    
            //cout << i << " " <<index << endl;
//             printf("P[i].total_free_loops=%lf\n",P[i].total_free_loops);
            term1x = (P[i].Pn[index].tot_2strand_bridges + P[i].Pn[index].tot_3strand_bridges)*
                ( (P[i].Pn[index].omega_ij_prime*ux) /P[i].Pn[index].omega_ij);
            
            term1y = (P[i].Pn[index].tot_2strand_bridges + P[i].Pn[index].tot_3strand_bridges)*
                ( (P[i].Pn[index].omega_ij_prime*uy) /P[i].Pn[index].omega_ij);  
            
            term1z = (P[i].Pn[index].tot_2strand_bridges + P[i].Pn[index].tot_3strand_bridges)*
                ( (P[i].Pn[index].omega_ij_prime*uz) /P[i].Pn[index].omega_ij);
                
            sum_term1x=sum_term1x+term1x;
            sum_term1y=sum_term1y+term1y;
            sum_term1z=sum_term1z+term1z;
                
            if((P[i].type==0 && P[index].type==0) || (P[i].type==1 && P[index].type==1))
            {
                term2x = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

                term2y = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

                term2z = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }

            if(P[i].type==0 && P[index].type==1)
            {
                term2x = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

                term2y = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

                term2z = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }
            if(P[i].type==1 && P[index].type==0)
            {
                term2x = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

                term2y = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

                term2z = (P[i].total_free_loops) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                    (P[index].total_free_loops) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }
//              printf("P[i].total_free_loops=%lf term2x=%lf P[i].Pn[index].e_ij_prime=%lf\n",P[i].total_free_loops,term2x,P[i].Pn[index].e_ij_prime);
                //if (term2==0.0)     cout << P[i].Pn[index].e_ij_prime << " " << P[i].Pn[index].e_ij_prime << endl;
                
            sum_term2x = sum_term2x + term2x;
            sum_term2y = sum_term2y + term2y;
            sum_term2z = sum_term2z + term2z;
               
            double omega_0_inert = 4*M_PI*SQ(Rc)*0.75;
            if (P[i].Pn[index].d < 2*Rc+0.75*L) {
                double rr=P[i].Pn[index].d;
                overlap(rr,10.75,0.75,&ov); //inner shpehere and corona  overlap 
                overlap_prime(rr,10.75,0.75,&ov_prime);
                    //term3 = 2 * 500.0 * log(1.0-P[i].Pn[index].e_ij/omega_0);
//                    term3x = beta *2 * 500.0 * ( (P[i].Pn[index].e_ij_prime*ux)/(omega_0-P[i].Pn[index].e_ij) );
//                    term3y = beta *2 * 500.0 * ( (P[i].Pn[index].e_ij_prime*uy)/(omega_0-P[i].Pn[index].e_ij) );
//                    term3z = beta *2 * 500.0 * ( (P[i].Pn[index].e_ij_prime*uz)/(omega_0-P[i].Pn[index].e_ij) );
                    
                term3x = beta * 2* 500.0 * ( (ov_prime*ux)/(omega_0_inert-ov) );
                term3y = beta * 2* 500.0 * ( (ov_prime*uy)/(omega_0_inert-ov) );
                term3z = beta * 2* 500.0 * ( (ov_prime*uz)/(omega_0_inert-ov) );
                sum_term3x = sum_term3x+term3x;
                sum_term3y = sum_term3y+term3y;
                sum_term3z = sum_term3z+term3z;
            }
        }
        P[i].Fx=(sum_term1x-sum_term2x-sum_term3x)/beta;
        P[i].Fy=(sum_term1y-sum_term2y-sum_term3y)/beta;            
        P[i].Fz=(sum_term1z-sum_term2z-sum_term3z)/beta;
    }
    for (q=0; q<ns; q++)
    {
        int index_cs=P[q].Sc;
        double dz_surface=P[index_cs].z-Sz;
        double uz=dz_surface/sqrt(dz_surface*dz_surface);
            
        term1z=  (P[index_cs].tot_2strand_bridges_surface + P[index_cs].tot_3strand_bridges_surface) * ( (P[index_cs].omega_is_prime*uz) / P[index_cs].omega_is );

        term2z=(P[index_cs].total_free_loops) *((P[index_cs].e_is_prime*uz)/P[index_cs].omega_i)+
            nfS_C*((P[index_cs].e_si_prime*uz)/omega_S);

        double omega_0_inert = 4*M_PI*SQ(Rc)*0.75*L;
        if (P[index_cs].Sd < Rc+0.75*L) {
            double rr=P[index_cs].Sd;
            overlap_surface(rr,Rc+0.75*L+0.0,Rc+0.75*L-0.0,0.0,&ov); // it gives e_is not e_si
            overlap_surface_prime(rr,Rc+0.75*L+0.0,Rc+0.75*R2s-0.0,0.0,&ov_prime);
                   
            term3z = beta * 2 * 500.0 * ( (ov_prime*uz)/(omega_0_inert-ov) );
        }
        else {term3x=0.0; term3y=0.0; term3z=0.0;}
        //P[index_cs].Fx=P[index_cs].Fx + term1x + term2x + term3x;
        //P[index_cs].Fy=P[index_cs].Fy + term1y + term2y + term3y;
        P[index_cs].Fz=P[index_cs].Fz + (term1z - term2z - term3z)/beta;
        //P[index_cs].Fz=P[index_cs].Fz + (term1z);
    }
    /****Repulsion from the upper plate*****/
    for (q=0; q<ns_upper; q++)
    {
        int index_cs=P[q].Sc_upper;
        double omega_0_inert = 4*M_PI*SQ(Rc)*0.75;
        double dz_surface=P[index_cs].z-Sz_upper;
        double uz=dz_surface/sqrt(dz_surface*dz_surface);

        double rr=Lz-P[index_cs].z;

        overlap_surface(rr,Rc+0.75*L+0.0,Rc+0.75*L-0.0,0.0,&ov); // it gives e_is not e_si
        overlap_surface_prime(rr,Rc+0.75*L+0.0,Rc+0.75*L-0.0,0.0,&ov_prime);
        term3z = beta * 2 * 500.0 * ( (ov_prime*uz)/(omega_0_inert-ov) );
        P[index_cs].Fz=P[index_cs].Fz + ( (- term3z)/beta);
    }
}

double Isolated_colloids(int NA, int NB, double DG0_L_star,double beta,double rho0, double omega_0)
{
    double nf_A_bulk, nf_B_bulk,PF_L,DG_ii_L,betaDG_ii_L,nf_A_bulk_new,nf_B_bulk_new,ql;
    double error=0.000001;
    bool flag;
    
    nf_A_bulk = rand()% NA+1;
    nf_B_bulk = rand()% NB+1;
                
    //cout << nf_A_bulk << " " << nf_B_bulk << endl;
    PF_L = 1/(omega_0*rho0);
                
    DG_ii_L = DG0_L_star/beta-log(PF_L)/beta;
                
    betaDG_ii_L=beta*DG_ii_L;
    ql=exp(-betaDG_ii_L);
    //ql=20.0;
    flag=true;
    while (flag==true)
    {            
        nf_A_bulk_new=NA/(1 + nf_B_bulk*ql);
        nf_B_bulk_new=NB/(1 + nf_A_bulk*ql);
        double diff_A=fabs(nf_A_bulk_new-nf_A_bulk);
        double diff_B=fabs(nf_B_bulk_new-nf_B_bulk);
                    
        //cout << diff_A << " " << diff_B << endl; 
                    
        if (diff_A<error && diff_B<error) flag=false;
        nf_A_bulk=nf_A_bulk_new;
        nf_B_bulk=nf_B_bulk_new;
    }
           
    //double betaFreeEnergy_bulk=NA*log((nf_A_bulk*nf_B_bulk)/(NA*NB))+(NA-nf_A_bulk);
    double FreeEnergy_bulk=(NA*log((nf_A_bulk*nf_B_bulk)/(NA*NB))+(NA-nf_A_bulk));
    cout << nf_A_bulk << " " << nf_B_bulk << " " << nf_A_bulk*nf_B_bulk*ql << endl;
    return (FreeEnergy_bulk);
    
    
}



void initialisation_after_removal(int Np, int ns, int ns_upper)
{
    int i,j;
    for (i=0; i<Np; i++)
    {

        //P[i].maxnn=0;
        //P[i].maxn=0;
        //P[i].nf=N;
        P[i].omega_i=0;
        P[i].omega_i_inert=0.0;
        P[i].count_regularization=0;
        //P[i].F=0.0;
//         P[i].Pn[i].aon_L=0.0;
//         P[i].Pn[i].pon_L=0.0;
//         P[i].Pn[i].aoff_L=0.0;
//         P[i].Pn[i].poff_L=0.0;
        P[i].Pn[i].aon_L_a1a2=0.0;
        P[i].Pn[i].aoff_L_a1a2=0.0;
        P[i].Pn[i].aon_L_b1b2=0.0;
        P[i].Pn[i].aoff_L_b1b2=0.0;

        P[i].A_on_a1a2_sum= 0.0;
        P[i].A_off_a1a2_sum= 0.0;
        P[i].A_on_a1b1_sum= 0.0;
        P[i].A_off_a1b1_sum= 0.0;
        P[i].A_on_a2b2_sum= 0.0;
        P[i].A_off_a2b2_sum= 0.0;
        P[i].A_on_b1b2_sum= 0.0;
        P[i].A_off_b1b2_sum= 0.0;
        P[i].A_on_a1a2_b1_sum= 0.0;
        P[i].A_off_a1a2_b1_sum= 0.0;
        P[i].A_on_a1_b1_a2_sum= 0.0;
        P[i].A_off_a1_b1_a2_sum= 0.0;
        P[i].A_on_a1a2_b2_sum= 0.0;
        P[i].A_off_a1a2_b2_sum= 0.0;
        P[i].A_on_a2_b2_a1_sum= 0.0;
        P[i].A_off_a2_b2_a1_sum= 0.0;
        P[i].A_on_b1b2_a1_sum= 0.0;
        P[i].A_off_b1b2_a1_sum= 0.0;
        P[i].A_on_b1_a1_b2_sum= 0.0;
        P[i].A_off_b1_a1_b2_sum= 0.0;
        P[i].A_on_b1b2_a2_sum= 0.0;
        P[i].A_off_b1b2_a2_sum= 0.0;
        P[i].A_on_b2_a2_b1_sum= 0.0;
        P[i].A_off_b2_a2_b1_sum= 0.0;
        P[i].A_on_sum_L_a1a2= 0.0;
        P[i].A_off_sum_L_a1a2= 0.0;
        P[i].A_on_sum_L_b1b2= 0.0;
        P[i].A_off_sum_L_b1b2= 0.0;

        P[i].A_on_a2a1_sum= 0.0;
        P[i].A_off_a2a1_sum= 0.0;
        P[i].A_on_b2b1_sum= 0.0;
        P[i].A_off_b2b1_sum= 0.0;


//         P[i].A_on_a1C_sum_S= 0.0;
//         P[i].A_off_a1C_sum_S= 0.0;
//         P[i].A_on_b2C_sum_S= 0.0;
//         P[i].A_off_b2C_sum_S= 0.0;
//         P[i].A_on_a1a2C_sum_S= 0.0;
//         P[i].A_off_a1a2C_sum_S= 0.0;
//         P[i].A_on_a1_C_a2_sum_S= 0.0;
//         P[i].A_off_a1_C_a2_sum_S= 0.0;
        P[i].A_tot_PP_2strand= 0.0;
        P[i].A_tot_PP_3strand= 0.0;
        P[i].A_tot_loop= 0.0;
        P[i].A_tot_PP= 0.0;
//         P[i].A_tot_sur_2strand= 0.0;
//         P[i].A_tot_sur_3strand= 0.0;
//         P[i].A_tot_PS= 0.0;
//         P[i].A_tot_PPS= 0.0;

        P[i].P_on_a1a2_sum= 0.0;
        P[i].P_off_a1a2_sum= 0.0;
        P[i].P_on_a1b1_sum= 0.0;
        P[i].P_off_a1b1_sum= 0.0;
        P[i].P_on_a2b2_sum= 0.0;
        P[i].P_off_a2b2_sum= 0.0;
        P[i].P_on_b1b2_sum= 0.0;
        P[i].P_off_b1b2_sum= 0.0;
        P[i].P_on_a1a2_b1_sum= 0.0;
        P[i].P_off_a1a2_b1_sum= 0.0;
        P[i].P_on_a1_b1_a2_sum= 0.0;
        P[i].P_off_a1_b1_a2_sum= 0.0;
        P[i].P_on_a1a2_b2_sum= 0.0;
        P[i].P_off_a1a2_b2_sum= 0.0;
        P[i].P_on_a2_b2_a1_sum= 0.0;
        P[i].P_off_a2_b2_a1_sum= 0.0;
        P[i].P_on_b1b2_a1_sum= 0.0;
        P[i].P_off_b1b2_a1_sum= 0.0;
        P[i].P_on_b1_a1_b2_sum= 0.0;
        P[i].P_off_b1_a1_b2_sum= 0.0;
        P[i].P_on_b1b2_a2_sum= 0.0;
        P[i].P_off_b1b2_a2_sum= 0.0;
        P[i].P_on_b2_a2_b1_sum= 0.0;
        P[i].P_off_b2_a2_b1_sum= 0.0;
        P[i].P_on_sum_L_a1a2= 0.0;
        P[i].P_off_sum_L_a1a2= 0.0;
        P[i].P_on_sum_L_b1b2= 0.0;
        P[i].P_off_sum_L_b1b2= 0.0;

        P[i].P_on_a2a1_sum= 0.0;
        P[i].P_off_a2a1_sum= 0.0;
        P[i].P_on_b2b1_sum= 0.0;
        P[i].P_off_b2b1_sum= 0.0;

        P[i].P_tot_PP_2strand= 0.0;
        P[i].P_tot_PP_3strand= 0.0;
        P[i].P_tot_loop= 0.0;




//         P[i].Pn[i].nl_a1a2 = 0.0;
//         P[i].Pn[i].nl_b1b2 = 0.0;
//         P[i].Pn[i].nl = 0.0;

        P[i].total_free_loops=0.0;
        P[i].total_bridge_WC=0.0;

        for (j=0; j<Np; j++)
        {
            //P[i].Pn[j].nb=0;
            //P[j].Pn[i].nb=0;
//             P[i].Pn[j].aoff_AB=0.0;
//             P[i].Pn[j].aoff_BA=0.0;
//             P[i].Pn[j].aon_AB=0.0;
//             P[i].Pn[j].aon_BA=0.0;
            P[i].Pn[j].aon_a1a2=0.0;
            P[i].Pn[j].aoff_a1a2=0.0;
            P[i].Pn[j].aon_a1b1=0.0;
            P[i].Pn[j].aoff_a1b1=0.0;
            P[i].Pn[j].aon_a2b2=0.0;
            P[i].Pn[j].aoff_a2b2=0.0;
            P[i].Pn[j].aon_b1b2=0.0;
            P[i].Pn[j].aoff_b1b2=0.0;

            P[i].Pn[j].aon_a2a1=0.0;
            P[i].Pn[j].aoff_a2a1=0.0;
            P[i].Pn[j].aon_b2b1=0.0;
            P[i].Pn[j].aoff_b2b1=0.0;

            P[i].Pn[j].aon_a1a2_b1=0.0;
            P[i].Pn[j].aoff_a1a2_b1=0.0;
            P[i].Pn[j].aon_a1_b1_a2=0.0;
            P[i].Pn[j].aoff_a1_b1_a2=0.0;
            P[i].Pn[j].aon_a1a2_b2=0.0;
            P[i].Pn[j].aoff_a1a2_b2=0.0;
            P[i].Pn[j].aon_a2_b2_a1=0.0;
            P[i].Pn[j].aoff_a2_b2_a1=0.0;
            P[i].Pn[j].aon_b1b2_a1=0.0;
            P[i].Pn[j].aoff_b1b2_a1=0.0;
            P[i].Pn[j].aon_b1_a1_b2=0.0;
            P[i].Pn[j].aoff_b1_a1_b2=0.0;
            P[i].Pn[j].aon_b1b2_a2=0.0;
            P[i].Pn[j].aoff_b1b2_a2=0.0;
            P[i].Pn[j].aon_b2_a2_b1=0.0;
            P[i].Pn[j].aoff_b2_a2_b1=0.0;

            P[i].Pn[j].pon_a1a2=0.0;
            P[i].Pn[j].poff_a1a2=0.0;
            P[i].Pn[j].pon_a1b1=0.0;
            P[i].Pn[j].poff_a1b1=0.0;
            P[i].Pn[j].pon_a2b2=0.0;
            P[i].Pn[j].poff_a2b2=0.0;
            P[i].Pn[j].pon_b1b2=0.0;
            P[i].Pn[j].poff_b1b2=0.0;

            P[i].Pn[j].pon_a2a1=0.0;
            P[i].Pn[j].poff_a2a1=0.0;
            P[i].Pn[j].pon_b2b1=0.0;
            P[i].Pn[j].poff_b2b1=0.0;

            P[i].Pn[j].pon_a1a2_b1=0.0;
            P[i].Pn[j].poff_a1a2_b1=0.0;
            P[i].Pn[j].pon_a1_b1_a2=0.0;
            P[i].Pn[j].poff_a1_b1_a2=0.0;
            P[i].Pn[j].pon_a1a2_b2=0.0;
            P[i].Pn[j].poff_a1a2_b2=0.0;
            P[i].Pn[j].pon_a2_b2_a1=0.0;
            P[i].Pn[j].poff_a2_b2_a1=0.0;
            P[i].Pn[j].pon_b1b2_a1=0.0;
            P[i].Pn[j].poff_b1b2_a1=0.0;
            P[i].Pn[j].pon_b1_a1_b2=0.0;
            P[i].Pn[j].poff_b1_a1_b2=0.0;
            P[i].Pn[j].pon_b1b2_a2=0.0;
            P[i].Pn[j].poff_b1b2_a2=0.0;
            P[i].Pn[j].pon_b2_a2_b1=0.0;
            P[i].Pn[j].poff_b2_a2_b1=0.0;

            P[i].Pn[j].tot_2strand_bridges=0.0;
            P[i].Pn[j].tot_3strand_bridges=0.0;

                
            P[i].Pn[j].e_ij=0.0;
            P[i].Pn[j].e_ij_prime=0.0;
            P[i].Pn[j].omega_ij=0.0;
            P[i].Pn[j].omega_ij_prime=0.0;
            P[i].Pn[j].regu_e_ij=0.0;
            P[i].Pn[j].regu_e_ij_prime=0.0;
            //P[i].Pn[j].nn_index=-1;
            //P[i].Pn[j].n_index=-1;
//             P[i].Pn[j].nn_index=0.0;
//             P[i].Pn[j].n_index=0.0;
        }
    }
    
    
    for (int q=0; q<ns; q++)
    {
        int index_cs=P[q].Sc;
//         P[i].aon_AC_S=0.0;
//         P[i].aoff_AC_S=0.0;
//         P[i].pon_AC_S=0.0;
//         P[i].poff_AC_S=0.0;

//         P[index_cs].nbS_a1C = 0.0;
//         P[index_cs].nbS_a1a2C = 0.0;
//         P[index_cs].nbS_b2C = 0.0;

        P[index_cs].aon_a1C_S=0.0;
        P[index_cs].aoff_a1C_S=0.0;
        P[index_cs].aon_b2C_S=0.0;
        P[index_cs].aoff_b2C_S=0.0;

        P[index_cs].aon_a1a2_C_S=0.0;
        P[index_cs].aoff_a1a2_C_S=0.0;
        P[index_cs].aon_a1_C_a2_S=0.0;
        P[index_cs].aoff_a1_C_a2_S=0.0;

        P[index_cs].A_on_a1C_sum_S=0.0;
        P[index_cs].A_off_a1C_sum_S=0.0;
        P[index_cs].A_on_b2C_sum_S=0.0;
        P[index_cs].A_off_b2C_sum_S=0.0;
        P[index_cs].A_on_a1a2C_sum_S=0.0;
        P[index_cs].A_off_a1a2C_sum_S=0.0;
        P[index_cs].A_on_a1_C_a2_sum_S=0.0;
        P[index_cs].A_off_a1_C_a2_sum_S=0.0;

        P[index_cs].A_tot_sur_2strand=0.0;
        P[index_cs].A_tot_sur_3strand=0.0;

        P[index_cs].A_tot_PS=0.0;
        P[index_cs].A_tot_PPS=0.0;

        P[index_cs].P_on_a1C_sum_S= 0.0;
        P[index_cs].P_off_a1C_sum_S= 0.0;
        P[index_cs].P_on_b2C_sum_S= 0.0;
        P[index_cs].P_off_b2C_sum_S= 0.0;
        P[index_cs].P_on_a1a2C_sum_S= 0.0;
        P[index_cs].P_off_a1a2C_sum_S= 0.0;
        P[index_cs].P_on_a1_C_a2_sum_S= 0.0;
        P[index_cs].P_off_a1_C_a2_sum_S= 0.0;

        P[index_cs].P_tot_sur_2strand= 0.0;
        P[index_cs].P_tot_sur_3strand= 0.0;
        
        P[index_cs].Sd=0.0;             //Distance between particles and surface
        P[index_cs].e_is=0.0;
        P[index_cs].e_is_prime=0.0;
        P[index_cs].e_si=0.0;
        P[index_cs].e_si_prime=0.0;
        P[index_cs].omega_is=0.0;
        P[index_cs].omega_is_prime=0.0;
        P[q].Sc=-1;
//         P[q].Sc=0;
    }
    for (int q_upper=0; q_upper<ns_upper; q_upper++)
    {
        int index_cs_upper=P[q_upper].Sc_upper;
        P[index_cs_upper].Sd_upper=0.0;
        P[q_upper].Sc_upper=-1;
    }
}

    

void initialisation(int Np, int ns, int ns_upper)
{
    int i,j;
    for (i=0; i<Np; i++)
    {

        P[i].maxnn=0;
        P[i].maxn=0;
        //P[i].nf=N;
        P[i].omega_i=0;
        P[i].omega_i_inert=0.0;
        P[i].count_regularization=0;
        //P[i].F=0.0;
//         P[i].Pn[i].aon_L=0.0;
//         P[i].Pn[i].pon_L=0.0;
//         P[i].Pn[i].aoff_L=0.0;
//         P[i].Pn[i].poff_L=0.0;
        P[i].Pn[i].aon_L_a1a2=0.0;
        P[i].Pn[i].aoff_L_a1a2=0.0;
        P[i].Pn[i].aon_L_b1b2=0.0;
        P[i].Pn[i].aoff_L_b1b2=0.0;

        P[i].A_on_a1a2_sum= 0.0;
        P[i].A_off_a1a2_sum= 0.0;
        P[i].A_on_a1b1_sum= 0.0;
        P[i].A_off_a1b1_sum= 0.0;
        P[i].A_on_a2b2_sum= 0.0;
        P[i].A_off_a2b2_sum= 0.0;
        P[i].A_on_b1b2_sum= 0.0;
        P[i].A_off_b1b2_sum= 0.0;
        P[i].A_on_a1a2_b1_sum= 0.0;
        P[i].A_off_a1a2_b1_sum= 0.0;
        P[i].A_on_a1_b1_a2_sum= 0.0;
        P[i].A_off_a1_b1_a2_sum= 0.0;
        P[i].A_on_a1a2_b2_sum= 0.0;
        P[i].A_off_a1a2_b2_sum= 0.0;
        P[i].A_on_a2_b2_a1_sum= 0.0;
        P[i].A_off_a2_b2_a1_sum= 0.0;
        P[i].A_on_b1b2_a1_sum= 0.0;
        P[i].A_off_b1b2_a1_sum= 0.0;
        P[i].A_on_b1_a1_b2_sum= 0.0;
        P[i].A_off_b1_a1_b2_sum= 0.0;
        P[i].A_on_b1b2_a2_sum= 0.0;
        P[i].A_off_b1b2_a2_sum= 0.0;
        P[i].A_on_b2_a2_b1_sum= 0.0;
        P[i].A_off_b2_a2_b1_sum= 0.0;
        P[i].A_on_sum_L_a1a2= 0.0;
        P[i].A_off_sum_L_a1a2= 0.0;
        P[i].A_on_sum_L_b1b2= 0.0;
        P[i].A_off_sum_L_b1b2= 0.0;

        P[i].A_on_a2a1_sum= 0.0;
        P[i].A_off_a2a1_sum= 0.0;
        P[i].A_on_b2b1_sum= 0.0;
        P[i].A_off_b2b1_sum= 0.0;


//         P[i].A_on_a1C_sum_S= 0.0;
//         P[i].A_off_a1C_sum_S= 0.0;
//         P[i].A_on_b2C_sum_S= 0.0;
//         P[i].A_off_b2C_sum_S= 0.0;
//         P[i].A_on_a1a2C_sum_S= 0.0;
//         P[i].A_off_a1a2C_sum_S= 0.0;
//         P[i].A_on_a1_C_a2_sum_S= 0.0;
//         P[i].A_off_a1_C_a2_sum_S= 0.0;
        P[i].A_tot_PP_2strand= 0.0;
        P[i].A_tot_PP_3strand= 0.0;
        P[i].A_tot_loop= 0.0;
        P[i].A_tot_PP= 0.0;
//         P[i].A_tot_sur_2strand= 0.0;
//         P[i].A_tot_sur_3strand= 0.0;
//         P[i].A_tot_PS= 0.0;
//         P[i].A_tot_PPS= 0.0;

        P[i].P_on_a1a2_sum= 0.0;
        P[i].P_off_a1a2_sum= 0.0;
        P[i].P_on_a1b1_sum= 0.0;
        P[i].P_off_a1b1_sum= 0.0;
        P[i].P_on_a2b2_sum= 0.0;
        P[i].P_off_a2b2_sum= 0.0;
        P[i].P_on_b1b2_sum= 0.0;
        P[i].P_off_b1b2_sum= 0.0;
        P[i].P_on_a1a2_b1_sum= 0.0;
        P[i].P_off_a1a2_b1_sum= 0.0;
        P[i].P_on_a1_b1_a2_sum= 0.0;
        P[i].P_off_a1_b1_a2_sum= 0.0;
        P[i].P_on_a1a2_b2_sum= 0.0;
        P[i].P_off_a1a2_b2_sum= 0.0;
        P[i].P_on_a2_b2_a1_sum= 0.0;
        P[i].P_off_a2_b2_a1_sum= 0.0;
        P[i].P_on_b1b2_a1_sum= 0.0;
        P[i].P_off_b1b2_a1_sum= 0.0;
        P[i].P_on_b1_a1_b2_sum= 0.0;
        P[i].P_off_b1_a1_b2_sum= 0.0;
        P[i].P_on_b1b2_a2_sum= 0.0;
        P[i].P_off_b1b2_a2_sum= 0.0;
        P[i].P_on_b2_a2_b1_sum= 0.0;
        P[i].P_off_b2_a2_b1_sum= 0.0;
        P[i].P_on_sum_L_a1a2= 0.0;
        P[i].P_off_sum_L_a1a2= 0.0;
        P[i].P_on_sum_L_b1b2= 0.0;
        P[i].P_off_sum_L_b1b2= 0.0;

        P[i].P_on_a2a1_sum= 0.0;
        P[i].P_off_a2a1_sum= 0.0;
        P[i].P_on_b2b1_sum= 0.0;
        P[i].P_off_b2b1_sum= 0.0;

        P[i].P_tot_PP_2strand= 0.0;
        P[i].P_tot_PP_3strand= 0.0;
        P[i].P_tot_loop= 0.0;




//         P[i].Pn[i].nl_a1a2 = 0.0;
//         P[i].Pn[i].nl_b1b2 = 0.0;
//         P[i].Pn[i].nl = 0.0;

        P[i].total_free_loops=0.0;
        P[i].total_bridge_WC=0.0;

        for (j=0; j<Np; j++)
        {
            //P[i].Pn[j].nb=0;
            //P[j].Pn[i].nb=0;
//             P[i].Pn[j].aoff_AB=0.0;
//             P[i].Pn[j].aoff_BA=0.0;
//             P[i].Pn[j].aon_AB=0.0;
//             P[i].Pn[j].aon_BA=0.0;
            P[i].Pn[j].aon_a1a2=0.0;
            P[i].Pn[j].aoff_a1a2=0.0;
            P[i].Pn[j].aon_a1b1=0.0;
            P[i].Pn[j].aoff_a1b1=0.0;
            P[i].Pn[j].aon_a2b2=0.0;
            P[i].Pn[j].aoff_a2b2=0.0;
            P[i].Pn[j].aon_b1b2=0.0;
            P[i].Pn[j].aoff_b1b2=0.0;

            P[i].Pn[j].aon_a2a1=0.0;
            P[i].Pn[j].aoff_a2a1=0.0;
            P[i].Pn[j].aon_b2b1=0.0;
            P[i].Pn[j].aoff_b2b1=0.0;

            P[i].Pn[j].aon_a1a2_b1=0.0;
            P[i].Pn[j].aoff_a1a2_b1=0.0;
            P[i].Pn[j].aon_a1_b1_a2=0.0;
            P[i].Pn[j].aoff_a1_b1_a2=0.0;
            P[i].Pn[j].aon_a1a2_b2=0.0;
            P[i].Pn[j].aoff_a1a2_b2=0.0;
            P[i].Pn[j].aon_a2_b2_a1=0.0;
            P[i].Pn[j].aoff_a2_b2_a1=0.0;
            P[i].Pn[j].aon_b1b2_a1=0.0;
            P[i].Pn[j].aoff_b1b2_a1=0.0;
            P[i].Pn[j].aon_b1_a1_b2=0.0;
            P[i].Pn[j].aoff_b1_a1_b2=0.0;
            P[i].Pn[j].aon_b1b2_a2=0.0;
            P[i].Pn[j].aoff_b1b2_a2=0.0;
            P[i].Pn[j].aon_b2_a2_b1=0.0;
            P[i].Pn[j].aoff_b2_a2_b1=0.0;

            P[i].Pn[j].pon_a1a2=0.0;
            P[i].Pn[j].poff_a1a2=0.0;
            P[i].Pn[j].pon_a1b1=0.0;
            P[i].Pn[j].poff_a1b1=0.0;
            P[i].Pn[j].pon_a2b2=0.0;
            P[i].Pn[j].poff_a2b2=0.0;
            P[i].Pn[j].pon_b1b2=0.0;
            P[i].Pn[j].poff_b1b2=0.0;

            P[i].Pn[j].pon_a2a1=0.0;
            P[i].Pn[j].poff_a2a1=0.0;
            P[i].Pn[j].pon_b2b1=0.0;
            P[i].Pn[j].poff_b2b1=0.0;

            P[i].Pn[j].pon_a1a2_b1=0.0;
            P[i].Pn[j].poff_a1a2_b1=0.0;
            P[i].Pn[j].pon_a1_b1_a2=0.0;
            P[i].Pn[j].poff_a1_b1_a2=0.0;
            P[i].Pn[j].pon_a1a2_b2=0.0;
            P[i].Pn[j].poff_a1a2_b2=0.0;
            P[i].Pn[j].pon_a2_b2_a1=0.0;
            P[i].Pn[j].poff_a2_b2_a1=0.0;
            P[i].Pn[j].pon_b1b2_a1=0.0;
            P[i].Pn[j].poff_b1b2_a1=0.0;
            P[i].Pn[j].pon_b1_a1_b2=0.0;
            P[i].Pn[j].poff_b1_a1_b2=0.0;
            P[i].Pn[j].pon_b1b2_a2=0.0;
            P[i].Pn[j].poff_b1b2_a2=0.0;
            P[i].Pn[j].pon_b2_a2_b1=0.0;
            P[i].Pn[j].poff_b2_a2_b1=0.0;

            P[i].Pn[j].tot_2strand_bridges=0.0;
            P[i].Pn[j].tot_3strand_bridges=0.0;

//             P[i].Pn[j].nb_a1a2 = 0.0;
//             P[i].Pn[j].nb_a1b1 = 0.0;
//             P[i].Pn[j].nb_a2b2 = 0.0;
//             P[i].Pn[j].nb_b1b2 = 0.0;
//
//             P[i].Pn[j].nb_a1a2b1 = 0.0;
//             P[i].Pn[j].nb_a1a2b2 = 0.0;
//             P[i].Pn[j].nb_b1b2a1 = 0.0;
//             P[i].Pn[j].nb_b1b2a2 = 0.0;
//
//             ////
//             P[j].Pn[i].nb_a1a2 = 0.0;
//             P[j].Pn[i].nb_a1b1 = 0.0;
//             P[j].Pn[i].nb_a2b2 = 0.0;
//             P[j].Pn[i].nb_b1b2 = 0.0;
//
//             P[j].Pn[i].nb_a1a2b1 = 0.0;
//             P[j].Pn[i].nb_a1a2b2 = 0.0;
//             P[j].Pn[i].nb_b1b2a1 = 0.0;
//             P[j].Pn[i].nb_b1b2a2 = 0.0;

                
//             P[i].Pn[j].poff_AB=0.0;
//             P[i].Pn[j].poff_BA=0.0;
//             P[i].Pn[j].pon_AB=0.0;
//             P[i].Pn[j].pon_BA=0.0;
                
            P[i].Pn[j].e_ij=0.0;
            P[i].Pn[j].e_ij_prime=0.0;
            P[i].Pn[j].omega_ij=0.0;
            P[i].Pn[j].omega_ij_prime=0.0;
            P[i].Pn[j].regu_e_ij=0.0;
            P[i].Pn[j].regu_e_ij_prime=0.0;
            P[i].Pn[j].nn_index=-1;
            P[i].Pn[j].n_index=-1;
//             P[i].Pn[j].nn_index=0.0;
//             P[i].Pn[j].n_index=0.0;
        }
    }
    
    
    for (int q=0; q<ns; q++)
    {
        int index_cs=P[q].Sc;
//         P[i].aon_AC_S=0.0;
//         P[i].aoff_AC_S=0.0;
//         P[i].pon_AC_S=0.0;
//         P[i].poff_AC_S=0.0;

//         P[index_cs].nbS_a1C = 0.0;
//         P[index_cs].nbS_a1a2C = 0.0;
//         P[index_cs].nbS_b2C = 0.0;

        P[index_cs].aon_a1C_S=0.0;
        P[index_cs].aoff_a1C_S=0.0;
        P[index_cs].aon_b2C_S=0.0;
        P[index_cs].aoff_b2C_S=0.0;

        P[index_cs].aon_a1a2_C_S=0.0;
        P[index_cs].aoff_a1a2_C_S=0.0;
        P[index_cs].aon_a1_C_a2_S=0.0;
        P[index_cs].aoff_a1_C_a2_S=0.0;

        P[index_cs].A_on_a1C_sum_S=0.0;
        P[index_cs].A_off_a1C_sum_S=0.0;
        P[index_cs].A_on_b2C_sum_S=0.0;
        P[index_cs].A_off_b2C_sum_S=0.0;
        P[index_cs].A_on_a1a2C_sum_S=0.0;
        P[index_cs].A_off_a1a2C_sum_S=0.0;
        P[index_cs].A_on_a1_C_a2_sum_S=0.0;
        P[index_cs].A_off_a1_C_a2_sum_S=0.0;

        P[index_cs].A_tot_sur_2strand=0.0;
        P[index_cs].A_tot_sur_3strand=0.0;

        P[index_cs].A_tot_PS=0.0;
        P[index_cs].A_tot_PPS=0.0;

        P[index_cs].P_on_a1C_sum_S= 0.0;
        P[index_cs].P_off_a1C_sum_S= 0.0;
        P[index_cs].P_on_b2C_sum_S= 0.0;
        P[index_cs].P_off_b2C_sum_S= 0.0;
        P[index_cs].P_on_a1a2C_sum_S= 0.0;
        P[index_cs].P_off_a1a2C_sum_S= 0.0;
        P[index_cs].P_on_a1_C_a2_sum_S= 0.0;
        P[index_cs].P_off_a1_C_a2_sum_S= 0.0;

        P[index_cs].P_tot_sur_2strand= 0.0;
        P[index_cs].P_tot_sur_3strand= 0.0;
        
        P[index_cs].Sd=0.0;             //Distance between particles and surface
        P[index_cs].e_is=0.0;
        P[index_cs].e_is_prime=0.0;
        P[index_cs].e_si=0.0;
        P[index_cs].e_si_prime=0.0;
        P[index_cs].omega_is=0.0;
        P[index_cs].omega_is_prime=0.0;
        P[q].Sc=-1;
//         P[q].Sc=0;
    }
    for (int q_upper=0; q_upper<ns_upper; q_upper++)
    {
        int index_cs_upper=P[q_upper].Sc_upper;
        P[index_cs_upper].Sd_upper=0.0;
        P[q_upper].Sc_upper=-1;
    }
}

//void initialisation_cpy(int Np, int ns, int ns_upper)
//{
//    int i,j;
//    for (i=0; i<Np; i++)
//    {
//        P_cpy[i].maxnn=0;
//        P[i].maxn=0;
//        //P[i].nf=N;
//        P_cpy[i].omega_i=0.0;
//        P_cpy[i].omega_i_inert=0.0;
//        //P[i].Fx=0.0;
//        //P[i].Fy=0.0;
//        //P[i].Fz=0.0;
//        P_cpy[i].Pn[i].nl_d=0.0;
//        P_cpy[i].count_regularization=0;
//        P_cpy[i].nf_A_old=0.0;
//        P_cpy[i].nf_B_old=0.0;
//        for (j=0; j<Np; j++)
//        {
//            P_cpy[i].Pn[j].e_ij=0.0;
//            P_cpy[i].Pn[j].e_ij_prime=0.0;
//            P_cpy[i].Pn[j].omega_ij=0.0;
//            P_cpy[i].Pn[j].omega_ij_prime=0.0;
//            P_cpy[i].Pn[j].regu_e_ij=0.0;
//            P_cpy[i].Pn[j].regu_e_ij_prime=0.0;
//            P_cpy[i].Pn[j].nn_index=-1;
//            //P[i].Pn[j].n_index=-1;
//            P_cpy[i].Pn[j].nb_d=0;
//            P_cpy[i].Pn[j].d=0.0;
//        }  
//    }
//    for (int q=0; q<ns; q++)
//    {
//        int index_cs=P_cpy[q].Sc;
//        
//        P_cpy[index_cs].nbS_d=0.0;
//        P[index_cs].Sd=0.0;
//        P_cpy[index_cs].e_is=0.0;
//        P_cpy[index_cs].e_is_prime=0.0;
//        P_cpy[index_cs].e_si=0.0;
//        P_cpy[index_cs].e_si_prime=0.0;
//        P_cpy[index_cs].omega_is=0.0;
//        P_cpy[index_cs].omega_is_prime=0.0;
//        P_cpy[q].Sc=-1;
//    }
//    for (int q_upper=0; q_upper<ns_upper; q_upper++)
//    {
//        int index_cs_upper=P_cpy[q_upper].Sc_upper;
//        
//        P_cpy[index_cs_upper].Sd_upper=0.0;
//        P_cpy[q_upper].Sc_upper=-1;
//    }
//}

