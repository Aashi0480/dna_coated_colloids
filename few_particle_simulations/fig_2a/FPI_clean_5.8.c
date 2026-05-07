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
# define maxN 50
# define SQ(x) ((x)*(x))
# define CUB(x) ((x)*(x)*(x))
using namespace std;
void distance_ads(double,double,double,int,double,double,double,double *,double *,double *,double *);
void distance(int,int,double,double,double,double *,double *,double *,double *);
double Isolated_colloids(int Na1, int Na2, double DG0_L_star,double beta, double rho0, double omega_0);
void neighborlist(int Np,int *ns,int *ns_upper,double Lx,double Ly,double Lz,double Sz, double Sz_upper, double Rc,double L,double R2s);
void overlapvolume(int,double,double,double,double);
void overlapvolume_surface(int ,double ,double ,double ,double , double *);
void overlap(double,double,double,double *);
void overlap_prime(double,double,double,double *);
void overlap_surface(double,double,double,double,double *);
void overlap_surface_prime(double,double,double,double,double *);
void normal(double *);
void fixpointiteration(int, int, int, int, int, int, int, double, double, double, double, double, double, double, double, double,double*,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double,double);
void Force(int Np, int ns, int ns_upper, double Rc, double L, double R2s,double Sz, double Sz_upper, double omega_S, double nfS_C_old, double beta,double Lx, double Ly, double Lz,double omega_0_inert);

struct Rates
{
	double cum;
	int IPi;			 //Maximum number of free linkers
	int IPj;		     //Maximum number of neighbors
	int type;
};
Rates R[maxNp];
struct Pairs
{
	int nn_index; 			              //index of neighbors
    int nb;
	int nb_AB,nb_BA;				    //number of bridges
	int n_index;
    int nl;
    double d;
    double e_ij;
    double e_ij_prime;
    double omega_ij;
    double omega_ij_prime;
    double regu_e_ij,regu_e_ij_prime;
    double nl_d,nb_d,nb_d_01,nb_d_10,nb_d_00,nb_d_11;
    double nb_d3S,nb_d3S_a1,nb_d3S_a2,nb_d3S_b1,nb_d3S_b2,nb_d3S_A,nb_d3S_B;
    double nl_d_A,nl_d_B; //no. of loops in A & B colloids separately
    double nb_d_a1; //no. of a1 involved in 2 strand bridges
    double nb_d_a2; //no. of a2 involved in 2 strand bridges
    double nb_d_b1; //no. of b1 involved in 2 strand bridges
    double nb_d_b2; //no. of b2 involved in 2 strand bridges
    double nb_d_a1_AA; //no. of a1 involved in 2 strand bridges
    double nb_d_a2_AA; //no. of a2 involved in 2 strand bridges
    double nb_d_b1_BB; //no. of b1 involved in 2 strand bridges
    double nb_d_b2_BB; //no. of b2 involved in 2 strand bridges
    int type;
};
Pairs Pn[maxNp];
struct Particles
{
    int nf;
	int nf_a1,nf_a2,nf_b1,nf_b2,maxn,maxnn;			//Maximum number of free linkers, Maximum number of neighbors
    int Sc,Sc_upper;                             //index of the neighbor particles of the surface
    int type;
    int count_regularization;
	double x,y,z;
    double nf_a1_old,nf_a2_old;
    double nf_a1_new,nf_a2_new;
    double nf_b1_old,nf_b2_old;
    double nf_b1_new,nf_b2_new;
    double nf_A_GC_before,nf_B_GC_before;
    double nbS_AC;
    double nbS_d, total_bridge_WC,total_free_loops,total_free_loops_A,total_free_loops_B;
    double nbS_d_A,nbS_d_B; //no. of bridges with surface for A & B colloid separately
    double Sd;                          //Distance between the neighbor and surface
	double omega_i,omega_i_inert;
    double Fx,Fy,Fz;
    double e_is;
    double e_is_prime;
    double e_si;
    double e_si_prime;
    double omega_is;
    double omega_is_prime;
    double lfb; //for three strand complex; first loop in the colloid, then bridge with the surface
    double bfb; //total no. of first bridges in between colloids , then bridge with surface
    Pairs Pn[maxNp];
};
Particles P[maxNp];
Particles P_cpy[maxNp];



int main()
{
    int Np;
    int N;
    int NsC;
    double Rc;
	double L;
    double R2s;
	double Lx;
	double Ly;
	double Lz;
    double beta;
	double Delta_tB;
    double t_tot;
    double rho0;
    double DG0_surface, DG0, DG0_L, DG0ta1a2b1, DG0ta1a2b2, DG0tb1b2a1, DG0tb1b2a2, DG0_2S_A_sur, DG0_2S_B_sur, DG0_3S_sur,DG0_a1b1,DG0_a2b2,DG0_L_a1a2,DG0_L_b1b2,DG0_br_a1a2,DG0_br_b1b2,DGT,n_alpha,n_T;
    double rho_ideal;
    char name1[50];
    FILE *in;
    in=fopen("INPUT_k0k1_blue","r");
    fscanf(in,"%s %d",name1,&Np);
    fscanf(in,"%s %d",name1,&N);
    fscanf(in,"%s %d",name1,&NsC);
    fscanf(in,"%s %lf",name1,&Rc);
    fscanf(in,"%s %lf",name1,&L);
    fscanf(in,"%s %lf",name1,&R2s);
    fscanf(in,"%s %lf",name1,&Lx);
    fscanf(in,"%s %lf",name1,&Ly);
    fscanf(in,"%s %lf",name1,&Lz);
    fscanf(in,"%s %lf",name1,&beta);
    fscanf(in,"%s %lf",name1,&Delta_tB);
    fscanf(in,"%s %lf",name1,&t_tot);
    fscanf(in,"%s %lf",name1,&DG0_surface);
    fscanf(in,"%s %lf",name1,&DG0);
    fscanf(in,"%s %lf",name1,&DGT);
    fscanf(in,"%s %lf",name1,&n_alpha);
    fscanf(in,"%s %lf",name1,&n_T);
    fscanf(in,"%s %lf",name1,&DG0_br_a1a2);
    fscanf(in,"%s %lf",name1,&DG0_br_b1b2);
    fscanf(in,"%s %lf",name1,&DG0_a1b1);
    fscanf(in,"%s %lf",name1,&DG0_a2b2);
    fscanf(in,"%s %lf",name1,&DG0_L);
    fscanf(in,"%s %lf",name1,&DG0_L_a1a2);
    fscanf(in,"%s %lf",name1,&DG0_L_b1b2);
    fscanf(in,"%s %lf",name1,&DG0_2S_A_sur);
    fscanf(in,"%s %lf",name1,&DG0_2S_B_sur);
    fscanf(in,"%s %lf",name1,&DG0_3S_sur);
    fscanf(in,"%s %lf",name1,&DG0ta1a2b1);
    fscanf(in,"%s %lf",name1,&DG0ta1a2b2);
    fscanf(in,"%s %lf",name1,&DG0tb1b2a1);
    fscanf(in,"%s %lf",name1,&DG0tb1b2a2);
    fscanf(in,"%s %lf",name1,&rho0);
    fscanf(in,"%s %lf",name1,&rho_ideal);

    fclose(in);
    int Na1=N/4;
    int Na2=N/4;
    int Nb1=N/4;
    int Nb2=N/4;
    double V_linkers=4*M_PI*SQ(Rc)*L;
    double Density_linkers=N/V_linkers;


    double V=Lx*Ly*Lz;
    double DG0_surface_star=beta*DG0_surface-log(1/rho0*CUB(L));
    double DG0_star=beta*DG0-log(1/rho0*CUB(L));
    double DG0_L_star=beta*DG0_L-log(1/rho0*CUB(L));

    int nfS_C;
    double nfS_C_old,nfS_C_new;

	double dx,dy,dz,d,X,Y,Z;
	double t,t_bar,tao;

    int n_deposited_particles, n_checked_particles;

	int count,count_nn,count_n;
	int seed =time(NULL)*getpid();

    srand (seed);
    cout << "seed=" << seed << endl;

	int i,j,k,E,q;
	bool flag;
    ofstream myfile;
    ofstream myfile1;
    ofstream myfile2;
    myfile1.open ("k0k1_blue_clean_5.8.dump");

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

    for (i=0; i<Np; i++)
    {
        if(P[i].type==0)
        {
            P[i].nf_a1=Na1; ///////////////////////
            P[i].nf_a2=Na2; //////////////////////
        }
        if(P[i].type==1)
        {
            P[i].nf_b1=Nb1; /////////////////////
            P[i].nf_b2=Nb2; /////////////////////
        }
    }
    nfS_C=NsC;         //number of free surface linker (initiation) C type
    double e_is,e_is_prime,e_si,e_si_prime,w_is,w_is_prime,e_si_sum,omega_S;
    double omega_S_0=Lx*Ly*R2s;
    double omega_0=4*M_PI*SQ(Rc)*L;
    double omega_0_inert=4*M_PI*SQ(Rc)*0.75*L;

    double e=2.71828;
    double nfS_C_cpy, omega_S_cpy, ns_cpy, ns_upper_cpy;
/***********************************Particle position *******************************************************************/

            P[0].x = 30;
            P[0].y = 30;
            P[0].z = 5.8;
            P[0].type=0;




    n_deposited_particles=Np;

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

    double Sx=Lx/2.0;
    double Sy=Ly/2.0;
    double Sz=0.0;
    double Sz_upper=Lz;
    double ov, ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
    double Dcs,dcs;
    double FreeEnergy_bulk;
    FreeEnergy_bulk = Isolated_colloids(Na1, Na2, DG0_L_star,beta,rho0,omega_0);
/**************************** neighbor calculations *********************************************************************/
    int count_BD_steps=0;
    t=0.0;
    int count_test=0;
	while (t < t_tot)  //For Brownian dynamics; for the time being it is commented out
	{
        count_test=count_test+1;
        int ns;
        int ns_upper;
        neighborlist(Np,&ns,&ns_upper,Lx,Ly,Lz,Sz, Sz_upper, Rc,L,R2s);
        overlapvolume(Np,Rc,L,omega_0,omega_0_inert);
        overlapvolume_surface(ns,Rc,L,R2s,omega_S_0,&omega_S);
        fixpointiteration(Np, Na1, Na2, Nb1, Nb2, ns, NsC, beta, rho0, Rc, L, R2s, DG0_star, DG0_L_star, DG0_surface_star, omega_S,&nfS_C_old,DG0_2S_A_sur,DG0_2S_B_sur,DG0_3S_sur, DG0ta1a2b1, DG0ta1a2b2, DG0tb1b2a1, DG0tb1b2a2,DG0,DG0_L,DG0_a1b1,DG0_a2b2,DG0_L_a1a2,DG0_L_b1b2,DG0_br_a1a2,DG0_br_b1b2,DGT,n_alpha,n_T);
        Force(Np,ns,ns_upper,Rc,L,R2s,Sz,Sz_upper,omega_S, nfS_C_old, beta,Lx,Ly,Lz,omega_0_inert);
/**************************************  Brownian dynamics  ****************************/
        int i_new,Np_new;
        double insetion_prob,remove_prob;
        double rand_BD_GC=1;
        if (rand_BD_GC<1/1000.0)
        {
            double Rx,Ry,Rz,Rn;
            double X1_new,X2_new,Y1_new,Y2_new,Z1_new,Z2_new;
            double D;
            for (i=0; i<Np; i++)
            {
                normal(&Rx);
                normal(&Ry);
                normal(&Rz);
                double X_new=P[i].x + beta*P[i].Fx*Delta_tB*L + sqrt(2*Delta_tB)*Rx;
                double Y_new=P[i].y + beta*P[i].Fy*Delta_tB*L + sqrt(2*Delta_tB)*Ry;
                double Z_new=P[i].z + beta*P[i].Fz*Delta_tB*L + sqrt(2*Delta_tB)*Rz;

                if      (X_new <  0)    X_new = X_new + Lx;
                else if (X_new >  Lx)   X_new = X_new - Lx;

                if      (Y_new <  0)    Y_new = Y_new + Ly;
                else if (Y_new >  Ly)   Y_new = Y_new - Ly;



                P[i].x=X_new;
                P[i].y=Y_new;
                P[i].z=Z_new;
            }
        }


/**************************************  Particle insertion or delete  ****************************/
        else if (rand_BD_GC<1/1000.0)
        {
            double Free_energy_old=0.0;
            double Free_energy_regularization=0.0;
            double Free_energy_regularization_surface=0.0;
            double omega_0_inert = 4*M_PI*SQ(Rc)*0.75;
            for (i=0; i<Np; i++)
            {
                if (P[i].z>Rc+L+R2s)    P[i].nbS_d=0.0;
                Free_energy_old = Free_energy_old + ( (Na1*log(P[i].nf_a1_old/(double)Na1) + Na2*log(P[i].nf_a2_old/(double)Na2) + P[i].total_bridge_WC/2.0 + P[i].Pn[i].nl_d + P[i].nbS_d/2.0) - ((Na1+Na2)*log(P[i].omega_i/omega_0)) );
                if (P[i].count_regularization>0)
                {
                   for(j=0; j<P[i].maxnn; j++)
                   {
                       int index=P[i].Pn[j].nn_index;
                       double r=P[i].Pn[index].d;
                       if (r < 2*Rc+0.75*L)
                       {
                           Free_energy_regularization= 500.0 * log(1.0-P[i].Pn[index].regu_e_ij/omega_0_inert); //log(omega_0_inert-P[i].Pn[index].regu_e_ij)
                           Free_energy_old=Free_energy_old-Free_energy_regularization;
                       }
                   }
                }
                if (P[i].z<Rc+0.75*L)
                {
                    double rr=P[i].z;
                    overlap_surface(rr,Rc+0.75*L+0.0,Rc+0.75*R2s-0.0,0.0,&ov); // it gives e_is not e_si
                    Free_energy_regularization_surface= 500.0 * log(1.0-ov/omega_0_inert);
                    Free_energy_old=Free_energy_old-Free_energy_regularization_surface;
                }
            }
            Free_energy_old = Free_energy_old + NsC*log(nfS_C_old/(double)NsC) + (NsC-nfS_C_old)/2.0 - NsC*log(omega_S/omega_S_0);  //Basically Free_energy_old is the betaFree_energy_old=beta*Free_energy_old

            double rand_insertion_delete=rand()/(double)RAND_MAX;

            /************insertion *****************************************************/
            if (rand_insertion_delete <0.5) //insertion
            {
                X =(rand()/(double)RAND_MAX)*Lx;
                Y =(rand()/(double)RAND_MAX)*Ly;
                Z =(rand()/(double)RAND_MAX)*Lz;
                int count_deposition_ne=0;
                for (i=0; i<Np; i++)
                {
                    distance_ads(X,Y,Z,i,Lx,Ly,Lz,&dx,&dy,&dz,&d);

                    if ( ( d>2*Rc) && Z>Rc && Z<(Lz-Rc))    count_deposition_ne=count_deposition_ne+1;
                    else break;
                }

                if (count_deposition_ne==Np)
                {
                    nfS_C_cpy=nfS_C_old; omega_S_cpy=omega_S; ns_cpy=ns; ns_upper_cpy=ns_upper;
                    for (i=0; i<Np; i++)
                    {
                        P_cpy[i]=P[i];
                        P[i].maxnn=0;
                        P[i].maxn=0;
                        P[i].omega_i=0.0;
                        P[i].omega_i_inert=0.0;
                        P[i].Pn[i].nl_d=0.0;

                        for (j=0; j<Np; j++)
                        {
                            P[i].Pn[j].e_ij=0.0;
                            P[i].Pn[j].e_ij_prime=0.0;
                            P[i].Pn[j].omega_ij=0.0;
                            P[i].Pn[j].omega_ij_prime=0.0;
                            P[i].Pn[j].regu_e_ij=0.0;
                            P[i].Pn[j].regu_e_ij_prime=0.0;
                            P[i].Pn[j].nn_index=0;
                            P[i].Pn[j].n_index=0;
                            P[i].Pn[j].nb_d=0;
                        }
                    }
                    for (q=0; q<ns; q++)
                    {
                        int index_cs=P[q].Sc;
                        P[index_cs].nbS_d=0.0;
                        P[index_cs].Sd=0.0;
                        P[index_cs].e_is=0.0;
                        P[index_cs].e_is_prime=0.0;
                        P[index_cs].e_si=0.0;
                        P[index_cs].e_si_prime=0.0;
                        P[index_cs].omega_is=0.0;
                        P[index_cs].omega_is_prime=0.0;
                    }
                    nfS_C_old=NsC; omega_S=omega_S_0; ns=0; ns_upper=0;

                    Np_new=Np+1;
                    i_new=Np_new-1;
                    P[i_new].x=X;P[i_new].y=Y;P[i_new].z=Z;

                    neighborlist(Np_new,&ns,&ns_upper,Lx,Ly,Lz,Sz, Sz_upper, Rc,L,R2s);
                    overlapvolume(Np_new,Rc,L,omega_0,omega_0_inert);
                    overlapvolume_surface(ns,Rc,L,R2s,omega_S_0,&omega_S);
                    fixpointiteration(Np_new, Na1, Na2, Nb1, Nb2, ns, NsC, beta, rho0, Rc, L, R2s, DG0_star, DG0_L_star, DG0_surface_star, omega_S,&nfS_C_old,DG0_2S_A_sur,DG0_2S_B_sur,DG0_3S_sur, DG0ta1a2b1, DG0ta1a2b2, DG0tb1b2a1, DG0tb1b2a2, DG0, DG0_L, DG0_a1b1, DG0_a2b2, DG0_L_a1a2, DG0_L_b1b2, DG0_br_a1a2, DG0_br_b1b2, DGT, n_alpha, n_T);

                    double Free_energy_new=0.0;
                    double Free_energy_regularization=0.0;
                    double omega_0_inert = 4*M_PI*SQ(Rc)*0.75;
                    for (i=0; i<Np_new; i++)  //Np_new
                    {
                        if (P[i].z>Rc+R2s+L)    P[i].nbS_d=0.0;
                        Free_energy_new = Free_energy_new + ( (Na1*log(P[i].nf_a1_old/(double)Na1) + Na2*log(P[i].nf_a2_old/(double)Na2) + P[i].total_bridge_WC/2.0 + P[i].Pn[i].nl_d + P[i].nbS_d/2.0) - ((Na1+Na2)*log(P[i].omega_i/omega_0)) );
                        if (P[i].count_regularization>0)
                        {
                            for(j=0; j<P[i].maxnn; j++)
                            {
                                int index=P[i].Pn[j].nn_index;
                                double r=P[i].Pn[index].d;
                                if (r < 2*Rc+0.75*L)
                                {
                                    Free_energy_regularization= 500.0 * log(1.0-P[i].Pn[index].regu_e_ij/omega_0_inert); //log(omega_0_inert-P[i].Pn[index].regu_e_ij)
                                    Free_energy_new=Free_energy_new-Free_energy_regularization;
                                }
                            }
                        }
                        if (P[i].z<Rc+0.75*L)
                        {
                            double rr=P[i].z;
                            overlap_surface(rr,Rc+0.75*L+0.0,Rc+0.75*L-0.0,0.0,&ov); // it gives e_is not e_si
                            Free_energy_regularization_surface= 500.0 * log(1.0-ov/omega_0_inert);
                            Free_energy_new=Free_energy_new-Free_energy_regularization_surface;
                        }
                    }
                    Free_energy_new=Free_energy_new + NsC*log(nfS_C_old/(double)NsC) + (NsC-nfS_C_old)/2.0 - NsC*log(omega_S/omega_S_0) ;


                    double Free_energy_diff= Free_energy_new - (Free_energy_old + FreeEnergy_bulk);

                    double insertion_prob=(V/(double) Np_new)*rho_ideal*e*exp(-Free_energy_diff);
                    double insertion_rand=rand()/(double)RAND_MAX;

                    for (i=0; i<Np_new; i++) myfile << P[i].maxnn << endl;
                    myfile << " t= " << t << " old_insertion= " << Free_energy_old/(double)Np << " new = " << Free_energy_new/(double)Np_new << " diff= " << Free_energy_diff << " " << Free_energy_new << " " << Free_energy_old+FreeEnergy_bulk << " " << Free_energy_new- (Free_energy_old+FreeEnergy_bulk) << " Prob= " << insertion_prob << endl;

                    if (insertion_prob > insertion_rand) Np=Np_new;
                    else
                    {
                        for (i=0; i<Np; i++) P[i]=P_cpy[i];
                        Np_new=Np;
                        nfS_C_old=nfS_C_cpy; omega_S=omega_S_cpy; ns=ns_cpy; ns_upper=ns_upper_cpy;
                    }
                }
            }
            /************remove ***********************************************/
            else if (rand_insertion_delete >=0.5)
            {
                int Rand_particle=rand() % Np;
                X=P[Rand_particle].x; Y=P[Rand_particle].y; Z=P[Rand_particle].z;

                int count_deposition_ne=0;
                for (i=0; i<Np; i++)
                {
                    if ( i!=Rand_particle ) {distance_ads(X,Y,Z,i,Lx,Ly,Lz,&dx,&dy,&dz,&d);
                    if (  d>2*Rc && Z>Rc )  count_deposition_ne=count_deposition_ne+1;
                    else break;}
                }

                if (count_deposition_ne==Np-1)
                {
                    nfS_C_cpy=nfS_C_old; omega_S_cpy=omega_S; ns_cpy=ns; ns_upper_cpy=ns_upper;
                    for (i=0; i<Np; i++)
                    {
                        P_cpy[i]=P[i];
                        P[i].maxnn=0;
                        P[i].maxn=0;
                        P[i].omega_i=0.0;
                        P[i].omega_i_inert=0.0;
                        P[i].Pn[i].nl_d=0.0;

                        for (j=0; j<Np; j++)
                        {
                            P[i].Pn[j].e_ij=0.0;
                            P[i].Pn[j].e_ij_prime=0.0;
                            P[i].Pn[j].omega_ij=0.0;
                            P[i].Pn[j].omega_ij_prime=0.0;
                            P[i].Pn[j].regu_e_ij=0.0;
                            P[i].Pn[j].regu_e_ij_prime=0.0;
                            P[i].Pn[j].nn_index=0;
                            P[i].Pn[j].n_index=0;
                            P[i].Pn[j].nb_d=0;
                        }
                    }
                    for (q=0; q<ns; q++)
                    {
                        int index_cs=P[q].Sc;
                        P[index_cs].nbS_d=0.0;
                        P[index_cs].Sd=0.0;
                        P[index_cs].e_is=0.0;
                        P[index_cs].e_is_prime=0.0;
                        P[index_cs].e_si=0.0;
                        P[index_cs].e_si_prime=0.0;
                        P[index_cs].omega_is=0.0;
                        P[index_cs].omega_is_prime=0.0;
                    }
                    nfS_C_old=NsC; omega_S=omega_S_0; ns=0; ns_upper=0;

                    Np_new=Np-1;
                    i_new=Np-1;
                    P[Rand_particle].x=P[i_new].x;P[Rand_particle].y=P[i_new].y; P[Rand_particle].z=P[i_new].z;

                    neighborlist(Np_new,&ns,&ns_upper,Lx,Ly,Lz,Sz, Sz_upper, Rc,L,R2s);
                    overlapvolume(Np_new,Rc,L,omega_0,omega_0_inert);
                    overlapvolume_surface(ns,Rc,L,R2s,omega_S_0,&omega_S);
                    fixpointiteration(Np_new, Na1, Na2, Nb1, Nb2, ns, NsC, beta, rho0, Rc, L, R2s, DG0_star, DG0_L_star, DG0_surface_star, omega_S,&nfS_C_old,DG0_2S_A_sur,DG0_2S_B_sur,DG0_3S_sur, DG0ta1a2b1, DG0ta1a2b2, DG0tb1b2a1, DG0tb1b2a2, DG0, DG0_L, DG0_a1b1, DG0_a2b2,DG0_L_a1a2,DG0_L_b1b2,DG0_br_a1a2,DG0_br_b1b2,DGT,n_alpha,n_T);

                    double Free_energy_new=0.0;
                    for (i=0; i<Np_new; i++) //
                    {
                        if (P[i].z>Rc+L+R2s)    P[i].nbS_d=0.0;
                        Free_energy_new = Free_energy_new + ( (Na1*log(P[i].nf_a1_old/(double)Na1) + Na2*log(P[i].nf_a2_old/(double)Na2) + P[i].total_bridge_WC/2.0 + P[i].Pn[i].nl_d + P[i].nbS_d/2.0) - ((Na1+Na2)*log(P[i].omega_i/omega_0)) );
                        if (P[i].count_regularization>0)
                        {
                            for(j=0; j<P[i].maxnn; j++)
                            {
                                int index=P[i].Pn[j].nn_index;
                                double r=P[i].Pn[index].d;
                                if (r < 2*Rc+0.75*L)
                                {
                                    Free_energy_regularization= 500.0 * log(1.0-P[i].Pn[index].regu_e_ij/omega_0_inert); //log(omega_0_inert-P[i].Pn[index].regu_e_ij)
                                    Free_energy_new=Free_energy_new-Free_energy_regularization;
                                }
                            }
                        }
                        if (P[i].z<Rc+0.75*L)
                        {
                            double rr=P[i].z;
                            overlap_surface(rr,Rc+0.75*L+0.0,Rc+0.75*L-0.0,0.0,&ov); // it gives e_is not e_si
                            Free_energy_regularization_surface= 500.0 * log(1.0-ov/omega_0_inert);
                            Free_energy_new=Free_energy_new-Free_energy_regularization_surface;
                        }
                    }
                    Free_energy_new=Free_energy_new + NsC*log(nfS_C_old/(double)NsC) + (NsC-nfS_C_old)/2.0 - NsC*log(omega_S/omega_S_0);

                    double Free_energy_diff= (Free_energy_new+FreeEnergy_bulk) - Free_energy_old;

                    double remove_prob = (Np/V)*(1/(rho_ideal*e))*exp(-Free_energy_diff);
                    double remove_rand=rand()/(double)RAND_MAX;

                    for (i=0; i<Np_new; i++) myfile << P[i].maxnn << endl;
                    myfile << " t= " << t << " old_removal= " << Free_energy_old/(double)Np << " new = " << Free_energy_new/(double)Np_new << " diff= " << Free_energy_diff << " " << Free_energy_new+FreeEnergy_bulk << " " << Free_energy_old << " " << Free_energy_new- (Free_energy_old+FreeEnergy_bulk) << " Prob= " << remove_prob << endl;

                    if (remove_prob>remove_rand) Np=Np_new;
                    else
                    {
                        for (i=0; i<Np; i++) P[i]=P_cpy[i];
                        Np_new=Np;
                        nfS_C_old=nfS_C_cpy; omega_S=omega_S_cpy; ns=ns_cpy; ns_upper=ns_upper_cpy;
                    }
                }
            }
        }


/**************************************  Data saved  ****************************/
        if (count_BD_steps % 1000==0)
        {

        myfile1 << "ITEM: TIMESTEP" <<  endl;
        myfile1 << t <<  endl;
        myfile1 << "ITEM: NUMBER OF ATOMS" <<  endl;
        myfile1 << Np <<  endl;
        myfile1 << "ITEM: BOX BOUNDS xy xz yz" <<  endl;
        myfile1 << "0 " << Lx << " 0" << endl;
        myfile1 << "0 " << Ly << " 0" << endl;
        myfile1 << "0 " << Lz << " 0" << endl;

        myfile1 << "ITEM: ATOMS type id x y z Fx Fy Fz nbS nl 2sa1 2sa2 3sa1 3sa2 2sb1 2sb2 3sb1 3sb2 nb3s nn nfa1 nfa2 nfb1 nfb2 lfbs nfS_C omega_i omega_S" <<  endl;


        double tot_bridges;

        for(i=0; i<Np; i++)
        {
            double total_3S_bridge=0.0;

            double tot_a1_in_Bridge = 0.0;
            double tot_a2_in_Bridge = 0.0;
            double tot_a1_in_3SBridge=0.0;
            double tot_a2_in_3SBridge=0.0;

            double tot_b1_in_Bridge=0.0;
            double tot_b2_in_Bridge=0.0;
            double tot_b1_in_3SBridge=0.0;
            double tot_b2_in_3SBridge=0.0;

            for (j=0; j<P[i].maxnn; j++)
            {
                int index=P[i].Pn[j].nn_index;
                if(P[i].type==0 && P[index].type==1)
                {
                    tot_a1_in_Bridge=tot_a1_in_Bridge+P[i].Pn[index].nb_d_a1;
                    tot_a2_in_Bridge=tot_a2_in_Bridge+P[i].Pn[index].nb_d_a2;
                    tot_a1_in_3SBridge=tot_a1_in_3SBridge+P[i].Pn[index].nb_d3S_a1;
                    tot_a2_in_3SBridge=tot_a2_in_3SBridge+P[i].Pn[index].nb_d3S_a2;
                }
                if(P[i].type==0 && P[index].type==0)
                {
                    tot_a1_in_Bridge=tot_a1_in_Bridge+P[i].Pn[index].nb_d_a1_AA;
                    tot_a2_in_Bridge=tot_a2_in_Bridge+P[i].Pn[index].nb_d_a2_AA;
                }
                if(P[i].type==1 && P[index].type==0)
                {
                    tot_b1_in_Bridge=tot_b1_in_Bridge+P[i].Pn[index].nb_d_b1;
                    tot_b2_in_Bridge=tot_b2_in_Bridge+P[i].Pn[index].nb_d_b2;
                    tot_b1_in_3SBridge=tot_b1_in_3SBridge+P[i].Pn[index].nb_d3S_b1;
                    tot_b2_in_3SBridge=tot_b2_in_3SBridge+P[i].Pn[index].nb_d3S_b2;
                }
                if(P[i].type==1 && P[index].type==1)
                {
                    tot_b1_in_Bridge=tot_b1_in_Bridge+P[i].Pn[index].nb_d_b1_BB;
                    tot_b2_in_Bridge=tot_b2_in_Bridge+P[i].Pn[index].nb_d_b2_BB;
                }

                if(P[i].type==0)
                {
                    total_3S_bridge=total_3S_bridge+P[i].Pn[index].nb_d3S_A;
                }
                if(P[i].type==1)
                {
                    total_3S_bridge=total_3S_bridge+P[i].Pn[index].nb_d3S_B;
                }

            }



            myfile1 << P[i].type << " " << i+1 << " " << P[i].x << " " << P[i].y << " " << P[i].z << " " << P[i].Fx << " " << P[i].Fy << " " << P[i].Fz  << " " << P[i].nbS_d << " " << P[i].Pn[i].nl_d << " " << tot_a1_in_Bridge << " " << tot_a2_in_Bridge << " " << tot_a1_in_3SBridge << " " << tot_a2_in_3SBridge << " " << tot_b1_in_Bridge << " " << tot_b2_in_Bridge << " " << tot_b1_in_3SBridge << " " << tot_b2_in_3SBridge  << " " << total_3S_bridge << " " << P[i].maxnn << " " <<P[i].nf_a1_old << " " << P[i].nf_a2_old << " " <<P[i].nf_b1_old << " " << P[i].nf_b2_old << " " << P[i].lfb << " " << nfS_C_old << " " << P[i].omega_i << " " << omega_S << endl;



        }
       }


        t = t + Delta_tB;
        count_BD_steps=count_BD_steps+1;
        cout << t << endl;
/***************************************************************************************************************************/

        for (i=0; i<Np; i++)
        {
            P[i].maxnn=0;
            P[i].maxn=0;
            P[i].omega_i=0.0;
            P[i].omega_i_inert=0.0;
            P[i].Fx=0.0;
            P[i].Fy=0.0;
            P[i].Fz=0.0;
            P[i].Pn[i].nl_d=0.0;
            P[i].Pn[i].nl_d_A=0.0;
            P[i].Pn[i].nl_d_B=0.0;

            for (j=0; j<Np; j++)
            {
                P[i].Pn[j].e_ij=0.0;
                P[i].Pn[j].e_ij_prime=0.0;
                P[i].Pn[j].omega_ij=0.0;
                P[i].Pn[j].omega_ij_prime=0.0;
                P[i].Pn[j].regu_e_ij=0.0;
                P[i].Pn[j].regu_e_ij_prime=0.0;
                P[i].Pn[j].nn_index=0;
                P[i].Pn[j].n_index=0;
                P[i].Pn[j].nb_d=0;

                P[i].Pn[j].nb_d_01=0.0;
                P[i].Pn[j].nb_d_10=0.0;
                P[i].Pn[j].nb_d_00=0.0;
                P[i].Pn[j].nb_d_11=0.0;
                P[i].Pn[j].nb_d3S_A=0.0;
                P[i].Pn[j].nb_d3S_B=0.0;
            }
        }
        for (q=0; q<ns; q++)
        {
            int index_cs=P[q].Sc;
            P[index_cs].nbS_d=0.0;
            P[index_cs].Sd=0.0;
            P[index_cs].e_is=0.0;
            P[index_cs].e_is_prime=0.0;
            P[index_cs].e_si=0.0;
            P[index_cs].e_si_prime=0.0;
            P[index_cs].omega_is=0.0;
            P[index_cs].omega_is_prime=0.0;

            P[index_cs].nbS_d_A=0.0;
            P[index_cs].nbS_d_B=0.0;
            P[index_cs].lfb=0.0;
        }


/***************************************************************************************************************************/

    }
    myfile.close();
    myfile1.close();
    return 0;
}






/******* Functions defination **************/
/******* Functions defination **************/
/******* Functions defination **************/
/******* Functions defination **************/

void neighborlist(int Np,int *ns, int *ns_upper,double Lx,double Ly,double Lz,double Sz, double Sz_upper,double Rc,double L,double R2s)
{
    *ns=0;
    *ns_upper=0;
    int i,j,count_n,count_nn,count_regularization;
    double dx,dy,dz,d;
    for(i=0; i<Np; i++)
    {
        count_regularization=0;
        count_n=0;
        for (j=i+1; j<Np; j++)
        {
            distance(i,j,Lx,Ly,Lz,&dx,&dy,&dz,&d);
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
                if (d < 2*Rc+0.75*L) count_regularization=count_regularization+1;
            }
            P[i].count_regularization=count_regularization;
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
        }
        else if (d_upper<Rc+0.75*L)
        {
            P[*ns_upper].Sc_upper=i;
            *ns_upper=*ns_upper+1;
        }
    }
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

void distance_ads(double X, double Y, double Z,int i, double Lx, double Ly, double Lz,double *dx,double *dy,double *dz, double *d)
{
    *dx = P[i].x-X;
    *dx = *dx - Lx * round(*dx/Lx);

    *dy = P[i].y-Y;
    *dy = *dy - Ly * round(*dy/Ly);

    *dz = P[i].z-Z;

    *d=sqrt( SQ(*dx) + SQ(*dy) + SQ(*dz));
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
            printf(" r=%lf",r);
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
            else if (r<(2*Rc+1.9*L) && r>=2*Rc+L) ////
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
            else if (r<2*Rc+L && r>2*Rc)
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

            else if (r<=2*Rc) //
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
        printf(" w_ij=%lf",w_ij);
        printf(" e_ij_sum=%lf\n",e_ij_sum);
    }

}
void overlapvolume_surface(int ns,double Rc,double L,double R2s,double omega_S_0, double *omega_S)
{
    int i;
    double e_is, e_si, w_is, e_is_prime, e_si_prime, w_is_prime;
    double ov, ov1, ov2, ov3, ov_prime, ov1_prime, ov2_prime, ov3_prime;
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
        else if (r<(Rc+L+R2s) && r >= (Rc+L+0.9*R2s)){
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
        else if (r < (Rc+L+0.9*R2s) && r >= (Rc+L)){
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
            overlap_surface(r,Dcs,dcs,R2,&ov3); //overlap between (Rc) and corona of surface
            overlap_surface_prime(r,Dcs,dcs,R2,&ov3_prime);

            w_is=ov1-ov2-ov3;
            w_is_prime=ov1_prime-ov2_prime-ov3_prime;

            e_is=ov3;       // colloids linkers are blocked by the surface
            e_is_prime=ov3_prime;


            e_si=ov2;       //surface linkes are blocked by the sphere
            e_si_prime=ov2_prime;
        }

        P[index_cs].e_is = e_is;
        P[index_cs].e_is_prime = e_is_prime;
        P[index_cs].e_si = e_si;
        P[index_cs].e_si_prime = e_si_prime;
        P[index_cs].omega_is = w_is;
        P[index_cs].omega_is_prime = w_is_prime;
        P[index_cs].omega_i = P[index_cs].omega_i-e_is;
        e_si_sum=e_si_sum+e_si;
        printf("%d P[index_cs].e_is %lf P[index_cs].e_is_prime %lf P[index_cs].e_si %lf P[index_cs].e_si_prime %lf P[index_cs].omega_is %lf P[index_cs].omega_is_prime %lf P[index_cs].omega_i %lf\n",index_cs,P[index_cs].e_is,P[index_cs].e_is_prime,P[index_cs].e_si,P[index_cs].e_si_prime,P[index_cs].omega_is,P[index_cs].omega_is_prime,P[index_cs].omega_i);

    }

    *omega_S=omega_S_0-e_si_sum;
}

void fixpointiteration(int Np, int Na1, int Na2,int Nb1, int Nb2, int ns, int NsC, double beta, double Rc, double L, double R2s, double rho0, double DG0_star, double DG0_L_star, double DG0_surface_star, double omega_S, double *nfS_C_old,double DG0_2S_A_sur,double DG0_2S_B_sur,double DG0_3S_sur,double DG0ta1a2b1,double DG0ta1a2b2,double DG0tb1b2a1,double DG0tb1b2a2,double DG0, double DG0_L, double DG0_a1b1, double DG0_a2b2, double DG0_L_a1a2, double DG0_L_b1b2, double DG0_br_a1a2, double DG0_br_b1b2, double DGT, double n_alpha, double n_T)
{
    double Sz=0.0;
    int i,j,q;
    double sum,sum_B;
    double sum_a1,sum_a2,sum_b1,sum_b2,sum_C;
    bool flag;
    double Np_roots;
    int index;
    int Niteration=0;
    double error_a1,error_a2,error_b1,error_b2,error_C;
    double Free_energy_regularization;
    double betaDG_ij;
    double betaDG_ij_surface;
    double betaDG_ii_L;
    double nfS_C_new;

    double K0 = exp(beta*(-DG0-DGT));
    double K1 = exp(beta*(-DGT));
    for (i=0; i<Np; i++)
    {
        if (P[i].type==0)
        {
        P[i].nf_a1_old=rand()% Na1+1; //Na1=na1
        P[i].nf_a2_old=rand()% Na2+1; //Na2=na2
        P[i].nf_a1_new=P[i].nf_a1_old;
        P[i].nf_a2_new=P[i].nf_a2_old;
        }

        else if (P[i].type==1)
        {
        P[i].nf_b1_old= rand()% Nb1+1; //Na1=na1
        P[i].nf_b2_old= rand()% Nb2+1; //Na2=na2
        P[i].nf_b1_new=P[i].nf_b1_old;
        P[i].nf_b2_new=P[i].nf_b2_old;
        }

     }
     *nfS_C_old=rand()% NsC+1;   //NsC=Nb1
    nfS_C_new=*nfS_C_old;

    flag=true;
    while (flag==true && Niteration < 10000000)
    {
        sum_C=0.0;
        for(q=0; q<ns; q++)  //ns= no. of particle that overlap with the surface, here ns=1
        {
            int index_cs=P[q].Sc; //Sc= index of particle close to the surface
            double PF = P[index_cs].omega_is/(P[index_cs].omega_i*omega_S*rho0); //eqn4(for suface) arg of log

            double PF_L = 1/(P[index_cs].omega_i*rho0);



            if(P[index_cs].type==0)
            {
                sum_C=sum_C+P[index_cs].nf_a1_old*PF*K0+P[index_cs].nf_a1_old*P[index_cs].nf_a2_old*PF_L*PF*K0*K1*(n_alpha+1);
            }
            if(P[index_cs].type==1)
            {
                sum_C=sum_C+P[index_cs].nf_b2_old*PF*K0;
            }

        }
        nfS_C_new=NsC/(1+sum_C); //eqn S12(for surface) (3rd eqn full) (no change required in the previous eqn)
        error_C=nfS_C_new-*nfS_C_old;
        *nfS_C_old=nfS_C_new;
        Np_roots=0.0;
        for (i=0; i<Np; i++)
        {
            sum_a1=0.0;
            sum_a2=0.0;
            sum_b1=0.0;
            sum_b2=0.0;
            for (j=0; j<P[i].maxnn; j++)
            {
                index=P[i].Pn[j].nn_index;
                double PF = P[i].Pn[index].omega_ij/(P[i].omega_i*P[index].omega_i*rho0);


                double PF_L = 1/(P[i].omega_i*rho0);
                double PF_Lj = 1/(P[index].omega_i*rho0);
                if(P[i].type==0 && P[index].type==1)
                {

                    sum_a1=sum_a1+P[index].nf_b1_old*K0*PF+P[index].nf_b1_old*P[index].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_a2_old*P[index].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_a2_old*P[index].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_L;

                    sum_a2=sum_a2+P[index].nf_b2_old*K0*PF+P[index].nf_b2_old*P[index].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_a1_old*P[index].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_a1_old*P[index].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_L;
                }

                if(P[i].type==0 && P[index].type==0)
                {
                    sum_a1=sum_a1+P[index].nf_a2_old*K0*PF;
                    sum_a2=sum_a2+P[index].nf_a1_old*K0*PF;
                }

                if(P[i].type==1 && P[index].type==0)
                {
                    sum_b1=sum_b1+P[index].nf_a1_old*K0*PF+P[index].nf_a1_old*P[index].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_b2_old*P[index].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_b2_old*P[index].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_L;

                    sum_b2=sum_b2+P[index].nf_a2_old*K0*PF+P[index].nf_a2_old*P[index].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_b1_old*P[index].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_b1_old*P[index].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_L;
                }

                if(P[i].type==1 && P[index].type==1)
                {
                    sum_b1=sum_b1+P[index].nf_b2_old*K0*PF;
                    sum_b2=sum_b2+P[index].nf_b1_old*K0*PF;
                }


            }
            double dz = P[i].z- Sz;
            double d = sqrt(SQ(dz));
            if (d<Rc+L+R2s && d>Rc)//to check if particles overlap with surface //R2s= width of surface //Rc= radius of colloid
            {                                   //L= length of DNA linkers
                double PF = P[i].omega_is/(P[i].omega_i*omega_S*rho0);

                double PF_L = 1/(P[i].omega_i*rho0);

                if(P[i].type==0)
                {
                sum_a1=sum_a1+ *nfS_C_old *PF*K0+ *nfS_C_old * P[i].nf_a2_old * PF_L*PF*K0*K1*(n_alpha+1); //(if overlap) modified eqn S12 1st eqn,denominator 2nd & 3rd term

                sum_a2=sum_a2+ *nfS_C_old * P[i].nf_a1_old * PF_L*PF*K0*K1*(n_alpha+1); //if overlap,modified eqn S12 2nd eqn
                }
                if(P[i].type==1)
                {
                sum_b2=sum_b2+ *nfS_C_old *PF*K0;
                }
            }
            double PF_L = 1/(P[i].omega_i*rho0);

            if(P[i].type==0)
            {
            sum_a1=sum_a1+P[i].nf_a2_old*K0*PF_L; //eqn S12 1st eqn summetion part
            sum_a2=sum_a2+P[i].nf_a1_old*K0*PF_L; //eqn S12 2nd eqn summetion part
            }
            if(P[i].type==1)
            {
            sum_b1=sum_b1+P[i].nf_b2_old*K0*PF_L;
            sum_b2=sum_b2+P[i].nf_b1_old*K0*PF_L;
            }

            if(P[i].type==0)
            {
            P[i].nf_a1_new=Na1/(1+sum_a1); //modified eqn S12 1st eqn as per our requirement
            error_a1=P[i].nf_a1_old-P[i].nf_a1_new;
            P[i].nf_a2_new=Na2/(1+sum_a2); //modified eqn S12 2nd eqn as per our requirement
            error_a2=P[i].nf_a2_old-P[i].nf_a2_new;
            }
            if(P[i].type==1)
            {
            P[i].nf_b1_new=Nb1/(1+sum_b1); //modified eqn S12 1st eqn as per our requirement
            error_b1=P[i].nf_b1_old-P[i].nf_b1_new;
            P[i].nf_b2_new=Nb2/(1+sum_b2); //modified eqn S12 1st eqn as per our requirement
            error_b2=P[i].nf_b2_old-P[i].nf_b2_new;
            }

            if (fabs(error_a1) < 0.0000001 && fabs(error_a2) < 0.0000001 && fabs(error_b1) < 0.0000001 && fabs(error_b2) < 0.0000001 && fabs(error_C) < 0.0000001) Np_roots=Np_roots+1;

            if(P[i].type==0)
            {
            P[i].nf_a1_old=P[i].nf_a1_new;
            P[i].nf_a2_old=P[i].nf_a2_new;
            }
            if(P[i].type==1)
            {
            P[i].nf_b1_old=P[i].nf_b1_new;
            P[i].nf_b2_old=P[i].nf_b2_new;
            }
        }

        if (Np_roots==Np)   flag=false;
        Niteration=Niteration+1;
    }


    for (i=0; i<Np; i++)
    {
        double PF_L = 1/(P[i].omega_i*rho0);


        if(P[i].type==0)
        {
            P[i].Pn[i].nl_d = P[i].nf_a1_old*P[i].nf_a2_old*K0*PF_L ;
        }
        if(P[i].type==1)
        {
            P[i].Pn[i].nl_d = P[i].nf_b1_old*P[i].nf_b2_old*K0*PF_L;
        }


        double total_bridge=0.0;
        double total_3S_bridge=0.0;

        double tot_a1_in_Bridge=0.0;
        double tot_a2_in_Bridge=0.0;
        double tot_b1_in_Bridge=0.0;
        double tot_b2_in_Bridge=0.0;

        double tot_a1_in_3SBridge=0.0;
        double tot_a2_in_3SBridge=0.0;
        double tot_b1_in_3SBridge=0.0;
        double tot_b2_in_3SBridge=0.0;

        for (j=0; j<P[i].maxnn; j++)
        {
            int index=P[i].Pn[j].nn_index;

            double PF = P[i].Pn[index].omega_ij/(P[i].omega_i*P[index].omega_i);
            double PF_Lj = 1/(P[index].omega_i*rho0);


            if(P[i].type==0 && P[index].type==1)
            {
            P[i].Pn[index].nb_d_a1=P[i].nf_a1_old*P[index].nf_b1_old*K0*PF;
            P[i].Pn[index].nb_d_a2=P[i].nf_a2_old*P[index].nf_b2_old*K0*PF;

            P[i].Pn[index].nb_d_01=P[i].nf_a1_old*P[index].nf_b1_old*K0*PF+P[i].nf_a2_old*P[index].nf_b2_old*K0*PF;
            }
            if(P[i].type==1 && P[index].type==0)
            {
            P[i].Pn[index].nb_d_b1=P[index].nf_a1_old*P[i].nf_b1_old*K0*PF;
            P[i].Pn[index].nb_d_b2=P[index].nf_a2_old*P[i].nf_b2_old*K0*PF;

            P[i].Pn[index].nb_d_10=P[index].nf_a1_old*P[i].nf_b1_old*K0*PF+P[index].nf_a2_old*P[i].nf_b2_old*K0*PF;
            }

            if(P[i].type==0 && P[index].type==0)
            {
            P[i].Pn[index].nb_d_a1_AA=P[i].nf_a1_old*P[index].nf_a2_old*K0*PF;
            P[i].Pn[index].nb_d_a2_AA=P[i].nf_a2_old*P[index].nf_a1_old*K0*PF;

            P[i].Pn[index].nb_d_00=(P[i].nf_a1_old*P[index].nf_a2_old+P[i].nf_a2_old*P[index].nf_a1_old)*K0*PF;
            }
            if(P[i].type==1 && P[index].type==1)
            {
            P[i].Pn[index].nb_d_b1_BB=P[index].nf_b2_old*P[i].nf_b1_old*K0*PF;
            P[i].Pn[index].nb_d_b2_BB=P[index].nf_b1_old*P[i].nf_b2_old*K0*PF;

            P[i].Pn[index].nb_d_11=(P[index].nf_b2_old*P[i].nf_b1_old+P[index].nf_b1_old*P[i].nf_b2_old)*K0*PF;
            }


            if(P[i].type==0 && P[index].type==1)
            {


                P[i].Pn[index].nb_d3S_a1=P[i].nf_a1_old*P[i].nf_a2_old*P[index].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_a1_old*P[i].nf_a2_old*P[index].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_L+P[index].nf_b1_old*P[index].nf_b2_old*P[i].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_Lj;

                P[i].Pn[index].nb_d3S_a2=P[i].nf_a1_old*P[i].nf_a2_old*P[index].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_a1_old*P[i].nf_a2_old*P[index].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_L+P[index].nf_b1_old*P[index].nf_b2_old*P[i].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_Lj;

            }

            if(P[i].type==1 && P[index].type==0)
            {


                P[i].Pn[index].nb_d3S_b1=P[index].nf_a1_old*P[index].nf_a2_old*P[i].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_b1_old*P[i].nf_b2_old*P[index].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_b1_old*P[i].nf_b2_old*P[index].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_L;

                P[i].Pn[index].nb_d3S_b2=P[index].nf_a1_old*P[index].nf_a2_old*P[i].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_b1_old*P[i].nf_b2_old*P[index].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_b1_old*P[i].nf_b2_old*P[index].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_L;
            }

            if(P[i].type==0 && P[index].type==1)
            {
            tot_a1_in_Bridge=tot_a1_in_Bridge+P[i].Pn[index].nb_d_a1;
            tot_a2_in_Bridge=tot_a2_in_Bridge+P[i].Pn[index].nb_d_a2;

            }

            if(P[i].type==0 && P[index].type==0)
            {
            tot_a1_in_Bridge=tot_a1_in_Bridge+P[i].Pn[index].nb_d_a1_AA;
            tot_a2_in_Bridge=tot_a2_in_Bridge+P[i].Pn[index].nb_d_a2_AA;

            }

            if(P[i].type==1 && P[index].type==0)
            {
            tot_b1_in_Bridge=tot_b1_in_Bridge+P[i].Pn[index].nb_d_b1;
            tot_b2_in_Bridge=tot_b2_in_Bridge+P[i].Pn[index].nb_d_b2;

            }

            if(P[i].type==1 && P[index].type==1)
            {
            tot_b1_in_Bridge=tot_b1_in_Bridge+P[i].Pn[index].nb_d_b1_BB;
            tot_b2_in_Bridge=tot_b2_in_Bridge+P[i].Pn[index].nb_d_b2_BB;

            }

            if(P[i].type==0 && P[index].type==1)
            {
            tot_a1_in_3SBridge=tot_a1_in_3SBridge+P[i].Pn[index].nb_d3S_a1;
            tot_a2_in_3SBridge=tot_a2_in_3SBridge+P[i].Pn[index].nb_d3S_a2;
            }
            if(P[i].type==1 && P[index].type==0)
            {
            tot_b1_in_3SBridge=tot_b1_in_3SBridge+P[i].Pn[index].nb_d3S_b1;
            tot_b2_in_3SBridge=tot_b2_in_3SBridge+P[i].Pn[index].nb_d3S_b2;
            }




            if(P[i].type==0 && P[index].type==1)
            {
                P[i].Pn[index].nb_d3S_A = P[i].nf_a1_old*P[i].nf_a2_old*P[index].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_a1_old*P[i].nf_a2_old*P[index].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_L+P[index].nf_b1_old*P[index].nf_b2_old*P[i].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[index].nf_b1_old*P[index].nf_b2_old*P[i].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_Lj;

            }
            if(P[i].type==1 && P[index].type==0)
            {
                P[i].Pn[index].nb_d3S_B = P[index].nf_a1_old*P[index].nf_a2_old*P[i].nf_b1_old*K0*K1*(n_alpha+1)*PF*PF_Lj+P[i].nf_b1_old*P[i].nf_b2_old*P[index].nf_a1_old*K0*K1*(n_alpha+1)*PF*PF_L+P[i].nf_b1_old*P[i].nf_b2_old*P[index].nf_a2_old*K0*K1*(n_alpha+1)*PF*PF_L+P[index].nf_a1_old*P[index].nf_a2_old*P[i].nf_b2_old*K0*K1*(n_alpha+1)*PF*PF_Lj;

            }
            else if(P[i].type==0 && P[index].type==0)
            {
                P[i].Pn[index].nb_d3S = 0.0;
            }
            else if(P[i].type==1 && P[index].type==1)
            {
                P[i].Pn[index].nb_d3S = 0.0;
            }

            if(P[i].type==0)
            {
                total_3S_bridge=total_3S_bridge+P[i].Pn[index].nb_d3S_A;
            }

            if(P[i].type==1)
            {
                total_3S_bridge=total_3S_bridge+P[i].Pn[index].nb_d3S_B;
            }


        }
        if(P[i].type==0)
        {
            P[i].total_free_loops_A=P[i].nf_a1_old+P[i].nf_a2_old+P[i].Pn[i].nl_d;
        }
        if(P[i].type==1)
        {
            P[i].total_free_loops_B=P[i].nf_b1_old+P[i].nf_b2_old+P[i].Pn[i].nl_d;
        }

        if(P[i].type==0)
        {
        total_bridge=tot_a1_in_Bridge+tot_a2_in_Bridge;
        }
        if(P[i].type==1)
        {
        total_bridge=tot_b1_in_Bridge+tot_b2_in_Bridge;
        }

    }


    for (q=0; q<ns; q++)
    {
        int index_cs=P[q].Sc;
        double PF = P[index_cs].omega_is/(P[index_cs].omega_i*omega_S*rho0);
        double PF_L = 1/(P[index_cs].omega_i*rho0);



        if(P[index_cs].type==0)
        {
            P[index_cs].nbS_d_A=P[index_cs].nf_a1_old * *nfS_C_old *PF*K0;
        }
        if(P[index_cs].type==1)
        {
            P[index_cs].nbS_d_B=P[index_cs].nf_b2_old * *nfS_C_old *PF*K0;
        }

        P[index_cs].nbS_d=P[index_cs].nbS_d_A+P[index_cs].nbS_d_B;


        if(P[index_cs].type==0)
        {
            P[index_cs].lfb=*nfS_C_old*P[index_cs].nf_a1_old*P[index_cs].nf_a2_old * PF_L*PF*K0*K1*(n_alpha+1);
        }
        else
        {
            P[index_cs].lfb = 0.0;
        }


    }
}


void Force(int Np, int ns, int ns_upper, double Rc, double L, double R2s,double Sz, double Sz_upper, double omega_S, double nfS_C_old, double beta,double Lx, double Ly, double Lz, double omega_0_inert)
{
    int i,j,q;
    double dx,dy,dz,d,ov,ov_prime;
    double term1x, term1y, term1z;
    double term2x, term2y, term2z;
    double term3x, term3y, term3z;
    double sum_term1x, sum_term1y, sum_term1z;
    double sum_term2x, sum_term2y, sum_term2z;
    double sum_term3x, sum_term3y, sum_term3z;

    double sum_term1x_10, sum_term1y_10, sum_term1z_10;
    double sum_term2x_10, sum_term2y_10, sum_term2z_10;

    double sum_term1x_01, sum_term1y_01, sum_term1z_01;
    double sum_term2x_01, sum_term2y_01, sum_term2z_01;

    double sum_term1x_00, sum_term1y_00, sum_term1z_00;
    double sum_term2x_00, sum_term2y_00, sum_term2z_00;

    double sum_term1x_11, sum_term1y_11, sum_term1z_11;
    double sum_term2x_11, sum_term2y_11, sum_term2z_11;

    double term1x_01,term1x_10,term1x_00,term1x_11;
    double term1y_01,term1y_10,term1y_00,term1y_11;
    double term1z_01,term1z_10,term1z_00,term1z_11;

    double term2x_01,term2x_10,term2x_00,term2x_11;
    double term2y_01,term2y_10,term2y_00,term2y_11;
    double term2z_01,term2z_10,term2z_00,term2z_11;


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

        double sum=0.0;

        for (j=0; j<P[i].maxnn; j++)
        {
            int index=P[i].Pn[j].nn_index;
            distance(i,index,Lx,Ly,Lz,&dx,&dy,&dz,&d);

            double ux=dx/P[i].Pn[index].d;
            double uy=dy/P[i].Pn[index].d;
            double uz=dz/P[i].Pn[index].d;



            if(P[i].type==0 && P[index].type==1)
            {
                printf("hii01\n");
                term1x_01 = (P[i].Pn[index].nb_d_01+P[i].Pn[index].nb_d3S_A) * ( (P[i].Pn[index].omega_ij_prime*ux) /P[i].Pn[index].omega_ij);
                term1y_01 = (P[i].Pn[index].nb_d_01+P[i].Pn[index].nb_d3S_A) * ( (P[i].Pn[index].omega_ij_prime*uy) /P[i].Pn[index].omega_ij);
                term1z_01 = (P[i].Pn[index].nb_d_01+P[i].Pn[index].nb_d3S_A) * ( (P[i].Pn[index].omega_ij_prime*uz) /P[i].Pn[index].omega_ij);
            }
            if(P[i].type==1 && P[index].type==0)
            {
                printf("hii10\n");
                term1x_10 = (P[i].Pn[index].nb_d_10+P[i].Pn[index].nb_d3S_B) * ( (P[i].Pn[index].omega_ij_prime*ux) /P[i].Pn[index].omega_ij);
                term1y_10 = (P[i].Pn[index].nb_d_10+P[i].Pn[index].nb_d3S_B) * ( (P[i].Pn[index].omega_ij_prime*uy) /P[i].Pn[index].omega_ij);
                term1z_10 = (P[i].Pn[index].nb_d_10+P[i].Pn[index].nb_d3S_B) * ( (P[i].Pn[index].omega_ij_prime*uz) /P[i].Pn[index].omega_ij);
            }
            if(P[i].type==0 && P[index].type==0)
            {
                printf("hii00\n");
                term1x_00 = (P[i].Pn[index].nb_d_00) * ( (P[i].Pn[index].omega_ij_prime*ux) /P[i].Pn[index].omega_ij);
                term1y_00 = (P[i].Pn[index].nb_d_00) * ( (P[i].Pn[index].omega_ij_prime*uy) /P[i].Pn[index].omega_ij);
                term1z_00 = (P[i].Pn[index].nb_d_00) * ( (P[i].Pn[index].omega_ij_prime*uz) /P[i].Pn[index].omega_ij);
            }
            if(P[i].type==1 && P[index].type==1)
            {
                printf("hii11\n");
                term1x_11 = (P[i].Pn[index].nb_d_11) * ( (P[i].Pn[index].omega_ij_prime*ux) /P[i].Pn[index].omega_ij);
                term1y_11 = (P[i].Pn[index].nb_d_11) * ( (P[i].Pn[index].omega_ij_prime*uy) /P[i].Pn[index].omega_ij);
                term1z_11 = (P[i].Pn[index].nb_d_11) * ( (P[i].Pn[index].omega_ij_prime*uz) /P[i].Pn[index].omega_ij);
            }



            if(P[i].type==0 && P[index].type==1)
            {
             sum_term1x=sum_term1x+term1x_01;
             sum_term1y=sum_term1y+term1y_01;
             sum_term1z=sum_term1z+term1z_01;
            }
            if(P[i].type==1 && P[index].type==0)
            {
             sum_term1x=sum_term1x+term1x_10;
             sum_term1y=sum_term1y+term1y_10;
             sum_term1z=sum_term1z+term1z_10;
            }
            if(P[i].type==0 && P[index].type==0)
            {
             sum_term1x=sum_term1x+term1x_00;
             sum_term1y=sum_term1y+term1y_00;
             sum_term1z=sum_term1z+term1z_00;
            }
            if(P[i].type==1 && P[index].type==1)
            {
             sum_term1x=sum_term1x+term1x_11;
             sum_term1y=sum_term1y+term1y_11;
             sum_term1z=sum_term1z+term1z_11;
            }

            printf("Ptype=%d Ntype=%d  sum_term1x_01=%lf\tsum_term1y_01=%lf\tsum_term1z_01=%lf\n",P[i].type,P[index].type,sum_term1x_01,sum_term1y_01,sum_term1z_01);


            if(P[i].type==0 && P[index].type==1)
            {
            term2x_01 = (P[i].total_free_loops_A) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                        (P[index].total_free_loops_B) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

            term2y_01 = (P[i].total_free_loops_A) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                         (P[index].total_free_loops_B) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

            term2z_01 = (P[i].total_free_loops_A) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                         (P[index].total_free_loops_B) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }

            if(P[i].type==1 && P[index].type==0)
            {
            term2x_10 = (P[i].total_free_loops_B) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                        (P[index].total_free_loops_A) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

            term2y_10 = (P[i].total_free_loops_B) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                         (P[index].total_free_loops_A) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

            term2z_10 = (P[i].total_free_loops_B) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                         (P[index].total_free_loops_A) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }

            if(P[i].type==0 && P[index].type==0)
            {
            term2x_00 = (P[i].total_free_loops_A) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                        (P[index].total_free_loops_A) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

            term2y_00 = (P[i].total_free_loops_A) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                         (P[index].total_free_loops_A) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

            term2z_00 = (P[i].total_free_loops_A) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                         (P[index].total_free_loops_A) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }

            if(P[i].type==1 && P[index].type==1)
            {
            term2x_11 = (P[i].total_free_loops_B) * ((P[i].Pn[index].e_ij_prime*ux)/P[i].omega_i) +
                        (P[index].total_free_loops_B) * ((P[index].Pn[i].e_ij_prime*ux)/P[index].omega_i);

            term2y_11 = (P[i].total_free_loops_B) * ((P[i].Pn[index].e_ij_prime*uy)/P[i].omega_i) +
                         (P[index].total_free_loops_B) * ((P[index].Pn[i].e_ij_prime*uy)/P[index].omega_i);

            term2z_11 = (P[i].total_free_loops_B) * ((P[i].Pn[index].e_ij_prime*uz)/P[i].omega_i) +
                         (P[index].total_free_loops_B) * ((P[index].Pn[i].e_ij_prime*uz)/P[index].omega_i);
            }



            if(P[i].type==0 && P[index].type==1)
            {
            sum_term2x = sum_term2x + term2x_01;
            sum_term2y = sum_term2y + term2y_01;
            sum_term2z = sum_term2z + term2z_01;
            }
            if(P[i].type==1 && P[index].type==0)
            {
            sum_term2x = sum_term2x + term2x_10;
            sum_term2y = sum_term2y + term2y_10;
            sum_term2z = sum_term2z + term2z_10;
            }
            if(P[i].type==0 && P[index].type==0)
            {
            sum_term2x = sum_term2x + term2x_00;
            sum_term2y = sum_term2y + term2y_00;
            sum_term2z = sum_term2z + term2z_00;
            }
            if(P[i].type==1 && P[index].type==1)
            {
            sum_term2x = sum_term2x + term2x_11;
            sum_term2y = sum_term2y + term2y_11;
            sum_term2z = sum_term2z + term2z_11;
            }
            if (P[i].Pn[index].d < 2*Rc+0.75*L) {
                term3x = beta * 2 * 500.0 * ( (P[i].Pn[index].regu_e_ij_prime*ux)/( omega_0_inert-P[i].Pn[index].regu_e_ij ) ); //
                term3y = beta * 2 * 500.0 * ( (P[i].Pn[index].regu_e_ij_prime*uy)/( omega_0_inert-P[i].Pn[index].regu_e_ij ) ); //P[i].omega_i_inert
                term3z = beta * 2 * 500.0 * ( (P[i].Pn[index].regu_e_ij_prime*uz)/( omega_0_inert-P[i].Pn[index].regu_e_ij ) );




                sum_term3x = sum_term3x+term3x;
                sum_term3y = sum_term3y+term3y;
                sum_term3z = sum_term3z+term3z;
            }
        }
         P[i].Fx=(sum_term1x-sum_term2x-sum_term3x)/beta;
         P[i].Fy=(sum_term1y-sum_term2y-sum_term3y)/beta;
         P[i].Fz=(sum_term1z-sum_term2z-sum_term3z)/beta;

    printf("Ptype=%d P[i].Fx=%lf P[i].Fy=%lf P[i].Fz=%lf\n",P[i].type,P[i].Fx,P[i].Fy,P[i].Fz);


    }

    for (q=0; q<ns; q++)
    {
        int index_cs=P[q].Sc;
        double dz_surface=P[index_cs].z-Sz;
        double uz=dz_surface/sqrt(dz_surface*dz_surface);


        if(P[index_cs].type==0)
        {
            term1z=  (P[index_cs].nbS_d_A+P[index_cs].lfb) * ( (P[index_cs].omega_is_prime*uz) / P[index_cs].omega_is );
            printf("P[index_cs].nbS_d_A=%lf\tP[index_cs].lfb=%lf\tterm1z=%lf\n",P[index_cs].nbS_d_A,P[index_cs].lfb,term1z);

        }
        if(P[index_cs].type==1)
        {
            term1z=  P[index_cs].nbS_d_B * ( (P[index_cs].omega_is_prime*uz) / P[index_cs].omega_is );
        }



        if(P[index_cs].type==0)
        {
        term2z=(P[index_cs].total_free_loops_A) *((P[index_cs].e_is_prime*uz)/P[index_cs].omega_i)+nfS_C_old*((P[index_cs].e_si_prime*uz)/omega_S);
        }

        if(P[index_cs].type==1)
        {
        term2z=(P[index_cs].total_free_loops_B) *((P[index_cs].e_is_prime*uz)/P[index_cs].omega_i)+nfS_C_old*((P[index_cs].e_si_prime*uz)/omega_S);
        }

        double omega_0_inert = 4*M_PI*SQ(Rc)*0.75*L;
        if (P[index_cs].Sd < Rc+0.75*L) {
            double rr=P[index_cs].Sd;
            overlap_surface(rr,Rc+0.75*L+0.0,Rc+0.75*L-0.0,0.0,&ov); // it gives e_is not e_si
            overlap_surface_prime(rr,Rc+0.75*L+0.0,Rc+0.75*R2s-0.0,0.0,&ov_prime);

            term3z = beta * 2 * 500.0 * ( (ov_prime*uz)/(omega_0_inert-ov) );
        }
        else {term3x=0.0; term3y=0.0; term3z=0.0;}
        P[index_cs].Fz=P[index_cs].Fz + (term1z - term2z - term3z)/beta;
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
        P[index_cs].Fz=P[index_cs].Fz + ((- term3z)/beta);
    }
}
double Isolated_colloids(int Na1, int Na2, double DG0_L_star,double beta,double rho0, double omega_0)
{
    double nf_A_bulk, nf_B_bulk,PF_L,DG_ii_L,betaDG_ii_L,nf_A_bulk_new,nf_B_bulk_new,ql;
    double error=0.000001;
    bool flag;

    nf_A_bulk = rand()% Na1+1;
    nf_B_bulk = rand()% Na2+1;

    PF_L = 1/(omega_0*rho0);

    DG_ii_L = DG0_L_star/beta-log(PF_L)/beta;

    betaDG_ii_L=beta*DG_ii_L;
    ql=exp(-betaDG_ii_L);
    flag=true;
    while (flag==true)
    {
        nf_A_bulk_new=Na1/(1 + nf_B_bulk*ql);
        nf_B_bulk_new=Na2/(1 + nf_A_bulk*ql);
        double diff_A=fabs(nf_A_bulk_new-nf_A_bulk);
        double diff_B=fabs(nf_B_bulk_new-nf_B_bulk);


        if (diff_A<error && diff_B<error) flag=false;
        nf_A_bulk=nf_A_bulk_new;
        nf_B_bulk=nf_B_bulk_new;
    }

    double FreeEnergy_bulk=(Na1*log((nf_A_bulk*nf_B_bulk)/(Na1*Na2))+(Na1-nf_A_bulk));
    return (FreeEnergy_bulk);
}




