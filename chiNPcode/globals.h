#include <complex>
#include <fstream>
#include <iostream>
//#include "malloc.h"
#include <cstdlib>
#include <stdio.h>
#include <string.h>

#ifdef PAR
#include "mpi.h"
#include "fftw3-mpi.h"
#else
#include "fftw3.h"
#endif

#define min(A,B) ((A)<(B) ? (A) : (B) )
#define max(A,B) ((A)>(B) ? (A) : (B) )

#define PI 3.14159265358979323846263383
using namespace std;

#define Dim 2 //Dimension of my system


// Global variables
#ifndef MAIN
extern
#endif
double dx[Dim], L[Dim], V, Vf, V_ab, V_1_fld_np,
       a_smear, a_squared, C, rho0, lam_mi, lam_pl,
       kappaN, chiABN, chiACN, chiBCN, fD, nD, nAH, nCH, Nch, phiAH, phiCH,
       nP, nFP,
       n_samples, ic_pre[2], ic_dir[2], ic_period[2],
       top_wall_lamA, top_wall_lamB, bot_wall_lamA, bot_wall_lamB,
       top_wall_lamC, bot_wall_lamC,
       wallT, wallXi, exp_nr_c[8][Dim],
       exp_nr_u[8][Dim], L_nr, R_nr, xi_nr,
       *theta, *phi, *theta_weights, *phi_weights, np_frac, smwp_min,
       exp_nr_chiAPN, exp_nr_chiBPN, exp_nr_chiCPN,
       L_low, L_high, L_step, brent_tol, error_tol,
       L_ideal, min_H_over_V, gsd_ds, gsd_tol,
       muA, muB, H;

#ifndef MAIN
extern
#endif
FILE *brent_otp, *dot;

#ifndef MAIN
extern
#endif
complex<double> I, *wpl, *wa, *wb, *wc, *wabp, *wabm,
                *wacp, *wacm, *wbcp, *wbcm,
                *smwa, *smwb, *smwc,
                **qd, **qddag, **qha, *qhc, **qhc_discrete,
                Qd, *rhoha, *rhoda, *rhodb, *rhohc, Qha, Qp, Qhc,
                *tmp, *tmp2, *gha, *ghc, *gaa, *gab, *gbb,
                *etap, *etam, 
                *avg_rhoda, *avg_rhodb, *avg_rhoha, *avg_rhohc,
                *avg_rho_fld_np, *avg_rho_fld_np_c,
                *rho_surf, *surfH, *rho_exp_nr, *exp_nrH,
                *rho_fld_np_c, *rho_fld_np, *fld_npH,
                Hcur, *poly_bond_fft,
                shift_wp , *hhat, **tmp_sph, ***tmp_aniso, *tmp_iso,
                ***Gamma_aniso, *Gamma_iso, ***smwp_aniso, *smwp_iso,
                *exp_neg_smwp_iso, ***exp_neg_smwp_aniso;

#ifndef MAIN
extern
#endif
int Nx[Dim], NxT[Dim], M, do_CL, iter, 
    print_freq, itermax, sample_freq, sample_wait, save_avg_freq,
    update_scheme, ML, NxL[Dim], zstart, size, myrank, nprocs,
    N, Nda, Ndb, Nah, Nch_discrete, c_type,
    ic_flag[3], do_brent, do_film, np_type, n_exp_nr, do_fld_np, Nu,
    keep_fields, first_sim, gsd_it_max;

#ifdef PAR
#ifndef MAIN
extern
#endif
MPI_Status status;

#endif 

#ifndef MAIN
extern
#endif
fftw_complex *fmin0, *fmot0 ;

#ifndef MAIN
extern
#endif
fftw_plan fwd0, fwd1, fbk0, fbk1;


complex<double> integ_trapPBC(complex<double>*);
double gasdev2(void);
double ran2(void);
double get_k( int , double* );
double get_k_alias( int , double* );
int stack(int*);
int stack2(int); 
int stack3(int, int);
void unstack(int, int*);
void unstack_local(int, int*);
int unstack_stack(int);
void fft_fwd_wrapper(complex<double>*, complex<double>*);
void fft_bck_wrapper(complex<double>*, complex<double>*);

void write_avg_rho(char* nm, complex<double>*);
void write_avg_dat(char* nm, complex<double>*);


double get_k(int, double*, int);

void calc_gd( complex<double>* , double );
complex<double> calc_wpart(void);

void initialize_1(void);
void initialize_2(void);
void allocate(void);
void calc_poly_density( void );

void zero_average(complex<double>*, int);
void read_input(void);
double pbc_mdr2( double*, double*, double*) ;
void write_data_bin( char* , complex<double>* ) ;
void field_gradient( complex<double>* , complex<double>* , int );
complex<double> calc_H( void ) ;
void write_outputs( void ) ;
