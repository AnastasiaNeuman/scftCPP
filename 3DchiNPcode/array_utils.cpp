#include "globals.h"
void calc_P_constants( int );
int stack_local( int* ) ;
void fft_fwd_wrapper(complex<double>*, complex<double>*, int);
void fft_bck_wrapper(complex<double>*, complex<double>*, int);
void zero_average(complex<double>*, int);
void initialize_averages( complex<double>* );
void accumulate_average_array( complex<double>* , complex<double>* ) ;
double get_k_global( int, double* ) ;

using namespace std;

void accumulate_all_averages() {

  if ( n_samples == 0.0 ) {
    initialize_averages( avg_rhoda );
    initialize_averages( avg_rhodb );
    initialize_averages( avg_rhoha );
    initialize_averages( avg_rho_fld_np );
    initialize_averages( avg_rho_fld_np_c );
  }

  accumulate_average_array( avg_rhoda , rhoda );
  accumulate_average_array( avg_rhodb , rhodb );
  accumulate_average_array( avg_rhoha , rhoha );
  accumulate_average_array( avg_rho_fld_np , rho_fld_np );
  accumulate_average_array( avg_rho_fld_np_c , rho_fld_np_c );

  n_samples += 1.0;
}

void initialize_averages( complex<double> *avg ) {
  for (int i=0; i<ML; i++)
    avg[i] = 0.0;
}

void accumulate_average_array( complex<double> *avg , complex<double>* dat ) {
  for (int i=0; i<ML; i++) 
    avg[i] += dat[i];
}

complex<double> get_avg(complex<double> *dat) {
  double real_sum_local = 0.0, real_sum_global = 0.0;
  double imag_sum_local = 0.0, imag_sum_global = 0.0;

  for (int i=0; i<ML; i++) {
    real_sum_local += real(dat[i]);
    imag_sum_local += imag(dat[i]);
  }

#ifdef PAR
  // Find the global sums
  MPI_Allreduce(&real_sum_local, &real_sum_global, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&imag_sum_local, &imag_sum_global, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#else
  // If no parallelization, the local sum is already the global sum
  real_sum_global = real_sum_local;
  imag_sum_global = imag_sum_local;
#endif

  return (real_sum_global + I * imag_sum_global) / double(M);
}

// Takes integer "id" in [0, ML] and finds the correct 
// unstacked integer value in [0, M]
int unstack_stack(int id) {

  int n[Dim];

  unstack_local(id, n );
  
  n[Dim-1] += zstart;

  return stack(n);
}

int stack_input(int x[Dim], int Nxx[Dim]) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nxx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nxx[1])*Nxx[0] );
}

// Stacks x using only local values
int stack_local(int x[Dim]) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*NxL[0]);
  else
    return  (x[0] + (x[1] + x[2]*NxL[1])*NxL[0] );
}

// Stacks vector x into 1D array index in [ 0, M ]
int stack( int x[Dim] ) {
  if (Dim==1)
    return x[0];
  else if (Dim==2)
    return (x[0] + x[1]*Nx[0]);
  else
    return  (x[0] + (x[1] + x[2]*Nx[1])*Nx[0] );
}

void unstack_input(int id, int nn[Dim], int Nxx[Dim]) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/Nxx[0];
    nn[0] = (id - nn[1]*Nxx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/Nxx[1]/Nxx[0];
    nn[1] = id/Nxx[0] - nn[2]*Nxx[1];
    nn[0] = id - (nn[1] + nn[2]*Nxx[1])*Nxx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}

// Receives index id in [ 0 , ML ] and turns it into 
// array nn[Dim]
void unstack_local(int id, int nn[Dim] ) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/NxL[0];
    nn[0] = (id - nn[1]*NxL[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/NxL[1]/NxL[0];
    nn[1] = id/NxL[0] - nn[2]*NxL[1];
    nn[0] = id - (nn[1] + nn[2]*NxL[1])*NxL[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}

// Receives index id in [0 , M ] and makes array
// nn[Dim] in [ 0 , Nx[Dim] ]
void unstack(int id, int nn[Dim] ) {

  if (Dim==1) {
    nn[0] = id;
    return;
  }
  else if (Dim==2) {
    nn[1] = id/Nx[0];
    nn[0] = (id - nn[1]*Nx[0]);
    return;
  }
  else if (Dim==3) {
    nn[2] = id/Nx[1]/Nx[0];
    nn[1] = id/Nx[0] - nn[2]*Nx[1];
    nn[0] = id - (nn[1] + nn[2]*Nx[1])*Nx[0];
  }
  else {
    cout << "Dim is goofy!" << endl;
    return;
  }
}

double get_r( int id , double r[Dim] ) {
  double r2 = 0.0;
  int i, id2, n[Dim];

  id2 = unstack_stack(id);

  unstack(id2, n);

  for ( i=0; i<Dim; i++) {
    r[i] = dx[i] * double( n[i] );
    
    if ( r[i] > L[i]/2.0 )
      r[i] -= L[i];
    else if ( r[i] <= -L[i]/2.0 )
      r[i] += L[i];
 
    r2 += r[i]*r[i];
  }

  return r2;
}

double get_k_alias( int id , double k[Dim] ) {

  double kmag = 0.0;
  int i, id2, n[Dim] , j , has_nyquist = 0;
  for ( i=0 ; i<Dim ; i++ )
    if ( Nx[i] % 2 == 0 )
      has_nyquist = 1;

  id2 = unstack_stack(id);

  unstack(id2, n);

  if ( Nx[0] % 2 == 0 && n[0] == Nx[0] / 2 )
    k[0] = 0.0 ;
  else if ( double(n[0]) < double(Nx[0]) / 2.)
   k[0] = 2*PI*double(n[0])/L[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/L[0];

  if (Dim>1) {
    if ( Nx[1] % 2 == 0 && n[1] == Nx[1] / 2 )
      k[1] = 0.0 ;
    else if ( double(n[1]) < double(Nx[1]) / 2.)
      k[1] = 2*PI*double(n[1])/L[1];
    else
      k[1] = 2*PI*double(n[1]-Nx[1])/L[1];
  }

  if (Dim==3) {
    if ( Nx[2] % 2 == 0 && n[2] == Nx[2] / 2 )
      k[2] = 0.0 ;
    else if ( double(n[2]) < double(Nx[2]) / 2.)
      k[2] = 2*PI*double(n[2])/L[2];
    else
      k[2] = 2*PI*double(n[2]-Nx[2])/L[2];
  }

  // Kills off the Nyquist modes
  if ( id2 != 0 && has_nyquist ) {
    for ( i=0 ; i<Dim ; i++ ) {
      if ( k[i] == 0.0 ) {
        for ( j=0 ; j<Dim ; j++ )
          k[j] = 0.0 ;
        kmag = 0.0;
        break;
      }
    }
  }
  
  for (i=0; i<Dim; i++)
    kmag += k[i]*k[i];

  return kmag;

}

// Receives index id in [ 0 , ML ] and returns 
// proper k-value, whether running in parallel or not
double get_k(int id, double k[Dim]) {

  double kmag = 0.0;
  int i, id2, n[Dim];

  id2 = unstack_stack(id);

  unstack(id2, n);

  if ( double(n[0]) < double(Nx[0]) / 2. )
   k[0] = 2*PI*double(n[0])/L[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/L[0];

  if (Dim>1) {
   if ( double(n[1]) < double(Nx[1]) / 2. )
    k[1] = 2*PI*double(n[1])/L[1];
   else
    k[1] = 2*PI*double(n[1]-Nx[1])/L[1];
  }

  if (Dim==3) {
   if ( double(n[2]) < double(Nx[2]) / 2. )
     k[2] = 2*PI*double(n[2])/L[2];
   else
     k[2] = 2*PI*double(n[2]-Nx[2])/L[2];
  }

  for (i=0; i<Dim; i++)
    kmag += k[i]*k[i];

  return kmag;

}

// Receives index id in [ 0 , ML ] and returns 
// proper k-value, whether running in parallel or not
double get_k_global(int id2, double k[Dim]) {

  double kmag = 0.0;
  int i, n[Dim];

  unstack(id2, n);

  if ( double(n[0]) < double(Nx[0]) / 2. )
   k[0] = 2*PI*double(n[0])/L[0];
  else
   k[0] = 2*PI*double(n[0]-Nx[0])/L[0];

  if (Dim>1) {
   if ( double(n[1]) < double(Nx[1]) / 2. )
    k[1] = 2*PI*double(n[1])/L[1];
   else
    k[1] = 2*PI*double(n[1]-Nx[1])/L[1];
  }

  if (Dim==3) {
   if ( double(n[2]) < double(Nx[2]) / 2. )
     k[2] = 2*PI*double(n[2])/L[2];
   else
     k[2] = 2*PI*double(n[2]-Nx[2])/L[2];
  }

  for (i=0; i<Dim; i++)
    kmag += k[i]*k[i];

  return kmag;

}

// Sets the average of tp to zero
void zero_average(complex<double>* tp) {

  int i;

  complex<double> integ;
  
  integ = integ_trapPBC(tp);

  integ *= (1.0 / V);

  for (i=0; i<M; i++)
    tp[i] -= integ;

}

// Allocates the memory for global variables
void allocate(void) {

  if (myrank == 0) printf("---Allocating memory---\n");

  int i;
  int Nf[Dim], alloc_size;

  for (i=0; i<Dim; i++) 
    Nf[i] = Nx[Dim-i-1];

  long long total_alloced = 0 ;

  ptrdiff_t Dm = Dim, Nfp[Dim], NxLtp, ztp;
  for (i=0; i<Dim; i++)
    Nfp[i] = Nf[i];

  size1 = fftw_mpi_local_size_many(Dm, Nfp, 1, 0, MPI_COMM_WORLD,
            &NxLtp, &ztp );

  NxL[Dim-1] = NxLtp;
  for (i=0; i<Dim-1; i++)
    NxL[i] = Nx[i];


  zstart = ztp;
  
  fmin0 = (fftw_complex*) fftw_malloc( size1 * sizeof(fftw_complex) );
  fmot0 = (fftw_complex*) fftw_malloc( size1 * sizeof(fftw_complex) );
  
  fwd0 = fftw_mpi_plan_dft(Dim, Nfp, fmin0, fmot0, 
      MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE );
  fbk0 = fftw_mpi_plan_dft(Dim, Nfp, fmin0, fmot0, 
      MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE );

  ML = 1;
  for (i=0; i<Dim; i++)
    ML *= NxL[i];
  
  total_alloced += size1*sizeof(fftw_complex)*2 ;

  // Set up the memory to allocate //
  alloc_size = NxL[0] ;

  for ( i=1 ; i<Dim ; i++ ) 
    alloc_size *= NxL[i];

  // Allocate the fields
  wpl = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wa = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wb = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wc = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  smwa = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  smwb = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  smwc = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wabp = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wabm = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wacp = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wacm = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wbcp = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  wbcm = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  total_alloced += alloc_size * sizeof(complex<double>) * 13;
  
  if (do_CL) {
    etap = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
    etam = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
    total_alloced += alloc_size * sizeof(complex<double>) * 2;
  }

  tmp  = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  tmp2 = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  total_alloced += alloc_size * sizeof(complex<double>) * 2;

  // Arrays required for spherical nanoparticles (no orientation dependence)
  if (do_fld_np && np_type == 1) {
    Gamma_iso = (complex<double>*)
                fftw_malloc(alloc_size*sizeof(complex<double>));
    smwp_iso = (complex<double>*)
               fftw_malloc(alloc_size*sizeof(complex<double>));
    exp_neg_smwp_iso = (complex<double>*)
                       fftw_malloc(alloc_size*sizeof(complex<double>));
    total_alloced += alloc_size * sizeof(complex<double>) * 3;
  }

  // Arrays required for nanorods (orientation dependent)
  if (do_fld_np && np_type == 2) {
    tmp_sph = (complex<double>**)
              fftw_malloc(Nu * sizeof(complex<double>*));
    tmp_aniso = (complex<double>***)
                fftw_malloc(Nu * sizeof(complex<double>**));
    Gamma_aniso = (complex<double>***)
                  fftw_malloc(Nu * sizeof(complex<double>**));
    smwp_aniso = (complex<double>***)
                 fftw_malloc(Nu * sizeof(complex<double>**));
    exp_neg_smwp_aniso = (complex<double>***)
                         fftw_malloc(Nu * sizeof(complex<double>**));
    for (i=0; i<Nu; i++) {
      tmp_sph[i] = (complex<double>*)
                   fftw_malloc(2*Nu * sizeof(complex<double>));
      tmp_aniso[i] = (complex<double>**)
                     fftw_malloc(2*Nu * sizeof(complex<double>*));
      Gamma_aniso[i] = (complex<double>**)
                       fftw_malloc(2*Nu * sizeof(complex<double>*));
      smwp_aniso[i] = (complex<double>**)
                      fftw_malloc(2*Nu * sizeof(complex<double>*));
      exp_neg_smwp_aniso[i] = (complex<double>**)
                              fftw_malloc(2*Nu * sizeof(complex<double>*));
      for (int j=0; j<2*Nu; j++) {
        tmp_aniso[i][j] = (complex<double>*)
                          fftw_malloc(alloc_size*sizeof(complex<double>));
        Gamma_aniso[i][j] = (complex<double>*)
                            fftw_malloc(alloc_size*sizeof(complex<double>));
        smwp_aniso[i][j] = (complex<double>*)
                           fftw_malloc(alloc_size*sizeof(complex<double>));
        exp_neg_smwp_aniso[i][j] = (complex<double>*)
                              fftw_malloc(alloc_size*sizeof(complex<double>));
      }
    }
    total_alloced += 2 * Nu * Nu * sizeof(complex<double>);
    total_alloced += 4 * 2 * Nu * Nu * alloc_size * sizeof(complex<double>);

    // Allocate theta and phi stuff
    theta = (double*) malloc(Nu * sizeof(double));
    theta_weights = (double*) malloc(Nu * sizeof(double));
    phi = (double*) malloc(2 * Nu * sizeof(double));
    phi_weights = (double*) malloc(2 * Nu * sizeof(double));
    total_alloced += 6 * Nu * sizeof(double);
  } // if (do_fld_np && np_type == 2)

  // Allocate the density operators
  rho_surf = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  rho_exp_nr = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  rho_fld_np_c = (complex<double>*)
                 fftw_malloc(alloc_size * sizeof(complex<double>));
  rho_fld_np = (complex<double>*)
               fftw_malloc(alloc_size * sizeof(complex<double>));
  surfH = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  exp_nrH = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  rhoda = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  rhodb = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  rhoha = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  rhohc = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  total_alloced += alloc_size * sizeof(complex<double>) * 10;

  if ( do_CL ) {
    avg_rhoda = (complex<double>*)
                fftw_malloc(alloc_size*sizeof(complex<double>));
    avg_rhodb = (complex<double>*)
                fftw_malloc(alloc_size*sizeof(complex<double>));
    avg_rhoha = (complex<double>*)
                fftw_malloc(alloc_size*sizeof(complex<double>));
    avg_rhohc = (complex<double>*)
                fftw_malloc(alloc_size*sizeof(complex<double>));
    avg_rho_fld_np = (complex<double>*)
                     fftw_malloc(alloc_size*sizeof(complex<double>));
    avg_rho_fld_np_c = (complex<double>*)
                       fftw_malloc(alloc_size*sizeof(complex<double>));
    total_alloced += alloc_size * sizeof(complex<double>) * 6;
  }

  // Debye functions
  gaa = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  gab = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  gbb = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  gha = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  ghc = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  total_alloced += alloc_size * sizeof(complex<double>) * 5;
  
  poly_bond_fft = (complex<double>*) fftw_malloc(alloc_size *
                                                 sizeof(complex<double>));
  hhat = (complex<double>*) fftw_malloc(alloc_size * sizeof(complex<double>));
  total_alloced += alloc_size * sizeof(complex<double>) * 2;

  qd = (complex<double>**) fftw_malloc(N * sizeof(complex<double>*));
  qddag = (complex<double>**) fftw_malloc(N * sizeof(complex<double>*));
  for (i=0; i<N; i++) {
    qd[i] = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
    qddag[i] = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  }
  total_alloced += alloc_size * sizeof(complex<double>) * N * 2;

  qha = (complex<double>**) fftw_malloc(Nah*sizeof(complex<double>*));
  for (i=0; i<Nah; i++) {
    qha[i] = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  }
  total_alloced += alloc_size * sizeof(complex<double>) *  Nah;

  qhc = (complex<double>*) fftw_malloc(alloc_size*sizeof(complex<double>));
  total_alloced += alloc_size * sizeof(complex<double>);

  if (c_type == 0) {
    qhc_discrete = (complex<double>**)
                   fftw_malloc(Nch_discrete*sizeof(complex<double>*));
    for (i=0; i<Nch_discrete; i++) {
      qhc_discrete[i] = (complex<double>*)
                        fftw_malloc(alloc_size*sizeof(complex<double>));
    }
    total_alloced += alloc_size * sizeof(complex<double>) *  Nch_discrete;
  }
 
  printf("Processor %d allocated %lf MB\n",
         myrank, double(total_alloced)/1.0E6);

#ifdef PAR
  // This just makes sure all processors finished before announcing allocation
  MPI_Barrier( MPI_COMM_WORLD );
#endif

  if (myrank == 0) {
    printf("---Memory allocation complete---\n\n"); 
    fflush(stdout);
  }
}

///////////////////////////////////////////////////////////////
// Calculates the gradient of a field in the "dir" direction //
// using spectral methods.  FFT, mult. by I*k[dir], iFFT     //
///////////////////////////////////////////////////////////////
void field_gradient( complex<double> *in , complex<double> *out , int dir ) {

  int i ;
  double kv[Dim] , k2 ;

  fft_fwd_wrapper( in , out );

  for ( i=0 ; i<ML ; i++ ) {
    k2 = get_k_alias( i , kv ) ;
    out[i] *= I * kv[ dir ] ;
  }

  fft_bck_wrapper( out , out ) ;

}

double pbc_mdr2( double x1[Dim] , double x2[Dim] , double dr[Dim] ) {
  int i;
  double mdr2 = 0.0 ;
  for ( i=0 ; i<Dim ; i++ ) {
    dr[i] = x1[i] - x2[i] ;

    if ( dr[i] >= 0.5*L[i] ) dr[i] -= L[i] ;
    else if ( dr[i] < -0.5*L[i] ) dr[i] += L[i] ;

    mdr2 += dr[i] * dr[i] ;
  }

  return mdr2 ;
}

// Returns u dot r where both u and r are in Cartesian coords 
double dot_prod( double u[Dim], double r[Dim] ) {

  int i;
  double sum = 0;

  for (i=0; i<Dim; i++)
    sum += u[i] * r[i];

  return sum;

}

// Returns magnitude of u x r where both u and r are in Cartesian coords 
double cross_prod( double u[Dim], double r[Dim] ) {

  double dum;
  double mag = 0;

  if (Dim==2) return fabs(u[0]*r[1] - u[1]*r[0]);
  else {
    // Add magnitude of i component squared
    dum = u[1]*r[2] - u[2]*r[1];
    mag += dum * dum;
    // Add magnitude of j component squared
    dum = u[2]*r[0] - u[0]*r[2];
    mag += dum * dum;
    // Add magnitude of k component squared
    dum = u[0]*r[1] - u[1]*r[0];
    mag += dum * dum;
    return sqrt(mag);
  }
}
