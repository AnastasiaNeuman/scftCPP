#include "globals.h"
void calc_P_constants(int);
void read_resume_files(void);
void init_particles(void);
void sphere_init(void);
complex<double> integ_sphere(complex<double>**);
void calc_gaa(complex<double>* , double);
void calc_gbb(complex<double>* , double);
void calc_gab(complex<double>* , double);
double get_r(int , double*) ;
double dot_prod(double[Dim], double[Dim]);
double cross_prod(double[Dim], double[Dim]);
void write_data(char*, complex<double>*); // Added for debugging
void write_outputs(void);
void init_negative_k(void) ;
void init_fields(void);
void init_Gamma_sphere(void);
void init_Gamma_rod(void);
void explicit_nanorod(double, double, double, double[Dim], double[Dim]);
void explicit_nanosphere(double, double, double[Dim]);

void initialize_1() {

  if (myrank == 0) printf("\n---Initialization 1---\n");

  I = complex<double>(0.0,1.0);
  N = Nda + Ndb;
  fD = ( N>0 ? double(Nda)/double(N) : 0.0 );
  if (a_smear == -1.0) a_squared = 1.0 / double(N-1);
  else a_squared = a_smear * a_smear;

  // Initialize M, V, dx
  M = 1; V = 1.0;
  for (int i=0; i<Dim; i++) {
    M *= Nx[i];
    V *= L[i];
    dx[i] = L[i] / double(Nx[i]);
  }

  if (myrank == 0) {
    printf("V = %lf\n", V);
    printf("N = %d\n", N);
    printf("fD = %lf\n", fD);
    printf("a^2 = %lf\n", a_squared);
    printf("---Initialization 1 complete---\n\n");
    fflush(stdout);
  }

  n_samples = 0.0;
}

// Initializes fields and other variables
void initialize_2() {
  int i, j, k;
  double k2 , kv[Dim];
  double mdr2, mdr , r1[Dim], r2[Dim], dr[Dim] , x[Dim];

  if (myrank == 0) printf("---Initialization 2---\n");

  // Initialize hhat
  for (i=0; i<ML; i++) {
    k2 = get_k(i , kv);
    hhat[i] = exp(-k2 * a_squared / 2.0) ;
  }

  // Zero all densities //
  for (i=0; i<ML; i++) {
    rhoda[i] = rhodb[i] = rhoha[i] = rho_fld_np_c[i] = rho_fld_np[i] = 0.0;
    rho_surf[i] = surfH[i] = rho_exp_nr[i] = exp_nrH[i] = 0.0;
  }

  // Initialize confining wall //
  if (do_film) {
    for (i=0; i<ML; i++) {

      // Distance from current position to z=0 plane
      r1[Dim-1] = 0.0;
      get_r(i, x);
      for (int j=0; j<Dim-1; j++)
        r1[j] = x[j];
      mdr2 = pbc_mdr2(x, r1, dr);
      mdr = sqrt(mdr2);

      rho_surf[i] = 0.5 * (1.0 - erf((fabs(mdr) - wallT) / wallXi));

    }
    write_data_bin("rho_surf", rho_surf);
    fft_fwd_wrapper(rho_surf, surfH);
  }

  // Loop over nanorods and add nanorod density to rho_exp_nr
  for (i=0; i<n_exp_nr; i++) {
    if (np_type == 1)
      explicit_nanosphere(R_nr, xi_nr, exp_nr_c[i]);
    else if (np_type == 2)
      explicit_nanorod(L_nr, R_nr, xi_nr, exp_nr_c[i], exp_nr_u[i]);
    else {
      printf("np_type must be 1 (spheres) or 2 (rods) if n_exp_nr>0 and/or"
             " do_fld_np is set to 1 (true). Exiting...\n");
      exit(1);
    }
  }
  // If there are nanorods (and/or nanospheres) write data and take FT
  if (n_exp_nr>0) {
    write_data_bin("rho_exp_nr",  rho_exp_nr);
    fft_fwd_wrapper(rho_exp_nr, exp_nrH);
  }

  ////////////////////////
  // Define the volumes //
  ////////////////////////

  double V_nps, V_1_exp_np, V_exp_nps, V_fld_nps ;
  // Total volume of all explicit nanoparticles
  V_exp_nps = real(integ_trapPBC(rho_exp_nr));
  // "Free" Volume (everything excluding walls)
  Vf = V - real(integ_trapPBC(rho_surf)) - V_exp_nps ;
  // Volume of A and B (exclude C)
  V_ab = Vf;
  if (do_fld_np) {
    // Error-check np_frac
    if (np_frac < 0.0) {
      if (myrank == 0) printf("Negative np_frac is not allowed\n");
      exit(1);
    }
    else if (np_frac == 0.0 && myrank == 0)
      printf("WARNING: You turned on field-based nanoparticles but set "
             "np_frac to 0. That's just plain silly.\n");
    // Total volume of all explicit and field-based nanoparticles
    V_nps = np_frac * V_ab;
    // Total volume of all field-based nanoparticles (total np vol minus
    // explicit np vol)
    V_fld_nps = V_nps - V_exp_nps;
  }
  else {
    // If not doing field-based nps, total volume of all nanoparticles is just
    // the total volume of all explicit nanoparticles
    V_nps = V_exp_nps;
    // Set np_frac = volume of explicit nanoparticles / free volume
    np_frac = V_nps /(V_ab + V_nps);
  }
  
  // Volume of just one explicit nanoparticle if applicable
  if (n_exp_nr > 0)
    V_1_exp_np = V_exp_nps / double(n_exp_nr);
  else
    V_1_exp_np = 0.0;

 // double phiAH = zAH / ( 1.0 + zAH + zCH ) ;
 // double phiCH = zCH / ( 1.0 + zAH + zCH ) ;
    
  // Number of molecules (or nanoparticles) of each component
  if (phiAH > 0.0)  {
    nD = C * Vf * (1-phiAH); // # of diblock chains
    nAH = C * Vf * phiAH*(double(N)/double(Nah)); // # of homopolymer chains
  }
    else
      nAH = 0.0;
  
  if (c_type == 0 && phiCH > 0.0) {
    Nch = Nch_discrete;
    nCH = C * Vf * phiCH * double(N) / double(Nch_discrete);
    if (myrank==0) {
      printf("Nch=%d, phiCH=%lf, nCH=%lf\n", Nch_discrete, phiCH, nCH);
      fflush(stdout);
    }
  }
  else if (c_type == 1 && phiCH > 0.0) {
    Nch_discrete = 0;
    Nch = C * Vf * phiCH * double(N) / nCH;
    if (myrank==0) {
      printf("Nch=%lf, phiCH=%lf, nCH=%lf\n", Nch, phiCH, nCH);
      fflush(stdout);
    }
  }
  else {
    Nch = nCH = 0.0;
    Nch_discrete = 0;
  }

  // Define rho0
  rho0 = N * C;

  // Initialize fields
  if (keep_fields && !first_sim) {
    // Do nothing if we're keeping the fields and it's not the first
    // simulation. The final fields from the previous simulation carry on to
    // the next one. (first_sim gets marked false at end of simulate() method)
    if (myrank == 0) printf("Keeping fields, not reinitializing\n");
  }
  else {
    init_fields();
    if (myrank == 0) {
      if (first_sim) printf("Initialized fields\n");
      else printf("Reinitialized fields\n");
    }
  }
  
  // Initialize Gamma (and V_1_fld_np) if doing field-based nps
  if (do_fld_np) {
    if (np_type == 1) {
      init_Gamma_sphere();
    }
    else if (np_type == 2) {
      init_Gamma_rod();
    }
    else {
      printf("np_type=%d is not an acceptable value. Use 1 for spheres and 2 "
             "for rods.\n", np_type);
      exit(1);
    }
  }
  
  // Number of nanoparticles
  if (n_exp_nr == 0 && !do_fld_np)
    nP = 0.0;
  else if (do_fld_np)
    // Define nP based on nP*Vp/V = rho0*np_frac
    nP = rho0 * np_frac * V_ab / V_1_fld_np;
  else
    nP = V_nps / V_1_exp_np;

  if (do_fld_np) {
    // Number of field-based nanoparticles = total - explicit
    nFP = nP - n_exp_nr;
    if (nFP < 0.0) {
      if (myrank == 0)
        printf("The nanoparticle volume fraction you're using is lower than "
               "the volume fraction taken up by the explicit nanoparticles\n");
      exit(1);
    }
  }
  else
    nFP = 0.0;

  if (myrank == 0) {
    printf("Total V_segment actual: %lf\n", nD*N + nAH*Nah + nCH*Nch);
    printf("Total V_segment theoretical: %lf\n", C * Vf * double(N));
    printf("Total V_segment actual db: %lf\n", nD*N);
    printf("Total V_segment theoretical db: %lf\n",
            C * Vf * double(N) * (1-phiAH) );
    printf("Total V_segment actual ha: %lf\n", nAH*Nah);
    printf("Total V_segment theoretical ha: %lf\n",
           C * Vf * double(N) * phiAH);
    printf("Total V_segment actual hc: %lf\n", nCH*Nch);
    printf("Total V_segment theoretical hc: %lf\n",
           C * Vf * double(N) * phiCH);
    printf("V - Vf: %lf\n", V - Vf);
    printf("rho0 = %lf\n", rho0);
    printf("V_1_exp_np = %lf\n", V_1_exp_np);
    printf("V_1_fld_np = %lf\n", V_1_fld_np);
    printf("nD = %lf\n", nD);
    printf("nAH = %lf\n", nAH);
    printf("nP = %lf\n", nP);
    printf("n_exp_nr = %d\n", n_exp_nr);
    printf("nFP = %lf\n", nFP);
    printf("nCH = %lf\n", nCH);
    printf("Nch = %lf\n", Nch);
  }

  // Initialize Debye functions
  calc_gaa(gaa, fD);
  calc_gbb(gbb, fD);
  calc_gab(gab, fD);
  
  if (Nah > 0) calc_gd( gha, double(Nah-1)/double(N-1) );
  if (nCH > 0.0) calc_gd( ghc, double(Nch-1)/double(N-1) );

  ///////////////////////////////
  // Set up bonding potentials //
  ///////////////////////////////
  if (myrank == 0) 
    printf("Setting up Gaussian bonds\n");
  for (i=0; i<ML; i++) {
    k2 = get_k(i, kv);
    poly_bond_fft[i] = exp( -k2/double(N-1) ) ;
  }

  iter = 0;

  if (myrank == 0) {
    printf("---Initialization 2 complete---\n\n");
    fflush(stdout);
  }
} // initialize_2

// Initializes fields based on pattern specified in input file, then
// overwrites with restart file values if present.
void init_fields() {
  int sincos_dir;
  double x[Dim];
  // Initialize fields based on custom pattern specified in input file
  for (int i=0; i<ML; i++) {
    // Initialize c fields
    int i_global = unstack_stack(i);
    int nn[Dim];

    unstack(i_global, nn);

    // Get x, the current position
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j]) * dx[j];

    // Initialize C layer on top of simulation box
    wacp[i] = wbcp[i] = 0.0;
    if (nCH > 0.0) {

      if ( x[Dim-1] > L[Dim-1] / 2.0 ) {
        wacm[i] = wbcm[i] = -abs(ic_pre[0]) ;
      }
      else {
        wacm[i] = wbcm[i] = abs(ic_pre[0]) ;
      }
//      if (real(rho_surf[i] + rho_exp_nr[i]) > 0.1)
//        wacm[i] = wbcm[i] = 0.0;
//      else {
//        double L_film = L[Dim-1];
//        double x_interface;
//        if (do_film) {
//          L_film -= 2.0 * wallT;
//          x_interface = wallT + L_film * (1.0 - phiCH);
//        }
//        else
//          x_interface = L_film * (1.0 - phiCH);
//        
//        if (x[Dim-1] > x_interface)
//          wacm[i] = wbcm[i] = -abs(ic_pre[0]);
//        else
//          wacm[i] = wbcm[i] = 0.0;
//      }
    }
    else
      wacm[i] = wbcm[i] = 0.0;

    wpl[i] = 10.0*I*C/double(N) ;
    wabp[i] = 0.0;//-I*chiABN;
   // wabm[i] =2.0*(ran2() - .5);

    // Sine and cosine pattern; trying this to force macrophase separation 
    if (ic_flag[0] < -1) {
      if (real(rho_surf[i] + rho_exp_nr[i]) > 0.5)
        // 0 if there's some wall or nanorod present
        wabm[i] = 0.0;
      else if (real(wacm[i]) < 0.0)
        wabm[i] = 0.0;
      else { // set custom wabm pattern
        wabm[i] = 1.0;
        for (int j=0; j<2; j++) { // initial condition flag sets
          sincos_dir = ic_dir[j];
          if (ic_flag[j] == -2) // cosine flag
            wabm[i] *= ic_pre[j]*tanh( cos(2.0 * PI * ic_period[j]
                  * x[sincos_dir] / L[sincos_dir])/0.2 );
          else if (ic_flag[j] == -3) // sin flag
            wabm[i] *= ic_pre[j]*tanh( sin(2.0 * PI * ic_period[j]
                  * x[sincos_dir] / L[sincos_dir])/0.2 );
          else if (ic_flag[j] == -4) // field = 1 everywhere
             wabm[i] = 1;
          else if (ic_flag[j] < -4) { // -4 flag is just a factor of 1
            printf("Initial condition flags smaller than -4 not supported\n");
            exit(1);
          }
        } // j (initial condition flag sets)

        //if (real(wabm[i]) > 0.0 && nCH > 0.0)
        //  wacm[i] = wabm[i];
        //else if (real(wabm[i]) < 0.0 && nCH > 0.0)
        //  wbcm[i] = -wabm[i];

      } // end setting custom wabm pattern
    } // if ic_flag[0] < -1
    // Random fields
    else if (ic_flag[0] == -1) {
      wabm[i] = ic_pre[0] * ran2();
     // wacm[i] = ic_pre[0] * ran2();
     // wbcm[i] = ic_pre[0] * ran2();
    }
//    // Zero fields
    else 
      wabm[i] = 0.0 ;
    
    get_r(i, x);
   // wabm[i] = 1.2*(ran2()-0.5) + sin(16.0*PI*x[0]/L[0]) - 2.0*tanh( (cos(4.0*PI*x[Dim-1]/L[Dim-1])-0.875)/0.2 ) ;
   // wbcm[i] = 2.0*(ran2()-0.5) ;
  } // i for loop 

  // Overwrite fields with restart files if available
  read_resume_files();
}

// Initialize 1 explicit nanorod by adding appropriate densities to the
// rho_exp_nr array
// 
// INPUTS:
//    len         Nanorod length
//    rad         Nanorod radius
//    xi          Thickness of interface between nanorod and surroundings
//    rel_center  Position vector for the nanosphere center. Relative to the
//                  length of each dimension (values from 0 to 1)
//    u           Nanorod orientation vector in cartesian coordinates. Points
//                  along long axis of nanorod. Must be unit vector.
void explicit_nanorod(double len, double rad, double xi,
                      double rel_center[Dim], double u[Dim]) {
  int i_global;
  int nn[Dim];
  double u_dot_r, u_cross_r;
  double center[Dim], dr[Dim], x[Dim];

  // Multiply center (values between 0 and 1) by L to get the absolute
  // center point of the particle
  for (int i=0; i<Dim; i++)
    center[i] = rel_center[i] * L[i];

  for (int i=0; i<ML; i++) {
    // Get nn, the int array showing which point in the each dimension we're at
    i_global = unstack_stack(i);
    unstack(i_global, nn);

    // Get x, the current position
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j]) * dx[j];

    // Get dr, the distance vector from nanorod center to current position
    pbc_mdr2(x, center, dr);

    // Compute explicit nanorod density
    u_dot_r = dot_prod(u, dr);
    u_cross_r = cross_prod(u, dr);
    rho_exp_nr[i] += 0.25 * erfc( (fabs(u_dot_r)-0.5*len)/xi )
                          * erfc( (fabs(u_cross_r)-rad)/xi );
  }
} // explicit_nanorod

// Initialize 1 explicit nanosphere by adding appropriate densities to the
// rho_exp_nr array
// 
// INPUTS:
//    rad         Nanosphere radius
//    xi          Thickness of interface between nanosphere and surroundings
//    rel_center  Position vector for the nanosphere center. Relative to the
//                  length of each dimension (values from 0 to 1)
void explicit_nanosphere(double rad, double xi, double rel_center[Dim]) {
  int i_global;
  int nn[Dim];
  double dr2, dr_abs;
  double center[Dim], dr[Dim], x[Dim];

  // Multiply rel_center (values between 0 and 1) by L to get the absolute
  // center point of the particle
  for (int i=0; i<Dim; i++)
    center[i] = rel_center[i] * L[i];

  for (int i=0; i<ML; i++) {
    // Get nn, the int array showing which point in the each dimension we're at
    i_global = unstack_stack(i);
    unstack(i_global, nn);

    // Get x, the current position
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j]) * dx[j];

    // Get dr_abs, the magnitude of the distance from nanorod center to current
    // position
    dr2 = pbc_mdr2(x, center, dr);
    dr_abs = sqrt(dr2);

    // Compute explicit nanosphere density (excluding rho0)
    rho_exp_nr[i] += 0.5 * erfc( (dr_abs-rad)/xi );
  }
} // explicit_nanosphere
// Initialize Gamma_iso for spherical (isotropic) particles
void init_Gamma_sphere() {
  if (R_nr <= 0.0 && myrank == 0) {
    printf("R_nr=%lf is invalid. Try again\n", R_nr);
    exit(1);
  }

  int nn[Dim];
  double dr[Dim], x[Dim], origin[Dim] = { 0.0 };

  for (int i=0; i<ML; i++) {
    // Get distance from nanorod center to current position (dr)
    int i_global = unstack_stack(i);
    unstack(i_global, nn);
    for (int j=0; j<Dim; j++)
      x[j] = double(nn[j])*dx[j];
    double dr2 = pbc_mdr2(x, origin, dr);
    double dr_abs = sqrt(dr2);

    // Compute nanosphere density (including rho0)
    Gamma_iso[i] = 0.5 * rho0 * erfc( (dr_abs-R_nr)/xi_nr );
    
    // Multiply in factor of V because this will speed up convolution
    // calculations later, which all need a factor of V
    Gamma_iso[i] *= V;
  }

  // Calculate volume of 1 nanosphere for use later. The 1/V factor is to
  // cancel out the extra factor of V multiplied in above.
  V_1_fld_np = real(integ_trapPBC(Gamma_iso)) / V;

  // Fourier transform Gamma_iso and leave it that way. It's only used for
  // convolutions which are all done in k-space anyway.
  fft_fwd_wrapper(Gamma_iso, Gamma_iso);
} // END init_Gamma_sphere

// Initialize Gamma_aniso for nanorods (anisotropic)
void init_Gamma_rod() {
  if (Dim < 3) {
    if (myrank == 0)
      printf("Field-based anisotropic nanoparticles currently only supported"
             " in 3D\n");
    exit(1);
  }

  if (L_nr <= 0.0 || R_nr <= 0.0) {
    if (myrank == 0) {
      printf("L_nr=%lf, R_nr=%lf is invalid. Try again\n", L_nr, R_nr);
    }
    exit(1);
  }

  // Read legendre abscissa and weights once with processor 0
  if (myrank == 0) {
    sphere_init();
  }     
  // Broadcast the values to all other processors
  MPI_Bcast(theta,           Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(theta_weights,   Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(phi,           2*Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(phi_weights,   2*Nu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // Compute nanorod density (including rho0)
  for (int i=0; i<Nu; i++) {
    for (int j=0; j<2*Nu; j++) {
      for (int k=0; k<ML; k++) {
        double r[Dim], u[Dim], u_dot_r, u_cross_r;
        get_r(k, r);
        // Convert from spherical coords to x, y, and z
        u[0] = sin(theta[i]) * cos(phi[j]);
        u[1] = sin(theta[i]) * sin(phi[j]);
        u[2] = cos(theta[i]);
        // Get absolute values of dot and cross products of u and r
        u_dot_r = abs(dot_prod(u, r));
        u_cross_r = abs(cross_prod(u, r));
        // Nanorod density equation
        Gamma_aniso[i][j][k] = 0.25 * rho0
          * erfc( (u_dot_r-0.5*L_nr) / xi_nr )
          * erfc( (u_cross_r-R_nr) / xi_nr );
        // Multiply by V so we don't have to do it more expensively during
        // convolution
        Gamma_aniso[i][j][k] *= V;
      } // k
      tmp_sph[i][j] = integ_trapPBC(Gamma_aniso[i][j]);
      // Fourier transform Gamma_aniso and leave it that way. It's only used
      // for convolutions which are all done in k-space anyway.
      fft_fwd_wrapper(Gamma_aniso[i][j], Gamma_aniso[i][j]);
    } // j
  } // i

  // Calculate average volume of 1 nanorod based on 
  // Vp = integral ( dr * 1/(4PI) * integral( du * Gamma ) )
  // The extra V in the denominator is to cancel out the extra V included in
  // Gamma to speed up convolution later
  V_1_fld_np = real( integ_sphere(tmp_sph) ) / (4 * PI * V);

} // init_Gamma_rod
