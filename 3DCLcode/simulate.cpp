#include "globals.h"
#include "r2.h"
void write_outputs(void);
void accumulate_all_averages(void);
void save_averages(void);
void update_1s(void);
void update_Euler(void);
void write_data(char*, complex<double>*); // for debugging

// This is the routine that is essentially the main routine in a code that
// doesn't use Brent's method. Int calculates the equilibrium structure and
// energy and returns the real part of the hamiltonian.
double simulate() {
  complex<double> Hcur, Ho, H;
  double error;
  FILE *otp;
  if ( myrank == 0 ) {
    otp = fopen("data.dat", "w");
    if (otp == NULL) {
      printf("Failed to open data.dat!\n");
      exit(1);
    }
  }

  if (myrank == 0) {
    printf("------Starting simulation for L[0]=%lf, L[1]=%lf------\n",
           L[0], L[1]);
  }

  // Initialize variables and fields
  initialize_1();
  initialize_2();
 
#ifdef PAR
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  calc_poly_density();


  if (do_fld_np) {
    // Field nanoparticle sanity check
    complex<double> np_check = integ_trapPBC(rho_fld_np_c)
                               + double(n_exp_nr);
    complex<double> npVp_check = integ_trapPBC(rho_fld_np)
                                 + integ_trapPBC(rho_exp_nr) * rho0;
    complex<double> C_check_easy = ( nD + nAH*Nah/double(N)
                                     + nP*V_1_fld_np/double(N) ) / V_ab;
    complex<double> C_check_hard = ( nD + nAH*Nah/double(N)
                                     + npVp_check/double(N) ) / V_ab;
    if (myrank == 0) {
      printf("Field-based nanoparticle calculations sanity check:\n");
      printf("nP=%lf, np_check=%lf\n", nP, real(np_check));
      printf("nFP = %lf, n_exp_nr = %d, nFP + n_exp_nr = %lf =? nP\n",
              nFP,       n_exp_nr,      nFP + n_exp_nr );
      printf("nP*Vp = %lf, npVp_check=%lf\n",
              nP*V_1_fld_np, real(npVp_check) );
      printf("C = %lf, C_check_easy=%lf, C_check_hard=%lf\n\n",
              C,  real(C_check_easy), real(C_check_hard) );
    }
  }

  if (nCH > 0.0) {
    // C homopolymer sanity check
    double ncNc_check = real( integ_trapPBC(rhohc) );
    if (myrank == 0)
      printf("ncNc_check = %lf, nc*Nc = %lf\n", ncNc_check, nCH * Nch); 
  }
  double smrhodadb = real( integ_trapPBC(rhoda)+integ_trapPBC(rhodb) );
  double smrhoha = real( integ_trapPBC(rhoha) );

  if (myrank == 0) {
    printf("Initial densities calculated!\n");
    printf("Segment counts:\n");
    printf("nD * N = %lf integ(rhoda + rhodb) = %lf\n",
           nD * N, smrhodadb);
    printf("nAH * Nah = %lf integ(rhoha) = %lf\n",
           nAH * Nah, smrhoha);

    printf("V*rho0: %lf total segs: %lf\n", V*rho0, nD*N + nAH*Nah ) ;
    fflush(stdout);
  }

  Ho = calc_H();

  if (myrank == 0) printf("Starting H: %lf\n\n", real(Ho));

  write_outputs();

  if (myrank == 0) printf("---Entering main loop---\n"); 


  ///////////////
  // MAIN LOOP //
  ///////////////
  FILE *zot ;
  if ( myrank == 0 ) {
    zot = fopen("z_phi.dat","w") ;
    fprintf(zot, "# chiABN=%lf\n", chiABN ) ;
  }
  printf("V: %lf Vf: %lf\n", V, Vf ) ;
 // while ( ( zah_factor > 1.0 && zAH < zah_final ) ||
  //      ( zah_factor < 1.0 && zAH > zah_final ) ) {

    for (iter=0; iter<=itermax; iter++) {
      if (update_scheme == 0)
        update_Euler();
      else
        update_1s();
 
      if (do_CL && iter >= sample_wait && iter % sample_freq == 0) 
        accumulate_all_averages();
 
      ////////////
      // OUTPUT //
      ////////////
      if (iter % print_freq == 0) {
        Ho = Hcur;
        H = Hcur = calc_H();
        if (Hcur != Hcur) { // If Hcur is NaN
          printf("Crashed! iteration %d\n", iter);
          printf("H: %lf Qd: %lf\n", real(H), real(Qd));
          write_data_bin("crashed.wpl", wpl);
          write_data_bin("crashed.wabp", wabp);
          write_data_bin("crashed.wabm", wabm);
          exit(1);
        }
        error = abs(H - Ho) / V / double(print_freq);

        double dtmp;

        if ( myrank == 0 ) 
          printf("Iter: %d, H=%lf", iter, real(H) );

        if (do_CL)
          printf(" + i%lf", imag(H) );
        
        if (nD > 0.0) {
          dtmp = real(integ_trapPBC(rhoda)+integ_trapPBC(rhodb))/Vf/rho0 ;
          if ( myrank == 0 ) printf(", -Qd=%lf phiAD=%lf", real(-Qd), dtmp );
        }

        if (nAH > 0.0) {
          dtmp = real(integ_trapPBC(rhoha))/Vf/rho0 ;
          if ( myrank == 0 ) printf(", -Qha=%lf phiAH=%lf", real(-Qha), dtmp ) ;
        }

        if (nCH > 0.0) {
          dtmp = real(integ_trapPBC(rhohc))/Vf/rho0 ;
          if ( myrank == 0 ) printf(", -Qhc=%lf phiCH=%lf", real(-Qhc), dtmp ) ;
        }

        if (do_fld_np)
          printf(", -log(Qp)=%lf", real(-log(Qp)+smwp_min) );
        
        if ( myrank == 0 ) printf(", err=%1.1e\n", error);

        dtmp = real(integ_trapPBC(rhoda)+integ_trapPBC(rhodb)
                   +integ_trapPBC(rhoha)+integ_trapPBC(rhohc)) ;

        if ( myrank == 0 ) printf("Segment total: %lf  V*rho0: %lf\n", dtmp, Vf*rho0) ;
        fflush(stdout);

        if (myrank == 0) {
          fprintf(otp, "%d %5.6lf %1.3e %5.6lf %5.6lf %5.6lf %1.3e", 
                  iter, real(H), imag(H), real(-log(Qd)), real(-log(Qha)),
                  real(-log(Qp)+smwp_min), error);
          fprintf(otp, "\n");
          fflush(otp);
        }
        write_outputs();
      } // Output
 
      if (do_CL && iter>sample_wait && (iter-sample_wait)%save_avg_freq == 0)
        save_averages();
 
      if (!do_CL && iter > 25 && error < error_tol) {     
        if (myrank == 0) {
          printf("Tolerance reached. Error = %.4e\n", error);
          printf("---Main loop complete---\n\n");
        }
        break;
      }
 
    }// Main Loop
    // Below is where phi is calculated?
    // Need to change phi to an input rather than n
    // below I calculate mu for A and B 
    double muA = real(-log(Qha));
    double muB = real(-log(Qd));
    double dtmp = real(integ_trapPBC(rhoha))/Vf/rho0 ;
    complex<double> dU;
    complex<double> TS;
    complex<double> term1();
    dU = -1.0*term1();
    TS = real(H) - dU;
    if ( myrank == 0 ) {
      fprintf(zot, "muA muB p H U S\n");
      fprintf(zot, "%lf %lf %lf %lf %lf %lf\n", muA, muB, real(H)/V, real(H),real(dU),real(TS)) ;
      fflush(zot) ;
    }
 // zAH *= zah_factor ;

    if ( myrank == 0 ) {
      printf("*****************\n");
   //   printf("** zAH = %1.3lf **\n", zAH);
      printf("*****************\n");
    }
  //zah_loop
  // Close output stream
  if ( myrank == 0 ) {
    fclose(otp);
    fclose(zot) ;
  }

  double H_over_V = real(H) / V;
  if (first_sim || H_over_V < min_H_over_V) {
    min_H_over_V = H_over_V;
    L_ideal = L[0];
  }

  if (myrank == 0) {
    // Output results to standard output
    printf("------Completed L[0]=%lf simulation. H=%lf and H/V=%lf------\n\n",
           L[0], real(H), H_over_V);
    if (do_brent) {
      // Open brent.dat
      brent_otp = fopen("brent.dat", "a");
      if (brent_otp == NULL) {
        printf("Failed to open brent.dat!\n");
        exit(1);
      }

      // Output length, H, and H/V data to brent.dat
      for (int i=0; i<Dim; i++)
        fprintf(brent_otp, "%5.6lf ", L[i]);
      fprintf(brent_otp, "%5.6lf %5.6lf\n", real(H), H_over_V);
      fclose(brent_otp);

      printf("Current global minimum: L=%lf, H/V=%lf\n\n", L_ideal, min_H_over_V);
    }
    fflush(stdout);
  }

  // The first simulation is complete by this point so set first_sim to 0
  // no matter what it's current value is
  first_sim = 0;

  return H_over_V;
}
complex<double> term1() {
      int i;
      for ( i=0; i<ML ; i++ )
        tmp[i] = (wabp[i] * wabp[i] + wabm[i] * wabm[i] ) * C/chiABN;
      return integ_trapPBC( tmp );
}
