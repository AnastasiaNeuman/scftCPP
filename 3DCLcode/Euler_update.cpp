#include "globals.h"
void generate_1s_noise( complex<double>* , double ) ;

void update_Euler( ) {

  int i;
  complex<double> evar , F, numer;

  if (nD > 0.0) {
    fft_fwd_wrapper(rhoda, rhoda);
    fft_fwd_wrapper(rhodb, rhodb);
  }
   
  if (nAH > 0.0) 
    fft_fwd_wrapper(rhoha, rhoha); 
 

  fft_fwd_wrapper(wpl, wpl);

  if (chiABN > 0.0) {
    fft_fwd_wrapper(wabm, wabm);
    fft_fwd_wrapper(wabp, wabp);
  }

  if (do_CL) 
    generate_1s_noise(etap, lam_pl);

  for (i=0; i<ML; i++) {
    // Update w+ field //
    if (i == 0 && myrank == 0) 
      evar = 1.0 ;
    else
      evar = 0.0 ;
    F = (kappaN <= 0.0 ? 0.0 : C * wpl[i] / kappaN) 
        + I * C * (surfH[i] + exp_nrH[i] - evar)
        + I * rho_fld_np[i] / double(N)
        + I * hhat[i] / double(N) * (rhoha[i] + rhoda[i] + rhodb[i]);
    numer = wpl[i] - lam_pl * F;
    if (do_CL) 
      numer += etap[i];
    wpl[i] = numer;
  }
  
  // Update AB fields //
  if (chiABN > 0.0) {
    if (do_CL) {
      generate_1s_noise(etap, lam_pl);
      generate_1s_noise(etam, lam_mi);
    }

    for (i=0; i<ML; i++) {
      // AB+ //
      F = 2.0 * C * wabp[i] / chiABN
          + I * rho_fld_np[i] / double(N)
          + I * hhat[i] / double(N) * (rhoda[i] + rhodb[i] + rhoha[i]);
      numer = wabp[i] - lam_pl * F;
      if ( do_CL )
        numer += etap[i];
      wabp[i] = numer;

      // AB- //
      F = 2.0 * C * wabm[i] / chiABN
          + rho_fld_np[i] / double(N)
          + hhat[i] / double(N) * (rhodb[i] - rhoda[i] - rhoha[i]);
      numer = wabm[i] - lam_mi * F;
      if (do_CL) 
        numer += etam[i];
      wabm[i] = numer;
    }
  }
    
  else {
    for (i=0; i<ML; i++)
      wabp[i] = wabm[i] = 0.0 ;
  }

  fft_bck_wrapper(wpl, wpl);

  if (chiABN > 0.0) {
    fft_bck_wrapper(wabm, wabm);
    fft_bck_wrapper(wabp, wabp);
  }

  calc_poly_density() ;
}
