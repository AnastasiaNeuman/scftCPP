#include "globals.h"
void generate_1s_noise( complex<double>* , double ) ;

void update_1s( ) {

  int i ;
  complex<double> evar , F, A , numer, denom ;

  if (nD != 0.0) {
    fft_fwd_wrapper(rhoda, rhoda);
    fft_fwd_wrapper(rhodb, rhodb);
  }
   
  if (nAH != 0.0) 
    fft_fwd_wrapper(rhoha, rhoha); 

  if (nCH != 0.0)
    fft_fwd_wrapper(rhohc, rhohc);

  if (do_fld_np)
    fft_fwd_wrapper(rho_fld_np, rho_fld_np);

  // Update w+ field //

  if (do_CL) 
    generate_1s_noise(etap, lam_pl);

  fft_fwd_wrapper(wpl, wpl);

  for (i=0; i<ML; i++) {
    // The 1 at the k=0 mode takes care of the delta function arising from
    // taking a fourier transform of a constant
    if (i == 0 && myrank == 0)
      evar = 1.0 ;
    else
      evar = 0.0 ;
    
    F = (kappaN <= 0.0 ? 0.0 : C*wpl[i]/kappaN) 
      + I * C * (surfH[i] + exp_nrH[i] - evar)
      + I * rho_fld_np[i] / double(N)
      + I*hhat[i]/double(N) * (rhoha[i] + rhoda[i] + rhodb[i] + rhohc[i]);
    A = (kappaN <= 0.0 ? 0.0 : C/kappaN)
      + nD * double(N) * hhat[i] * hhat[i] * (gaa[i] + 2.0 * gab[i] + gbb[i]) / V
      + nAH * double(Nah * Nah) / double(N) * hhat[i] * hhat[i] * gha[i] / V
      + nCH * double(Nch * Nch) / double(N) * hhat[i] * hhat[i] * ghc[i] / V;
    numer = wpl[i] - lam_pl * (F - A * wpl[i]);
    if (do_CL) 
      numer += etap[i];
    denom = 1.0 + lam_pl * A;
    wpl[i] = numer / denom; 
  }

  fft_bck_wrapper( wpl , wpl ) ;
  
  // Update AB fields //
  if (chiABN > 0.0) {
    fft_fwd_wrapper(wabm, wabm); 
    fft_fwd_wrapper(wabp, wabp); 

    if (do_CL) {
      generate_1s_noise(etap, lam_pl);
      generate_1s_noise(etam, lam_mi);
    }

    for ( i=0 ; i<ML ; i++ ) {
      // AB+ //
      F = 2.0 * C * wabp[i] / chiABN
        + I * rho_fld_np[i] / double(N)
        + I * hhat[i] / double(N) * (rhoda[i] + rhodb[i] + rhoha[i]);
      A = 2.0 * C / chiABN 
        + nD*double(N)*hhat[i]*hhat[i] * (gaa[i] + 2.0*gab[i] + gbb[i]) / V 
        + nAH * Nah * Nah / double(N) * hhat[i] * hhat[i] * gha[i] / V ;
      numer = wabp[i] - lam_pl * ( F - A * wabp[i] ) ;
      if ( do_CL )
        numer += etap[i] ;
      denom = 1.0 + lam_pl * A ;
      wabp[i] = numer / denom ;

      // AB- //
      F = 2.0 * C * wabm[i] / chiABN
          + rho_fld_np[i] / double(N)
          + hhat[i] / double(N) * ( rhodb[i] - rhoda[i] - rhoha[i] );
      A = 2.0 * C / chiABN ;
      numer = wabm[i] - lam_mi * ( F - A * wabm[i] ) ;
      if ( do_CL ) 
        numer += etam[i] ;
      denom = 1.0 + lam_mi * A ;
      wabm[i] = numer / denom ;
    }

    fft_bck_wrapper( wabm , wabm ) ;
    fft_bck_wrapper( wabp , wabp ) ;
  } // if chiABN > 0.0
  else {
    for ( i=0 ; i<ML ; i++ )
      wabp[i] = wabm[i] = 0.0 ;
  }

  // Update AC fields //
  if (chiACN > 0.0) {
    fft_fwd_wrapper(wacm, wacm); 
    fft_fwd_wrapper(wacp, wacp); 

    if (do_CL) {
      generate_1s_noise(etap, lam_pl);
      generate_1s_noise(etam, lam_mi);
    }

    for (i=0; i<ML; i++) {
      // AC+ //
      F = 2.0 * C * wacp[i] / chiACN
          + I * rho_fld_np[i] / double(N)
          + I * hhat[i] / double(N) * (rhoda[i] + rhoha[i] + rhohc[i]);
      A = 2.0 * C / chiACN
          + nAH * Nah * Nah / double(N) * hhat[i] * hhat[i] * gha[i] / V
          + nCH * Nch * Nch / double(N) * hhat[i] * hhat[i] * ghc[i] / V;
      numer = wacp[i] - lam_pl * (F - A * wacp[i]);
      if (do_CL)
        numer += etap[i];
      denom = 1.0 + lam_pl * A;
      wacp[i] = numer / denom;

      // AC- //
      F = 2.0 * C * wacm[i] / chiACN
          + rho_fld_np[i] / double(N)
          + hhat[i] / double(N) * (rhohc[i] - rhoda[i] - rhoha[i]);
      A = 2.0 * C / chiACN;
      numer = wacm[i] - lam_mi * (F - A * wacm[i]);
      if (do_CL)
        numer += etam[i];
      denom = 1.0 + lam_mi * A;
      wacm[i] = numer / denom;
    }

    fft_bck_wrapper(wacm, wacm);
    fft_bck_wrapper(wacp, wacp);
  } // if chiACN > 0.0
  else {
    for (i=0; i<ML; i++)
      wacp[i] = wacm[i] = 0.0;
  }

  if (chiBCN > 0.0) {
    fft_fwd_wrapper(wbcm, wbcm); 
    fft_fwd_wrapper(wbcp, wbcp); 

    if (do_CL) {
      generate_1s_noise(etap, lam_pl);
      generate_1s_noise(etam, lam_mi);
    }

    for (i=0; i<ML; i++) {
      // BC+ //
      F = 2.0 * C * wbcp[i] / chiBCN
          + I * hhat[i] / double(N) * (rhodb[i] + rhohc[i]);
      A = 2.0 * C / chiBCN
          + nCH * Nch * Nch / double(N) * hhat[i] * hhat[i] * ghc[i] / V;
      numer = wbcp[i] - lam_pl * (F - A * wbcp[i]);
      if (do_CL)
        numer += etap[i];
      denom = 1.0 + lam_pl * A;
      wbcp[i] = numer / denom;

      // BC- //
      F = 2.0 * C * wbcm[i] / chiBCN
          + hhat[i] / double(N) * (rhohc[i] - rhodb[i]);
      A = 2.0 * C / chiBCN;
      numer = wbcm[i] - lam_mi * (F - A * wbcm[i]);
      if (do_CL)
        numer += etam[i];
      denom = 1.0 + lam_mi * A;
      wbcm[i] = numer / denom;
    }

    fft_bck_wrapper(wbcm, wbcm);
    fft_bck_wrapper(wbcp, wbcp);
  } // if chiBCN > 0.0
  else {
    for (i=0; i<ML; i++)
      wbcp[i] = wbcm[i] = 0.0;
  }

  calc_poly_density() ;
}


// Generates Gaussian noise in k-space with appropriate statistics //
// for the 1s updating scheme. 
void generate_1s_noise( complex<double> *et, double lambda ) {

  int i ;
  double scale = sqrt( 2.0 * lambda ) ;
  for ( i=0 ; i<Dim ; i++ )
    scale /= sqrt( dx[i] ) ;

  for ( i=0 ; i<ML ; i++ ) 
    et[i] = scale * gasdev2() ;

  fft_fwd_wrapper( et , et ) ;

}
