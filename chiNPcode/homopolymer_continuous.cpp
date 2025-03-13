#include "globals.h"

complex<double> get_avg(complex<double>*);

complex<double> homopolymer_continuous_gsd(complex<double> *WA,
    complex<double> *q, complex<double> *rho, double nH, double Nh, double vol,
    double ds, double tol, int it_max) {
  
  int it, i;
  double error = 1.0;
  complex<double> *q_old, *dum_arr, *lambda;

  q_old = (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);
  dum_arr = (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);
  lambda = (complex<double>*) fftw_malloc(sizeof(complex<double>)*size);

  for (i=0; i<ML; i++) {
    if (iter == 0) {
      q_old[i] = 1.0;
      //if (myrank==0)
      //  printf("iter=%d, initializing q_old with 1.0's\n", iter);
    }
    else
      q_old[i] = q[i];
  }

  for (it=0; it<it_max; it++) {
    for (i=0; i<ML; i++)
      q[i] = q_old[i] * exp(-0.5*WA[i]*ds);

    fft_fwd_wrapper(q, q);
    for (i=0; i<ML; i++)
      q[i] *= poly_bond_fft[i];
    fft_bck_wrapper(q, q);

    for (i=0; i<ML; i++) {
      q[i] *= exp(-0.5*WA[i]*ds);
      dum_arr[i] = q[i] * q[i];
    }

    complex<double> dum = integ_trapPBC(dum_arr);
    for (i=0; i<ML; i++) {
      lambda[i] = q_old[i] / q[i];
      q[i] *= sqrt(vol / dum);
    }

    if (error < tol && it > 3) {
      break;
    }
    else {
      // if (myrank == 0) printf("sqrt(vol/dum)=%lf\n", sqrt(vol/dum));
      for (i=0; i<ML; i++)
        dum_arr[i] = abs(q[i]-q_old[i]);

      error = real(integ_trapPBC(dum_arr));
      // if (myrank == 0) printf("it %d, qc error = %lf\n", it, error);
      for (i=0; i<ML; i++)
        q_old[i] = q[i];
    }
  }

  complex<double> avg_lambda = get_avg(lambda);
  if (myrank == 0 && iter % print_freq == 0) printf("avg_lambda=%lf\n",
                                                    real(avg_lambda));
  // if (myrank == 0) printf("avg_lambda = %lf + %lfi\n", real(avg_lambda),
  //                        imag(avg_lambda) );

  complex<double> Q = exp( -avg_lambda * Nh / double(N) );

  for (i=0; i<ML; i++) {
    rho[i] = nH * Nh / vol * q[i] * q[i];
  }
  if (nCH > 0.0 && iter % print_freq == 0) {
    // C homopolymer sanity check
    double ncNc_check = real( integ_trapPBC(rhohc) );
    if (myrank == 0)
      printf("ncNc_check = %lf, nc*Nc = %lf\n", ncNc_check, nCH * Nch); 
  }

  if (iter % print_freq == 0 && myrank == 0) {
      printf("%d GSD Iterations\n", it);
  }

  free(q_old);
  free(dum_arr);

  return Q;
}
