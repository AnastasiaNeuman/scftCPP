#include "globals.h"
#include "r2.h"
double get_slope(double);
void change_L0(double);
double simulate(void);

// This implementation of Brent's method optimizes box size in the x-direction
// (L[0]) and maintains the same x-y ratio. The z dimension doesn't change.
double brent_method(double lowerLimit, double upperLimit, double errorTol) {
  double a = lowerLimit;
  double b = upperLimit;
  double c = 0;
  double d = upperLimit * 1.25;

  // Create blank file that will store final results for each simulation run.
  // One line is added to this file at the end of each pass through the
  // simulate method
  brent_otp = fopen("brent.dat", "w");
  if (brent_otp == NULL) {
    printf("Failed to open brent.dat!\n");
    exit(1);
  }
  fclose(brent_otp);

  double fa = get_slope(a);
  double fb = get_slope(b);

  double fc = 0;
  double s = 0;
  double fs = 0;

  // if f(a) f(b) >= 0 then error-exit
  if (fa * fb >= 0 ) {
    if (myrank == 0) {
      printf("Bad length range for Brent's method! Change range and retry\n");
    }
    exit(1);
  }

  // if |f(a)| < |f(b)| then swap (a,b) end if
  if (abs(fa) < abs(fb)) {
    double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp;
  }

  c = a;
  fc = fa;
  bool mflag = true;

  // Iterate max 5 times so this doesn't go on forever
  for (int j=0; j<5; j++) {
    if (myrank == 0)
      printf("---Beginning of iteration %d of Brent's method---\n", j);
    if ((fa != fc) && (fb != fc)) {
      // Inverse quadratic interpolation
      s = a*fb*fc / (fa-fb) / (fa-fc) +
          b*fa*fc / (fb-fa) / (fb-fc) +
          c*fa*fb / (fc-fa) / (fc-fb);
      if (myrank == 0) {
        printf("Used inverse quadratic interpolation with:\n");
        printf("a=%lf, b=%lf, c=%lf, fa=%lf, fb=%lf, fc=%lf\n",
                a,     b,     c,     fa,     fb,     fc);
        printf("to get s=%lf\n\n", s);
      }
    }
    else {
      // Secant method (linear interpolation)
      s = b - fb * (b-a) / (fb-fa);
      if (myrank == 0) {
        printf("Used secant method with\n");
        printf("a=%lf, b=%lf, fa=%lf, fb=%lf\n", a, b, fa, fb);
        printf("to get s=%lf\n\n", s);
      }
    }

    // Use bisection method if any of the following conditions is true:
    // Condition 1: s is not between (3a+b)/4 and b
    double tmp2 = (3 * a + b) / 4;
    bool cond1 = (s<tmp2 && s<b) || (s>tmp2 && s>b);
    // Condition 2: mflag is set and |s-b|>=|b-c|/2
    bool cond2 = mflag && (abs(s-b) >= abs(b-c)/2);
    // Condition 3: mflag is cleared and |s-b|>=|c-d|/2
    bool cond3 = !mflag && (abs(s-b) >= abs(c-d)/2);
    // Condition 4: mflag is set and |b-c|<|delta|
    bool cond4 = mflag && (abs(b-c) < errorTol);
    // Condition 5: mflag is cleared and |c-d|<|delta|
    bool cond5 = !mflag && (abs(c-d) < errorTol);

    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      s = (a + b) / 2;
      mflag = true;
      if (myrank == 0) {
        printf("Used bisection method with\n");
        printf("a=%lf and b=%lf to get s=%lf\n\n", a, b, s);
      }
    }
    else
      mflag = false;

    fs = get_slope(s);
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0) { b = s; fb = fs; }
    else { a = s; fa = fs; }

    // if |f(a)| < |f(b)| then swap (a,b) end if
    if (abs(fa) < abs(fb)) {
      double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp;
    }

    if (myrank == 0) {
      printf("a=%lf, b=%lf, c=%lf, fa=%lf, fb=%lf, fc=%lf, s=%lf\n\n",
              a,     b,     c,     fa,     fb,     fc,     s);
      printf("---End of iteration %d of Brent's method---\n\n", j);
    }

    // exit loop if f(s)=0 or |b-a| is small enough
    if (fs==0.0 || abs(b-a) < errorTol) break;
  } // j for loop

  return b;
}

// Returns the slope of the H/V vs L graph at the L-value input to get_slope.
// This is a part of Brent's method where H/V is calculated at 2 nearby L 
// values by running separate simulations, then the slope is calculated
double get_slope(double l) {
  double L0_1 = l - 0.5 * L_step;
  double L0_2 = l + 0.5 * L_step;
  if (myrank == 0) {
    printf("Calculating slope at L[0]=%lf using simulations at %lf and %lf\n\n",
           l, L0_1, L0_2);
  }

  // Change L[0] to l-L_step/2 and keep same L[1] / L[0] ratio for 3D
  change_L0(L0_1);
  double H_over_V_1 = simulate();

  // Add L_step to L[0] and keep same L[1] / L[0] ratio for 3D
  change_L0(L0_2);
  double H_over_V_2 = simulate();

  // Return slope of H/V vs L graph
  return (H_over_V_2 - H_over_V_1) / L_step;
}

// Little function to change L[0] to len passed in and keep same Lx/Ly ratio
// for 3D simulations
void change_L0(double len) {
  if (Dim==3) L[1] *= len / L[0];
  L[0] = len;
}
