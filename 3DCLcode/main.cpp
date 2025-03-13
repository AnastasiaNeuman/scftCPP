#define MAIN
#include "globals.h"
#include "r2.h"
void write_outputs(void);
void accumulate_all_averages(void);
void update_1s(void);
void update_Euler(void);
void write_data(char*, complex<double>*); // for debugging
double brent_method(double, double, double);
void change_L0(double);
double simulate(void);
double sweep_z(void);

// The main routine is essentially a wrapper that allocates memory, then either
// runs a single simuation or implements Brent's method to find the optimum box
// size before calculating the structure and energy at that box size.
int main(int argc, char** argv) {

#ifdef PAR
  // MPI initialization stuff
  MPI_Init( &argc , &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
  fftw_mpi_init();
#else
  myrank = 0;
#endif

  // Other initialization stuff
#ifdef SETIDUM 
  idum = 100;
#else
  idum = -long(time(0)) * (myrank + 1);
#endif
  read_input();
  initialize_1();
  allocate();

  // Flag that's only marked true for first simulation init
  first_sim = 1;

  // Check if the do_brent flag is true, then either run Brent's method to
  // optimize box size or just run a single simulation.
  if (do_brent) {
    // Brent's method to optimize box size
    double L_brent = brent_method(L_low, L_high, brent_tol);
    if (myrank == 0) {
      printf("L_brent = %lf\n\n", L_brent);
    }
    // Simulate at L_brent if that wasn't just done. There was not necessarily
    // a simulation performed at the L value returned by Brent's method,
    // especially with the current (as of 6/10/15) implementation of get_slope.
    if (L[0] != L_brent) {
      change_L0(L_brent);
      simulate();
    }
    if (myrank == 0) printf("L_ideal = %lf\n\n", L_ideal);
    // Run a final simulation at the global optimum L. After running a
    // simulation at L_brent, it might turn out that a different point close
    // to L_brent had a slightly lower energy. If that is the case, run a final
    // simulation at the length that gives the global minimum energy.
    if (L[0] != L_ideal) {
      change_L0(L_ideal);
      simulate();
    }
  } // if do_brent
  else
    // Single simulation if not doing Brent's method
    simulate();

#ifdef PAR
  MPI_Finalize();
#endif

  return 0;
}
