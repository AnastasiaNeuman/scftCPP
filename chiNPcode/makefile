# To use this makefile, set FFTWHOME as an environment variable containing the
# path to the FFTW3 software. The contents of this directory should be the
# usual "bin include lib share". In the Comet cluster, FFTWHOME is
# automatically set if you have "module load fftw/3.3.4" in ~/.bashrc

#CC = g++
FFTWHOME = /opt/seas/pkg/gcc/fftw3/mpi/double/3.3.7
#FFTWHOME=${HOME}/Install/fftw3-intel-mpi
#FFTWHOME=${HOME}/Install/fftw3
CC = mpic++
DFLAGS = -DPAR
CFLAGS = -O3 -w -Wall
INCLUDE = -I${FFTWHOME}/include 
LDFLAGS = -lm -lfftw3_mpi -lfftw3 -L${FFTWHOME}/lib

SRCS = main.cpp calc_poly_density.cpp homopolymer_discrete.cpp \
       io_utils.cpp integrate_utils.cpp ran2.cpp array_utils.cpp calc_h.cpp \
       fft_wrappers.cpp calc_debye.cpp initialize.cpp \
       read_input.cpp diblock_discrete.cpp 1s_update.cpp \
       Euler_update.cpp brent.cpp simulate.cpp sphere_utils.cpp \
       homopolymer_continuous.cpp \
			 

OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c ${INCLUDE} $<

a.out:	${OBJS}
	${CC} ${CFLAGS} -o $@ ${OBJS} ${LDFLAGS}

clean:
	rm -f *.o a.out
