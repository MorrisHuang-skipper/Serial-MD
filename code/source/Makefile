FC=gfortran
FFLAGS=-O3 -fopenmp -std=gnu -g -Wall -fbounds-check
SRC= math.f95 inp.f95 mic.f95 main.f95
OBJ=${SRC:.f95=.o}

%.o: %.f95
	${FC} ${FFLAGS} -o $@ -c $<

morris: ${OBJ}
	${FC} ${FFLAGS} -o $@ ${OBJ}

clean:
	@rm -f *.mod *.o morris