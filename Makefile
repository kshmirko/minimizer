FC=gfortran
BASEPATH=/usr/local/Cellar
FFLAGS=-Ofast -flto -fexternal-blas -fcheck= -I/$(BASEPATH)/netcdf/4.5.0/include -I$(BASEPATH)/nlopt/include -I$(BASEPATH)/openblas/0.2.20/include -Bshared
LDFLAGS=-L$(BASEPATH)/netcdf/4.5.0/lib -lnetcdf -lnetcdff -L$(BASEPATH)/nlopt/lib -lnlopt -lm -L$(BASEPATH)/openblas/0.2.20/lib -lpthread -lgfortran -lopenblas

OBJ=bin/libutils.a \
		src/main.o

LIBOBJ=src/types.o \
				src/ncdfapi.o \
				src/optprops.o \
				src/polarization.o \
				src/mathutils.o \
				src/distrtypes.o \
				src/dataio.o\
				src/minimize_funct.o

.PHONY: all clean
	
	
%.o: %.f03
	echo $@
	$(FC) -c -J src/ -o $@ $< $(FFLAGS)
	
all: libutils.a minimizer
	
libutils.a: $(LIBOBJ)
	ar rcs bin/$@ $(LIBOBJ)
	
minimizer: $(OBJ)
	$(FC) -o bin/$@ $^ $(FFLAGS) $(LDFLAGS)
	
clean:
	rm src/*.mod
	rm src/*.o
	rm bin/*.a
	rm bin/minimizer


	