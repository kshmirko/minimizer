FC=gfortran
FFLAGS=-Ofast -flto -fcheck= -I/usr/local/Cellar/netcdf/4.5.0/include -I/usr/local/Cellar/nlopt/include -Bshared
LDFLAGS=-L/usr/local/Cellar/netcdf/4.5.0/lib -lnetcdf -lnetcdff -L/usr/local/Cellar/nlopt/lib -lnlopt -lm

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


	