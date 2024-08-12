FC = gfortran
#FC = g95

#FFLAGS = -O3 -fopenmp
FFLAGS = -O3 -mcmodel=medium
#FFLAGS = -O0 -ggdb

#LFLAGS = -O3 -fopenmp
LFLAGS = -O3 -mcmodel=medium
#LFLAGS = -O0 -ggdb

#LIBS = -lgomp /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libblas.a
LIBS = /usr/lib/x86_64-linux-gnu/liblapack.a /usr/lib/x86_64-linux-gnu/libblas.a

#DEBUG = -fsanitize=address
DEBUG = 

#OBJECTS = interface.o silic2ddos.o tightb.o
OBJECTS = interface.o run2.o 

MODULES = interface.mod

DATA = 

.PHONY: clean

run2.dat: run2.exe
	./run2.exe > saida.txt
#	gdb run2.exe

run2.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) $(DEBUG) -o run2.exe $(LIBS)

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(MODULES) $(DATA) $(FIGURES) run2.exe

help:
	@echo "Valid targets:"
	@echo "  run2.exe"
	@echo "  run2.o"
	@echo "  clean:  removes .o, .txt, .ps, and .exe files"
