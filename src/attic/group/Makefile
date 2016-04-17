LAPACKLIB=/usr/lib/liblapack.a /usr/lib/libblas.a
F90=pgf90
BIN=/home/morse/bin
FAST=-fast -Mfree=f90
NOPT=-Mfree=f90

####################################
# Utilities
####################################
.PHONY: clean
.PHONY: veryclean

clean:
	rm -f *.o *.mod 

new_group:	new_group.f \
        ../../unit_cell_mod.o ../../group_mod.o \
	../../io_mod.o ../../const_mod.o
	$(F90) $(FAST) -o $(BIN)/new_group new_group.f \
        ../../unit_cell_mod.o ../../group_mod.o \
	../../io_mod.o ../../const_mod.o

