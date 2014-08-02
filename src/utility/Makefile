BIN=/home/qin/bin
F90=pgf90
FAST=-fast -Mfree=f90
NOPT=-Mfree=f90

####################################
# Utilities
####################################
.PHONY: clean
.PHONY: veryclean

clean:
	rm -f *.o *.mod 

switch:	switch.f \
	../build/const_mod.o
	$(F90) $(FAST) -o $(BIN)/switch switch.f \
	../build/const_mod.o

nto2:	nto2.f \
	../build/const_mod.o
	$(F90) $(FAST) -o $(BIN)/nto2 nto2.f \
	../build/const_mod.o

nto3:	nto3.f \
	../build/const_mod.o
	$(F90) $(FAST) -o $(BIN)/nto3 nto3.f \
	../build/const_mod.o

