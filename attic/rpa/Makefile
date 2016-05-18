include config.mk

# ============ Root Directories ================================
# REPO = Root directory for pscf source code repository
# SRC = $(REPO)/src = directory for source files
# BLD = $(REPO)/make = build directory for building with make
REPO=..
SRC=$(REPO)/src
BLD=$(REPO)/make

# ======  Forpedo Fortran preprocessor =========================
# FORPEDO = path to Forpedo preprocessor

FORPEDO=$(REPO)/tools/python/preprocess-0.6.1/preprocess.py

# ====== Main Program ==========================================
.PHONY: all clean 

all:
	make rpa_mod.o

clean:
	rm -f  *.o  *.mod

rpa_mod.o: $(SRC)/rpa_mod.f\
	chemistry_mod.o const_mod.o
	$(F90) $(FAST) -c $(SRC)/rpa_mod.f

