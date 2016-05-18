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
.PHONY: all clean install uninstall

all:
	make pscf

clean:
	rm -f  *.o  *.mod  *.f  pscf

deform_mod.o: $(SRC)/deform_mod.f\
	const_mod.o string_mod.o io_mod.o\
	group_mod.o space_groups_mod.o unit_cell_mod.o\
	grid_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/deform_mod.f

