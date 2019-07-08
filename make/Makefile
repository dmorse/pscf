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

# ====== Modules =================================================
include modules.mk

# ====== Main Program ==========================================
.PHONY: all clean install uninstall

all:
	make pscf

clean:
	rm -f  *.o  *.mod  *.f90  pscf

pscf:	pscf.o\
	const_mod.o string_mod.o io_mod.o version_mod.o field_io_mod.o\
	unit_cell_mod.o group_mod.o space_groups_mod.o basis_mod.o\
	chemistry_mod.o scf_mod.o iterate_mod.o sweep_mod.o\
	step_mod.o grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o\
	response_pd_mod.o spinodal_mod.o\
	response_mod.o response_step_mod.o extrapolate_mod.o
	$(F90) $(LINK) -o pscf pscf.o \
	const_mod.o string_mod.o io_mod.o version_mod.o field_io_mod.o \
	unit_cell_mod.o group_mod.o space_groups_mod.o basis_mod.o \
	chemistry_mod.o scf_mod.o iterate_mod.o sweep_mod.o \
	step_mod.o grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o \
	response_pd_mod.o spinodal_mod.o \
	response_mod.o response_step_mod.o extrapolate_mod.o \
	$(LIBDIRS) $(FFTWLIB) $(LAPACKLIB) -lm

pscf.o: pscf.f90 const_mod.o io_mod.o version_mod.o\
	unit_cell_mod.o chemistry_mod.o field_io_mod.o\
	sweep_mod.o iterate_mod.o scf_mod.o basis_mod.o\
	spinodal_mod.o\
	response_mod.o $(FFT_FILE).o response_step_mod.o extrapolate_mod.o
	$(F90) $(FAST) -c pscf.f90

pscf.f90: $(SRC)/pscf.fp.f90
	$(FORPEDO) $(DEVEL) $(SRC)/pscf.fp.f90 > pscf.f90

# ====== Install =================================================

BIN_DIR=$(INSTALL)/bin
PYTHON_DIR=$(INSTALL)/lib/python2.7/site-packages

install: pscf
	mkdir -p $(BIN_DIR)
	cp pscf $(BIN_DIR)/pscf
	cp $(REPO)/tools/bin/pscf-env $(BIN_DIR)
	cp $(REPO)/tools/bin/pscf-read-sweep $(BIN_DIR)
	mkdir -p $(PYTHON_DIR)/pscf
	cp $(REPO)/tools/python/pscf/*.py $(PYTHON_DIR)/pscf

uninstall:
	rm -rf $(BIN_DIR)/pscf*
	rm -rf $(PYTHON_DIR)/pscf/*

