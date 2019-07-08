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

# ====== Modules =================================================

const_mod.o: $(SRC)/const_mod.f90
	$(F90) $(FAST) -c $(SRC)/const_mod.f90

# io Directory
string_mod.o: $(SRC)/io/string_mod.f90
	$(F90) $(FAST) -c $(SRC)/io/string_mod.f90

io_mod.o: $(SRC)/io/io_mod.f90 const_mod.o string_mod.o
	$(F90) $(FAST) -c $(SRC)/io/io_mod.f90

version_mod.o: $(SRC)/io/version_mod.f90
	$(F90) $(FAST) -c $(SRC)/io/version_mod.f90

field_io_mod.o: $(SRC)/io/field_io_mod.f90 const_mod.o io_mod.o\
	string_mod.o unit_cell_mod.o chemistry_mod.o basis_mod.o\
	$(FFT_FILE).o grid_basis_mod.o
	$(F90) $(FAST) -c $(SRC)/io/field_io_mod.f90

# grid directory
grid_mod.o: $(SRC)/grid/grid_mod.f90 const_mod.o\
	group_mod.o unit_cell_mod.o
	$(F90) $(FAST) -c $(SRC)/grid/grid_mod.f90

$(FFT_FILE).o: $(SRC)/grid/$(FFT_FILE).f90 const_mod.o
	$(F90) $(FAST) -c $(SRC)/grid/$(FFT_FILE).f90

# crystal directory
group_mod.o: $(SRC)/crystal/group_mod.f90 const_mod.o version_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/group_mod.f90

unit_cell_mod.o: $(SRC)/crystal/unit_cell_mod.f90 const_mod.o\
	io_mod.o group_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/unit_cell_mod.f90

space_groups_mod.o: $(SRC)/crystal/space_groups_mod.f90 const_mod.o\
	group_mod.o
	$(F90) $(NOPT) -c $(SRC)/crystal/space_groups_mod.f90

basis_mod.o: $(SRC)/crystal/basis_mod.f90\
	const_mod.o string_mod.o io_mod.o\
	group_mod.o space_groups_mod.o unit_cell_mod.o\
	grid_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/basis_mod.f90

grid_basis_mod.o: $(SRC)/crystal/grid_basis_mod.f90 const_mod.o\
	grid_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/grid_basis_mod.f90

deform_mod.o: $(SRC)/crystal/deform_mod.f90\
	const_mod.o string_mod.o io_mod.o\
	group_mod.o space_groups_mod.o unit_cell_mod.o\
	grid_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/deform_mod.f90

# scf directory
chemistry_mod.o: $(SRC)/scf/chemistry_mod.f90 const_mod.o io_mod.o
	$(F90) $(FAST) -c $(SRC)/scf/chemistry_mod.f90

chain_mod.o: $(SRC)/scf/chain_mod.f90 const_mod.o\
	chemistry_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/scf/chain_mod.f90

step_mod.o: $(SRC)/scf/step_mod.f90 const_mod.o\
	$(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/scf/step_mod.f90

scf_mod.f90: $(SRC)/scf/scf_mod.fp.f90
	$(FORPEDO) $(DEVEL) $(SRC)/scf/scf_mod.fp.f90 > scf_mod.f90

scf_mod.o: scf_mod.f90 const_mod.o io_mod.o\
	basis_mod.o chemistry_mod.o step_mod.o\
	grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o
	$(F90) $(FAST) -c scf_mod.f90

#rpa directory
rpa_mod.o: $(SRC)/rpa/rpa_mod.f90\
	chemistry_mod.o const_mod.o
	$(F90) $(FAST) -c $(SRC)/rpa/rpa_mod.f90

spinodal_mod.o: $(SRC)/rpa/spinodal_mod.f90 const_mod.o io_mod.o\
	response_pd_mod.o chemistry_mod.o
	$(F90) $(FAST) -c $(SRC)/rpa/spinodal_mod.f90

# iterate directory
iterate_mod.f90: $(SRC)/iterate/iterate_mod.fp.f90
	$(FORPEDO) $(DEVEL) $(SRC)/iterate/iterate_mod.fp.f90 > iterate_mod.f90

iterate_mod.o: iterate_mod.f90 const_mod.o\
	scf_mod.o basis_mod.o chemistry_mod.o unit_cell_mod.o\
	response_pd_mod.o
	$(F90) $(FAST) -c iterate_mod.f90

response_pd_mod.o: $(SRC)/iterate/response_pd_mod.f90 const_mod.o io_mod.o\
	chemistry_mod.o basis_mod.o unit_cell_mod.o scf_mod.o\
	grid_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/iterate/response_pd_mod.f90

sweep_mod.o: $(SRC)/iterate/sweep_mod.f90 const_mod.o io_mod.o\
	chemistry_mod.o unit_cell_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/iterate/sweep_mod.f90

# response directory
response_step_mod.o: $(SRC)/response/response_step_mod.f90\
	chemistry_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/response/response_step_mod.f90

extrapolate_mod.o: $(SRC)/response/extrapolate_mod.f90\
	const_mod.o
	$(F90) $(FAST) -c $(SRC)/response/extrapolate_mod.f90

response_mod.o:	$(SRC)/response/response_mod.f90\
	chemistry_mod.o const_mod.o chain_mod.o\
	grid_mod.o $(FFT_FILE).o group_mod.o response_step_mod.o\
	extrapolate_mod.o field_io_mod.o spinodal_mod.o
	$(F90) $(FAST) -c $(SRC)/response/response_mod.f90

# ========================================================
# tests subdirectory 

rpa_test.o: $(SRC)/tests/rpa/test.f90 rpa_mod.o\
	chemistry_mod.o io_mod.o string_mod.o const_mod.o
	$(F90) $(FAST) -o rpa_test.o -c $(SRC)/tests/rpa/test.f90

rpa_test: rpa_test.o rpa_mod.o\
	chemistry_mod.o io_mod.o string_mod.o const_mod.o
	$(F90) $(LINK) -lm -o $(SRC)/tests/rpa/exc\
	rpa_test.o rpa_mod.o\
	chemistry_mod.o io_mod.o string_mod.o const_mod.o\
	$(LAPACKLIB)

# ========================================================
# Code to convert 2D groups to Fortran
 
new_2dgroups: $(SRC)/tests/group/2dgroups/new_2dgroups.f90\
	const_mod.o io_mod.o version_mod.o string_mod.o\
	group_mod.o unit_cell_mod.o
	$(F90) $(LINK) -o $(SRC)/tests/group/2dgroups/new_2dgroups\
	$(SRC)/tests/group/2dgroups/new_2dgroups.f90\
	const_mod.o io_mod.o string_mod.o version_mod.o\
	group_mod.o unit_cell_mod.o
>>>>>>> master

