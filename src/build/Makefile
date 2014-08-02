# ============ Root Directories ================================
# SCF = Root directory for pscf project
# SCF = $(SCF)/src = directory for source files
SCF=../..
SRC=..

# ====== Executable - Directory and File Name ==================
# BIN = path to directory for executable file
# EXE = name of executable file
BIN=$(SCF)
EXE=pscf

# ============ Compiler and Options ============================
# Fortran compiler and settings for high and low optimization

# Gnu Compiler Collection (gcc) - gfortran
# F90=gfortran
# FAST=-O3 -ffree-form 
# NOPT=-ffree-form
# LINK=-O3 -ffree-form -lgfortran
 
# Intel compiler
# F90=ifort
# FAST=-O3 -free
# NOPT=-free
# LINK=-O3 -free

# Portland group compiler
  F90=pgf90
  FAST=-O3 -Mfree=f90
  NOPT=-Mfree=f90
  LINK=-O3 -Mfree=90

# ============ FFTW Fast Fourier Transform Library =============
# FFTWLIB  = path(s) to directory(ies) containing FFTW library
# FFT_FILE = fft2_mod for FFTW version 2 
#          = fft3_mod for FFTW version 3

# Use fftw3 on fedora core 7 (fftw3 installed by yum)
# FFTWLIB=-L/usr/lib -lfftw3
# FFT_FILE=fft3_mod

# Use fftw2 on flory.cems.umn.edu
  FFTWLIB=-L/usr/local/lib -lrfftw -lfftw
  FFT_FILE=fft2_mod

# Use fftw3 on nf.msi.umn.edu
# FFTWLIB=-L/usr/local/fftw-3.0.1/lib -lfftw3
# fft_file=fft3_mod

# use fftw3 on altix.msi.umn.edu
# FFTWLIB=-L/soft/local/fftw/3.0.1-fma/lib -lfftw3
# fft_file=fft3_mod

# Use fftw3 on calhoun.msi.umn.edu
# FFTWLIB=-L/soft/fftw/lib -lfftw3
# FFT_FILE=fft3_mod

# ============ LAPACK Linear Algebra Library ==================
# LAPACKLIB = path(s) to directories containing lapack and blas

# On fedora core 7 (lapack installed by yum)
# LAPACKLIB=-L/usr/lib -llapack -lblas
   
# On flory.cems.umn.edu or nf.msi.umn.edu
  LAPACKLIB=-llapack -lblas

# on altix.msi.umn.edu
# LAPACKLIB=-L/usr/local/LAPACK-8.0b -llapack -L/usr/local/LAPACK -lblas

# on calhoun.msi.umn.edu
# LAPACKLIB=/cluster/apps/intel/cmkl/9.1.021/lib/em64t/libmkl_lapack.a /cluster/apps/intel/cmkl/9.1.021/lib/em64t/libmkl_lapack.so /cluster/apps/intel/cmkl/9.1.021/lib/em64t/libmkl_em64t.a -lpthread

# ======  Forpedo Fortran preprocessor =========================
# FORPEDO = path to Forpedo preprocessor
# DEVEL   = '' (empty string) for stable version, 
#         = '-D DEVEL' for development version

FORPEDO=$(SCF)/tools/python/preprocess-0.6.1/preprocess.py
#DEVEL=-D DEVEL
DEVEL=

# ======   Dummy targets to clean up debris ===================

.PHONY: clean

clean:
	rm -f  *.o  *.mod  *~\
	$(SRC)/*~  $(SRC)/*/*~  $(SRC)/*/*/*~  $(SRC)/tests/*/exc

very_clean:
	rm -f  *.o  *.mod  *.f  *~\
	$(SRC)/*~  $(SRC)/*/*~  $(SRC)/*/*/*~  $(SRC)/tests/*/exc\

# ====== Main Program ==========================================

pscf_pd.f: $(SRC)/pscf_pd.fp.f
	$(FORPEDO) $(DEVEL) $(SRC)/pscf_pd.fp.f > pscf_pd.f

pscf_pd.o: pscf_pd.f const_mod.o io_mod.o version_mod.o\
	unit_cell_mod.o chemistry_mod.o field_io_mod.o\
	sweep_mod.o iterate_mod.o scf_mod.o basis_mod.o\
	spinodal_mod.o\
	response_mod.o $(FFT_FILE).o response_step_mod.o extrapolate_mod.o
	$(F90) $(FAST) -c pscf_pd.f

pscf:	pscf_pd.o\
	const_mod.o string_mod.o io_mod.o version_mod.o field_io_mod.o\
	unit_cell_mod.o group_mod.o space_groups_mod.o basis_mod.o\
	chemistry_mod.o scf_mod.o iterate_mod.o sweep_mod.o\
	step_mod.o grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o\
	response_pd_mod.o spinodal_mod.o\
	response_mod.o response_step_mod.o extrapolate_mod.o
	$(F90) $(LINK) $(FFTWLIB) $(LAPACKLIB) -lm -o $(BIN)/$(EXE) pscf_pd.o \
	const_mod.o string_mod.o io_mod.o version_mod.o field_io_mod.o \
	unit_cell_mod.o group_mod.o space_groups_mod.o basis_mod.o \
	chemistry_mod.o scf_mod.o iterate_mod.o sweep_mod.o \
	step_mod.o grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o \
	response_pd_mod.o spinodal_mod.o \
	response_mod.o response_step_mod.o extrapolate_mod.o

# ====== Modules =================================================

const_mod.o: $(SRC)/const_mod.f
	$(F90) $(FAST) -c $(SRC)/const_mod.f

# io Directory
string_mod.o: $(SRC)/io/string_mod.f
	$(F90) $(FAST) -c $(SRC)/io/string_mod.f

io_mod.o: $(SRC)/io/io_mod.f const_mod.o string_mod.o
	$(F90) $(FAST) -c $(SRC)/io/io_mod.f

version_mod.o: $(SRC)/io/version_mod.f
	$(F90) $(FAST) -c $(SRC)/io/version_mod.f

field_io_mod.o: $(SRC)/io/field_io_mod.f const_mod.o io_mod.o\
	string_mod.o unit_cell_mod.o chemistry_mod.o basis_mod.o\
	$(FFT_FILE).o grid_basis_mod.o
	$(F90) $(FAST) -c $(SRC)/io/field_io_mod.f

#grid directory
grid_mod.o: $(SRC)/grid/grid_mod.f const_mod.o\
	group_mod.o unit_cell_mod.o
	$(F90) $(FAST) -c $(SRC)/grid/grid_mod.f

$(FFT_FILE).o: $(SRC)/grid/$(FFT_FILE).f const_mod.o
	$(F90) $(FAST) -c $(SRC)/grid/$(FFT_FILE).f

# crystal directory
group_mod.o: $(SRC)/crystal/group_mod.f const_mod.o version_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/group_mod.f

unit_cell_mod.o: $(SRC)/crystal/unit_cell_mod.f const_mod.o\
	io_mod.o group_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/unit_cell_mod.f

space_groups_mod.o: $(SRC)/crystal/space_groups_mod.f const_mod.o\
	group_mod.o
	$(F90) $(NOPT) -c $(SRC)/crystal/space_groups_mod.f

basis_mod.o: $(SRC)/crystal/basis_mod.f\
	const_mod.o string_mod.o io_mod.o\
	group_mod.o space_groups_mod.o unit_cell_mod.o\
	grid_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/basis_mod.f

grid_basis_mod.o: $(SRC)/crystal/grid_basis_mod.f const_mod.o\
	grid_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/grid_basis_mod.f

deform_mod.o: $(SRC)/crystal/deform_mod.f\
	const_mod.o string_mod.o io_mod.o\
	group_mod.o space_groups_mod.o unit_cell_mod.o\
	grid_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/crystal/deform_mod.f

# scf directory
chemistry_mod.o: $(SRC)/scf/chemistry_mod.f const_mod.o io_mod.o
	$(F90) $(FAST) -c $(SRC)/scf/chemistry_mod.f

chain_mod.o: $(SRC)/scf/chain_mod.f const_mod.o\
	chemistry_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/scf/chain_mod.f

step_mod.o: $(SRC)/scf/step_mod.f const_mod.o\
	$(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/scf/step_mod.f

scf_mod.f: $(SRC)/scf/scf_mod.fp.f
	$(FORPEDO) $(DEVEL) $(SRC)/scf/scf_mod.fp.f > scf_mod.f

scf_mod.o: scf_mod.f const_mod.o io_mod.o\
	basis_mod.o chemistry_mod.o step_mod.o\
	grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o
	$(F90) $(FAST) -c scf_mod.f

#rpa directory
rpa_mod.o: $(SRC)/rpa/rpa_mod.f\
	chemistry_mod.o const_mod.o
	$(F90) $(FAST) -c $(SRC)/rpa/rpa_mod.f

spinodal_mod.o: $(SRC)/rpa/spinodal_mod.f const_mod.o io_mod.o\
	response_pd_mod.o chemistry_mod.o
	$(F90) $(FAST) -c $(SRC)/rpa/spinodal_mod.f

# iterate directory
iterate_mod.f: $(SRC)/iterate/iterate_mod.fp.f
	$(FORPEDO) $(DEVEL) $(SRC)/iterate/iterate_mod.fp.f > iterate_mod.f

iterate_mod.o: iterate_mod.f const_mod.o\
	scf_mod.o basis_mod.o chemistry_mod.o unit_cell_mod.o\
	response_pd_mod.o
	$(F90) $(FAST) -c iterate_mod.f

response_pd_mod.o: $(SRC)/iterate/response_pd_mod.f const_mod.o io_mod.o\
	chemistry_mod.o basis_mod.o unit_cell_mod.o scf_mod.o\
	grid_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/iterate/response_pd_mod.f

sweep_mod.o: $(SRC)/iterate/sweep_mod.f const_mod.o io_mod.o\
	chemistry_mod.o unit_cell_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/iterate/sweep_mod.f


# response directory
response_step_mod.o: $(SRC)/response/response_step_mod.f\
	$(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/response/response_step_mod.f

extrapolate_mod.o: $(SRC)/response/extrapolate_mod.f\
	const_mod.o
	$(F90) $(FAST) -c $(SRC)/response/extrapolate_mod.f

response_mod.o:	$(SRC)/response/response_mod.f\
	chemistry_mod.o const_mod.o chain_mod.o\
	grid_mod.o $(FFT_FILE).o group_mod.o response_step_mod.o\
	extrapolate_mod.o field_io_mod.o spinodal_mod.o
	$(F90) $(FAST) -c $(SRC)/response/response_mod.f

# ========================================================
# tests subdirectory 

rpa_test.o: $(SRC)/tests/rpa/test.f rpa_mod.o\
	chemistry_mod.o io_mod.o string_mod.o const_mod.o
	$(F90) $(FAST) -o rpa_test.o -c $(SRC)/tests/rpa/test.f

rpa_test: rpa_test.o rpa_mod.o\
	chemistry_mod.o io_mod.o string_mod.o const_mod.o
	$(F90) $(LINK) -lm -o $(SRC)/tests/rpa/exc\
	rpa_test.o rpa_mod.o\
	chemistry_mod.o io_mod.o string_mod.o const_mod.o\
	$(LAPACKLIB)

# ========================================================
# Code to convert 2D groups to Fortran
 
new_2dgroups: $(SRC)/tests/group/2dgroups/new_2dgroups.f\
	const_mod.o io_mod.o version_mod.o string_mod.o\
	group_mod.o unit_cell_mod.o
	$(F90) $(LINK) -o $(SRC)/tests/group/2dgroups/new_2dgroups\
	$(SRC)/tests/group/2dgroups/new_2dgroups.f\
	const_mod.o io_mod.o string_mod.o version_mod.o\
	group_mod.o unit_cell_mod.o

