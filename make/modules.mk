
# Utility modules

const_mod.o: $(SRC)/const_mod.f90
	$(F90) $(FAST) -c $(SRC)/const_mod.f90

string_mod.o: $(SRC)/string_mod.f90
	$(F90) $(FAST) -c $(SRC)/string_mod.f90

io_mod.o: $(SRC)/io_mod.f90 const_mod.o string_mod.o
	$(F90) $(FAST) -c $(SRC)/io_mod.f90

version_mod.o: $(SRC)/version_mod.f90
	$(F90) $(FAST) -c $(SRC)/version_mod.f90

# FFT grid modules

grid_mod.o: $(SRC)/grid_mod.f90 const_mod.o\
	group_mod.o unit_cell_mod.o
	$(F90) $(FAST) -c $(SRC)/grid_mod.f90

$(FFT_FILE).o: $(SRC)/$(FFT_FILE).f90 const_mod.o
	$(F90) $(FAST) -c $(SRC)/$(FFT_FILE).f90

# Crystallography modules

group_mod.o: $(SRC)/group_mod.f90 const_mod.o version_mod.o
	$(F90) $(FAST) -c $(SRC)/group_mod.f90

unit_cell_mod.o: $(SRC)/unit_cell_mod.f90 const_mod.o\
	io_mod.o group_mod.o
	$(F90) $(FAST) -c $(SRC)/unit_cell_mod.f90

space_groups_mod.o: $(SRC)/space_groups_mod.f90 const_mod.o\
	group_mod.o
	$(F90) $(NOPT) -c $(SRC)/space_groups_mod.f90

basis_mod.o: $(SRC)/basis_mod.f90\
	const_mod.o string_mod.o io_mod.o\
	group_mod.o space_groups_mod.o unit_cell_mod.o\
	grid_mod.o
	$(F90) $(FAST) -c $(SRC)/basis_mod.f90

grid_basis_mod.o: $(SRC)/grid_basis_mod.f90 const_mod.o\
	grid_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/grid_basis_mod.f90

field_io_mod.o: $(SRC)/field_io_mod.f90 const_mod.o io_mod.o\
	string_mod.o unit_cell_mod.o chemistry_mod.o basis_mod.o\
	$(FFT_FILE).o grid_basis_mod.o
	$(F90) $(FAST) -c $(SRC)/field_io_mod.f90

# SCFT modules

chemistry_mod.o: $(SRC)/chemistry_mod.f90 const_mod.o io_mod.o
	$(F90) $(FAST) -c $(SRC)/chemistry_mod.f90

chain_mod.o: $(SRC)/chain_mod.f90 const_mod.o\
	chemistry_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/chain_mod.f90

step_mod.o: $(SRC)/step_mod.f90 const_mod.o\
	$(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/step_mod.f90

scf_mod.f90: $(SRC)/scf_mod.fp.f90
	$(FORPEDO) $(DEVEL) $(SRC)/scf_mod.fp.f90 > scf_mod.f90

scf_mod.o: scf_mod.f90 const_mod.o io_mod.o\
	basis_mod.o chemistry_mod.o step_mod.o\
	grid_mod.o chain_mod.o $(FFT_FILE).o grid_basis_mod.o
	$(F90) $(FAST) -c scf_mod.f90

spinodal_mod.o: $(SRC)/spinodal_mod.f90 const_mod.o io_mod.o\
	response_pd_mod.o chemistry_mod.o
	$(F90) $(FAST) -c $(SRC)/spinodal_mod.f90

# Iteration modules

iterate_mod.f90: $(SRC)/iterate_mod.fp.f90
	$(FORPEDO) $(DEVEL) $(SRC)/iterate_mod.fp.f90 > iterate_mod.f90

iterate_mod.o: iterate_mod.f90 const_mod.o\
	scf_mod.o basis_mod.o chemistry_mod.o unit_cell_mod.o\
	response_pd_mod.o
	$(F90) $(FAST) -c iterate_mod.f90

response_pd_mod.o: $(SRC)/response_pd_mod.f90 const_mod.o io_mod.o\
	chemistry_mod.o basis_mod.o unit_cell_mod.o scf_mod.o\
	grid_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/response_pd_mod.f90

sweep_mod.o: $(SRC)/sweep_mod.f90 const_mod.o io_mod.o\
	chemistry_mod.o unit_cell_mod.o basis_mod.o
	$(F90) $(FAST) -c $(SRC)/sweep_mod.f90

# Linear response modules

response_step_mod.o: $(SRC)/response_step_mod.f90\
	chemistry_mod.o $(FFT_FILE).o
	$(F90) $(FAST) -c $(SRC)/response_step_mod.f90

extrapolate_mod.o: $(SRC)/extrapolate_mod.f90\
	const_mod.o
	$(F90) $(FAST) -c $(SRC)/extrapolate_mod.f90

response_mod.o:	$(SRC)/response_mod.f90\
	chemistry_mod.o const_mod.o chain_mod.o\
	grid_mod.o $(FFT_FILE).o group_mod.o response_step_mod.o\
	extrapolate_mod.o field_io_mod.o spinodal_mod.o
	$(F90) $(FAST) -c $(SRC)/response_mod.f90

