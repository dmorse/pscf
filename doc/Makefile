# = Paths to Script and Python Interpreter (modify as needed)
# PYTHON = python interpreter command (must be in $PATH)
# GDOC   = path to GDoc.py python documentation generation script

SRC=../src
BLD=../make
PYTHON=python
GDOC=../tools/python/GDoc/GDoc.py

# ====== Target to make html user manual  =======================
# Note: This requires the python2.7-sphinx python package

user-man:
	cd user-man; make html

# ====== Target to make html developer/api documentation ========
#
.PHONY: devel-man 

devel-man:	
	cp $(BLD)/pscf.f90 $(SRC)/pscf.f90
	cp $(BLD)/scf_mod.f90 $(SRC)/scf/scf_mod.f90
	cp $(BLD)/iterate_mod.f90 $(SRC)/iterate/iterate_mod.f90
	$(PYTHON) $(GDOC)
	rm $(SRC)/pscf.f90
	rm $(SRC)/scf/scf_mod.f90
	rm $(SRC)/iterate/iterate_mod.f90

.PHONY: clean

clean:
	rm -rf  devel-man/*
	rm -rf  user-man/_build

