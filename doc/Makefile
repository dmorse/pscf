# = Paths to Script and Python Interpreter (modify as needed)
# PYTHON = python interpreter command (must be in $PATH)
# GDOC   = path to GDoc.py python documentation generation script

SRC=../src
BLD=../bld
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
	cp $(BLD)/pscf_pd.f $(SRC)/pscf_pd.f
	cp $(BLD)/scf_mod.f $(SRC)/scf/scf_mod.f
	cp $(BLD)/iterate_mod.f $(SRC)/iterate/iterate_mod.f
	$(PYTHON) $(GDOC)
	rm $(SRC)/pscf_pd.f
	rm $(SRC)/scf/scf_mod.f
	rm $(SRC)/iterate/iterate_mod.f

.PHONY: clean

clean:
	rm -rf  devel-man/*
	rm -rf  user-man/_build

