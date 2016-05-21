# = Paths to Script and Python Interpreter (modify as needed)
# PYTHON = python interpreter command (must be in $PATH)
# GDOC   = path to GDoc.py python documentation generation script

SRC=../src
BLD=../make
PYTHON=python
GDOC=../tools/python/GDoc/GDoc.py

# ====== Target to make html user manual  =======================
# Note: This requires the python2.7-sphinx python package

user-man: user-man/*.rst
	cd user-man; make html

# ====== Target to make html developer/api documentation ========
#
.PHONY: devel-man 

devel-man:	
	cp $(BLD)/pscf.f $(SRC)/pscf.f
	cp $(BLD)/scf_mod.f $(SRC)/scf_mod.f
	cp $(BLD)/iterate_mod.f $(SRC)/iterate_mod.f
	$(PYTHON) $(GDOC)
	rm $(SRC)/pscf.f
	rm $(SRC)/scf_mod.f
	rm $(SRC)/iterate_mod.f

.PHONY: clean

clean:
	rm -rf  devel-man/*
	rm -rf  user-man/_build

