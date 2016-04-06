# = Paths to Script and Python Interpreter (modify as needed)
# PYTHON = python interpreter command (must be in $PATH)
# GDOC   = path to GDoc.py python documentation generation script

SRC=../src
BUILD=$(SRC)/build
PYTHON=python
GDOC=../tools/python/GDoc/GDoc.py

# ====== Target to make html api documentation =================
#
.PHONY: api

api:	
	cp $(BUILD)/pscf_pd.f $(SRC)/pscf_pd.f
	cp $(BUILD)/scf_mod.f $(SRC)/scf/scf_mod.f
	cp $(BUILD)/iterate_mod.f $(SRC)/iterate/iterate_mod.f
	$(PYTHON) $(GDOC)
	rm $(SRC)/pscf_pd.f
	rm $(SRC)/scf/scf_mod.f
	rm $(SRC)/iterate/iterate_mod.f

# ======   clean - dummy target to clean up debris =============

.PHONY: clean

clean:
	rm -rf  *~ api

