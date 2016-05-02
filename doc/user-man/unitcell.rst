
.. _unitcell-page:

**********
Unit Cell
**********

The variables in the UNIT_CELL section of the parameter file specify the Bravais
lattice of the crystal. The variables required in this section, in the required
order, are:

  ===============  ===================================================
  Variable         Meaning
  ===============  ===================================================
  dim              number of periodic directions (1,2, or 3)
  crystal_system   lattice system string identifier
  N_cell_param     the number of required unit cell parameters 
  cell_param       an array of unit cell parameter values
  ===============  ===================================================

The "crystal_system" for a three dimensional crystal string could, for 
example, be "cubic", "tetragonal", "orthorhombic" or several other 
possibilities.  All allowed values of the "crystal_system" string are 
given below 1, 2 and 3 dimensionally periodic systems. The number
N_cell_param of parameters required to specify the unit cell is also
listed for each crystal_system, along with the meaning and order of
appearance of different parameters.  The value of N_cell_param, is, 
for example, 1 for a 1D lamellar phase, a 2D hexagonal phase or a 
3D cubic phase, 3 for an orthorhombic 3D phase or 6 for a triclinic 
3D phase.

Below, we discuss all possible crystal systems of 1, 2, and 3 
dimensional crystals separately. Lists of the names of all possible
space groups for each crystal system are given on the following page.

1D Crystal Systems
==================

The only allowed crystal system name for a one-dimensional crystal is 
s "lamellar". Only one unit cell parameter is required to specify a 
lamellar unit cell (i.e., N_cell_param = 1). The value of that parameter
is equal to the layer spacing. 


2D Crystal Systems
==================

For two dimensional crystals (dim=2), let the parameters a and b 
denote the lengths of two independent Bravais lattice basis vectors. 
For oblique crystals, gamma denotes the angle between these two 
basis vectors, in radians. 

============  ============ ============
systems       N_cell_param cell_param
============  ============ ============
square        1            a

hexagonal     1            a

rectangular   2            a, b

oblique       3            a, b, gamma
============  ============ ============


3D Crystal Systems
===================

For three dimensional crystals (dim=3), the parameters a, b, and c 
are the lengths of the three Bravais lattice vectors, alpha is the 
angle between b and c, beta is the angle between c and a, and gamma 
is the angle between a and b. 

============= ============ ============================
3D systems    N_cell_param cell_param
============= ============ ============================
cubic         1            a
tetragonal    2            a, c
orthorhombic  3            a, b, c
monoclinic    4            a, b, c, beta
hexagonal     2            a, c
trigonal      2            a, alpha
triclinic     6            a, b, c alpha, beta, gamma
============= ============ ============================

