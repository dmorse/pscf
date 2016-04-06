
***************
Bravais Lattice
***************

The variables in the UNIT_CELL section of the input script specify the Bravais
lattice of the crystal. The variables required in this section, in the required
order, are:

  ===============  ===================================================
  Variable         Meaning
  ===============  ===================================================
  dim              the dimensionality (i.e., number of periodic dimensions)
                   of the crystal, which must be 1, 2, 3.
  crystal_system   the crystal lattice system, which is a string that 
                   can be cubic, tetragonal, etc. for 3D crystals.
  N_cell_param     the number of parameters required to specify the 
                   dimensions and shape of the unit cell 
  cell_param       an array of the unit cell parameters
  ===============  ===================================================

The value of the number of parameters required to describe the unit cell,
N_cell_param, is, for example, 1 for a 1D lamellar phase, a 2D hexagonal 
phase or a 3D cubic phase, 3 for an orthorhombic 3D phase or 6 for a 
triclinic 3D phase.

Below, we discuss each of the possible crystals systems for 1, 2 and 3
dimensional crystals, and specify the meaning and order of the elements
of the parameter array cell_param

1D Crystal Systems
==================

The only allowed crystal system name for a one-dimensional crystal is 
'lamellar', for which N_cell_param = 1. The value of the single element
of cell_param denotes the layer spacing.


2D Crystal Systems
==================

For two dimensional crystals (dim=2), the parameters a and b are
the lengths of the two Bravais lattice basis vectors. For an oblique 
crystal, gamma is the angle between them, in radians. 

============  ============ ============
systems       N_cell_param cell_param
============  ============ ============
square        1            a

rectangular   2            a, b

hexagonal     1            a

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

