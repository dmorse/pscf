
*******************
Input Script Format
*******************

The main program reads an input script containing the parameters and 
instructions for a calculation. The script is divided into sections, 
each of which contains a different type of information.  Each section
is preceded by a blank line and starts with a line containing a
section title in all capital letters (i.e., 'CHEMISTRY', 'UNIT_CELL', 
etc.) Each block may contain values of a sequence of variables. The 
name of each variable appears on a line by itself, followed by the 
value or (for arrays) values on one more subsequent lines.  The
order in which variables must appear within a section is fixed. The
program stops when it encounters the block title 'FINISH'. 

.. example-sec:

Example
=======

An example of a complete script is shown below. This example is for 
a system containing a triblock copolymer containing three chemically 
distinct blocks in a solvent that is chemically identical to one of 
the blocks. The first line identifies the version of the file format 
(version 1.0, for the format documented here).  The remainder of the 
script is divided into sections, each of which may contain input 
variables, as discussed above. 

::

   format  1  0
   
   MONOMERS
   N_monomer           
                 2
   kuhn                
     0.4000000E+00  0.6000000E+00 
   
   CHAINS
   N_chain              
                 1
   N_block             
                 3
   block_monomer  
                 1              2              3
   block_length   
     1.2000000E+02  0.7000000E+02  0.6000000E+02
   
   SOLVENTS
   N_solvent              
                 1
   solvent_monomer
                 2
   solvent_size
               1.0
   
   COMPOSITION
   ensemble            
                 0
   phi_chain      
     0.8000000E+00
   phi_solvent      
     0.2000000E+00
   
   INTERACTIONS
   interaction_type
               'chi'
   chi                 
     1.2000000E-02
   
   UNIT_CELL
   dim                 
                 1
   crystal_system      
        'lamellar'
   N_cell_param        
                 1
   cell_param          
     1.3200000E+01
   
   DISCRETIZATION
   ngrid
                32
   ds
              1.00
   
   FILE_PREFIXES
   input_prefix        
             'in.'
   output_prefix       
                ''
   
   BASIS
   group_name          
              '-1'
   
   ITERATE
   max_itr             
                20
   error_max           
     1.0000000E-08
   domain              
                 T
   itr_algo
              'NR'
   N_cut
               100
   
   SWEEP
   s_max               
      10.00000E+00
   d_chi
     1.0000000E+00
   end_increments
   
   FINISH

.. _script-sections-sec:

Overview 
========
 
Primary Sections
----------------

The following list shows the titles of the blocks required to calculate 
a 'sweep' of solutions for a sequence of incrementally different different 
parameters, in the order in which they appear in the above example script. 
Subsequent sections describe each of the corresponding blocks of the input 
file in detail. To solve the SCF problem for a single set of parameters,
leave out the penulimate SWEEP section.

  ================================  ==============================================================
  Section                           Description
  ================================  ==============================================================
  :ref:`script-monomers-sub`        # of monomers and kuhn lengths
  :ref:`script-chains-sub`          Chain species, block sequences and lengths, etc.
  :ref:`script-solvents-sub`        Solvent species, chemical identities, volumes
  :ref:`script-solvents-sub`        Statistical ensemble and mixture composition
  :ref:`script-unitcell-sub`        Dimensionality (1,2 or 3), lattice, and unit cell parameters
  :ref:`script-discretization-sub`  Numbers of spatial grid points and 'time' step ds.
  :ref:`script-prefixes-sub`        Prefixes for paths to input and output files
  :ref:`script-basis-sub`           Read space group and construct symmetrized basis functions
  :ref:`script-iterate-sub`         Solve SCFT for one set of parameters
  :ref:`script-sweep-sub`           Solve SCFT for a sweep of consecutive parameters
  :ref:`script-finish-sub`          Stop program
  ================================  ==============================================================
 
Linear Response
----------------

To calculate the self-consistent-field linear susceptibility of a
periodic microstructure, introduce a RESPONSE section after ITERATE 
and before FINISH, and leave out the SWEEP section.



 ===============  =========================================================
 Section          Description
 ===============  =========================================================
 :ref:`response`  Calculate linear response matrix at one or more k vectors
 ===============  =========================================================

Utilities
---------

The following sections invoke actions that are not essential to the operation 
of the program:

  :ref:`rescale_omega`   Output symmetry elements of space group to file
  :ref:`output_waves`    Output plane waves and coefficients used in basis functions
  :ref:`field_to_grid`   Convert field from basis function coefficients to values on a grid
  :ref:`rescale_omega`   Redefine monomer reference volume, rescale parameters and omega

Details about the contents and purpose of each section are given below.


Format and Variable Conventions
===============================

The following sections contains a table listing the parameters that must be input 
in a corresponding section of the input file. A few comments about how to read 
these tables:

Array Parameters and Indices
----------------------------

Some required input parameters are one or two dimensional arrays. 
Array parameters are indicated by writing the name of the parameter 
with an index: For kuhn(im) denotes a one dimensional array of statistical 
segment lengths for different monomer types.  The symbol 'im' indicates 
an index for monomer type. Two dimensional arrays are shown with two
indices. The meaning and range of each type of index is summarized 
in the following table:

Meaning of Array Indices
------------------------

  ========= =====================  ================
  Indices   Meaning                Range   
  ========= =====================  ================
  im, in    monomer types          1,...N_monomer
  ic        chain/polymer species  1,...N_chain
  ib        blocks within a chain  1,...N_block(ic)
  is        solvent species        1,...N_solvent
  ========= =====================  ================
 
Array Parameter Formats
-----------------------

For array parameters, the input script expects the elements of the array 
to appear in the input script in a specific format. Generally, arrays that 
contain a polymer or solvent molecular species index are input with the 
required information about each molecule on a separate line, while values 
associated with different monomer types or with different blocks within a 
molecule are listed sequentially on a single line.  The expected format for 
each array parameter in specified by a code in the 'Format' column of the 
table. The meaning of each code is described below:

  =======  ==================================================
  Format   Meaning   
  =======  ==================================================
  R        1D array, row format (all values in a single line) 
  C        1D array, column format (one value per line) 
  MR       2D array, multiple rows of different length 
  LT       2D array, lower triangular 
  =======  ==================================================

Within each line, values may be separated by any amount of whitespace.
In the row (R) format for 1D arrays, all values appear on a single line 
separated by whitespace. In the column format (C), each value appears on 
a separate line. In the multiple row (MR) format, which is used for the
arrays block_monomer(ib,ic) and block_length(ib,ic), each line of data 
contains the values for all of the blocks of one chain molecule, with 
N_block(ic) values in the line for molecule number ic.
The lower triangular (LT) format for square 2D arrays is used for the
array chi(im,in) of Flory-Huggin interaction parameters. In this format,
a symmetric array with zero diagonal elements is input in the form:

   chi(2,1)
   chi(3,1) chi(3,2)
   .....

in which line i contains elements chi(i+1,j) for j< i. For a 
system with only two monomer types (e.g., a diblock copolymer melt
or a binary homopolymer blend) only a single value on a single line 
is required. 



Conditionally Required Parameters
---------------------------------
Some variables may be present or absent depending on the value
of a previous variable.  These conditions, if any, are given in 
the 'Absent if' column.

Reference Volume
----------------
Values of the parameters block_length, solvent_size, kuhn, and
chi all depend on the choice of a value for a reference volume
used to define an effective repeat unit.  Each element of the 
variable block_length represents the number of "monomers" in a 
block of a block copolymer, defined to be the ratio of the block 
volume to the chosen reference volume.  Similarly, the variable 
solvent_size is given by ratio of the solvent volume to the 
reference volume. The values of the chi parameters are proportional
to the reference volume, while kuhn lengths are proportional to
the square root of the reference volume. The program does not
require a value for the reference volume as an input - the
choice only effects the values required for other quantities.

Length Units
------------

Any units of length can be used for the kuhn lengths and the
unit cell dimensions, as long as the same units are used for 
all quantities with units of length. One can use either a
physical unit, such as nanometers or Angstroms, or 
dimensionless units in which one or more of the statistical 
segment lengths is set to unity. 

Discussion of Example
---------------------

First, let us discuss the 'CHEMISTRY' block of the above example.
This block into sets of variables that contain information about 
monomer properties, polymer properties, solvent properties, and 
mixture composition.

The first set of variables (N_monomer, names, and kuhn) contain 
information about the monomers: N_monomer=3 is the total number 
of monomer types. We associate labels name(1)='A', name(2)='B', 
and name(3)='C' with these monomer types.  The monomers have 
statistical segment lengths of kuhn(1)=0.6 nm, kuhn(2)=0.5 nm, 
and kuhn(3)=0.7 nm.  

The next set of variables specifies the chi parameters. These can
be specified either by giving "bare" values for each interaction
parameter, in which case chi_flag = 'B', as in this example, or
as values of the form chi(i,j) = chiA(i,j)/T + chiB(i,j), in which
case we would set chi_flag = 'T'. Because the matrix chi(i,j) is 
symmetric with zeros on the diagonal, chi parameters are input as 
a lower triangular matrix of the form:

  chi(2,1)
  chi(3,1)  chi(3,2) 
  ....

Only a single value on a single line would be required for 
N_monomer=2.  For systems with chi_flag = 'T', the matrices chi_A 
and chi_B would be input in the same form. For this system, 
chi(1,2) = 0.15, chi(2,3)=0.20, chi(1,3)=0.10. 


The next set of variables describes the structure of the polymers.
The system contains two types of polymer (N_chain=2). Properties
of individual molecules are then given with one molecular per line. 
The two molecule types are a triblock, with with N_block(1)=3, and 
a homopolymer with N_block(2)=1.  The array block_monomer 
specifies the chemical identity of each block in each molecule,
with one molecule per line in the input format.  In the example, 
the copolymer is an ACB copolymer and the homopolymer contains C 
monomers. Array block_length specifies the number of monomers in
each block, in the same format.

Individual Script Sections
==========================

.. _script-monomers-sub:

MONOMERS
--------

Chemistry Parameters

  ============  ========  =========================================   =============
  Variable      Type      Description                                 Format
  ============  ========  =========================================   =============
  N_monomer     integer   Number of monomer types
  kuhn(im)      real      statistical segment length of monomer im    R
  ============  ========  =========================================   =============

.. _script-interaction-sub:

INTERACTION
-----------

Interaction Parameters

  ============ ======= ==================================  ======  ============
  Variable     Type    Description                         Format  Required if
  ============ ======= ==================================  ======  ============
  chi_flag     char(1) 'B' for 'bare' chi values
                       'T' for chi=chi_A/T + chi_B
  chi(im,in)   real    Flory-Huggins parameter ('bare')    LT      chi_flag='B'
  chi_A(im,in) real    Enthalpic coefficient for chi(T)    LT      chi_flag='T'
  chi_B(im,in) real    Entropic contribution to chi(T)     LT      chi_flag='T'
  Temperature  real    Absolute temperature                        chi_flag='T'
  ============ ======= ==================================  ======  ============

.. _script-chains-sub:

CHAINS
------

Chain Parameters

Variable             Type     Description                                  Format
N_chain              integer  Number of chain species
N_block(ic)          integer  Number of blocks in species ic               C
block_monomer(ib,ic) integer  Monomer type for block ib of species ic      MR
block_length(ib,ic)  real     Number of monomers in block ib of species ic MR

.. _script-solvents-sub:

SOLVENTS
--------

Solvent Parameters

Variable
Type
Description
Format


  N_solvent
  integer
  Number of solvent species



  solvent_monomer(is)
  integer
  Monomer type for solvent is
  C


  solvent_size(is)
  real
  Volume of solvent is
  C



.. _script-composition-sub:

COMPOSITION
-----------

Composition Parameters

  =============== ======== ========================================= ======= ============================
  Variable        Type     Description                               Format  Required if:
  =============== ======== ========================================= ======= ============================
  ensemble        integer
  phi_chain(ic)   real     volume fraction of chain species ic       C       ensemble=0 and N_chain > 0
  phi_solvent(is) real     volume fraction of solvent species is     C       ensemble=0 and N_solvent > 0
  mu_chain(ic)    real     chemical potential of chain species is    C       ensemble=1 and N_chain > 0
  mu_solvent(ic)  real     chemical potential of solvent species ic  C       ensemble=1 and N_solvent > 0
  =============== ======== ========================================= ======= ============================

.. _script-unitcell-sub:

UNIT_CELL
---------

The variables in the UNIT_CELL section contain the information
necessary to define the unit cell type, size, and shape.

  ================ ============== ============================================
  Variable         Type           Description
  ================ ============== ============================================
  dim              integer        dimensionality =1, 2, or 3
  crystal_system   character(60)  unit cell type
                                  (cubic, tetragonal, orthorhombic, etc.)
  N_cell_param     integer        # parameters required to describe unit cell
  cell_param(i)    real           N_cell_param unit cell parameters
  ================ ============== ============================================

The array cell_param contains N_cell_param elements, which are input in
row format, with all elements in a single line.

.. _script-discretization-sub:

DISCRETIZATION
--------------

The discretization section defines the grid used to spatially discretize
the modified diffusion equaiton and the size ds of the "step" ds in the
time-like contour length variable used to integral this equation.

DISCRETIZATION 
--------------

Parameters

  ========= ========  ====================================== ====
  Variable    Type    Description                            Form
  ========= ========  ====================================== ====
  ngrid(id) integer  # grid points in direction id=1,..,dim  R
  ds        real      contour length step size
  ========= ========  ====================================== ====

The integer array ngrid(id) is input in row format, with dim
(i.e., 1,2 or 3) values on a line.

.. _script-prefixes-sub:

FILE_PREFIXES
-------------

The FILE_PREFIXES section inputs the prefixes that are used to construct
the names of the input and output files. The input prefix is concatenated
with 'omega' to construct the name of the input file. The output prefix
is concatenated with the suffixes 'out', 'rho', and 'omega' to create
the name of the output summary, output monomer concentration field, and
output omega field files. This, to specify file name 'in.omega', 'out',
'rho', and 'omega' in the current directory, you would set in_prefix to
'in.', and the output prefix to the blank string ''. To specify an input
file from another directory, you would set in_prefix to the path to that
directory, followed by a trailing '/' directory separator.  Both string
variables are required, and must appear in the order listed below.

  ==========  ============= ===========================================
  Variable    Type          Description
  ==========  ============= ===========================================
  in_prefix   character(60) prefix to *omega input file
  out_prefix  character(60)  prefix to *rho, *omega, *out output files
  ==========  ============= ===========================================

.. _script-basis-sub:

BASIS
-----

The BASIS block instructs the code to construct symmetrized
basis functions that are invariant under the operations of
a specified space group.  It contains only one variable,
named "group", which is a string containing either the name
of one of the standard space groups (which are hard coded
into the program) or the path to a file that contains the
elements of the group. After reading this string from file,
basis functions are constructed by the make_basis routine
of module basis_mod.

  ======== =============  ==========================================
  Variable Type           Description
  ======== =============  ==========================================
  group    character(60)  name of group, or file that contains group
  ======== =============  ==========================================

The file format for a group file is determined by the input_group
routine in module group_mod. Some simple 2D examples of the format
are provided in src/tests/group.

.. _script-iterate-sub:

ITERATE
-------

The ITERATE block reads in variables required by our iteration
algorithm, and attempts to iteratively solve the SCFT equations
for one set of input parameters.

  ========= ============= =====================================================
  Variable  Type          Description
  ========= ============= =====================================================
  max_itr   integer       maximum allowed number of iterations
  max_error real          tolerance - maximum error after convergence
  domain    logical       variable unit cell if true, fixed unit cell if false
  itr_algo  character(10) character code for iteration algorithm
  N_cut     integer       dimension of cutoff Jacobian
  ========= ============= =====================================================

For now, the value of the 'itr_algo' variable must be 'NR', for Newton-Raphson.
The variable is included in order to allow us to add other iteration algorithms
in the future.

.. _script-sweep-sub:

SWEEP
-----

The presence of a SWEEP section instructs the program to solve the SCFT for
a sequence of nearby values of parameters along a path through parameter
space (a 'sweep'). We define a sweep contour variable s that varies from 0
up to a maximum value s_max, in increments of 1. For each integer step in the
sweep parameter, each of the relevant parameters in CHEMISTRY section (i.e.,
any parameter for which a floating point value or values are specified in the
input script) may be incremented by a user specified amount. For simulations
with a fixed unit cell (domain=1), the elements of the unit_cell_param array
may also be incremented. The desired increment for any variable <;name&gt;
is specified by the value or (for an array) values of a corresponding
increment variable named d_<;name>. Any number of increments may be specified.
Variables that are not incremented do not need to be referred to explicitly -
increments of zero are assigned default. When an array variable is incremented,
however, increment values must be specified for all of the elements of the
array.  The reading of increment variables ends when the program encounters
the line 'end_increments'.

  ============= =============== =======================================
  Variable      Type            Description
  ============= =============== =======================================
  s_max         real            maximum value of sweep contour variable
  s_<name>      type of <name>  increment in variable <name>
  end_increment none            indicates end of the list of increments
  ============= =============== =======================================

.. _script-response-sub:

RESPONSE
--------

The presence of a RESPONSE section instructs the program to
calculate the linear response matrix for a converged ordered
structure at one or more k-vectors in the first Brillouin
zone. If the linear response is calculated for more than one
k-vector, they must lie along a line in k-space, separated by
a user defined vector increment.

  ========= ===========  =====================================
  Variable  Type         Description
  ========= ===========  =====================================
  pertbasis char         'PW' => plane wave basis
                         'SYM' => symmetrized basis functions
  k_group   character    Group used to construct symmetrized
                         basis functions
  kdim      int          # dimensions in k-vector (kdim >= dim)
  kvec0(i)  real         initial k-vector, i=1,...,kdim
  dkvec(i)  real         increment in k-vector
  nkstep    integer      # of k-vectors
  ========= ===========  =====================================

.. _script-finish-sub:

FINISH
------

The FINISH string causes the program to terminate.

