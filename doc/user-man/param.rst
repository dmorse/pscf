
.. _param-page:

**************
Parameter File
**************

The main program reads an parameter file containing the parameters and
instructions for a calculation. This file  is divided into sections,
each of which contains a different type of information.  Each section
is preceded by a blank line and starts with a line containing a
section title in all capital letters (i.e., 'CHEMISTRY', 'UNIT_CELL',
etc.) Each block may contain values of a sequence of variables. The
name of each variable appears on a line by itself, followed by the
value or (for arrays) values on one more subsequent lines.  The
order in which variables must appear within a section is fixed. The
program stops when it encounters the block title 'FINISH'.

.. _example-sec:

Example
=======

An example of a complete parameter file is shown below. This example is
for a system containing a triblock copolymer containing three chemically
distinct blocks in a solvent that is chemically identical to one of
the blocks. The first line identifies the version of the file format
(in this case, version 1.0).  The remainder of the file is divided into
sections, each of which begins with a line containing a capitalized label,
such as MONOMERS, CHAINS, ec. The first few sections in this example
simply provide blocks input data. The ITERATE and SWEEP sections instead
contain the instructions required to initiate a computation. Execution
stops when a FINISH block is encountered.

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

   BASIS
   group_name
              '-1'

   ITERATE
   input_filename
        'in.omega'
   output_prefix
            'out/'
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


The MONOMERS block contains information about the monomers used in this
calculation, including the number N_monomer of monomer types and the
statistical segment length of each type, given as elemetns of the
one-dimensional array named "kuhn".

The CHAINS block describes the structure and composition of all polymer
chains, which must linear block polymers or hompolymers.

.. _param-overview-sec:

Overview of Sections
====================

**Primary Sections**

The following list shows the titles of the blocks required to complete most
standard computations, in the order in which they normally appear.
Subsequent sections describe each of the corresponding blocks of the input
file in detail. To solve the SCF problem for a single set of parameters,
leave out the penulimate SWEEP section.

  ===============================  ====================================================
  Section                          Description
  ===============================  ====================================================
  :ref:`param-monomers-sub`        # of monomers and kuhn lengths
  :ref:`param-chains-sub`          Chain species, block sequences and lengths, etc.
  :ref:`param-solvents-sub`        Solvent species, chemical identities, volumes
  :ref:`param-composition-sub`     Statistical ensemble and mixture composition
  :ref:`param-interaction-sub`     Interaction parameters (excess free energy)
  :ref:`param-unitcell-sub`        Unit cell dimension, lattice type, and parameters
  :ref:`param-discretization-sub`  Spatial grid dimensions and 'time' step ds.
  :ref:`param-basis-sub`           Construct symmetry adapted basis 
  :ref:`param-iterate-sub`         Solve SCFT for one set of parameters
  :ref:`param-sweep-sub`           Solve SCFT for multiple sets of parameters
  :ref:`param-response-sub`        Compute linear susceptibility of ordered phase
  :ref:`param-finish-sub`          Stop program
  ===============================  ====================================================

Several standard types of computation are possible using the blocks listed above:

   - Iterate: To solve solve SCF equations for a single state point, include
     all of the listed below sections except the SWEEP and RESPONSE sections.

   - Sweep: To compute a sequence of different states along a line in parameter
     space, include both an ITERATE and SWEEP function, but not a RESPONSE
     section. The ITERATE section must precede the SWEEP section, and is used
     to obtain a solution for the initial choice of parameters.

   - Response: To compute the self-consistent-field or RPA linear susceptibility of a
     periodic microstructure, include ITERATE and RESPONSE sections, but do not include
     a SWEEP section.

**Miscellaneous Utilities**

The following sections are used to invoke a variety of data processing operations or
transformations on fields or parameters, or to output additional information.

  ============================== ===============================================
  Section                        Description
  ============================== ===============================================
  :ref:`param-fieldtorgrid-sub`  Read field file in symmetry-adapated format
                                 and output file in coordinate grid format
  :ref:`param-rgridtofield-sub`  Read field in coordinate grid file format
                                 and output in symmetry-adapated format
  :ref:`param-kgridtorgrid-sub`  Read field in k-space and output in r-space
  :ref:`param-rhotoomega-sub`    Read rho field, compute and output omega field
  :ref:`param-rescale-sub`       Redefine monomer reference volume 
  :ref:`param-waves-sub`         Output map of waves to basis functions
  :ref:`param-group-sub`         Output all elements of space group
  ============================== ===============================================

Further details about the contents and purpose of each section are given below.

.. _param-conventions-sec:

Parameter Conventions
======================

PSCF does not impose the use of a particular system of units
for lengths. Any system of units can be used for entering values
of the monomer statistical segment lengths and the unit cell
dimensions, as long as the same unit of length are used for all
relevant quantities.  One can use either a physical unit, such
as nanometers or Angstroms, or dimensionless units in which one
or more of the statistical segment lengths is set to unity.


SCFT also leaves the user some freedom to redefine what he or
she means by a "monomer", which need not correspond to a chemical
repeat unit.  The choice of values of the parameters block_length,
solvent_size, kuhn, and chi to represent a particular experimental
system all depend on the choice of a value for a reference volume
used to define an effective repeat unit.  Each element of the
variable block_length represents the number of "monomers" in a
block of a block copolymer, defined to be the ratio of the block
volume to the chosen reference volume.  Similarly, the variable
solvent_size is given by ratio of the solvent volume to the
reference volume. The values of the chi parameters are proportional
to the reference volume, while kuhn lengths are proportional to
the square root of the reference volume.  Note that PSCF does not
require the user to input a value for the monomer reference volume
- the choice only effects the values required for other quantities.

All parameters that are represented internally as characters or
character strings must appear in the parameter file with single
quotes, e.g., as 'chi' or 'out.'.

.. _param-array-sec:

Array Formats
==============

Some input parameters are one or two-dimensional array. Here, we discuss how
the dimension and format of these parameters is indicated in subsequent sections
that describe the parameters required in different sections of the input
script.

Below, the discussion of possible section of an parameter file contains a table
listing the required parameters and meaning. One or two-dimensional parameters
are indicated in these tables by displaying the name of each array variable
with an appropriate number of indices.  One dimensional parameters are thus
indicated by writing the name of the parameter with one index: For example,
in the description of the MONOMERS section, kuhn(im) denotes a one dimensional
array of statistical segment lengths for different monomer types.  Two
dimensional arrays are shown with two indices.

The meaning and range of each such array index is indicated by using a set of
standard variable names to indicate different types of indices, with different
ranges of allowed values. For example, in the remainder of this page, the
symbol 'im' is always used to indicates an index for a monomer type.  The
meaning and range of every index symbol is summarized in the following table:

Meaning of Array Indices:

  ========= =====================  ================
  Indices   Meaning                Range
  ========= =====================  ================
  im, in    monomer types          1,...,N_monomer
  ic        chain/polymer species  1,...,N_chain
  ib        blocks within a chain  1,...,N_block(ic)
  is        solvent species        1,...,N_solvent
  id        Cartesian direction    1,...,dim
  ========= =====================  ================

For each array parameter, the elements of the array are expected to appear
in the parameter file in a specific format. Generally, arrays that contain
a polymer or solvent molecular species index are input with the required
information about each molecule on a separate line, while values
associated with different monomer types or with different blocks within
a molecule are listed sequentially on a single line. The expected format
for each array parameter in specified by a code labeled "Format" in each
the table of parameters for each section. The meaning of each array format
code is specified below:

Array Format Codes:

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
a symmetric array with zero diagonal elements is input in the form::

   chi(2,1)
   chi(3,1) chi(3,2)
   .....

in which line i contains elements chi(i+1,j) for j< i. For a
system with only two monomer types (e.g., a diblock copolymer melt
or a binary homopolymer blend), only the single value chi(2,1) on
a single line is required.

.. _param-sections-sec:

Individual Sections
====================

Each of the following subsections describes the format of one possible
section of the parameter file. Array-valued parameters are indicated using
the conventions described above.  Some variables may be present or absent
depending on the value of a previous variable.  These conditions, if any,
are given in a column entitled 'Required if' or 'Absent if'.


.. _param-monomers-sub:

MONOMERS
--------

Chemistry Parameters

  ===========  ========  =========================================   ==========
  Variable     Type      Description                                 Format
  ===========  ========  =========================================   ==========
  N_monomer    integer   Number of monomer types
  kuhn(im)     real      statistical segment length of monomer im    R
  ===========  ========  =========================================   ==========

.. _param-chains-sub:

CHAINS
------

Chain Parameters

  ==================== ======== ============================================ ======
  Variable             Type     Description                                  Format
  ==================== ======== ============================================ ======
  N_chain              integer  Number of chain species
  N_block(ic)          integer  Number of blocks in species ic               C
  block_monomer(ib,ic) integer  Monomer type for block ib of species ic      MR
  block_length(ib,ic)  real     Number of monomers in block ib of species ic MR
  ==================== ======== ============================================ ======

The block_monomer and block_length arrays are entered in a format in which each
line contains the data with one polymer species, so that the number of entries
in line ic must equal to the value of N_block(ic), i.e., to the number of blocks
in chain species ic. The length of each block in an incompressible mixture is
equal to the volume occupied by that block (computed using the density of the
corresponding hompolymer) divided by the monomer reference volume.

.. _param-solvents-sub:

SOLVENTS
--------

Solvent Parameters

  ==================== ======== ============================= ======
  Variable             Type     Description                   Format
  ==================== ======== ============================= ======
  N_solvent            integer  Number of solvent species
  solvent_monomer(is)  integer  Monomer type for solvent is   C
  solvent_size(is)     real     Volume of solvent is          C
  ==================== ======== ============================= ======

The parameter solvent_size is given by the ratio of the actual volume
occupied by a particular solvent to the monomer reference volume.

.. _param-composition-sub:

COMPOSITION
-----------

Composition Parameters:

  =============== ======== ========================================= ======
  Variable        Type     Description                               Format  
  =============== ======== ========================================= ======
  ensemble        integer  0 if canonical, 1 if grand
  phi_chain(ic)   real     volume fraction of chain species ic       C       
  phi_solvent(is) real     volume fraction of solvent species is     C       
  mu_chain(ic)    real     chemical potential of chain species is    C       
  mu_solvent(ic)  real     chemical potential of solvent species ic  C       
  =============== ======== ========================================= ======

The integer parameter "ensemble" determines the choice of statistical ensemble, 
and should be set to 0 for canonical (NVT) ensemble and to 1 for grand-canonical
ensemble. The remainder of the section then contains only the input parameters
required in the specified ensemble: If canonical ensemble is specified (ensemble=0), 
then the rest of the section must contain values for the parameters phi_chain and 
(if N_solvent > 0) phi_solvent that specify the volume fractions of all species.
The example parameter file shows this for a canonical ensemble simulations of a
single-component polymer melt.  If grand canonical ensemble is specified (ensemble=1), 
then the rest of the section must contain values for the parameters mu_chain and 
(if N_solvent > 0) mu_solvent that specify values for the chemical potentials of 
all species. Chemical potentials are specified as free energies per molecule in 
units with :math:`k_{B}T=1`. Values of phi_solvent (in canonical ensemble) or
mu_solvent (in grand-canonical ensemble) should be given if and only if there
are solvent species present, i.e., if N_solvent > 0.

.. _param-interaction-sub:

INTERACTION
-----------

Interaction Parameters

  ============ ======= ================================= ======  
  Variable     Type    Description                       Format  
  ============ ======= ================================= ======  
  chi_flag     char(1) 'B' => bare chi,
                       'T' => chi=chi_A/T + chi_B
  chi(im,in)   real    Flory-Huggins parameter ('bare')  LT      
  chi_A(im,in) real    Enthalpic coefficient for chi(T)  LT      
  chi_B(im,in) real    Entropic contribution to chi(T)   LT      
  Temperature  real    Absolute temperature                       
  ============ ======= ================================= ======

The parameter "chi_flag" determines whether the Flory-Huggins interation 
parameters should be input by specifying values, if chi_flag = 'B', or by
specifying a temperature dependence of the form A/T + B, if chi_flag = 'T'.
The array chi should be present if and only if chi_flag = 'B', while the
parameters chi_A and chi_B should be present if and only if chi_flag = 'T'.

.. _param-unitcell-sub:

UNIT_CELL
---------

The variables in the UNIT_CELL section contain the information necessary to define
the unit cell type, and the unit cell dimensions and shape.


  ================ ============== ============================================ ======
  Variable         Type           Description                                  Format
  ================ ============== ============================================ ======
  dim              integer        dimensionality =1, 2, or 3
  crystal_system   character(60)  unit cell type (cubic, tetragonal, etc.)
  N_cell_param     integer        # parameters required to describe unit cell
  cell_param(i)    real           N_cell_param unit cell parameters            R
  ================ ============== ============================================ ======

The array cell_param contains N_cell_param elements, which are input in row format,
with all elements in a single line. Further information about the allowed values of
the crystal_system string and the number and type of parameters required by each
type of lattice is given in the :ref:`unitcell-page`  page.


.. _param-discretization-sub:

DISCRETIZATION
--------------

The discretization section defines the grid used to spatially discretize
the modified diffusion equation and the size ds of the "step" ds in the
time-like contour length variable used to integral this equation.

Parameters

  ========= ========  ====================================== ====
  Variable  Type      Description                            Form
  ========= ========  ====================================== ====
  ngrid(id) integer   # grid points in direction id=1,..,dim  R
  ds        real      contour length step size
  ========= ========  ====================================== ====

The integer array ngrid(id) is input in row format, with dim (i.e., 1,2 or 3)
values on a line, where dim is the dimensionality of space.

.. _param-basis-sub:

BASIS
-----

The BASIS block instructs the code to construct symmetrized basis 
functions that are invariant under the operations of a specified space 
group.  The file format for this block contains only one variable, 
named "group", which is a string identifier for the space group. 
The value of the "group" string can be either a standard name of 
one of the possible space groups or the path to a file that 
contains the elements of the group. The names of all possible space 
groups, in the form expected by PSCF, are presented in the page on 
:ref:`group-page`.

  ======== =============  ==========================
  Variable Type           Description
  ======== =============  ==========================
  group    character(60)  group name, or file name
  ======== =============  ==========================

.. _param-iterate-sub:

ITERATE
-------

The ITERATE command causes the program to read in an input omega file and 
then attempt to iteratively solve the SCFT equations for one set of input 
parameters. This is the workhorse of a SCFT computation. An ITERATE
section must immediately precede any SWEEP or RESPONSE section. 

If an ITERATE section is immediately preceded by a RESCALE section, it 
uses the rescaled version of the omega field that was read by the RESCALE 
command.  In this case the parameter file should not contain an 
input_filename parameter.

Parameters:

  ============== ============= =================================================
  Variable       Type          Description
  ============== ============= =================================================
  input_filename character(60) input omega file name
  output_prefix  character(60) prefix to all output files
  max_itr        integer       maximum allowed number of iterations
  max_error      real          tolerance - max. norm of residual
  domain         logical       unit cell is variable if true, rigid if false
  itr_algo       character(10) code for iteration algorithm
  N_cut          integer       dimension of cutoff Jacobian in NR algorithm
                               (required iff itr_algo = 'NR')
  N_hist         integer       Number of histories used in AM algorithm
                               (required iff itr_algo = 'AM')
  ============== ============= =================================================

Discussion:

The string "output_prefix" is concatenated with the suffixes 'out', 'rho', 
and 'omega' to create paths (file names) for the output summary, output 
monomer concentration (rho) field, and output chemical potential (omega) 
field files.  The output prefix string should usually be either the name 
of a subdirectory followed by a "/" directory separator string, such as 
'out/', in order to place these files in a separate directory, or a string 
that ends with a period, such as 'out.', to obtain files with file 
extensions '.out', '.rho' and '.omega'. In all of the examples, we set 
output_prefix = 'out/' to place all output files in a subdirectory.

The value of the "domain" logical parameter determines whether PSCF 
attempts to solve the self-consistent field equations in a fixed unit 
cell (if domain == F) or whether it adjusts the parameters of the unit 
cell so as to find a state of vanishing stress, and thus minimum free
energy (if domain == "T").

The value of the string "itr_algo" determines the choice of iteration
algorithm. The only valid values (thus var) are "NR" or "AM". 

If "itr_algo" is "NR", PSCF uses a quasi-Newton-Raphson iteration 
algorithm that is unique to this program. This algorithm constructs 
a physically motivated initial approximation for the Jacobian matrix
in which elements associated with long wavelength components of the
:\math:`\omega` field are computed numerically and shorter wavelength 
components are estimated. After construction and inversion of this 
initial estimate, Broyden updates of the inverse Jacobian are used to 
refine the estimate of the inverse Jacobian. This method requires a 
parameter "N_cut", which determines how many rows and columns of the 
Jacobian matrix are to be computed numerically. The time required to
construct the initial estimate of the Jacobian, which can become quite 
long for 3D problems that require many basis functions, increases 
linearly with "N_cut" . For problems involving relatively simple 3D 
unit cells of block copolymer melts, values of N_cut of order 100 
often provide a reasonable balance between accuracy and cost. One
important disadvantage of the "NR" algorithm is that it requires
storage of the full Jacobian matrix, which can become impossible for 
problems with more than about 10,000 basis functions.

If "iter_algo" is set to "AM", PSCF using an Anderson mixing algorithm
that uses much less memory. This algorithm requires an integer parameter 
"N_history" that determines how many previous iterations are stored and 
used to estimate each update. We often set N_history = 30.

.. _param-sweep-sub:

SWEEP
-----

The presence of a SWEEP section instructs the program to solve the SCFT 
for a sequence of nearby values of parameters along a line through 
parameter space (a 'sweep'). We define a sweep contour variable s that 
varies from 0 up to a maximum value s_max, in increments of 1. For each 
integer step in the sweep parameter, the user may specify a fixed increment
per step for any of the real parameters that are relevant to the problem.
The parameters that can be incremented include all of the real parameters 
in the MONOMERS, CHAINS, SOLVENTS, COMPOSITION, and INTERACTION section 
(i.e., all parameter in these sections for which a floating point value 
or an array of floating point values is given in the parameter file). For 
simulations with a fixed unit cell (domain=1), the elements of the 
unit_cell_param array may also be incremented. 

The desired increment per step for any variable <name> is specified by 
the value or (for an array) array of values of a corresponding increment 
variable named d_<name>. Any number of increments may be specified.
Variables that are not incremented do not need to be referred to explicitly -
increments of zero are assigned default. When an array-valued variable is 
incremented, however, increment values must be specified for all of the 
elements of the array.  The reading of increment variables ends when the 
program encounters the line containing the string "end_increments".

  ============== =============== =======================================
  Variable       Type            Description
  ============== =============== =======================================
  s_max          real            maximum value of sweep contour variable
  s_<name>       type of <name>  increment in variable <name>
  end_increments none            indicates end of the list of increments
  ============== =============== =======================================

.. _param-response-sub:

RESPONSE
--------

The presence of a RESPONSE section instructs the program to calculate 
the linear response matrix for a converged ordered structure at one or 
more k-vectors in the first Brillouin zone. If the linear response is 
calculated for more than one k-vector, they must lie along a line in 
k-space, separated by a user defined vector increment.

  ========= ===========  =====================================
  Variable  Type         Description
  ========= ===========  =====================================
  pertbasis char         If 'PW' => plane wave basis.
                         If 'SYM' => symmetrized basis functions
  k_group   character    Group used to construct symmetrized
                         basis functions
  kdim      int          # dimensions in k-vector (kdim >= dim)
  kvec0(i)  real         initial k-vector, i=1,...,kdim
  dkvec(i)  real         increment in k-vector
  nkstep    integer      # of k-vectors
  ========= ===========  =====================================

.. _param-fieldtorgrid-sub:

FIELD_TO_RGRID
--------------

This command reads a file containing a field in the symmetry-adapted
Fourier expansion format and outputs a representation containing
values of the field on a coordinate space grid. This and the other
commands to transform representation can be applied to either a rho
or omega field.

  ================  ============= ============================
  Variable          Type          Description
  ================  ============= ============================
  input_filename    character(60) input file name
                                  (symmetry-adapted format)
  output_filename   character(60) output file name
                                  (coordinate grid format)
  ================  ============= ============================

.. _param-rgridtofield-sub:

RGRID_TO_FIELD
--------------

This command performs the inverse of the transformation performed
by FIELD_TO_RGRID: It reads a file containing values of a field on
the nodes of a coordinate grid and outputs a file containing a
representationo as an symmetry-adapted Fourier expansion.

  ================ ============= ========================================
  Variable         Type          Description
  ================ ============= ========================================
  input_filename   character(60) input file name
                                 (coordinate grid) 
  output_filename  character(60) output file name
                                 (symmetry-adapted)
  ================ ============= ========================================

.. _param-kgridtorgrid-sub:

KGRID_TO_RGRID
--------------

This command inverts the operation applied by RGRID_TO_KGRID: It reads
a file containing values Fourier components of a field on wavevectors
on a k-space FFT grid, performs an inverse Fourier transform, and
outputs values of the field on a coordinate r-space grid.

  ================ ============= ============================
  Variable         Type          Description
  ================ ============= ============================
  input_filename   character(60) input file name
                                 (wavevector grid)
  output_filename  character(60) output file name
                                 (coordinate grid)
  ================ ============= ============================

.. _param-rhotoomega-sub:

RHO_TO_OMEGA
--------------

This command reads a file containing a monomer concetnration field
and outputs a corresponding initial guess for the omega field. Both
input and ouput files use the symmetry-adapted Fourier expansion
format. The omega field is computed by simply setting the Lagrange
multiplier pressure field to zero, giving a field that only contains
the contributions that arise from the excess interaction free
energy, e.g., terms that explicitly involve the Flory-Huggins chi
parameter. This command is intended to be used to generate an initial
guess for $\omega$ from an approximate structural model for the
volume fraction fields in a particular structure.

  ================  ============= ============================
  Variable          Type          Description
  ================  ============= ============================
  input_filename    character(60) input rho file name
                                  (symmetry-adapted)
  output_filename   character(60) output omega file name
                                  (symmetry-adapted)
  ================  ============= ============================

.. _param-rescale-sub:

RESCALE
-------

This command reads in an omega file, then applies a change in the
convention for the monomer reference volume to the omega field and 
to all parameters whose value depend upon an implicit choice of 
monomer reference volume. This command may only be called (if at all)
immediately prior to an ITERATE commands, in order to read in an
omega field and then change the convention for the monomer reference 
volume prior to solving the SCFT equations.

This command applies a change in the omega field and various 
properties that corresponds to a change of the monomer reference 
volume :math:`v` by a factor :math:`v \rightarrow v/\lambda`. The
scale factor :math:`\lambda` is given in the parameter file by 
the input variable "vref_scale".

  ================  ============= ============================
  Variable          Type          Description
  ================  ============= ============================
  input_filename    character(60) input omega file name
  vref_scale        real          scale factor
  ================  ============= ============================

This command applies the following set of transformations to each
block length :math:`N`, solvent size :math:`S`, statistical segment 
length :math:`b`, Flory-Huggins interaction parameter :math:`\chi`,
and monomer chemical potential field :math:`\omega`:

   ==================  ==============  ========================
   Variable type       Symbol          New value
   ==================  ==============  ========================
   block length        :math:`N`       :math:`N \lambda`
   solvent size        :math:`S`       :math:`S \lambda`
   monomer length      :math:`b`       :math:`b/\sqrt{\lambda}`
   interaction         :math:`\chi`    :math:`\chi/\lambda`
   field               :math:`\omega`  :math:`\omega/\lambda`
   ==================  ==============  ========================

The SCFT equations can be shown to be invariant under such a change 
in convention for the definition of a "monomer". Also note that 
this transformation leaves invariant any product :math:`\chi N` of 
a interaction parameter and a block or a chain length or any product 
:math:`\omega N` of a chemical potential field per monomer and the 
number of monomers in a block, both of which correspond to measures 
of the free energy of interaction of a block with its surroundings. 
The transformation also leaves invariant any product 
:math:`\sqrt{N} b` that corresponds to a random-walk coil size.

Applying this rescaling to an omega field that already solves the
SCFT equations for the choice of parameters given in the parameter
file simply generates an equivalent solution corresping to a 
rescaled choice of parameter values. Using the RESCALE command to 
read in a file containing such a converged solution should thus 
cause the subsequent ITERATE command to terminate immediately,
since the error should be less than the numerical threshhold on 
and output the new parameters to an output summary file and the 
rescaled omega field to an output omega file. 

.. _param-waves-sub:

OUTPUT_WAVES
------------

Output the relationship between plane waves and symmetry-adapted
basis functions, by outputting a file containing showing which
star each wavevector belongs to and the coefficients of the 
plane-wave within a symmetry adapted basis function assocated 
with that star.

  ================  ============= ============================
  Variable          Type          Description
  ================  ============= ============================
  output_filename   character(60) output file name
  ================  ============= ============================


.. _param-group-sub:

OUTPUT_GROUP
------------

Output all symmetry elements of the current space group to a file.
See the discussion of space group :ref:`group-symmetry-sec` for a
discussion of the internal representation of space groups and the
output file format.

Parameters:

  ================  ============= ============================
  Variable          Type          Description
  ================  ============= ============================
  output_filename   character(60) output file name
  ================  ============= ============================


.. _param-finish-sub:

FINISH
------

The FINISH string is the last section of any parameter file, and
causes program execution to terminate.

