.. _field-page:

***********
Field Files
***********

PSCF uses several file formats to describe fields, and uses the same 
set of formats for both "omega" and "rho" fields. In each of the
available file formats, for a system with with N_monomer monomer types, 
a single file contains a description of N_monomer "omega" or "rho" 
fields, each of which is associated with a specific monomer type. 
We refer to such a set of N_monomer fields in what follows as a 
multi-component field. 

PSCF can read, write and intercovert three different file formats for 
multi-component fields, which are based on different mathematical 
representations of a field. These are:

    * :ref:`field-symmetry-sec`
    * :ref:`field-grid-sec`
    * :ref:`field-fourier-sec`

Each of these three format is discussed in a separate section below.  
The symmetry-adaped format is based on generalized Fourier expansion 
using a basis of symmetry-adapted basis functions. This is the default 
file format that is used by the ITERATE command for the input omega 
field and the output omega and rho fields.  The coordinate space grid 
format contains the values of all fields on regular grid of points 
within one unit cell, and is useful as an input to external programs 
that can visualize a structure. 

.. _field-symmetry-sec:

Symmetry-Adapted Format
=======================

This default file format used by the ITERATE command is based on an
expansion of each field as an expansion in terms of basis functions
that exhibit the space group symmetry of the crystal.  In what follows,
we first discuss the underlying mathematical expansion, and then the
file format


**Mathematical Basis**

Consider a system with N_monomer monomer types indexed by an 
integer :math:`\alpha = 1, \ldots,` N_monomer. Let 
:math:`\phi_{\alpha}(\textbf{r})` denote a field (e.g., a volume fraction) 
associated with monomer type :math:`\alpha` . We approximate the field 
:math:`\phi_{\alpha}` as an expansion of the form

.. math::

    \phi_{\alpha}(\textbf{r}) = 
    \sum_{i=1}^{\texttt{N_star}} \phi_{i\alpha} f_{i}(\textbf{r})

in which each function :math:`f_{i}(\textbf{r})` is a real basis 
function, :math:`\phi_{i\alpha}` is an associated real coefficient, and 
:math:`\texttt{N_star}` is the number of basis functions used to
approximate the field. In a symmetry-adapted Fourier expansion of a 
field with a specified space group symmetry, each basis function 
:math:`f_{i}(\textbf{r})` is a real function that is invariant under 
all symmetry elements of the chosen space group, and is also an 
eigenfunction of the Laplacian, such that

.. math::

   -\nabla^{2}f_{i}(\textbf{r}) = \lambda_{i} f_{i}(\textbf{r})

for some :math:`\lambda_{i} \geq 0`. The basis functions form 
an orthogonal basis, which are normalized such that

.. math::

   \frac{1}{V} \int \! d^{D}r \; f_{i}(\textbf{r}) f_{j}(\textbf{r}) 
   = \delta_{ij}

where the integral is taken over one unit cell of a periodic structure 
in D-dimensional space and :math:`V` is the generalized volume (length, 
area or volume) of the unit cell. Here :math:`\delta_{ij}` denotes the
Kronecker delta function, which is defined to be :math:`\delta_{ij} = 1` 
for :math:`i=j` and 0 otherwise

Each symmetry-adapted basis function can be expressed as a superposition 
of plane waves with wavevectors 
:math:`\textbf{k}_{1}, \ldots, \textbf{k}_{M}` 
that are reciprocal lattice vectors of equal equal magnitude, as a sum 
of the form

.. math::

   f_{i}(\textbf{r}) = 
   \sum_{i=1}^{M} a_{i} e^{i\textbf{k}_{i}\cdots\textbf{r}_{i}}

where :math:`M` is the number of associated plane waves, and where the
the coefficients :math:`a_{1}, \ldots ,a_{M}` are complex numbers that 
all have equal magnitude 
:math:`|a_{1}| = |a_{2}| = \cdots |a_{M}| = 1/\sqrt{M}`. 

For any space group with inversion symmetry (i.e., a centrosymmetric 
group), the set of wavevectors associated with a basis function is a 
set of vectors that are related to one another by the point group 
symmetries of the crystal.  For non-centrosymmetric groups, the 
wavevectors associated with a real basis function are related by 
point group symmetries or by inversion (i.e., by the relation 
:math:`\textbf{k} \rightarrow - \textbf{k}`).  A group of reciprocal lattice 
wavevectors that are related by symmetries of the crystal is referred 
to here and in the PSCF source code as a "star".  The number of stars,
denoted by N_star in the corresponding file format, is equal to the
number of basis functions. In a cubic crystal, for example, each star
of wavevectors (and thus each basis functions) correspond to a family 
of wavevectors with integer indices {ijk} that are related to one another 
by permutations and/or sign changes. For example, the basis function 
associated with the {321} star has 48 associated wavevectors that 
include wavevectors with integer indices (1,2,3), (3,2,1), (-3,2,1),
etc, whereas the {200} star in a cubic crystal has only six wavevectors 
:math:`(\pm2,0,0), (0,\pm 2, 0), (0, 0, \pm 2)`.

Each basis function for d-dimensional crystal is uniquely identified 
in the symmetry-adapted file format by a label consisting of d integer 
indices. These indices correspond to the Miller indices for one of the 
wavevectors associated with the basis function.  Thus for example, we 
identity the basis function associated with the {321} family of wavevectors 
in a cubic crystal by a label "3 1 1".  The conventions for choosing 
which plane wave to use to identify the basis function is discussed 
in the comments provided in the source code of crystal_mod (See the html
developers manual for browseable version of this. The set of Miller 
indices output to file corresponds to the value of the variable 
wave_of_star.) 

**Example: 3D Gyroid Phase**

Below is an example of a "rho" file output from a simulation of a 
gyroid phase for a diblock copolymer melt. A section of the middle
of this file has been removed, as indicated by the vertical dots.

::

    format  1  0
   dim                 
                      3
   crystal_system      
                'cubic'
   N_cell_param        
                      1
   cell_param          
       3.6735414146E+00
   group_name          
             'I a -3 d'
   N_monomer           
                      2
   N_star              
                    235
     3.000000000000E-01  7.000000000000E-01       0   0   0     1
    -2.942897932802E-01  2.942897932848E-01       2   1   1    24
    -9.425546329793E-02  9.425546327223E-02       2   2   0    12
    -3.864399409689E-03  3.864399436086E-03       3   2   1    48
    -1.483047814338E-02  1.483047815806E-02       4   0   0     6
    -3.546446264855E-02  3.546446265383E-02       4   2   0    24
     3.138519869858E-02 -3.138519870524E-02       3   3   2    24
     2.003121375277E-02 -2.003121374994E-02       4   2   2    24
     1.572048423239E-02 -1.572048424396E-02       4   3   1    48
    -1.376822797257E-02  1.376822798292E-02       5   2   1    48
    -1.063353913450E-02  1.063353913935E-02       4   4   0    12
            .                   .                 .   .   .    .
            .                   .                 .   .   .    .
            .                   .                 .   .   .    .
    -7.575067702553E-05  7.575067344206E-05      13  13  10    24
    -2.570604494615E-05  2.570604263390E-05      14  12  10    24
    -5.627606758688E-05  5.627606408758E-05      14  14   8     6
     5.879116047898E-05 -5.879115755266E-05      14  14  12     6


**Example: 1D Lamellar Phase**

Below is an example the "rho" file output for a small simulation of a
lamellar phase of a diblock copolymer melt:

::

   format  1  0
   dim                 
                      1
   crystal_system      
             'lamellar'
   N_cell_param        
                      1
   cell_param          
       1.3835952906E+00
   group_name          
                   '-1'
   N_monomer           
                      2
   N_star              
                     21
     5.600000000000E-01  4.400000000000E-01       0     1
     2.179734275940E-01 -2.179734275841E-01       1     2
    -1.523969262415E-02  1.523969262143E-02       2     2
    -5.575240954520E-03  5.575240954490E-03       3     2
     1.108470498335E-03 -1.108470498556E-03       4     2
     1.455449531056E-04 -1.455449530934E-04       5     2
    -6.218980135235E-05  6.218980146350E-05       6     2
    -8.059872486808E-07  8.059872753625E-07       7     2
     2.826732709838E-06 -2.826732713547E-06       8     2
    -2.194238294935E-07  2.194238338772E-07       9     2
    -1.060764766149E-07  1.060764782164E-07      10     2
     1.946388906884E-08 -1.946388995126E-08      11     2
     3.010764186682E-09 -3.010764203812E-09      12     2
    -1.161872573075E-09  1.161872692383E-09      13     2
    -3.137859071779E-11  3.137865228352E-11      14     2
     5.685537948359E-11 -5.685537190418E-11      15     2
    -3.817653721188E-12  3.817577312625E-12      16     2
    -2.332684668702E-12  2.332625641218E-12      17     2
     4.053664853576E-13 -4.051318636739E-13      18     2
     3.071545504276E-14 -3.077687877704E-14      19     2
    -1.475930488937E-13 -4.916067553040E-14      20     1


**Description of Format**

This first part of such a field file is a header that ends with the
parameter N_star, which is the number of basis functions. This is 
followed by a data section that that is N_star rows long. Each row
in the data section contains the coefficients associated with one 
basis function in the symmetry-adapted Fourier expansion described 
above, along with some additional information that identifies the
basis function.

The structure of the header is similar to that of the parameter file.
The first line specifies a file format version number (file format v1.0). 
The rest of the header contains information that is required to interpret 
the field file, including the dimensionality of the structure (1,2, or 3) , 
the crystal system, the unit cell parameters, the space group, the number
of monomer types, and the number of basis functions, denoted here by
N_star. The second example above is for a lamellar structure with 
inversion symmetry, for which the space group symbol is -1. 

The data section contains N_star rows, each of which contains the
coefficients associated with one basis function, along with an identifier
for the basis function. The first N_monomer columns of row i (e.g., the
first two columns, in both of the above examples) contain the 
coefficients associated with different monomer types.  Specifically,
a coefficient :math:`\phi_{i\alpha}` associated with basis function 
:math:`i` and monomer type :math:`\alpha` is given in column 
:math:`\alpha` of row i of this data section.

In the file format for a crystal with dimension d (e.g., d=1 for a 
lamellar phase or d=3 for a gyroid phase) the next d columns, after
the columns containing the expansion coefficients, contain a set of 
d integers that identify each basis function. As discussed above,
these are integer indices for one of the wavevectors in the basis 
function. The last column is the number of wavevectors in an 
associated star of wavevectors, which we will refer to as the
multiplicity.

The first basis function in the symmetry adapted Fourier expansion, 
which is given in the first row of the data section, is always the
spatially homogeneous function :math:`f_{1}(\textbf{r}) = 1`. This 
constant function is associated with the single wavevector 
:math:`\textbf{k} = 0`, and identified in a 3D crystal by a label 
"0 0 0", with multiplicity 1.

The second row in the gyroid example contains the coefficients 
for the basis function associated with the {211} family of
wavevectors, which is identified in columns 3-5 by the label 
"2 1 1".  Because this family contains 24 wavevectors, the last
column lists a multiplicity of 24. The {211} family is the first 
star of non-zero wavevectors from which it possible to construct 
a nonzero basis function that is invariant under all of the 
symmetries of space group :math:`Ia\bar{3}d` of the gyroid 
structure. The stars that can be used to construct a basis 
function are precisely those that satisfy the reflection rules
for allowed reflections in scattering from a particular space
group symmetry, for which the {211} family gives the first 
allowed family of reflections in scattering from a gyroid
crystal. 

Consider the second example, which is 1D lamellar phase with 
inversion symmetry. The first basis function is the constant 
:math:`f_{1}(r)=1`, with a label "0" and a multiplicity 1. All
subsequent basis functions are cosine functions of the form 
:math:`f_{n+1}(r) = \sqrt{2}\cos(k_{n}r)` with 
:math:`k_{n} = 2\pi n/L` for a crystal with period :math:`L`, 
for which we see an integer label n. The multiplicity of each
cosine basis function is 2, as indicated in the last column,
since each such function can be expressed as a superposition 
of two plane waves of wavenumbers :math:`\pm k_{n}`. 

The rules for constructing real basis functions for 
non-centrosymmetric space groups is somewhat more complicated 
than for centrosymmetric groups.  When the group has no 
inversion symmetry, a basis function that is constructed 
by superposing plane waves that are related by symmetry 
elements of the space group will generally not be proportional
to a real function. The simplest example of this is a one 
dimensional crystal with no inversion symmetry (group 1), 
and thus no symmetry elements other than the identity. In
this case, no plane wave is related to any other by symmetry. 
The natural basis functions, from the point of view of symmetry
alone, are single complex exponential plane waves, but these
are complex functions of position.  In order to construct basis 
functions that are real, in this example, one must construct
two real superpositions of each pair of plane waves that are 
related by inversion (which is not a symmetry of the crystal).
The required basis functions in this case are both cosine
and sine functions. More generally, to form real basis 
functions in crystals with no inversion symmetry, we use 
generalizations of the cosine and sign functions that are 
construction by constructing two different superpositions 
of "stars" that are related to one another by inversion. 
Conventions used for doing this are described best in the 
comments in the source code for the basis_mod module.
 
.. _field-grid-sec:

Coordinate Grid Format
=======================

PSCF can also output the values of set of fields (one per 
monomer type) evaluated on all of the grid points of the FFT 
grid that is used to solve the modified diffusion equation.

**Example: 2D Hex Phase of Diblock Copolymer Melt**

Here is example of a converged omega field for a hex phase::

    format  1  0
   dim                 
                      2
   crystal_system      
            'hexagonal'
   N_cell_param        
                      1
   cell_param          
       1.7703537313E+00
   group_name          
              'P 6 m m'
   N_monomer           
                      2
   ngrid               
                     24                  24
          0.340581085      19.518839883
          0.570887775      19.658020087
          1.199229419      19.984609517
          2.070864605      20.233012735
          2.929754416      19.853514300
               .                 .
               .                 .
               .                 .
          0.999219800      19.890258066
          0.570887775      19.658020087


**Description of Format**

Like the others, this file format contains a header section with
crystallographic information followed by a data section. The header
section is similar that for the symmetry adapted format, except that
the last variable is an array "ngrid" of integers giving the number
of grid points in each direction. In this example, because it is a
two-dimensional crystal (dim = 2), this array contains two numbers,
both 24, indicating a grid in which there are 24 grid points along
each axis. To describe a hexagonal phase, we use a non-orthogonal
coordinate system in which each axis is parallel to one of the 
Bravais lattice vectors, which in a hexagonal phase have an angle
of 60 degrees between. 

The data section contains the values of fields associated 
with N_monomer monomer types at grid points given by

.. math::

    \textbf{r}(n_1, \ldots, n_{D}) = \sum_{i=1}^{\textrm{D}}
    \frac{n_{i}}{N_{i}}\textbf{a}_{1}

where $D$ is the dimensionality of the crystal (denoted by "dim" in 
the header and the parameter file), :math:`\textbf{a}_{i}` is a Bravais 
lattice vector, :math:`N_{i}` is the number of grid points along 
direction :math:`i`, :math:`\textbf{a}_{i}`, and $n_{i}$ is an integer 
in the range :math:`0 \leq n_{i} < N_{i}`.  The number of rows in the 
data section is equal to the total number of grid points, and each row 
contains values of all fields at a single grid point. The number of 
columns is equal to the number of monomer types, so that data in column 
:math:`alpha` contains the values of the field associated with monomer 
type :math:`\alpha`. 

Grid points are listed in order using index :math:`n_{1}` as the most 
rapidly varying (innermost) loop index. This is implemented in the 
field_io_mod module, in subroutines output_field_grid and 
input_field_grid as a fortran loop of the form::

   do n3 = 0, ngrid(3) - 1
     do n2 = 0, ngrid(2) - 1
       do n1 = 0, ngrid(1) - 1
          [Read or write data at grid point (n1, n2, n3)]
       enddo
     enddo
   enddo


.. _field-fourier-sec:

Wavevector Grid Format
=======================

Finally, PSCF can read and write the unsymmetrized discrete 
Fourier transform of a multi-component field, which is related
to the values on a grid a by discrete Fourier transform.  The 
required file format is very similar to that used for the 
coordinate space grid. The file consists of a header and a
data section. The format of the header is identical to that
used for the coordinate grid format, and includes a list of
the number of grid points used in each direction, denoted by
ngrid.

The data section contains the Fourier coefficients obtained by a 
discrete Fourier transform of each field at wavevectors given by

.. math::

    \textbf{k}(n_1, \ldots, n_{D}) = \sum_{i=1}^{\textrm{D}}
    n_{i}\textbf{b}_{i}

where :math:`D` is the dimensionality of the crystal (i.e., dim
in the header file), :math:`\textbf{b}_{i}` is a reciprocal lattice 
basis vector, :math:`N_{i}` is the number of grid points along 
direction :math:`i`, and :math:`n_{i}` is an integer in the 
range :math:`0 \leq n_{1} \leq N_{1}/2` for the first index and 
:math:`0 \leq n_{i} \leq N_{i} - 1` for indices :math:`i > 1`. 
The number of rows in the data section is equal to the total 
number of such wavevectors, and each row contains values of 
Fourier coefficients associated with a single wavevector, 
with coefficients for fields associated with different monomer 
types in different columnns. 

Coefficients for different wavevectors are output in sequential
order, using the last index (e.g., :math:`n_{3}` for a 3D crystal) 
as the most rapidly varying (inner-most) loop index. This is 
implemented by a fortran loop of the form::

   do n1 = 0, ngrid(1)/2
     do n2 = 0, ngrid(2) - 1
       do n3 = 0, ngrid(3) - 1
         [Read or write coefficients for (n1, n2, n3)]
        enddo
     enddo
   enddo

