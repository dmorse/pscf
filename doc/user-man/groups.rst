
.. _group-page:

*************
Space Groups
*************

The symbol for a space group may be entered as the value of "space_group" in 
the BASIS section of the parameter file. The tables below list the allowed 
space group symbols. 

1D Space Groups
===============

The only possible nontrivial symmetry for a one-dimensional lamellar 
phase is inversion symmetry. There are thus only two possible groups: 
The centrosymmetric group (group -1) and the non-centrosymmetric group
(group 1). Fields for a centrosymmetric lamellar phase are expanded
using a basis of cosine waves, while fields for a non-centrosymmetric 
phase are expanded using a basis that contains both cosines and 
sine waves.

======== ======  =================
Number   Symbol  Comments
======== ======  =================
1        -1      Inversion symmetry
2         1      No symmetry
======== ======  =================


2D Space Groups
===============

The names of all 17 possible 2D plane groups are given below in the text format 
expected by PSCF. The format used in PSCF for both 2D and 3D space group names 
is based on the names used in the international tables of crystallography, but 
allows space group names to be written as simple ascii text strings, which
contain spaces between elements of the space group name.

 ====== ======== ==============
 Number Symbol   Lattice System
 ====== ======== ==============
 1      p 1      oblique
 2      p 2      oblique
 3      p m      rectangular
 4      p g      rectangular
 5      c m      rectangular
 6      p 2 m m  rectangular
 7      p 2 m g  rectangular
 8      p 2 g g  rectangular
 9      c 2 m m  rectangular
 10     p 4      square
 11     p 4 m m  square
 12     p 4 g m  square
 13     p 3      hexagonal
 14     p 3 m 1  hexagonal
 15     p 3 1 m  hexagonal
 16     p 6      hexagonal
 17     p 6 m m  hexagonal
 ====== ======== ==============

3D Space Groups
===============

The names of all possible 3D space groups are given below in the text format 
expected by PSCF. These names are based on the names given in Hermann-Mauguin
or "international" notation used the international tables of crystallography, 
but are given in a format that allows space group names to be written as simple 
ascii text strings, with no special symbols or subscripts. In this format, for
example, the space group :math:`Ia\overline{3}d: of the gyroid phase (space 
group 230) is written as "I a -3 d". 

The following rules are applied to convert Hermann-Mauguin symbols into text 
strings:

   * A single space is introduced between different elements of the space 
     group name, with a few exceptions described below. 

   * Integers with overbars in the Hermann-Mauguin symbol, which indicate
     inversion (:math:`\overline{1}`) or a 3-, 4- or 6-fold rotoinversion 
     axis (:math:`\overline{3}`, :math:`\overline{4}`, or :math:`\overline{6}`), 
     are indicated in the PSCF text string by placing a "-" sign before 
     the overbarred integer. Thus for, example, "-3" represents the symbol
     :math:`\overline{3}` in the text space group name "I a -3 d"

   * Integers with subscripts, such as :math:`4_2`, which indicate screw 
     axes, are indicated in the text representation by placing the two 
     integers directly adjacent, with no intervening white space. Thus, 
     for example, :math:`4_2` is replaced by "42".

   * Symbols that are separated by a slash appear with no space on either 
     side of the slash. 

   * Different "settings" of the same space group, which correspond to 
     different definitions of the origin of space in the definition of
     the symmetry elements, are indicated by a colon followed by an 
     integer label at the end of the space group. 


 ========  =================
  Number   Symbol 
 ========  =================
    1      P 1 
    2      P -1 
    3      P 1 2 1 
    4      P 1 21 1 
    5      C 1 2 1 
    6      P 1 m 1 
    7      P 1 c 1 
    8      C 1 m 1 
    9      C 1 c 1 
   10      P 1 2/m 1 
   11      P 1 21/m 1 
   12      C 1 2/m 1 
   13      P 1 2/c 1 
   14      P 1 21/c 1 
   15      C 1 2/c 1 
   16      P 2 2 2 
   17      P 2 2 21 
   18      P 21 21 2 
   19      P 21 21 21 
   20      C 2 2 21 
   21      C 2 2 2 
   22      F 2 2 2 
   23      I 2 2 2
   24      I 21 21 21 
   25      P m m 2 
   26      P m c 21 
   27      P c c 2 
   28      P m a 2 
   29      P c a 21 
   30      P n c 2 
   31      P m n 21 
   32      P b a 2 
   33      P n a 21 
   34      P n n 2 
   35      C m m 2 
   36      C m c 21 
   37      C c c 2 
   38      A m m 2 
   39      A b m 2 
   40      A m a 2 
   41      A b a 2 
   42      F m m 2 
   43      F d d 2 
   44      I m m 2 
   45      I b a 2 
   46      I m a 2 
   47      P m m m 
   48      P n n n : 2 
   48      P n n n : 1 
   49      P c c m 
   50      P b a n : 2 
   50      P b a n : 1 
   51      P m m a 
   52      P n n a 
   53      P m n a 
   54      P c c a 
   55      P b a m 
   56      P c c n 
   57      P b c m 
   58      P n n m 
   59      P m m n : 2 
   59      P m m n : 1 
   60      P b c n 
   61      P b c a 
   62      P n m a 
   63      C m c m 
   64      C m c a 
   65      C m m m 
   66      C c c m 
   67      C m m a 
   68      C c c a : 2 
   68      C c c a : 1 
   69      F m m m 
   70      F d d d : 2 
   70      F d d d : 1 
   71      I m m m 
   72      I b a m 
   73      I b c a 
   74      I m m a 
   75      P 4 
   76      P 41 
   77      P 42 
   78      P 43 
   79      I 4 
   80      I 41 
   81      P -4 
   82      I -4 
   83      P 4/m 
   84      P 42/m 
   85      P 4/n : 2 
   85      P 4/n : 1 
   86      P 42/n : 2 
   86      P 42/n : 1 
   87      I 4/m 
   88      I 41/a : 2 
   88      I 41/a : 1 
   89      P 4 2 2 
   90      P 4 21 2 
   91      P 41 2 2 
   92      P 41 21 2 
   93      P 42 2 2 
   94      P 42 21 2 
   95      P 43 2 2 
   96      P 43 21 2 
   97      I 4 2 2 
   98      I 41 2 2 
   99      P 4 m m 
  100      P 4 b m 
  101      P 42 c m 
  102      P 42 n m 
  103      P 4 c c 
  104      P 4 n c 
  105      P 42 m c 
  106      P 42 b c 
  107      I 4 m m 
  108      I 4 c m 
  109      I 41 m d 
  110      I 41 c d 
  111      P -4 2 m 
  112      P -4 2 c 
  113      P -4 21 m 
  114      P -4 21 c 
  115      P -4 m 2 
  116      P -4 c 2 
  117      P -4 b 2 
  118      P -4 n 2 
  119      I -4 m 2 
  120      I -4 c 2 
  121      I -4 2 m 
  122      I -4 2 d 
  123      P 4/m m m 
  124      P 4/m c c 
  125      P 4/n b m : 2 
  125      P 4/n b m : 1 
  126      P 4/n n c : 2 
  126      P 4/n n c : 1 
  127      P 4/m b m 
  128      P 4/m n c 
  129      P 4/n m m : 2 
  129      P 4/n m m : 1 
  130      P 4/n c c : 2 
  130      P 4/n c c : 1 
  131      P 42/m m c 
  132      P 42/m c m 
  133      P 42/n b c : 2 
  133      P 42/n b c : 1 
  134      P 42/n n m : 2 
  134      P 42/n n m : 1 
  135      P 42/m b c 
  136      P 42/m n m 
  137      P 42/n m c : 2 
  137      P 42/n m c : 1 
  138      P 42/n c m : 2 
  138      P 42/n c m : 1 
  139      I 4/m m m 
  140      I 4/m c m 
  141      I 41/a m d : 2 
  141      I 41/a m d : 1 
  142      I 41/a c d : 2 
  142      I 41/a c d : 1 
  143      P 3 
  144      P 31 
  145      P 32 
  146      R 3 : H 
  146      R 3 : R 
  147      P -3 
  148      R -3 : H 
  148      R -3 : R 
  149      P 3 1 2 
  150      P 3 2 1 
  151      P 31 1 2 
  152      P 31 2 1 
  153      P 32 1 2 
  154      P 32 2 1 
  155      R 3 2 : H 
  155      R 3 2 : R 
  156      P 3 m 1 
  157      P 3 1 m 
  158      P 3 c 1 
  159      P 3 1 c 
  160      R 3 m : H 
  160      R 3 m : R 
  161      R 3 c : H 
  161      R 3 c : R 
  162      P -3 1 m 
  163      P -3 1 c 
  164      P -3 m 1 
  165      P -3 c 1 
  166      R -3 m : H 
  166      R -3 m : R 
  167      R -3 c : H 
  167      R -3 c : R 
  168      P 6 
  169      P 61 
  170      P 65 
  171      P 62 
  172      P 64 
  173      P 63 
  174      P -6 
  175      P 6/m 
  176      P 63/m 
  177      P 6 2 2 
  178      P 61 2 2 
  179      P 65 2 2 
  180      P 62 2 2 
  181      P 64 2 2 
  182      P 63 2 2 
  183      P 6 m m 
  184      P 6 c c 
  185      P 63 c m 
  186      P 63 m c 
  187      P -6 m 2 
  188      P -6 c 2 
  189      P -6 2 m 
  190      P -6 2 c 
  191      P 6/m m m 
  192      P 6/m c c 
  193      P 63/m c m 
  194      P 63/m m c 
  195      P 2 3 
  196      F 2 3 
  197      I 2 3 
  198      P 21 3 
  199      I 21 3 
  200      P m -3 
  201      P n -3 : 2 
  201      P n -3 : 1 
  202      F m -3 
  203      F d -3 : 2 
  203      F d -3 : 1 
  204      I m -3 
  205      P a -3 
  206      I a -3 
  207      P 4 3 2 
  208      P 42 3 2 
  209      F 4 3 2 
  210      F 41 3 2 
  211      I 4 3 2 
  212      P 43 3 2 
  213      P 41 3 2 
  214      I 41 3 2 
  215      P -4 3 m 
  216      F -4 3 m 
  217      I -4 3 m 
  218      P -4 3 n 
  219      F -4 3 c 
  220      I -4 3 d 
  221      P m -3 m 
  222      P n -3 n : 2 
  222      P n -3 n : 1 
  223      P m -3 n 
  224      P n -3 m : 2 
  224      P n -3 m : 1 
  225      F m -3 m 
  226      F m -3 c 
  227      F d -3 m : 2 
  227      F d -3 m : 1 
  228      F d -3 c : 2 
  228      F d -3 c : 1 
  229      I m -3 m 
  230      I a -3 d 
 ========  =================

.. _groups-symmetry-sec:

Symmetry Elements
=================

A list of all of the symmetry elements of any space group can be output to file by placing a "OUTPUT_GROUP" command in the parameter file at any point after the "BASIS" section.

Every space group symmetry can be expressed mathematically as an operation

.. math::

   \textbf{r} \rightarrow \textbf{A}\textbf{r} 
                    + \textbf{t}

Here, :math:`\textbf{r} = [r_{1}, \ldots, r_{D}]^{T}` is a dimensionless 
D-element column vector containing the components of a position within 
the unit cell in a basis of Bravais lattice vectors, :math:`\textbf{A}` 
is a :math:`D \times D` matrix that represents a point group symmetry 
operation (e.g., identity, inversion, rotation about an axis, or 
reflection through a plane), and :math:`\textbf{t}` is a dimenionless
D-element colummn vector that (if not zero) represents a translation 
by a fraction of a unit cell. Every group contains an identity element in 
which :math:`\textbf{A}` is the identity matrix and :math:`\textbf{t}=0`. 

The elements of the column vectors :math:`\textbf{r}` and :math:`\textbf{t}` 
in the above are dimensionless components defined using a basis of Bravais 
basis vectors. The position :math:`\textbf{r} = [1/2, 1/2, 1/2]^{T}` thus
always represents the center of a 3D unit cell. The Cartesian representation 
of a position vector is instead given by a sum

.. math::

   \sum_{i=1}^{D} r_{i}\textbf{a}_{i}


in which :math:`\textbf{a}_{i}` is the Cartesian representation of 
Bravais lattice vector number i. The elements of the dimensionless 
translation vector :math:`\textbf{t}` are always either zero or 
simple fractions such as 1/2, 1/4, or 1/3. For example, a symmetry 
element in a 3D BCC lattice in which :math:`\textbf{A}` is the identity 
matrix and :math:`\textbf{t} = [1/2, 1/2, 1/2]^{T}` represents the 
purely translational symmetry that relates the two equivalent positions 
per cubic unit cell in a BCC lattice. Similarly, a glide plane in 
a 3D crystal is represented by a diagonal :math:`\textbf{A}` matrix 
with values of :math:`\pm 1` on the diagonal that represents 
inversion through a plane and a translation vector that represents 
a translation by half a unit cell within that plane.

The OUTPUT_GROUP command outputs a list of symmetry elements in 
which each element is displayed by showing the elements of the 
matrix :math:`\textbf{A}` followed by elements of the associated 
column vector :math:`\textbf{t}`.

The Bravais lattice vectors used internally by PSCF for cubic, tetragonal, 
and orthorhombic 3D systems are orthogonal basis vectors for the simple 
cubic, tetragonal, or orthorhombic unit cells, which are aligned along 
the x, y, and z axes of a Cartesian coordinate system. Similarly, the 
basis vectors used for the 2D square and rectangular space groups are 
orthogonal vectors which form a basis for a cubic or rectangular
unit cell. The grid used to solve the modified diffusion equation is
based on the same choice of unit cell and, thus for example, uses a
regular grid within a cubic unit cell to represent fields in a BCC or 
FCC lattice.  For body-centered and space-centered lattice systems, 
it is worth nothing that this unit cell not a primitive (minimum 
size) unit cell of the crystal: For example, a cubic unit cell actually 
contains 2 equivalent primitive unit cells of a BCC lattice or 4 
primitive cells of an FCC lattice. 
 
One consequence of the fact that PSCF does not always use a primitive 
unit cell is that, in the Fourier expansion of the omega and rho fields,
the Fourier coefficients associated with some sets of symmetry-related 
wavevectors (some "stars") are required to vanish in order to satisfy 
the requirement that the field be invariant under all elements of the 
specified space group. The rules regarding which stars must have 
vanishing Fourier coefficients are the same as the rules for systematic 
cancellations of Bragg reflections in X-ray or neutron scattering from 
a crystal of the specified space group. The procedure used by PSCF to 
construct symmetry adapted basis functions automatially identifies and 
accounts for these systematic cancellations.

