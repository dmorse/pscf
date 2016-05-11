
.. _python-page:

*************
Python Tools
*************

The PSCF is distributed with a set of Python modules that are can 
be used to simplify many common data analysis and job preparation 
tasks.  All of these modules are part of a Python package named 
pscf. These Python modules are located within the source code 
repository in a subdirectory named tools/python/pscf subdirectory 
of the directory that contains the git repository (e.g., of pscf/git). 

Overview
=========

The most important of the Python modules distributed with PSCF are 
a set of modules that can read and write several of the file formats 
used by PSCF. Each of these modules is a file that contains the
definition of a single class, in which the name of the class is a
a capitalized version of the name of the file or module. The names
of these main modules and associated classes are as follows:

  ========== ========= =====================================================
  Module     Class     Description
  ========== ========= =====================================================
  paramfile  ParamFile Reads and stores contents of a parameter file
  fieldfile  FieldFile Reads and stores a symmetry-adapted field file
  outfile    OutFile   Reads and stores a output summary file
  sweep      Sweep     Reads and stores summary files produced by a sweep
  ========== ========= =====================================================

The constructor for each of these classes takes the name of a file or 
(for the Sweep class) a prefix common to a set of files as an argument, 
and immediately reads the specified file or set of files. The constructor
for each class thereby constructs a Python object that contains the entire
contents of the specified file or set of files. The values of each 
named variable or (for field files) array of values that are read 
from file is stored in a corresponding attribute (member data variables) 
of the associated Python object. If the name of a parameter or variable 
is given by a label in the file format, the name of the object 
attribute that holds the corresponding value is always the same as 
the label used in the file format.

The ParamFile, OutFile and FieldFile classes each also provide method a 
method named "write" that can write the contents of the object to an 
output file in a format that can be read by PSCF. This provides a 
basis for programmatically modifying files by reading in one file,
using Python to change values of one or more parameters or variables,
and then writing the modified object to a (generally) different file.
The FieldFile class also provides some specialized methods to carry out 
common manipulations on entire lists of Fourier coefficients.

The Sweep class reads in all of the numbered output summary files 
that are produced by a PSCF "SWEEP" command and stores the contents 
of these files as a Python list of OutFile objects. This provides a 
convenient basis for analyzing the results of a sweep, by looping 
over output files produced by different simulations and outputting 
values of selected variables (e.g., free energies) in a compact 
form suitable for plotting and/or further analysis.

To use any of these modules, one must first open a python interpreter, 
import the relevant class, and then construct an object of the desired 
type by passing the constructor the name of the file or directory of 
interest. This pattern is demonstrated in each of the examples given 
below.

Python Environment
==================

To use the pscf python package, Python obviously must be installed
on your computer. Because Python is also required to install PSCF 
from source, instructions for how to check whether Python is
installed (it usually is) or install it if necessary are given in
the documentation regarding :ref:`install-compile-cmake-sec`.

After PSCF is installed, whether it is installed using a binary 
installer or after compiling from source, copies of all associated 
python modules are placed in a directory::

   $(INSTALL_DIR)/lib/python2.7/site-packages/pscf

where $(INSTALL_DIR) is a place holder for the actual absolute path
to the root of the installation directory tree. This pscf/ subdirectory
must be included in the PYTHONPATH environment variable in order for 
the Python interpreter to find and use these modules. This can be 
accomplished by running script $(INSTALL)/bin/pscf-env, as discussed 
in more detail in the discussion of :ref:`install-compile-cmake-paths-sub`. 

After running the pscf-env script, to confirm that the $PYTHONPATH 
is set up correctly, type::

   > echo $PYTHONPATH

This will cause a string to be printed to the terminal that should
include the absolute path of the relevant directory, either by itself 
or as part of a colon separated list of directories.

ParamFile
==========

In the following example, we use class ParamFile to construct
an object that contains the contents of a single parameter 
input file named 'param' in the current directory, and print 
values of the parameter N_monomer (the number of monomer 
types) and the volume fraction of the first chain species.
To do this, one could type the following::

    > python
    python> from pscf.outfile import ParamFile
    python> param = ParamFile('param')
    python> print out.N_monomer
    python> print out.phi_chain[0]
    python> quit()

It is safe to try running this example on any of the examples 
sets of input files that are provided with PSCF - this example 
will not modify any input files. The import statement will 
produce an error if the python interpreter cannot find module 
pscf.outfile, which generally indicates that the PYTHONPATH is 
not set correctly.

The first line in the above example opens an interactive python 
interpreter, which is closed by the last line. The second line 
imports the class OutFile from module (or file) outfile of 
package pscf, making it available for use. The third line 
creates an object named out of type OutFile that contains the 
contents of the file whose name ('out') is passed to the OutFile 
constructor function. In python, like in C++ or Java, a 
constructor is invoked by using the name of the class as a 
function, which in this case takes a file name as a single 
argument. 

After the constructor is called, object out contains the 
values of all of the variables stored in the file 'out'.
This is demonstrated by the 4th and 5th line of the above
example, which print the values of two variables to the
terminal screen.  The value of any labelled parameter in 
the original file is stored in an attribute (or member 
variable) whose name is the same as the name of the 
label associated with the parameter in the original file. 

Parameters such as phi_chain whose values are stored in
PSCF in one-dimensional arrays in which different elements 
refer to e.g., different molecular species, are stored in 
the associated python object as python lists. Individual 
elements of a python list can be accessed using a subscript 
notation similar to that used for arrays in C, using indices
that are numbered consecutively from 0. This is demonstated 
by the 5th line of the above example, in which we use the 
symbol 'out.phi_chain[0]' to access the value of the volume 
fraction of the first (index 0) polymer species. 

Note: All pscf python modules use the C/Python convention 
in which C array and python list indices are numbered 
consecutively from zero, while PSCF uses the Fortran 
convention in which indices start from 1.  One consequence 
of this is that, for example, data associated with the 
second of two or more monomer types is associated with a 
list index of 1 in all python objects, but is labelled by 
an integer "2" in block_monomer array in the PSCF input 
parameter file, which uses the Fortran convention to 
assign monomer type label values. Users need to be aware 
of this difference and correct for it as necessary when 
using the python modules.

OutFile 
=======

Output summary files can be parsed, modified and output using a 
syntax essentially identical to that used for parameter files. 
In the following simple example, we read an output summary file 
in the working directory named 'out', and then print out the 
values of f_Helmholtz, the free energy per monomer, and 
mu_chain[0], the chemical potential of the first chain species::

    > python
    python> from pscf.outfile import ParamFile
    python> out = ParamFile('out')
    python> print out.f_Helmholtz
    python> print out.mu_chain[0]
    python> quit()

Because the first part of an output summary file has the same 
syntax as an input parameter file, an output summary file from
one simulation can be used as a starting point for creating a 
parameter file for a related system. This can be done either 
by manually editing and copying the output file, or by using
python to read the file, modify the values of a few parameters
and 


FieldFile
==========

A FieldFile object holds all of the information stored in the
symmetry-adapated field file format, including the values of 
the coefficients of all basis functions for the fields 
associated with each monomer type. A Field object is 
constructed using a syntax similiar to that for a ParamFile 
or OutFile object, by passing the name of an associated file
to the constructor.

**Attributes**

A symmetry-adapted field file contains a header with labelled
parameters and data section containing columns of number that
give coefficients of different basis functions in a symmetry
adapted Fourier expansion. The value of each of the parameters
that appears in the header is stored in an attribute with a
name given by the parameter label that appears in the file. 
The contents of the data section are stored in three attributes
named "fields", "waves" and "counts", as discussed below.

The attribute "fields" is a list of lists of Fourier 
coefficients.  Each element of list fields is a list that 
contains of the Fourier coefficients for one monomer type. 
Thus, for example, fields[1] is a list that contains the 
list of coefficients given in the second column (with indices 
numbered from 0) of the data section of the associated file. 
The item fields[1][13] is a real number that is equal to the
coefficient of basis function 13 (the 14th basis function,
with indices numbered from 0) of the field associated with
monomer type number 1 (i.e., the 2nd monomer type).

The attribute "waves" is a list in which element contains a
list of 1, 2, or 3 integers indices for a wavevector that
is contained in the associated basis function, which acts
as an identifier for the basis function.  Thus, for example, 
for the gyroid phase, the second basis function, with index 
1, is associated with the {211} family of plane waves. In 
this case, the value of waves[1] is a list of integers, 
waves[1] == [2, 1, 1], that identifies this basis function.

The attribute "counts" is containing the integers given in the
last column of the data section of a field file. Each of these
integers gives the number of wavevectors in a "star" of symmetry
related wavevectors that is associated with the corresponding 
basis function. Thus for example, in a file for a gyroid phase,
with space group "I a -3 d", for which waves[1] = [2, 1, 1],
count[1] == 24, because there are 24 wavevectors in the {211}
family of wavevectors of a cubic crystal. 

The syntax for creating a FieldFile object, which is illustrated
in the following example,  is essentially the same as that for 
creating a ParamFile or OutFile object::


    > python
    python> from pscf.fieldfile import FieldFile
    python> omega = FieldFile('omega')
    python> print omega.waves[1]
    python> print omega.fields[0][1]

In this example, we open and read a chemical potential field
file named 'omega' in the current directory, print the short
list of Miller indices that identifies basis function number 1 
(the second basis function), and print the value of the coefficient
of this basis function in the expansion of the chemical potential
field for monomer type number 0. 

Sweep
======

The Sweep class is a container that holds all of the data given
in the set of number output summary files produced by a PSCF 
SWEEP command.  

Recall that the SWEEP command performs a sequence of SCFT 
calculations along a line in parameter space, and produces a set 
of output files for each calculation, with file names that 
contain an integer index. The resulting output summary files 
have names of the form <output_prefix>i.out, where <output_prefix> 
denotes the output_prefix string parameter given in the input 
parameter file, and where i is an integer index. The index i 
has values in the range [0, s_max], where s_max is the maximum 
value given in the parameter file.  Typically, output_prefix 
is taken to be the name of a directory including a trailing 
backslash (/) directory separator, such as "out/".  In this 
case, the SWEEP produces a series of output summary files in 
the specified directory with names 0.out, 1.out, 2.out, etc.

The constructor for a Sweep object assumes that the SWEEP
command was run using a directory name with a trailing slash
as an output_prefix, and that the output directory thus 
contains a sequence of files with names 0.out, 1.out etc.
The constructor takes the name of the directory (with no
trailing slash) as an argument, and reads any sequence of
such output summary files that it finds, from which it 
creates a python list of OutFile objects, each of which
contains the contents of one output summary file. Each of 
the resulting Outfile objects can be accessed by applying 
the subscript [] operator directly to the Sweep object, 
thus emulating the syntax of a Python list. Thus if x is 
a Sweep object, x[8] is an OutFile object containing the 
contents of file named 8.out in the directory that was 
named in the Sweep constructor. The number of OutFile
objects in Sweep object named x is returned by the operator
len(x), as for a list.

The following example illustrates the syntax for creating
a Sweep object and accessing a particular variable in a
particular simulation::

    > python
    python> from pscf.sweep import Sweep
    python> x = Sweep('.')
    python> print len(x)
    python> print x[8].f_Helmholtz

In this example, we assume that the python interpreter 
was run from the directory containing a set of output 
summary files named 0.out, 1.out etc. The third line of 
this example thus reads all of the output files in the 
working directory, indicated here by the unix shorthand 
'.' that is passed to the constructor.  The fourth line 
prints the number of output files found to the screen.
The fifth line prints the value of the variable 
f_Helmholtz read from the file 8.out.

Users can aalso iterate over the list of OutFile 
objects contained in a Sweep object in order to output
or manipulate lists showing how selected variables change
within the sequence of calculations. This is shown in the 
following example::

    > python
    python> from pscf.sweep import Sweep
    python> x = Sweep('.')
    python> print len(x)
    python> file = open('free_energy','w')
    python> for outfile in x:
        ***     line = str(outfile.block_length[0][1]) + '  '
        ***     line += str(outfile.f_Helmholtz) 
        ***     print line
        ***     file.write(line + "\n")
    python>
    python> file.close()

The fifth line of this example uses the Python open()
function to open a new file named 'free_energy' for writing 
(mode = 'w'). The for loop produces a sequence of text 
lines containing two columns of numbers, in which the first
column contains values of the length block_length[0][1] 
of the second (index 1) block of the first (index 0) 
chain species, while the second column contains the 
value of f_Helmholtz, which is the Helmholtz free energy 
per monomer normalized by kT, and each line contains values
from a different output summary file. In this example, each 
line of this output is both printed to the screen and written 
to a file named free_energy. The penultimate line closes the
file before closing the python interpreter.

The type of operation given above, which produces 
a string of containing two columns of numbers, is
commonly needed to summarize information about a
sweep.  The Sweep class thus provides a method named 
write() that is designed to simplify this operation.
The write function takes two arguments, named expr1 
and expr2, each of which is literal string containing 
a mathematical expression written using the names of 
attributes as variable names. It returns a string 
containing two columns of numbers, in which each 
value in the first column is obtaining by evaluating
expression expr1 and each value is obtained by 
evaluating expr2, and in which each row represents
a pair of values obtained from a different simulation.
The above example could also be expressed using the 
write method as::

    > python
    python> from pscf.sweep import Sweep
    python> sweep = Sweep('.')
    python> text = sweep.write('block_length[0][1]','f_Helmholtz')
    python> print text
    python> file = open('free_energy','w')
    python> file.write(text)
    python> file.close()
    python> quit()

Each of the arguments of the write functions are text
strings, and need not simply be attribute names or items 
within attributes that are lists, as in the above example.
These arguments can be any text strings that can contain 
a valid python expression in which attribute names are 
used as variable names. The statement::

    text = sweep.write('2.0*block_length[0][1]','f_Helmholtz')

would thus also be valid, and would produce a string 
in which the first column contained values equal to 
twice the length of the 2nd block of the first chain 
species. 
