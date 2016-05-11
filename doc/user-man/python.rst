
.. _python-page:

**************
Python Modules
**************

The PSCF is distributed with a set of python modules that are can 
be used to simplify many common data analysis and job preparation 
tasks.  All of these modules are part of a python package named 
pscf. These python modules are located within the source code 
repository in a subdirectory named tools/python/pscf subdirectory 
of the directory that contains the git repository (e.g., of pscf/git). 

When PSCF is installed, using either a binary installer or after
compiling from source, copies of all associated python modules are 
placed in a subdirectory named lib/python2.7/site-packages/pscf of 
the installation directory. This directory must be added to the 
PYTHONPATH environment variable in order for the python interpreter 
to find and use these modules. This can be accomplished by using 
the 'source' command to run the pscf-env script that is installed
in the bin/ subdirectory of the installation directory, as described 
in the discussion of :ref:`compile-cmake-paths-sub`. 

To check if the $PYTHONPATH is set up correctly, type::

   > echo $PYTHONPATH

This should cause a string to be printed to the terminal that includes 
the absolute path of the relevant directory, either by itself or as 
part of a colon separated list of directories.

Overview
=========

The most important of the python modules distributed with PSCF are 
several modules that can read and write several of the file formats 
used by PSCF. Each of these modules is a file that contains the
definition of a single class, in which the name of the class is a
a capitalized version of the name of the file or module. The names
of these main modules and associated classes are as follows:

  ========== ========= =====================================================
  Module     Class     Description
  ========== ========= =====================================================
  paramfile  ParamFile Reads and store contents of a parameter file
  fieldfile  FieldFile Reads and contain a symmetry-adapted field file
  outfile    OutFile   Reads and contain a output summary file
  sweep      Sweep     Reads and stores summary files produced by a sweep
  ========== ========= =====================================================

The constructor for each of these classes takes the name of a file 
or (for the Sweep class) a directory as an argument, and immediately 
reads the specified file or set of files, thereby constructing an 
object that contains the contents of the specified file or set of 
files. The values of each variable or (for field files) array of
field coefficients that are read from file is stored in an attribute 
(member data variables) of the resulting object. If the name of a 
parameter or variable is given by a label in the file format, the 
name of the object attribute that holds a parameter value is always
the same as the label used in the original file.

To use any of these modules, one must first open a python interpreter, 
import the relevant class, and then construct an object of the desired 
type by passing the constructor the name of the file or directory of 
interest. This pattern is demonstrated in each of the examples given 
below.

The ParamFile, OutFile and Field classes each provide method a method
named "write" that can write the contents of the object to an output
file in a format that can be read by PSCF. One simple way that these
modules can be used is to modify files, by reading in a file, changing
the value of one variable and then writing the modified object to a
file with a different name.  The Field class also provides methods 
to carry out common manipulations on coefficients in the Fourier 
expansion of a set fields, such as adding a new monomer type, removing
a monomer type, switching fields associated with two monomer types, or
multiplying all the coefficients of a particular field by a constant.
The ability to read, modify and output parameter and field files 
provides a convenient framework for carrying out operations in which
the output of a simulation of one system is used as a starting point
for preparing input files for simulation of a related system. 

The Sweep class reads in all of the numbered output summary files 
that are produced by a PSCF "SWEEP" command. The SWEEP command 
performs a sequence of SCFT calculations along a line in parameter
space, and produces a sequence of output summary files whose names
are distinguished by an integer index. The names of the resulting
output files are of the form <output_prefix>i.out, in which 
<output_prefix> represents the output_prefix string parameter 
given in the input parameter file and i is an integer index,
numbered from zero. Typically, output_prefix is the name of a 
directory including a trailing backslash (/) directory separator,
such as "out/", in which case the SWEEP produces a series of
output summary files in the output directory with names 
0.out, 1.out, ...., n.out, where n is an integer that is 
specified by the value of the parameter s_max in the input 
parameter file. The constructor of a Sweep object reads all 
of these files and creates a python list of OutFile objects, 
in which the ith element of the list is an outfile object 
containing the contents of file i.out. This provides a very
convenient basis for framework for analyzing changes in free
energies or unit cell dimensions with changes in parameters, 
by allowing the user to loop over the output files produced by 
a sweep and output selected sets of parameter or output values 
for each simulation.

==================
ParamFile Examples
==================

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
convention in which indices start from 1.  One 
consequence of this is that, for example, data associated 
with the second of two or more monomer types is associated 
with a list index of 1 in all python objects, but is 
labelled by an integer "2" in block_monomer array in 
the PSCF input parameter file, which uses the Fortran 
convention to assign monomer type label values. Users 
need to be aware of this difference and correct for 
it as necessary when using the python modules.

================
OutFile Examples
================

Output summary files can be parsed, modified and output
using a syntax essentially identical to that used for
parameter files. In the following example, we read an
output summary file in the working directory named 'out',
and then print out the values of f_Helmholtz, the free
energy per monomer, and mu_chain[0], the chemical potential
of the first chain species::

    > python
    python> from pscf.outfile import ParamFile
    python> out = ParamFile('out')
    python> print out.f_Helmholtz
    python> print out.mu_chain[0]
    python> quit()

==============
Field Examples
==============

A Field object holds all of the information stored in the
symmetry-adapated field file format, including the values
of the coefficients of all basis functions for the fields 
associated with each monomer type. A Field object is 
constructed using a syntax similiar to that for a ParamFile
or OutFile object, by passing the name of an associated file
to the constructor.

A symmetry-adapted field file contains a header with labelled
parameters and data section containing columns of number that
give coefficients of different basis functions in a symmetry
adapted Fourier expansion. The value of each of the parameters
that appears in the header is stored in an attribute with a
name given by the parameter label that appears in the file. 
The contents of the data section are stored in three attributes
named "fields", "waves" and "counts".

The attribute fields is a list of lists of coefficients. 
Each element of lists fields is a list that contains a list 
of coefficient for one monomer type. Thus, for example, 
fields[1] is a list that contains the list of coefficients 
given in the second column (with indices numbered from 0)
of the data section of the associated file. The element
fields[1][13] contains the coefficient of basis function
number 13 (the 14th basis function) for the second monomer 
type. 

The attribute waves is a list in which element contains a
list of 1, 2, or 3 integers indices for a wavevector that
is contained in the associated basis function, which acts
as an identifier for the basis function.  Thus, for
example, for the gyroid phase, the second basis function,
with index 1, is associated with the {211} family of plane
waves. In this case, the value of waves[1] is a list of
integers, waves[1] == [2, 1, 1], that identifies this 
basis function.



