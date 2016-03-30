************
Installation
************

The src directory contains the source files for scf code for
periodic structures of incompressible blend of linear block
copolymers, homopolymers, and small molecule solvent. The code
is compiled from the src/build directory, using the Makefile
in the directory.

To compile the code, you will have to cd to src/build, edit the
Makefile, and issue the command 'make pscf' from within src/build.
These steps are described in more detail below

Customize the Makefile
----------------------

In Makefile in the src/build directory, you will need to set
values for a set of macro variables to values appropriate to your
system. The variables you may need to reset are:
 
==========  ===============================================
 SCF        root of scf directory tree.
 SRC        source file directory. Default: $(SCF)/src
 BIN        directory to which executable should be written
 EXE        name of executable file
 F90        path to executable for Fortran 90 compiler
 FAST       compiler options for high optimization
 NOPT       compiler options for no optimization
 LAPACKLIB  directory with Lapack libraries
 FFTWLIB    directory with FFTW library
==========  ===============================================

Compile and Link
----------------

To compile and link, from the src/build directory, issue the
command::

   > make pscf

This should fill the src/build directory with *.o and *.mod files,
and create an executable $(BIN)/$(EXE)

Cleaning Up
------------
	
To remove all of the *.o amd *.mod files from the src/build
directory, as well as all *~ files from src and its subdirectory,
if desired, enter the command::

   > make clean

