
.. _install-compile-make-sec:

Compiling from source, using make
=================================

It is also possible to compile using the unix make utility using a simple
Makefile that is located in the pscf/make directory. The instructions for
using make to compile from source are same on any unix-like operating system,
including Max OS X. The main difference among different unix environments is 
the locations of the required libraries. 

To compile the code in this way, proceed as follows:

   * Change the working directory (cd) to the pscf/make directory

   * Make a copy named makefile of the file Makefile_r, by entering::

        cd Makefile_r Makefile

   * Examine and edit the Makefile to reflect your environment (see below)

   * To compile, enter::

        > make pscf

     from within src/make.

Some of these steps are discussed in more detail below

Customize the Makefile
-----------------------

In the Makefile in the src/make directory (which you must create by
copying Makefile_r), you will need to set values for a set of macro 
variables to values appropriate to your system. Makefile variables 
you may need to reset are:
 
 =========  ========================================================
 BIN        directory in which executable should be installed
 EXE        name of executable file
 F90        path to executable for Fortran 90 compiler
 FAST       compiler options for high optimization
 NOPT       compiler options for no optimization
 LIBDIR     option setting any nonstandard directories for libraries
 LAPACKLIB  option setting name of lapack library (e.g., "-l liblapack")
 FFTWLIB    option setting name of the FFTW library (e.g., "-l fftw3")
 =========  ========================================================

The makefile contains values appropriate for a number of different common 
environments, most of which are commented out. Modify the Makefile file
as needed, but make sure you give exactly one uncommented definition for 
each variable, and comment out any unused definitions.

Compile and Link
-----------------

To compile and link, from the src/make directory, simply enter::

   > make pscf

This should fill the src/make directory with .o and .mod files, and 
create an executable $(BIN)/$(EXE). By default, this will create a program 
named pscf in the pscf/bin directory. The executable file can be relocated 
to somewhere else if you desire.

Setting your Path and Running 
------------------------------

To invoke the program, you will either need to:

   * Invoke the program using an absolute path name

   * Add the directory containing your executable to your command search
     PATH variable. To do so, enter:

         PATH=$PATH:$(BIN)
         export path

     where $(BIN) denotes the path in which you installed the executable.
     You may want to add this to your .bashrc or .profile file so that 
     this directory is added to your path when automatically when you 
     log in.

   * Move pscf to a directory such as /usr/local/bin that is already in 
     your $PATH. 

Cleaning Up
-----------
	
To remove all of the .o amd .mod files from the src/make directory::

   > make clean

If you have moved the executable, you will need to remove this manually
to fully clean up.
