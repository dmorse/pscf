
.. _install-compile-make-sec:

Compiling from source, using make
=================================

It is also possible to compile using a Makefile in the src/build directory. 
This does an "in source" build, in which all of the files generated during 
compilation are placed in the pscf/src/build/ directory. The instructions 
for doing this are the same on any unix-like operating system. The main 
difference among different unix environments is the locations of the 
required libraries. 

To compile the code in this way, you should:

   * Change the working directory to pscf/src/build.

   * Make a copy named makefile of the file Makefile_r, by entering::

        cd Makefile_r Makefile

   * Examine and edit the Makefile to reflect your environment

   * To compile, enter::

        > make pscf

     from within src/build.

Some of these steps are discussed in more detail below

Customize the Makefile
-----------------------

In the Makefile in the src/build directory (which you must create by
copying Makefile_r), you will need to set values for a set of macro 
variables to values appropriate to your system. Makefile variables 
you may need to reset are:
 
 =========  ========================================================
 BIN        directory in which executable should be installed
 EXE        name of executable file
 F90        path to executable for Fortran 90 compiler
 FAST       compiler options for high optimization
 NOPT       compiler options for no optimization
 LAPACKLIB  name of lapack library, and directory if not standard
 FFTWLIB    name of the FFTW library, and directory if not standard
 =========  ========================================================

The makefile contains values appropriate for a number of different common 
environments, most of which are commented out. Make sure you give exactly
one uncommented definition for each variable, and that you comment out any 
definitions you are not using.

Compile and Link
-----------------

To compile and link, from the src/build directory, simply enter::

   > make pscf

This should fill the src/build directory with .o and .mod files, and 
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
	
To remove all of the .o amd .mod files from the src/build directory, as 
well as any editor buffer files with a ~ suffix from src tree, enter::

   > make clean

If you have moved the executable, you will need to remove this manually
to fully clean up.
