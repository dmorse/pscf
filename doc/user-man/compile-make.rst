
.. _install-compile-make-sec:

Compiling from source, using make
=================================

It is also possible to compile using the unix make utility using a simple
Makefile that is provided in the make/ directory of the git repository. The 
instructions for using make to compile from source are same on any unix-like 
operating system, including Max OS X. The main difference among different 
unix environments is the locations of the required libraries. 

To compile the code in this way, proceed as follows:

   * Change the working directory (cd) to the pscf/git/make directory
     (i.e., to make/ subdirectory of the directory tree containing the
     pscf git repository).

   * Make a copy named config.mk of the file config.mk_r, by entering::

        cd config.mk_r config.mk

   * Examine and edit the config.mk file to reflect your environment, 
     and to choose an installation directory (see below).

   * To compile, enter::

        > make -j4 all

     from within pscf/git/make.

   * To install in the location specified by the $(INSTALL) makefile 
     variable (defined in config.mk), enter::

        > make install

   * Modify the PATH and PYTHONPATH environment variables, if needed, by 
     following the instructions given for compiling using cmake, in the 
     subsection :ref:`install-compile-cmake-paths-sub`.

Some of these steps are discussed in more detail below.

**Editing the config.mk configfuration file**

In the config.mk file in the src/make directory (which you must create by 
copying config.mk_r), you will need to set values for a set of macro 
variables to values appropriate to your system. Makefile variables you 
may need to reset are:
 
 =========  ========================================================
 F90        path to executable for Fortran 90 compiler
 FAST       compiler options for high optimization
 NOPT       compiler options for no optimization
 LIBDIR     option setting any nonstandard directories for libraries
 LAPACKLIB  option setting name of lapack library (e.g., "-l liblapack")
 FFTWLIB    option setting name of the FFTW library (e.g., "-l fftw3")
 INSTALL    root installation directory 
 =========  ========================================================

The INSTALL makefile variable in this makefile is equivalent to the 
MAKE_INSTALL_PREFIX variable that can be passed to cmake when compiling
using cmake.

The default config.mk file contains values appropriate for a number of 
different common environments, most of which are commented out. It also
contains some comments about appropriate choices for common environments. 
Modify your copy of config.mk as needed, but make sure you give exactly 
one uncommented definition for each variable, and comment out any unused 
definitions.

**Compile and Link**

To compile and link, from the git/make directory, simply enter::

   > make -j4 

This should create several .o object files and .mod module files in
the git/make directory, and then create an executable file named pscf
in the same directory. 

**Install**

To install a copy of the executable file in your chosen location, after
compiling, simply enter::

   > make install

This will install:

   * pscf and other executable files in directory $(INSTALL)/bin/

   * python modules in $(INSTALL)/lib/python2.7/site-packages/pscf/

where $(INSTALL) denotes the value of the makefile variable defined in 
the config.mk file.

**Cleaning Up**
	
To remove all generated files from the git/make directory, if desired, 
enter::

   > make clean

To remove all files installed in the INSTALL directory by the 
"make install" command, enter::

   > make uninstall

