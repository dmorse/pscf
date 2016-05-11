
.. _install-compile-make-sec:

Compiling from source, using make
=================================

It is also possible to compile using the unix make utility using a simple
Makefile that is provided in the make/ directory of the git repository. The 
instructions for using make to compile from source are same on any unix-like 
operating system, including Max OS X. The main difference among different 
unix environments is the locations of the required libraries. 

To compile the code in this way, proceed as follows:

   * Follow the directions given in the subsection of the previous
     page entitled :ref:`install-compile-cmake-dependencies-sub` in order
     to install dependencies on your operating systems. Install all
     required dependencies except cmake.

   * Follow the directions given in the subsection of the previous
     page entitled :ref:`install-compile-cmake-getsource-sub` in order
     to obtain the source code and create a directory structure.

   * Change the working directory (cd) to the pscf/git/make directory.
     Note that this is an existing subdirectory of the pscf/git directory, 
     and is different from the initially empty pscf/cmake directory from
     which we suggested that you should invokve cmake, if compiling with
     cmake. 

   * This directory should contain files named config.mk_r and Makefile.
     Make a copy of the file config.mk_r, by entering::

        cd config.mk_r config.mk

   * Examine and edit the new config.mk file to reflect your environment, 
     and to specify an installation directory. The need to manually edit 
     this configuration file is the main difference between using cmake 
     to and using make alone. See below for further instructions about
     this step.

   * To compile, enter::

        > make -j4 all

     from within pscf/git/make.

   * To install in the directory specified by the $(INSTALL) makefile 
     variable (as defined in config.mk), enter::

        > make install

   * Follow the directions given in the subsection of the previous
     page entitled :ref:`install-compile-cmake-paths-sub` in order
     to modify environment variables that define search paths.

Some of these steps are discussed in more detail below.

**Editing the config.mk configfuration file**

In the config.mk file in the src/make directory (which you must create 
by copying config.mk_r), you will need to set values for a set of macro 
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

The file config.mk_r is a default makefile that is stored in the 
repository. We require you to make an operational copy of this, named 
config.mk, so that you can modify your config.mk file as needed without 
changing the default copy. The config.mk_r file contains values for
Makefile variables appropriate for a several different common operating
system environments, most of which are commented out. It also contains 
some comments about appropriate choices for specific environments. 
Modify your copy of config.mk as needed, but avoid modifying config.mk_r,
and make sure you give exactly one uncommented definition for each 
variable in the config.mk file, and comment out any unused definitions.

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

