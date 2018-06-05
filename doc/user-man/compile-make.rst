
.. _install-compile-make-sec:

Compiling from source, using make
=================================

It is also possible to compile using the unix make utility alone, using 
a simple Makefile that is provided in the make/ subdirectory of the git 
repository. The instructions for using make to compile from source are 
the same on any unix-like operating system, including Max OS X. The main 
differences among different unix environments are the locations of the 
required libraries.

To compile the code in this way, proceed as follows:

   * Follow the directions given in the discussion of 
     :ref:`install-compile-cmake-dependencies-sub` on the previous 
     page. You will need to install all dependencies listed there
     except cmake.

   * Follow the directions given in the discussion of
     :ref:`install-compile-cmake-getsource-sub` on the previous page
     to create an appropriate directory structure and obtain the 
     source code. After this step, you should have a directory named
     pscf/ with a subdirectory named git/ that contains the contents 
     of the git repository. If you are not using cmake, then you do 
     not need to create a subdirectory of pscf/ named cmake/.

   * Change working directory (cd) to the directory pscf/git/make .
     Note that this is an existing subdirectory of the pscf/git 
     directory, and is different from the initially empty directory
     pscf/cmake from which we recommended cmake be invoked when 
     using cmake to compile.

   * The pscf/git/make directory will contain files named config.mk_r 
     and Makefile. Make a copy of the file config.mk_r, by entering::

        cd config.mk_r config.mk

   * Examine and edit the new config.mk file to reflect your environment, 
     and to specify an installation directory.  See below for further 
     instructions about this step.  The need to manually edit this 
     configuration file is the main difference between using cmake to 
     generate makefiles and using the simple makefile distributed with
     the source code.  

   * To compile, enter::

        > make -j4 all

     from within pscf/git/make. The "-j4" option is not necessary, and
     simply instructs the make utility to try to use up to 4 CPU cores 
     during compilation, if multiple cores are available.

   * To install in the directory specified by the $(INSTALL) makefile 
     variable (as defined in config.mk), enter::

        > make install

   * Follow the directions given in the subsection of the previous
     page entitled :ref:`install-compile-cmake-paths-sub` in order
     to modify environment variables that define search paths.

Several of these steps are discussed in more detail below.

**Editing the config.mk configfuration file**

In the config.mk file in the src/make directory (which you should have
created by making a copy of config.mk_r), you will need to set values for 
several macro variables to values appropriate to your system. Makefile 
variables you may need to reset are:
 
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
MAKE_INSTALL_PREFIX variable that can be passed to cmake when using
cmake to create a makefile.

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

