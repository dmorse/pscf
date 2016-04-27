
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

   * Make a copy named config.mk of the file config.mk_r, by entering::

        cd config.mk_r config.mk

   * Examine and edit the config.mk file to reflect your environment, 
     and to choose an installation directory (see below).

   * To compile, enter::

        > make all

     from within src/make.

   * To install in the location specified by $(INSTALL) makefile variable, 
     enter::

        > make install

    * To modify your $PATH and $PYTHONPATH environment variables to include
      the directories in which you have installed executables and python
      scripts, enter::

           > source $(INSTALL)/bin/pscf-env

      where $(INSTALL) denotes the base installation directory. This method 
      will only modify these environment variables temporarily until you 
      log out (on linux) or close the terminal (on a Mac). 

    * To setup your environment to add the appropriate paths to your 
      $PATH and $PYTHONPATH whenever you log in, add the above command,
      "source $(INSTALL)/bin/pscf-env", to the .profile configuration file
      in your home directory.
      
Some of these steps are discussed in more detail below

Customize the Makefile
-----------------------

In the config.mk file in the src/make directory (which you must create 
by copying config.mk_r), you will need to set values for a set of macro 
variables to values appropriate to your system. Makefile variables 
you may need to reset are:
 
 =========  ========================================================
 F90        path to executable for Fortran 90 compiler
 FAST       compiler options for high optimization
 NOPT       compiler options for no optimization
 LIBDIR     option setting any nonstandard directories for libraries
 LAPACKLIB  option setting name of lapack library (e.g., "-l liblapack")
 FFTWLIB    option setting name of the FFTW library (e.g., "-l fftw3")
 INSTALL    root directory for installation of executables and scripts
 =========  ========================================================

The default config.mk file contains values appropriate for a number of 
different common environments, most of which are commented out. Modify 
this configuration file as needed, but make sure you give exactly one 
uncommented definition for each variable, and comment out any unused 
definitions.

Compile and Link
-----------------

To compile and link, from the src/make directory, simply enter::

   > make -j4 

This should fill the src/make directory with .o and .mod files, and 
create an executable named pscf in the same directory. 

Installation
------------

To install the executable in your chosen location, simply enter::

   > make install

This will install:

   * pscf and other executable files in the $(INSTALL)/bin directory

   * python modules in $(INSTALL)/lib/python2.7/site-packages/pscf/

   * matlab scripts in $(INSTALL)/lib/matlab/

where $(INSTALL) denotes the makefile variable defined in the config.mk
file, which specifies the root of the installation directory tree.


Setting your Path and Running 
------------------------------

To invoke the pscf program, you will need to do one of the following:

   * Invoke pscf using an absolute path name

   * Install pscf in a directory such as /usr/local/bin that is already 
     in your command search $PATH. 

   * Add the directory containing your executable to your command search
     PATH variable. To do so manually, enter:

         PATH=$PATH:$(INSTALL)/bin
         export path

     where $(INSTALL)/bin denotes the directory in which you installed 
     the pscf executable. 

   * As discussed above, you can set both the $PATH and $PYTHONPATH 
     variables by using the "source" command to execute the pscf-env 
     script.

Cleaning Up
-----------
	
To remove all generated files from the pscf/make directory::

   > make clean

To remove all files installed in the INSTALL directory and start over, enter::

   > make uninstall

If you have manually moved the executable or other files, you may need to 
remove these files manually to fully clean up.

