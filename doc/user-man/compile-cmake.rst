
.. _install-compile-cmake-sec:

Compiling from source, using cmake
==================================

Overview
--------

In this page, we discuss instructions for how to compile PSCF from source
using the cmake build system on either Mac OS X or linux operating systems.
It is also possible to compile the code with just the unix make utility, 
using a makefile that is provided with the source code. Compilation using 
make alone is described on a separate page on :ref:`install-compile-make-sec`.
The advantage of using cmake is that cmake can generally find the paths to
the Lapack and FFTW libraries upon which PSCF depends, whereas you need to
figure out the locations of these libraries yourself if you use make.

Compiling with cmake involves the following steps:

    * Install other software on which PSCF depends
    * Create an appropriate directory structure in which to build
    * Obtain the source code
    * Compile and install PSCF 

Only the first step, of installing external dependencies, is substantially
different for different operating systems. We thus give separate instructions 
for this step separately for different operating systems at the end of this 
page.

In order to obtain the source code, one can either:

    * Download a compressed archive (.tar) package
    * Use the git version control software to clone the repository

We recommend using git if possible, since this makes it simple to update 
the code later, but both procedures are described in what follows.

The following software packages must be installed before using cmake to
compile PSCF, if you plan on using git to obtain the source code:

   * git (in order to clone the source code)
   * cmake (to build a makefile)
   * a fortran 90 compiler 
   * LAPACK linear algebra library
   * FFTW version 3.x fast fourier transform library

You do not need to install git if you plan to simply download the 
source code rather than using git. In what follows, we will assume
that you plan to use free gfortran Fortran 90 compiler, which is 
part of the Gnu Compiler Collection (gcc) suite of compilers. 

Installing Dependencies
-----------------------

Mac OS X
~~~~~~~~~

**Installing XCode**

To create an environment in which you can compile from source on OSX, you 
will generally first need to install the apple XCode development environment.
XCode is available gratis from the app store, but is a large package that can
take a long time to install (do it with a good internet connection).  The 
XCode package contains git, so it is not necessary to install git separately
The OXS operating system also appears to come with a version of LAPACK, and 
the BLAS library upon which it depends.

**Package Managers: HomeBrew vs. MacPorts**

The remaining dependencies (cmake, gfortran and fftw) can be most easily 
installed using either the MacPorts or Homebrew package manager systems.  
These are both systems for managing open-source unix software on the unix 
subsystem of the Mac OSX.  The choice of package managers is up to you, 
but you should avoid using both on the same machine.  If either Homebrew 
or MacPorts is already installed and in use on your Mac, use the existing 
system and do not install the alternative, because they do not play well 
together.  If neither Homebrew or MacPorts is installed, we have slight 
preference for Homebrew, which made it slightly easier to install the
dependencies required by PSCF. We have succeeded in building PSCF using 
both package managers on different machines that were running the latest
version of Mac OS X (El Capitan, X 10.11) Instructions for both package
managers are given separately below.

**Installing dependencies via Homebrew**

To install from a command line terminal using homebrew::

   > brew install cmake
   > brew install gcc --with-fortran
   > brew install fftw --with-fortran

**Installing dependencies via Macports**

After MacPorts is installed, to install the required dependencies 
using the most recent version of the gnu compiler collection (gcc), 
which is gcc 5.X at the time of writing, enter::

   > sudo port install cmake
   > sudo port install gcc5
   > sudo port install fftw-3 +gfortran

Note that MacPorts (unlike homebrew) requires you to use "sudo"
to execute installation with superuser/administrator privileges, 
and thus will ask for a password after each of the above commands.

The gcc5 MacPorts package installs the gfortran Fortran 90 compiler 
executable at /opt/local/bin/gfortran-mp-5 . Versions compiled with 
earlier versions of gcc (e.g., 4.9) seem to be placed in the same 
directory with a different numerical suffix, e.g., gfortran-mp-49.  
CMake appears to be unable to find this compiler executable without 
help.  To remedy this, you should set the FC environment variable 
(which indicates the path to a Fortran compiler) to point to the 
absolute path to the gfortran executable before attempting to 
compile, by entering, for example::

   > FC=/opt/local/bin/gfortran-mp-5
   > export FC

If expect to compile this and other fortran programs repeatedly, 
you may want to put this in your .profile or .bashrc bash 
configuration file.

Ubuntu Linux
~~~~~~~~~~~~

Use the Ubuntu software manager or the command line apt-get utility to 
install the following packages:

   * git
   * cmake
   * gfortran
   * libfftw3-dev
   * liblapack3

To use the apt-get utility from the command line, enter:

   > sudo apt-get cmake
   > sudo apt-get gfortran
   > sudo apt-get libfftw3-dev
   > sudo apt-get liblapack3

Developers: To build .deb package for installation of binary executables 
on other Ubuntu and debian systems, as well as .tar and .zip source code 
archives, after installing on your machine, simply enter::

   > make package

To check the .deb file for semi-detailed information::

    # This extracts multiple files
    ar -vx pscf-1.0.0-Linux.deb
    # See the files that would be installed
    tar tvfz data.tar.gz 

Fedora Linux
~~~~~~~~~~~~

Instructions for Fedora are similar to those for Ubuntu, except that one 
should use the native yum command line package manager or the Fedora 
graphical software manager to install dependencies. The required Fedora 
packages are:

   * cmake
   * gcc-gfortran
   * lapack-devel
   * fftw-devel

To install these packages from the command line, enter::

   > sudo yum install cmake
   > sudo yum install gcc-gfortran
   > sudo yum install lapack-devel
   > sudo yum install fftw-devel

Instructions for obtaining source code, compiling and installing are the same 
as for Max OSX and Ubuntu.

Developers: On a Fedora machine, you can build a .rpm package and .tar 
and .zip archives by entering::

   > make package

from within the build directory.

To check the RPM for detailed information (Metadata, Dependencies, and 
File Contents), enter::

   > rpm --info -qpR -qlvp pscf-1.0.0-Linux.rpm 

Systems with Linux Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~

The following instructions describe how to build PSCF in a user directory 
at the Minnesota Computer Institute (MSI) Mesabi computer, using linux 
modules and the Intel compiler.  Similar instructions should apply to 
other large supercomputer clusters that use linux modules.

To load the required modules, enter::

   > module load cmake
   > module load intel mkl
   > module load fftw

The remaining instruction for how to obtain and compile the source code 
are generally similar to thos given for OSX or Linux. The only difference 
is that, to use the Intel compiler, one must tell cmake to use the Intel 
compiler by adding the option "-DUSE_INTEL=1" to the cmake command. The 
required command is thus::

   > cmake -DUSE_INTEL=1 -DCMAKE_INSTALL_PREFIX=/path/to/install ../pscf

Creating a Directory Tree
-------------------------
We assume in what follows that you will use cmake to implement on "out-of-source" build, in which all of the files generated during compilation are placed in a different directory than the source code. To do this, we suggest that you create a directory named pscf/ with a subdirectory named build/, by entering::

     mkdir pscf
     cd pscf
     mkdir build

Obtaining the Source Code
-------------------------
The source code for pscf is stored in a repository on the github.com server, at:

    https://github.com/dmorse/pscf

A copy of the source code may be obtained either, by:

    * Downloading a zip file, or 
    * Using git to clone the source code.  

To download a zip file:

    * Point a browser at the pscf github repository

    * Click the button labelled "Download ZIP" near the upper right corner. 
      On Mac OS X and most linux systems, this will download a directory 
      named pscf-master into the users Downloads directory.

    * Move the pscf-master/ directory into the pscf/ directory, making it
      a subdirectory of pscf/

    * Rename the pscf/pscf-master/ directory as repo/, by changing directory
      to pscf and then entering::

         mv pscf-master repo

To use git to clone the repository, after git is installed on your machine:

    * Change directory to the pscf directory.

    * Clone the repository, by entering::

          git clone https://github.com/dmorse/simpatico.git

    * This should create a subdirectory of pscf/ that is also named pscf/. 
      To avoid confusion, we recommend that you change the subdirectory 
      name to pscf/repo/, exactly as described above for the case of a 
      directory created from a zip file. 

Compile and Install
-------------------

Before attempting to compile, you must install all required dependencies, 
by following instructions given below for each operating system, create
a pscf/ directory tree with a build/ subdirectory, and obtain the source 
code. At this point you should have a pscf/ directory structure::

    pscf/
       build/
       repo/

in which the build/ subdirectory is empty and the repo/ subdirectory 
contains the pscf source code, as obtained from the github repository.

To compile and install, cd to the build/ directory and, from there,
enter::

   > cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../repo
   > make -j 4
   > make install 

In the "cmake" command, the string "/path/to/install" is the root 
of path used for installation.  The last argument "../pscf". If 
you use "-DCMAKE_INSTALL_PREFIX=.", the executable and other files 
that you generate will be installed in tree rooted at the build 
directory (e.g., pscf-build). The final pscf executable is 
self-contained and can be copied to wherever you want after it is 
created.

Wherever you install the executable file, you will need to make sure that 
directory containing the executable (or a symlink to the executable) is 
in the bash PATH variable, so that the operating system can find the 
executable when it is invoked by name.

