
.. _install-compile-sec:

Compiling from Source, using cmake
==================================

The most recent version of the source code for PSCF can be obtained either by
cloning the code or by downloading a tar archive. We recommend using git if
possible, since it makes it simple to update the code, but both procedures 
for obtaining the source code are described in what follows.

The source code can be compiled using either cmake build system or, without
cmake, using a simple makefile that is provided in the src/build directory. 
One advantage of using cmake is that cmake can generally find the paths to
the Lapack and FFTW libraries upon which PSCF depends, whereas you need to
figure out the locations of the correct libraries yourself if you use the
simple makefile. The following instructions focus on compiling and installing
PSCF using cmake. A separate section at the end describes how to use a simple
makefile to compile from source without using cmake.

The following software packages must be installed and accessible before 
attempting to compile PSCF if you plan on using git to obtain the source
code and cmake to build the software:

   * git (in order to clone the source code)
   * cmake (to build a makefile)
   * a fortran 90 compiler 
   * LAPACK linear algebra library
   * FFTW version 3.x fast fourier transform library

In what follows we give instructions for how to build pscf on different
opeating systems. Instructions for most systems assume that you will use 
the free gfortran Fortran 90 compiler, which is part of the Gnu Compiler 
Collection (gcc) suite of compilers. 

Mac OS X
--------

Installing XCode
^^^^^^^^^^^^^^^^

To create an environment in which you can compile from source on OSX, you 
will generally first need to install the apple XCode development environment.
XCode is available gratis from the app store, but is a large package that can
take a long time to install (do it with a good internet connection).  The 
XCode package contains git, so it is not necessary to install git separately
The OXS operating system also appears to come with a version of LAPACK, and 
the BLAS library upon which it depends.

Package Managers: HomeBrew vs. MacPorts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Installing dependencies via Homebrew
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install from a command line terminal using homebrew::

   > brew install cmake
   > brew install gcc --with-fortran
   > brew install fftw --with-fortran

Installing dependencies via Macports
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Getting the source code
^^^^^^^^^^^^^^^^^^^^^^^

To obtain the most recent PSCF source code from github::

   > git clone https://github.com/dmorse/pscf.git

Compile and Install
^^^^^^^^^^^^^^^^^^^
Before compiling, you should make a new directory in which 
the program will be built "out-of-source". This build directory
should not be subdirectory of the pscf/ directory. The following 
assumes that the build directory is called pscf-build, and that 
it and pscf/ are subdirectories of the same parent directory.

Starting from the common parent directory of pscf/ and pscf-build/,
enter::

   > mkdir pscf-build
   > cd pscf-build
   > cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../pscf
   > make -j 4
   > make install 

In the "cmake" command, the string "/path/to/install" is the root 
of path used for installation. 
The last argument "../pscf" If you 
use "-DCMAKE_INSTALL_PREFIX=.", the executable and other
files that you generate will be installed in tree rooted
at the build directory (e.g., pscf-build). The final
pscf executable is self-contained and can be copied to 
wherever you want after it is created.

For developers: To build a Mac OSX .dmg binary installer,
as well as .tar and .zip source code archive files, when
working on a Mac, after completing compilation and 
installation, enter::

   > make package

Ubuntu or Debian Linux
----------------------

Use the Ubuntu software manager or the command line apt-get 
utility to install the following packages:

   * git
   * cmake
   * gfortran
   * libfftw3-dev
   * liblapack3

To obtain the PSCF source code from github, as for OS X, enter::

   > git clone git@github.com/dmorse/pscf.git

The steps to compile and install are also the same as for Mac OSX::

   > mkdir pscf-build
   > cd pscf-build 
   > cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../pscf
   > make -j 4
   > make install 

On linux, an executable file will be installed in the bin directory of the
directory "/path/to/install" that is passed to cmake.  The executable is 
movable, so you can place the executable in the build directory by entering

   > cmake -DCMAKE_INSTALL_PREFIX=.  ../pscf

(where the "." after the = sign represents the current directory), and then 
move the file to wherever you want. 

Wherever you install the executable file, you will need to make sure that 
directory containing the executable (or a symlink to the executable) is 
in the bash PATH variable, so that the operating system can find the 
executable when it is invoked by name.

Developers: To build .deb package for installation of binary executables 
on other Ubuntu and debian systems, as well as .tar and .zip source code 
archives, after installing on your machine, simply enter::

   > make package

To check the .deb file for semi-detailed information::

    # This extracts multiple files
    ar -vx pscf-1.0.0-Linux.deb
    # See the files that would be installed
    tar tvfz data.tar.gz 

Fedora / Redhat Linux
---------------------

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

Linux Modules and Intel Compiler
--------------------------------

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

