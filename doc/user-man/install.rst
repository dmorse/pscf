************
Installation
************

PSCF was designed to be used from the command line in a unix-like environment, 
and can be installed on either the linux operating system or in the Mac OS X.
There are two different ways to install a working version of PSCF in either 
operating system:

   * Install a precompiled binary
   * Compile and install from source code

Installation of a precompiled binary is somewhat simpler, but does not allow
the user to modify and recompile the program. The required steps for both
procedures are described separately below, including separate discussions of
installation in different operating system environments.

Binary Installation
===================

The steps to install a precompiled binary are different for different operating
systems, and somewhat different for different distributions of linux. 

Mac OSX
-------

The simplest way to install PSCF on a Mac is to use the Mac installer. The 
procedure is similar to that for installing any other application on a Mac:

  * Download the Mac pscf<version>.dmg installer from the PSCF home page

  * Double click on the installer file

  * Drag the file pscf_terminal application into the applications folder.

Double clicking the pscf_terminal application will open up a yellow terminal
window from which you can use standard unix commands to navigate within the
directory structure of your Mac, and from which you can invoke the pscf 
command.

Ubuntu or Debian Linux
----------------------

The Ubuntu and Debian distributions of the linux operating systems both use 
package management systems that use deb package file format, with the .deb 
extension. To install on Ubuntu or Debian:

  * Download the pscf<version>.deb package from the PSCF home page

  * Install the package by running::

       dpkg -i pscf<version>.deb

Fedora / Redhat Linux
---------------------

Fedora and Redhat distributions of the linux operating systems use package 
management systems that use .rpm package files. Instructions are similar to
those for Ubuntu/Debian, except for the use of a different package file 
format and package manager. In this case:

  * Download the pscf<version>.rpm package from the PSCF home page

  * To install, enter::

        sudo rpm -Uvh pscf<version>.rpm

Compiling from Source
=====================

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

We found that the gcc-5 port installs the gfortran Fortran 90 
compiler at /usr/opt/local/bin/gfortran-mp-5 . Versions compiled 
with other versions of gcc (e.g., 4.9) seem to be placed in the 
same directory with a different numerical suffix, e.g., 
gfortran-mp-49.  CMake was unable to find this compiler 
executable without help.  To remedy this, you should set the 
FC environment variable (the path to a Fortran compiler) to 
point to the absolute path to the gfortran executable before
attempting to compile, by entering, for example::

   > FC=/usr/opt/local/bin/gfortran-mp-5
   > export FC

If expect to compile this and other fortran programs repeatedly, 
you may want to put this in your .profile or .bashrc bash 
configuration file.

Getting the source code
^^^^^^^^^^^^^^^^^^^^^^^

To obtain the most recent PSCF source code from github::

   > git clone git@github.com/dmorse/pscf.git

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
   * fftw3-dev
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

Compiling via make, without cmake
---------------------------------

It is also possible to compile using a Makefile in the src/build directory. 
This does an "in source" build, in which all of the files generated during 
compilation are placed in the pscf/src/build/ directory. The instructions 
for doing this are the same on any unix-like operating system. The main 
difference among different unix environments is the locations of the 
required libraries. 

To compile the code in this way, you should:

   * cd to the pscf/src/build directory
   * Examine and edit the Makefile (as discussed below)
   * Enter 'make pscf' from within src/build.

These steps are described in more detail below

Customize the Makefile:
^^^^^^^^^^^^^^^^^^^^^^

In Makefile in the src/build directory, you will need to set values for a 
set of macro variables to values appropriate to your system. Makefile 
variables you may need to reset are:
 
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

The makefile contains values appropriate for a number of different common 
environments, most of which are commented out. Make sure you give only one
definition for each variable, and that you comment out any definitions you 
are not using.

Compile and Link
^^^^^^^^^^^^^^^^

To compile and link, from the src/build directory, issue the
command::

   > make pscf

This should fill the src/build directory with .o and .mod files, and 
create an executable $(BIN)/$(EXE). By default, this will create a program 
named pscf in the pscf/bin directory. The executable file can be relocated 
to somewhere else if you desire.

To invoke the program, you will either need to:

   * Invoke the program using an absolute path name

   * Add the directory containing your executable to your command search
     PATH variable. To do so, enter:

         PATH=$PATH:~$(SCF)/bin
         export path

     where $(SCF) should be replaced by the actual absolute path to the
     pscf/ directory. You may want to add this to your .bashrc or .profile 
     file so that this directory is added to your path when automatically 
     when you log in.

   * Move pscf to a directory such as /usr/local/bin that is already in 
     your $PATH. 

Cleaning Up
^^^^^^^^^^^
	
To remove all of the .o amd .mod files from the src/build directory, as 
well as any editor buffer files with a ~ suffix from src tree, enter::

   > make clean

