
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

    * Install dependencies (other software on which PSCF depends)
    * Obtain the source code
    * Compile and install PSCF 

Each of these steps is explained in greater detail below. 

Only the first step, of installing external dependencies, is 
substantially different for different operating systems, because
different operating systems use different systems to manage external
software packages and contain different software pre-installed. We 
thus give separate instructions for Mac OS and different Linux 
distributions for this process.

In order to obtain the source code, on any operating system, one 
can either:

    * Download a zip file
    * Use the git version control software to clone the repository

We recommend using git if possible, since this makes it simple to 
update the code later, but both procedures are described below.

The following software packages must be installed before using cmake 
to compile PSCF, if you plan on using git to obtain the source code:

   * git (in order to clone the source code)
   * cmake (to build a makefile)
   * a Fortran 90 compiler 
   * LAPACK linear algebra library
   * FFTW version 3.x fast fourier transform library

You do not need to install git if you plan to simply download the 
source code rather than using git. On a Mac, some of thes packages
come bundled with the XCode development environment, which must in
any case be installed before you try to compile from source. In 
what follows, we will assume that you plan to use free gfortran 
Fortran compiler, which is part of the Gnu Compiler Collection 
(gcc) suite of compilers. 

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
the BLAS library upon which it depends. It does not include cmake, gfortran,
or FFTW.

**Package Managers: HomeBrew vs. MacPorts**

The remaining dependencies (cmake, gfortran and FFTW) can be most easily 
installed using either the MacPorts or Homebrew package manager systems.  
These are both systems for managing open-source unix software on the unix 
subsystem of the Mac OSX.  The choice between these package managers is 
up to you, but you should avoid using both on the same machine.  If either 
Homebrew or MacPorts is already installed and in use on your Mac, use the 
existing system and do not install the alternative, because they do not 
play well together.  If neither Homebrew or MacPorts is installed, we have 
slight preference for Homebrew, which made it slightly easier to install 
the dependencies required by PSCF. We have succeeded in building PSCF using 
both package managers on different machines that are running the latest
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

Obtaining the Source Code
-------------------------
We assume in what follows that you will use cmake to implement on "out-of-source" build, in which all of the files generated during compilation are placed in a different directory than the source code. To begin, we suggest that you create a directory named pscf/ with a subdirectory named cmake/, by entering::

     mkdir pscf
     cd pscf
     mkdir build

The source code will be placed in another subdirectory of pscf/ which we will call git/ (for repository).

The source code for pscf is stored in a repository on the github.com server, at:

      https://github.com/dmorse/pscf

A copy of the source code may be obtained either, by:

    * Downloading a zip file, or 
    * Using git to clone the source code.  

To download a zip file:

    * Point your browser at the pscf github repository.

    * Click the button labelled "Download ZIP" near the upper right corner 
      of that web page. On Mac OS X and most linux systems, this will 
      create a subdirectory named pscf-master with your Downloads 
      directory.

    * Move the pscf-master/ directory into the pscf/ directory that you
      just created.

    * Rename the pscf/pscf-master/ directory as git/, by changing directory
      to pscf and then entering::

         mv pscf-master git

To use git to clone the repository, after git is installed on your machine:

    * Change directory to the pscf directory.

    * Clone the repository, by entering::

          git clone https://github.com/dmorse/simpatico.git

    * This should create a subdirectory of pscf/ that is also named pscf/. 
      To avoid confusion, we recommend that you change the subdirectory 
      name to pscf/git/, exactly as described above for the case of a 
      directory created from a zip file. 

At this point, by either method, you should have pscf/ directory structure::

    pscf/
       cmake/
       git/

in which the cmake/ subdirectory is empty and the git/ subdirectory contains 
the contents of github repository, including the source code.

Choosing an Install Directory
-----------------------------

After installing all dependencies and obtaining the source code, you are
ready to compile PSCF. 

Before compiling the code, you need to decide where you would like to install 
the pscf executable, along with several executable scripts, python modules, 
and matlab files. The build system created by cmake will install these files 
in subdirectories of a directory that we will refer to as the install directory.  
Specifically, it will install the pscf executable and several executable scripts 
in the bin/ subdirectory of the install directory, install python modules and
matlab scripts in different subdirectories of the lib/ subdirectory, and install 
several text files in the share/ subdirectory.  After installation, the install 
directory, denoted below by install/, will thus contain three subdirectories::

    install/
       bin/
       lib/
       share/

The build system will create these three subdirectories if they do not 
already exist. The default choice for the install directory is the system
/usr/local directory, which is a standard location on linux for a system 
administrator to install 'local' software on linux that is not part of the 
linux distribution.

We suggest that you consider the following three possible locations for the
install directory for pscf:

   * The pscf/ directory, which also contains the source code. 

   * A standard installation directory within your user directory.

   * The system /usr/local directory (the default).

The advantage of the first two options is that both of them install all 
of the software within your user directory tree, and thus do not require
adminstrative privileges. The further logistical advantage of the first 
option (installing within pscf/) is that it keeps all of the files in a 
single directory tree within your user directory that only contains files 
associated with pscf/, which makes it particularly easy to erase everything 
and start over if desired. 

The disadvantage of the first and second options, which both install files 
within your user directory, is that both of them will require that you
modify some operating system environment variables in order to conveniently
use pscf. Specifically, if you install files in non-standard locations, you
you will need to modify the PATH and PYTHONPATH environment variable thats
allow the operating system and the python interpreter, respectively, to 
find executable files and python modules when referrerd to by file name.
Conversely, the advantage of installing in the /usr/local directory is
that doing so causes executable files and python modules to be placed in
standard locations where they will be found automatically.

The only advantage of the second option (installation in a standard 
location within your user directory tree) relative to the first is that, 
if you plan to install multiple software packages from source and install 
all of them in this location, you can configure your enviroment to always 
look in appropriate subdirectories of this user directory for files, so
that you do not need to further modify environment variables every time 
you install new software within your user directory tree.  If you choose 
this option, it is conventional on some versions of linux to install 
software in a hidden subdirectory of you home directory named ".local".
We recommend this practice whether you are using linux or using the unix
command line interface of Mac OS X.  Note the dot in the beginning of 
the name ".local", which makes it a hidden directory that will not show 
up when you use "ls" to list files and directories in home directory, 
unless you add the "-a" option, as "ls -a", to show hidden files and 
directories.  Installing in this location will cause the creation of 
a tree of subdirectories of your private ${HOME}/.local directory that 
is analogous to the structure of the /usr/local directory.

Compiling and Installing
------------------------

As the first step of compiling and installing, change directory to the 
pscf/cmake/ directory. Then make sure the cmake/ directory is empty 
(removing all contents if necessary) and, from there, enter::

   > cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../git

In this command, the string "/path/to/install" denotes the path to
the root of the install directory.  The last argument "../git" is the
relative path to your copy of the source code repository from the 
pscf/cmake directory. 

To install in the pscf/ directory tree, you would thus enter::

   > cmake -DCMAKE_INSTALL_PREFIX=..  ../git

where ".." represents the parent pscf/ directory. This will cause the 
creation of bin/, lib/ and share/ subdirectories of the pscf/ directory, 
alongside cmake/ and git/. 

To install in the .local subdirectory of your home directory, instead
enter::

   > cmake -DCMAKE_INSTALL_PREFIX=~/.local  ../git

in which the tilde (~) is linux shortand for the users home directory.

Finally, to install in the /usr/local directory, you need adminstrator
privileges on your machine, and would enter::

   > sudo cmake ../git

In this case, you must use the "sudo" command to apply the command 
with "super-user" or administrator privileges, and you will be prompted 
for your password. No -DCMAKE_INSTALL_PREFIX=" option is required in 
this case /usr/local is the default installation location.

The cmake command described above should create several subdirectories of 
the pscf/cmake/ directory, which will contain files with instructions for 
building pscf. To finish compiling and installing, simply enter::

   > make -j 4
   > make install 

from the pscf/cmake directory. 


After the "make install" finishes execution, you can check that your chosen 
install directory contains subdirectories named bin/, lib/ and share/, and 
that the the bin/ subdirectory contains an executable file named pscf, along 
with several executable scripts whose names begin with the suffix "pscf-...".

