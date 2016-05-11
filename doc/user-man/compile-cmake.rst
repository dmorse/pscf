
.. _install-compile-cmake-sec:

Compiling from source, using cmake
==================================

.. _install-compile-cmake-overview-sub:

Overview
--------

In this page, we discuss instructions for how to compile PSCF from source
using the cmake build system, on either Mac OS X or linux operating systems.
It is also possible to compile the code using the unix make utility alone,
using a makefile that is provided with the source code. Compilation using 
make alone is described on a separate page on :ref:`install-compile-make-sec`.
The advantage of using cmake is that cmake can generally find the paths to
the Lapack and FFTW libraries upon which PSCF depends, whereas you are more
likely to need to figure out the locations of these libraries yourself if 
you use make.

Compiling with cmake involves the following steps:

    * :ref:`install-compile-cmake-dependencies-sub`
    * :ref:`install-compile-cmake-getsource-sub`
    * :ref:`install-compile-cmake-compile-sub`
    * :ref:`install-compile-cmake-paths-sub`

Each of these steps is explained in greater detail below. 

Only the first step, installing external dependencies, is substantially 
different for different operating systems. We thus give separate 
instructions for Mac OS and different Linux distributions for this part 
of the process.

To obtain the source code from the github repository in which it is 
stored one can either:

    * Download a zip file
    * Use the git version control software to clone the repository

We recommend using git, since this makes it simple to update the code 
later, but both procedures are described below.

.. _install-compile-cmake-dependencies-sub:

Installing Dependencies
-----------------------

The following software packages must be installed before using cmake 
to compile PSCF, if you plan on using git to obtain the source code:

   * git (in order to clone the source code)
   * cmake (to build a makefile)
   * a Fortran 90 compiler (to compile the source code)
   * Python (used by the build system)
   * LAPACK linear algebra library
   * FFTW version 3.x fast fourier transform library

You do not need to install git if you plan to simply download the 
source code rather than using git. On a Mac, some of these packages
come bundled with the XCode development environment, which must in
any case be installed before you try to compile software from source
on Mac OS X. Python is included as part of the most common linux 
distributions, and is also bundled with recent versions of Mac OS X. 

In what follows, we will assume that you plan to use free gfortran 
Fortran compiler, which is part of the Gnu Compiler Collection (gcc) 
suite of compilers, and give instructions for installing this
compiler. 

Mac OS X
~~~~~~~~~

**Installing XCode**

To create an environment in which you can compile from source on OSX, you 
will generally first need to install the apple XCode development environment.
XCode is available gratis from the app store, but is a large package that can
take a long time to install (do this with a good internet connection).  The 
XCode package contains git, so it is not necessary to install git separately.
The Mac OS X operating system also appears to come with a version of LAPACK, 
and the BLAS library upon which it depends. Neither the operating system nor
XCode provide cmake, gfortran, or FFTW.

**Checking for Python**

A Python interpreter is included in recent versions of Mac OS X. To check if Python is already installed, enter the command::

   > which python

or::

   > which python2.7

from a unix terminal window. If the system responds to either of these
commands by writing a location for a python executable file to the
terminal, then Python is already installed. If nothing is written to 
the terminal in response, then either Python is not installed or the 
operating system doesn't know where to find it. Instructions for
installing python, if necessary, are given below.

**Package Managers: HomeBrew vs. MacPorts**

The remaining dependencies (cmake, gfortran and FFTW) can be most easily 
installed using either the MacPorts or Homebrew package manager systems.  
These are both systems for managing open-source unix software on the unix 
subsystem of the Mac OSX.  The choice between these package managers is 
up to you, but you should avoid using both on the same machine.  If either 
Homebrew or MacPorts is already installed and in use on your Mac, use the 
existing system, and do not install the other, because they do not play 
well together.  If neither Homebrew or MacPorts is installed, we have a
slight preference for Homebrew, which we find makes it slightly easier to 
install the dependencies required by PSCF. We have succeeded in building 
PSCF using both package managers on different machines that are running 
the latest version of Mac OS X (El Capitan, X 10.11) Instructions for 
both package managers are given separately below.

**Installing dependencies via Homebrew**

To install from a command line terminal using homebrew::

   > brew install cmake
   > brew install gcc --with-fortran
   > brew install fftw --with-fortran

If python is required, enter::

   > brew install python

**Installing dependencies via Macports**

After MacPorts is installed, to install the required dependencies 
using the most recent version of the gnu compiler collection (gcc), 
which is gcc 5.X at the time of writing, enter::

   > sudo port install cmake
   > sudo port install gcc5
   > sudo port install fftw-3 +gfortran

If python is required, enter::

   > sudo port install python27

Note that MacPorts (unlike homebrew) requires you to use "sudo"
to execute installation with superuser/administrator privileges, 
and so will ask for a password after each of the above commands.

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

Use the Ubuntu software manager or the command line apt-get utility 
to install the following packages:

   * git
   * cmake
   * gfortran
   * libfftw3-dev
   * liblapack3

To use apt-get from the command line, enter::

   > sudo apt-get update
   > sudo apt-get install git
   > sudo apt-get install cmake
   > sudo apt-get install gfortran
   > sudo apt-get install libfftw3-dev
   > sudo apt-get install liblapack3

Fedora Linux
~~~~~~~~~~~~

Instructions for Fedora are similar to those for Ubuntu, except that one 
should use the native yum command line package manager or the Fedora 
graphical software manager to install dependencies. The required Fedora 
packages are:

   * git
   * cmake
   * gcc-gfortran
   * lapack-devel
   * fftw-devel

To install these packages from the command line, enter::

   > sudo yum install git-all
   > sudo yum install cmake
   > sudo yum install gcc-gfortran
   > sudo yum install lapack-devel
   > sudo yum install fftw-devel

For Fedora 22 and later, you may use the command "dnf" rather than "yum" 
to use the an updated version of the yum package manager. Instructions for 
obtaining source code, compiling and installing are the same as for Mac
OS X and Ubuntu operating systems.

Using Linux Modules
~~~~~~~~~~~~~~~~~~~~

Many large multi-user computer clusters use linux modules to allow users
to load software packages that they require, chosen from among a list of
available modules. The following instructions describe how to load the
required modules to build PSCF in a user directory on the Minnesota 
Supercomputer Institute (MSI) Mesabi computer, using linux modules and 
the Intel compiler.  Similar instructions should apply to other large 
clusters that use linux modules.

To load the required modules on Mesabi at MSI, and also choose the Intel
compiler, enter::

   > module load cmake
   > module load intel mkl
   > module load fftw

The remaining instruction for how to obtain and compile the source code 
are generally similar to thos given for OSX or Linux. The only difference 
is that, to use the Intel compiler, one must tell cmake to use the Intel 
compiler by adding the option "-DUSE_INTEL=1" to the cmake command. The 
required command is thus::

   > cmake -DUSE_INTEL=1 -DCMAKE_INSTALL_PREFIX=/path/to/install ../pscf

More generally, using the "-D" to define USE_INTEL=1 to search for an 
Intel compiler rather than using gnu fortran, on any operating system.

.. _install-compile-cmake-getsource-sub:

Obtaining the Source Code
-------------------------

We assume in what follows that you will use cmake to perform an
"out-of-source" build, in which all of the files generated during 
compilation are placed in a directory tree outside the source code tree. 
To begin, we recommend that you create a directory named pscf/ with a 
subdirectory named cmake/, by entering::

     mkdir pscf
     cd pscf
     mkdir build

The directory named cmake/ will be used as the build directory. The source 
code will be placed in another subdirectory of pscf/, which we will call 
git/ in this example, since it contains the contents of the git repository.

The source code for pscf is stored in a repository on the github.com 
server, at: 

      https://github.com/dmorse/pscf

A copy of the source code may be obtained either, by:

    * Downloading a zip file, or 
    * Using git to clone the source code.  

To download a zip file:

    * Point your browser at the pscf github repository.

    * Click the "Download ZIP" button near the upper right corner 
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

          git clone https://github.com/dmorse/pscf.git

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

.. _install-compile-cmake-compile-sub:

Compiling and Installing
------------------------

**Choose an Install Directory**

After installing all dependencies and obtaining the source code, you are ready 
to compile PSCF. Before compiling the code, you need to decide where you would 
like to install the pscf executable, along with several other executable 
scripts and python files.  The build system created by cmake will install 
these files in subdirectories of a directory that we will refer to as the 
install directory that you can specify on the command line of the "cmake"
command. After installation, the install directory (denoted by install/
below) will contain the following three subdirectories::

    install/
       bin/
       lib/
       share/

After installation, the bin/ subdirectory will contain the pscf executable 
and other executable files, the lib/ subdirectory will contain python 
modules and matlabe files and the share/ directory will contain several
text files containing information about the program.

We recommend that you choose one of the three following three possible 
locations for the install directory for pscf:

   * The pscf/ directory that contains the cmake/ and git/ subdirectories.

   * A standard location for installing software within your user directory.
     within your user directory. 

   * The system-wide /usr/local directory.

If you choose to install software within a standard location within your
user directory, one common choice for this is a hidden directory of your 
home directory named .local.

One advantage of the first two options listed above is that both install 
all of the software within your user directory, and thus do not require 
adminstrative privileges. This also makes it somewhat easier for you to
see what you have installed and remove it if ever desired. The further 
logistical advantage of the first option, of installing within the pscf/ 
directory that also contains the source code, is that it keeps all of the
files associated with PSCF in a single directory tree within the user 
directory.

The main disadvantage of both the first and second options is that, 
because both install files within your user directory, they both guarantee 
that you will have to modify some operating system environment variables 
in order to allow the operating system to find the PSCF executable and to 
allow the python intepreter to find python modules that are provided
to faciliitate data analysis. Conversely, the advantage of installing 
in /usr/local is that, because this puts the executable in a standard 
location, the operating system should be able to automatically find 
the pscf executable.  Instructions for modifying the relevant environment 
variables, if necessary, are given below.

**Invoke cmake**

The first step of compiling with cmake is to invoke the cmake command
in order to construct a set of makefiles that contain instructions for
building the system. To begin, change directory (cd) to the pscf/cmake/ 
directory. Then make sure the cmake/ directory is empty, and remove any 
contents if necessary. From there, enter::

   > cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ../git

In this command, the string "/path/to/install" denotes the path to the 
root of the install directory.  The last argument, "../git", is the
relative path to your copy of the source code repository, in pscf/git, 
from the pscf/cmake directory. 

To install within in the pscf/ directory tree, you would enter::

   > cmake -DCMAKE_INSTALL_PREFIX=..  ../git

where ".." represents the pscf/ directory, which is the parent of the
pscf/cmake directory from which the command is issued. This will cause 
the later creation of bin/, lib/ and share/ subdirectories of the 
pscf/ directory, alongside the cmake/ and git/ subdirectories.

To install in the .local subdirectory of your home directory, instead
enter::

   > cmake -DCMAKE_INSTALL_PREFIX=~/.local  ../git

in which the tilde (~) is linux shortand for the users home directory.

Finally, to install in the /usr/local directory, you need adminstrator
privileges on your machine, and would enter::

   > sudo cmake ../git

In this case, you must use the "sudo" command to apply the command 
with "super-user" or administrator privileges, and you will be prompted 
for your password. No "-DCMAKE_INSTALL_PREFIX=" option is required in 
this case, however, because /usr/local is the default installation 
that will be use by cmake if no alternative is specified.

**Invoke make**

The cmake command described above should create several subdirectories 
of the pscf/cmake/ directory, which contain makefiles with instructions 
for building pscf. To actually compile and install the program, simply
enter::

   > make -j4
   > make install 

from the pscf/cmake directory.  The "-j4" option simply instructs the
make utility to use up to 4 processor cores to compile, if available,
to speed up compilation. It is not required. The first "make" command
compiles the code and places all the files generated by compilation 
in the pscf/cmake directory. The "make install" command installs files 
in the chosen installation directory.

After the "make install" finishes execution, check that your chosen 
install directory contains subdirectories named bin/, lib/ and share/, 
and that the the bin/ subdirectory contains an executable file named pscf, 
along with several executable scripts whose names begin with the suffix 
"pscf-...". One of these should be a bash script named "pscf-env".

.. _install-compile-cmake-paths-sub:

Modifying Search Paths
-----------------------

If you install pscf in a directory within your home directory tree, 
you may need to modify a few environment variables to allow the
operating system to find the pscf program when it is invoked from 
the command line by name, and to allow the python interpreter to find 
some associated python modules that are useful for data analysis. 

**Changing Paths**

The simplest way to make the required changes to your user environment
is to cd to bin/ subdirectory of the root install directory and, from
there, enter::

    source ./pscf-env

This will run a script that is installed by PSCF, which adds the 
appropriate paths to your PATH and PYTHONPATH environment variables.

Alternatively, to make the required changes manually, you could simply 
enter the commands::

    PATH=$PATH:install/bin
    PYTHONPATH=$PYTHONPATH:install/lib/python2.7/site-packages

where "install" denotes an absolute path to your chosen installation
directory.

**Making Changes Permanent**

The above procedures (running pscf-env script or manually setting the
relevant environment variables) only modifies the $PATH and $PYTHONPATH
variables temporarily, until you close the terminal window or log out.
To have the appropriate directories added to these environment variables 
automatically whenever you log in or open a terminal, simply add the 
command::

   source install/bin/pscf-env 

to the .bashrc file or (on Mac OS X) .profile configurtion file
in your home directory. Here, the string "install/" is a placeholder 
for the absolute path to the pscf install directory.

**Configuration files: Linux vs. Mac OS X**

On linux, after a user logs in, the operating system looks for a file 
in the user directory named .profile or .bash_profile (in that order)
and executes the first of these files that finds, if any. When you 
open a new interactive shell that is not a login shell, e.g., by
opening a new termiinal, it instead looks for and (if it exists)
executes a file named .bashrc in the users home directory. To make
sure that the modifications of the environment are applied to both 
login and non-login terminals, the .bashrc file is normally executed 
by the .profile or .bash_profile file, by a command such as::

    if [ -f "${HOME}/.bashrc" ]; then
	source "${HOME}/.bashrc"
    fi

This part of the .profile or .bash_profile file checks if there is 
a .bashrc file in the users home directory and, if one is found, 
executes that file. With this configuration, commands that set up
environment variables should be added to the .bashrc file.

On Mac OS X, the Mac Terminal program instead executes the .profile
script whenever you open a terminal, rather than using different 
files for login and non-login terminals. The Mac Terminal program 
thus thus does not ever directly execute the .bashrc file. A Mac 
user that always uses the Mac Terminal program could thus either 
use the procedure described above (which would still work correctly), 
or simply place all commands that customize the user environment 
into the .profile script.
