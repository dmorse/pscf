
.. _package-sec:

Creating Binary Installers
===========================

This page for how to use cmake to create a binary .dmg installer file for Mac OS
X and a .deb or .rpm package for different variants of linux.  This information 
is not relevant to most users, and is provided primarily for core developers. 

The first step in creating a package, for any operating system, is to follow the 
instruction given on the page about :ref:`compile-cmake-sec` for installing 
dependencies and obtaining the source code. The remaining instructions given 
here assume that the dependencies are already installed, and that a copy of the
source code has been installed within the directory structure described in the
instructions for compiling from source.

Mac OS X
--------

On Mac OS X, after installing all dependencies and installing a copy of the
source code in a repository named git/, one must::

    * Change directory to the pscf/cmake directory.

    * From the directory, enter::

          > cmake -DBUILD_DMG=1 -DCMAKE_INSTALL_PREFIX=. ../git

    * Then enter::

          > make -j 4
          > make package

This should create an installer file named pscf<version>.dmg in the pscf/cmake 
directory.

Linux (Fedora or Ubuntu)
------------------------

On a linux system, one must::

    * Change directory to the pscf/cmake directory.

    * From the directory, enter::

          > cmake -DCMAKE_INSTALL_PREFIX=.. ../git

    * Then enter::

          > make -j 4
          > make install
          > make package

The make install command should install the software in subdirectories of
pscf named bin, lib and share. The package command should create either a
file named pscf<version>.rpm, if built on a system that uses the redhat 
packaging manager, or a file named pscf<version>.deb if built on a system
that uses .deb files.

