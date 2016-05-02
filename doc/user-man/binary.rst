
.. _install-binary-sec:

Binary Installation
===================

The steps to install a precompiled binary are different for different 
operating systems, and somewhat different for different distributions of 
linux.  The installer for Mac OS X installs a completely self-contained 
package that includes copies of all required external libraries. Instructions 
for installing binary .rpm and .deb packages for Redhat and Debian linux 
systems instead require the user to install packages containing the FFTW 
fast Fourier transform and Lapack linear algebra libraries before installing 
PSCF.

Mac OSX
-------

The procedure for installing PSCF on a Mac using a binary installer is 
similar to that for installing any application on a Mac:

  * Download the Mac .dmg installer from the PSCF home page. This is a
    file with a name of the form pscf-<version>-Darwin.dmg.

  * Open the .dmg file, and drag and drop the pscf_terminal icon file 
    into the Applications folder.

To run the program, then simply double click the pscf_terminal application.  
This will open up a yellow terminal window from which you can use standard 
unix commands to navigate within the directory structure of your Mac, and 
from which you can invoke the pscf command.

The first time you attempt to run pscf_terminal, Mac OSX security settings 
will prevent the application from starting because it is from an "Unknown 
Developer".  You will need to add an exception for this software. The 
following instructions are provided by Apple on 
https://support.apple.com/kb/PH18657?locale=en_US: 

1. In the Finder, locate the app you want to open. Don’t use Launchpad to 
   do this. Launchpad doesn’t allow you to access the shortcut menu.

2. Press the Control key, then click the app icon, then choose Open from 
   the shortcut menu.

3. Click Open.

The app is saved as an exception to your security settings, and you will
be able to open it in the future by double-clicking it just as you can 
any registered app.

Ubuntu or Debian Linux
----------------------

Ubuntu and Debian distributions of the linux operating systems both use 
variants of the debian package management system, which uses .deb package 
files.  To install from binary on Ubuntu:

  * First use the Ubuntu software center graphical installer or the 
    apt-get command line utility to install the following packages:
   
        - libfftw3-3
        - liblapack3

    To install these packages using apt-get, enter::

        sudo apt-get libfftw3-3
        sudo apt-get liblapack3

  * Download the .deb package from the PSCF home page. This is a file
    with a name of the form pscf-<version>-Linux.deb.

  * Install the pscf package by running::

        dpkg -i pscf-<version>-Linux.deb

  * If the above command fails because of a missing dependence, try 
    running::

        apt-get install pscf-<version>-Linux.deb

    This sometimes allows apt-get to attempt to fetch missing dependencies.


Fedora / Redhat Linux
---------------------

Redhat distributions of the linux operating systems, including Fedora 
and CentOS, use package management systems that use .rpm package files. 
Instructions are similar to those for Ubuntu/Debian, except for the use 
of a different package file format and package manager. In this case:

  * Use the yum command line utility or the Fedora graphical package 
    manager to install the following packages from the appropriate
    Fedora repository:
   
        - fftw
        - lapack

    To install these packages using yum from the command line, enter::

        > sudo yum install lapack
        > sudo yum install fftw

  * Download the pscf-<version>-Linux.rpm package from the PSCF home 
    page.

  * To install, enter::

        sudo rpm -Uvh pscf-<version>-Linux.rpm

