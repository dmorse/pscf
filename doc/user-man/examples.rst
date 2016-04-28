
.. _install-examples-sec:

Installing the examples
=======================

A library of examples of input files for pscf simulations is provided as
a separate git repository. The best way to become familiar with the format
of the pscf input and output files is to look at and run some of these
examples.  The following instructions assume that you will install the 
directory containing these examples as a directory named pscf/examples, 
alongside a pscf/git directory that contains the source code repository.

**Obtaining examples**

The repository containing the examples is located on the github.com server 
at: 

      https://github.com/dmorse/pscf-examples

As for the source code repository, a copy of the examples may be obtained 
either by downloading a zip file from the repository web page or using the 
git version control system to clone the repository. The instruction given 
below for both methods are very similar to those given for obtaining the 
source code:
 
To download a zip file:

    * Point your browser at the pscf-examples github repository.

    * Click the "Download ZIP" button near the upper right corner 
      of that web page. 

    * Move the pscf-examples-master/ directory into the pscf/ directory. 

    * Rename the pscf/pscf-examples-master/ directory as pscf/examples.

To use git to clone the repository, after git is installed on your machine:

    * Change directory to the pscf directory.

    * Clone the repository, by entering::

          git clone https://github.com/dmorse/pscf-examples.git

    * Change the name of the pscf/pscf-examples directory to pscf/examples/

At this point, by either method, you should have pscf/ directory structure::

    pscf/
       git/
       examples/

in which the pscf/git/ subdirectory contains the source code and the 
pscf/examples/ directory contains the examples.

**Directory structure**

The resulting directory has several subdirectories containing different
classes of examples involving diblock copolymer melts, triblock copolymer
melts, and mixtures of diblock copolymer and solvent. Lower level directories
that contain individual examples each have a file named param (the input
parameter file), a file named in.omega (the input chemical potential file),
and an initially empty directory named "out/" in which output files will
be placed. The instructions for running pscf given in the next section
use the same conventions for input file names. 

