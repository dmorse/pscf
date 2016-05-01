
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

A copy of this repository may be obtained (as for the source code 
repository) either by downloading a zip file from the repository web 
page or using git to clone the repository. The instruction given 
below for both methods are very similar to those given for obtaining 
the source code:
 
To download a zip file:

    * Point your browser at the pscf-examples github repository.

    * Click the "Download ZIP" button near the upper right corner 
      of that web page. 

    * Move the pscf-examples-master/ directory into your pscf/ directory,
      which should also contain a subdirectory that holds the repository 
      containing the source code.

    * Rename the resulting directory pscf/pscf-examples-master/ to 
      pscf/examples.

To clone the repository of examples, if git is installed on your machine:

    * Change directory to the parent pscf/ directory.

    * Clone the pscf-examples repository, by entering::

          git clone https://github.com/dmorse/pscf-examples.git

    * Change the name (i.e., mv) the pscf/pscf-examples directory to 
      pscf/examples/

Either method should give you a pscf/ directory with the structure::

    pscf/
       git/
       examples/

in which the pscf/git/ subdirectory contains the source code and the 
pscf/examples/ directory contains the examples.

The structure and contents of the examples directory is explained in 
the README file within that directory.

