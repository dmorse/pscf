
************
Introduction
************

PSCF is a numerical implementation of the Edwards-Helfand 
self-consistent field theory (SCFT) for periodic phases formed by 
liquids containing block copolymers. It is designed to describe 
spatially periodic structures of incompessible liquids that may 
contain a mixture of block copolymers, homopolymers, and small 
molecule solvents. It can describe structures with 1, 2, 3 or
dimensional periodicity with any type of unit cell and any
specified space group symmetry. The modified diffusion equation 
(MDE) is solved with an efficient pseudo-spectral method. The 
code was designed primarily for calculating phase diagrams of
block copolymer melts by comparing free energies of competing 
ordered phases, and contains several features that facilitate 
such calculations. 

Features
========

*  Arbitrary mixtures of block copolymers, homopolymers, and solvents 
*  Ordered phases with 1, 2, or 3 dimensional periodicity
*  Canonical and grand-canonical ensembles
*  Arbitrary 2 or 3D unit cell, with non-orthogonal axes
*  Imposition of any specified space-group symmetry
*  Hard-coded "database" of all 230 3D space groups and 17 2D plane groups
*  Variable unit cell algorithm, in which unit cell parameters adjust to
   minimize free energy
*  Continuation of solutions along lines in parameter space
*  Linear response calculations for periodic structures 

Documentation
=============

*  User manual (this file)
*  Developer/API manual, generated from comments in the source code
*  A library of examples, in a separate git repository located at

   http://github.com/dmorse/pscf-examples.git

The developer manual is available online via a link from the program
home page, and can also can be regenerated from the source code (see 
the file doc/README for instructions).

Dependencies
============
 
The program is written in Fortran 90. It depends the open source FFTW Fast 
Fourier Transform library and the LAPACK linear algebra library. To compile 
the code from source, you will thus need a Fortran 90 compiler and these 
two libraries.

Contributors
============

* David Morse
* Chris Tyler
* Jian Qin
* Amit Ranjan
* Raghuram Thiagarajan
* Akash Arora

