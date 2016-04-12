
************
Introduction
************

PSCF is a numerical implementation of the Edwards-Helfand 
self-consistent field theory (SCFT) for periodic phases formed by 
liquids containing block copolymers. It is designed to describe 
any periodic structures with 1, 2, or 3 dimensional periodicity of 
incompessible liquids that may contain a mixture of block copolymers, 
homopolymers, and small molecule solvents. It can describe periodic
structure structures with any unit cell (e.g., cubic, orthorhombic, 
triclinic, etc.) with any specified space group symmetry. The 
modified diffusion equation (MDE) is solved with an efficient 
pseudo-spectral method. The code was designed for calculating
phase diagrams by comparing free energies of competing ordered 
phases, and contains features that facilitate such calculations. 

In addition to its primary purpose of solving the SCF equations, 
PSCF can also calculate the linear susceptibility of ordered 
phases. This is an efficient pseudo-spectral implementation of 
the linear response theory of An-Chang Shi and coworkers. 

Features
========

*  Arbitrary mixtures of block copolymers, homopolymers, and solvents 
*  Ordered phases with 1, 2, or 3 dimensional periodicity
*  Canonical and grand-canonical ensembles
*  Arbitrary 2 or 3D unit cell, with non-orthogonal axes
*  Imposition of any specified space-group symmetry
*  Hard-coded "database" of all 230 3D space groups 
*  Variable unit cell algorithm, in which the unit cell parameters 
   are adjusted to minimize free energy. 
*  Efficient treatment of sequences of solutions in which input 
   parameters are changed by small amounts ("sweeps"), using continuation
   of solutions
*  RPA response and spinodal calculations for homogeneous phase 
*  Linear response calculations for periodic structures 


Documentation
=============

*  User manual (this file)
*  Developer manual, generated from comments in the source code
*  A library of examples, in a separate git repository

   http://github.com/dmorse/pscf-examples.git

See the file doc/README for instruction regarding how to regenerate the
developer manual.

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

