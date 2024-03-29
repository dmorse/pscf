# PSCF - Polymer Self-Consistent Field Theory (Fortran)

Copyright (2002-2017) Regents of the University of Minnesota

PSCF is a program for numerically solving the polymer self-consistent 
field theory (SCFT) for spatially periodic structures formed by block 
copolymer melts and mixtures of block copolymers with linear 
homopolymers and/or small molecule solvents.

## History and Versions

The version of PSCF provided here is a legacy version that was written 
in Fortran 90. The current version of PSCF is a rewritten C++/CUDA 
package that is available in a separate github repository at 
https://github.com/dmorse/pscfpp

Almost all of the features of this legacy Fortran version have been 
ported to the current C++/CUDA version of PSCF, which also provides
several important new features. The home page for PSCF, which provides 
a discussion of the differences between the two versions and links 
to both is available at https://pscf-home.cems.umn.edu

The Fortran version provided in this repository is no longer under 
development, and is not being patched to correct any newly reported 
bugs. New users are thus strongly encouraged to use the current 
C++/Cuda version.  Users of this legacy Fortran version are encouraged 
to switch to the current version when they can. 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation. A copy of this license is included in
the LICENSE file in the top-level PSCF directory.

## Dependencies

This version of PSCF depends upon the FFTW fast Fourier transform library 
and the LAPACK linear algebra library.  These packages must be installed
before attempting to compile the program from source.

## User Documentation

A web user manual is available at https://pscf.readthedocs.io

Instructions for compiling the program from source, as well how to
install a precompiled executable on some systems, are given in the 
user maual.

The source files for the user manual are text files that are stored in
the doc/user-man directory. The relevant files have file extension .rst.

We recommend that users also refer to the 2016 reference article from the 
journal Macromolecules, cited below. This article provides a brief 
introduction to SCFT, discusses some of the algorithms and file formats 
used in the PSCF program, and discusses some issues and questions faced 
by first time users.

## Developer Documentation

A local copy of the developer API documentation, which includes
documentation for all modules, subroutines and public variables, may
be regenerated by following instructions given in the file doc/README.
The resulting html pages are installed in the doc/devel-man directory.

## Reference Article 

If you use either version of PSCF in published work, we request that you 
cite the paper:

"Broadly Accessible Self-Consistent Field Theory for Block Copolymer
Materials Discovery", A. Arora, J. Qin, D.C. Morse, K.T. Delaney,
G.H. Fredrickson, F.S. Bates and K.D. Dorfman, 
*Macromolecules* **49**, 4675 (2016)

available electronically at 
http://pubs.acs.org/doi/10.1021/acs.macromol.6b00107

## Examples

A library of examples is provided in a separate github repository, 
located at https://github.com/dmorse/pscf-examples

## Links

  - Home Page:    https://pscf-home.cems.umn.edu/
  - Source Code:  https://github.com/dmorse/pscf
  - CI Server:    https://travis-ci.org/dmorse/pscf
  - User Manual:  https://pscf.readthedocs.io
  - API Manual:   https://dmorse.github.io/pscf-dev-man/toc.html

[buildstatus_image_travis]: https://travis-ci.org/dmorse/pscf.svg?branch=master
[travisci]: https://travis-ci.org/dmorse/pscf

[![Travis][buildstatus_image_travis]][travisci]

## Repository Directory Structure

    src/                          - Fortran 90 source files
    doc/                          - documentation files
    doc/user-man                  - - user manual source files
    tools/                        - Tools for processing output and source
    tools/matlab                  - - matlab scripts for visualization
    tools/python                  - - python modules
    make/                         - build directory for make

An annotated list of source files is given in the file src/SRC_FILES.
Before modifying any fortran files, see the note at the end of that
file regarding the use of preprocessor to generate some files.

## Contributors

- David Morse
- Chris Tyler
- Jian Qin
- Amit Ranjan
- Raghuram Thiagarajan
- Akash Arora

## Support 

Development of this Fortran version of PSCF has been supported by 
National Science Foundation grants NSF-DMR097338 and NSF-DMR130436.
