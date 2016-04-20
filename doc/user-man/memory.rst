

===========================
Appendix: Memory Counting
===========================

The majority of memory in a typical SCF calculation is consumed by the 
storage of mono-constrained partition functions (forward and backward) 
and, if solved using the quasi Newton-Raphson, by the storage of the
inverse Jacobian matrix needed in Broyden's method.

**Chain Partition Functions**

The memory (in bytes) needed to store the two mono-constrained partition 
function is determined by the following factors:

  ===============================   ==============================
  Factor                            Value
  ===============================   ==============================
  number of grid points             ngrid(1)*ngrid(2)*ngrid(3)        
  number of contour segments        2*chain_length/chain_step         
  bytes per real (double) number    8                              
  ===============================   ==============================

The product of these numbers yield:

   mem(Q) = (ngrid(1)*ngrid(2)*ngrid(3)) * (chain_length/chain_step) * 16


**Jacobian Matrix**

The dimension of Jacobian is "N_monomer * N_basis" as long as the N_cell_param
(<= 6) is ignored, thus the memory requirement, in bytes, is:

   mem(J) = (N_monomer*N_basis)^2 * 8

The number of basis function (N_basis) may be roughfully estimated, if 
needed, to be the number of the grid points divided by the number of 
space group symmetry elements.


**Numerical Example**

Imagine simulating a gyroid phase for an ABC triblock copolymer (N_monomer=3)
using normalized chain_length=1, chain_step=0.01, a grid with dimensions
ngrid(1)=ngrid(2)=ngrid(3)=56 The I a -3 d space group of the gyroid phase 
has 96 symmetry elements, so our simple estimate suggests 
N_basis = 56*56*56/96=1829.3. In fact, we obtain a slightly higher value of 
N_basis=1856. mem(Q) and mem(J) are found to be 0.28GB and 0.25GB, 
respectively.

