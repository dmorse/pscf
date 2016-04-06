<h2>Appendix: Memory Counting</h2>
<p>
The majority of memory in a typical SCF calculation is consumed by the storage
of mono-constrained partition functions (forward and backward) and by the storage
of te Jacobian matrix needed in the Broyden's iteration.
</p>

<p>
<b>Mono-constrained Partition Function:</b>
</p>
<p>
The memory (in bytes) needed to store the mono-constrained partition function
is determined by the following terms:
<table>
  <tr>
    <td>ngrid(1)*ngrid(2)*ngrid(3) </td>
    <td>number of grid points</td>
  </tr>
  <tr>
    <td>2*chain_length/chain_step  </td>
    <td>number of contour segments</td>
  </tr>
  <tr>
    <td>8                          </td>
    <td>number of bytes for a double precision variable</td>
  </tr>
</table>
The product of these numbers yield:
<pre>
   mem(Q) = (ngrid(1)*ngrid(2)*ngrid(3)) * (chain_length/chain_step) * 16
</pre>
</p>

<p>
<b>Jacobian Matrix:</b>
</p>
The dimension of Jacobian is "N_monomer * N_basis" as long as the N_cell_param
(<= 6) is ignored, thus the memory cost:
<pre>
   mem(J) = (N_monomer*N_basis)^2 * 8
</pre>
The number of basis function (N_basis) may be roughfully estimated to be the
number of the grid points divided by the number of space group symmetry
elements.

<p>
<b>Numerical Example:</b>
</p>
Imagine simulating a Gyroid phase for an ABC triblock copolymer (N_monomer=3)
using normalized chain_length=1, chain_step=0.01,
ngrid(1)=ngrid(2)=ngrid(3)=56, and N_basis=1856. mem(Q) and mem(J) are found to
be 0.28GB and 0.25GB, respectively.

