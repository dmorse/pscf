
.. _theory-page:

*********************
Theoretical Formalism
*********************

PSCF solves self-consistent field theory (SCFT) equations for an 
incompressible mixture of any number of linear block copolymer species and 
point-like solvent molecular species. 

Let :math:`C` be the number of distinct monomer or solvent types in the 
system. Let :math:`P` denote the number of polymer species and :math:`S` 
be the number of solvent species.  In the remainder of this discussion, 
se use a convention in which integer variables 
:math:`\alpha, \beta = 1, \ldots, C` are consistently used to indicate
monomer types, and variables :math:`i, j` to denote molecular species.
We assume that the species are ordered with polymer species before any
solvent species, so that species index values in the range 
:math:`i, j = 1, \ldots, P` denote polymeric species, and values in the
range :math:`P+1,\ldots, :math:`P+S` indicate point-like solvent species.

SCFT for a liquid of flexible polymers is based on a mean-field 
approximation that allows us to predict properties of an interacting 
liquid by considering considering the behavior of a corresponding gas of 
noninteracting molecules in a spatially inhomogeneous chemical potential 
landscape. In what follows, let :math:`\omega_{\alpha}(\textbf{r})` 
denote the chemical potential field for monomers of type :math:`\alpha`, 
which gives the free energy cost of placing such a monomer at location
:math:`\textbf{r}`. Let :math:`\rho_{\alpha}(\textbf{r})` denote the
corresponding average volume fraction of monomers of type :math:`\alpha`.

PSCF implements a version of SCFT for incompressible liquids in which
each species of molecule occupies a well-defined volume, independent 
of composition and (modest) changes in pressure. In what follows, 
:math:`N_{i}` is dimensionless measure of the volume per molecule (or size)
of species :math:`i`, which is defined for both polymeric and point-like 
molecules as the ratio of the volume occupied by one molecule of that 
species to an monomer reference volume. The length of each block 
within a block copolymer is specified similarly, as a ratio of the
block volume to a monomer reference volume. 

**Polymer Species**

For each polymeric species :math:`i` of overall chain lengtt :math:`N_{i}`, we 
define a pair of constrained partition functions :math:`q_{i}(\textbf{r}, s)` 
and :math:`q^{\dagger}_{i}(\textbf{r}, s)`. These are normalized partition 
functions for segments of chain containing monomers :math:`[0,s]` and 
:math:`[s,N_{i}]`, respectively, when the monomer at contour position 
:math:`s` is constrained to position :math:`\textbf{r}`. These obey a
pair of modified diffusion equations

.. math:

  \begin{eqnarray}
  \frac{\partial q_{i} (\textbf{r},s)}{\partial s} 
  & = & -H_{\alpha(s)}q_{i} (\textbf{r},s) 
  \\
  \frac{\partial q_{i}^{\dagger}(\textbf{r},s)}{\partial s} 
  & = & +H_{\alpha(s)}q_{i}(\textbf{r},s) 
  \end{eqnarray}

in which :math:`\alpha(s)` is the monomer type of the block containing 
monomer :math:`s` of polymer species :math:`i`, and in which :math:`H_{\alpha}` 
is a linear diferential operator

.. math:

  H_{\alpha} = -\frac{b_{\alpha}^{2}}{6}\nabla^{2} + \omega_{\alpha}  

in which :math:`b_{\alpha}` is the statistical segment length for monomers of
type :math:`\alpha`, and :math:`\omega_{\alpha}` is the corresponding chemical
potential field. These equations must be solved for :math:`0 < s < N_{i}` 
subject to an initial condition

.. math:

   q(\textbf{r},s=0) = q^{\dagger}(\textbf{r},s=N) = 1

for all :math:`\textbf{r}`, and boundary condition requiring that :math:`q` 
and :math:`q^{\dagger}` be periodic functions of :math:`\textbf{r}` with the 
periodicity of some specified Bravais lattice. 

The quantity :math:`Q_{i}` is a normalized overall partition function
for chains of species :math:`i`, given by an integral

.. math:

    Q_{i} = \frac{1}{V}\int d^{D}\textbf{r} q(\textbf{r},s=N)

in which the integral is taken over one unit cell of a periodid structure,
and :math:`V` denotes the generalized volume per unit cell of structure
with is periodic in :math:`D` dimensions (i.e., volume per unit cell for a 3D
crystal, area per unit cell for a 2D structure such as the hexagonal 
cylinder phase, and length per unit cell in a 1D lamellar phase).

The probability of finding a specific monomer :math:`s` of species 
:math:`i` at position :math:`\textbf{r}` is proportional to the product 
:math:`q(\textbf{r},s) q^{\dagger}(\textbf{r},s)` of the constrained
partition functions for the two chain segments that meet at that monomer.
Let :math:`\rho_{\alpha}^{(i)}` denote the contribution of a polymer
species :math:`i` that contains at least one block with monomers of 
type :math:`\alpha`. This quantity is given by a product

.. math:

   \rho_{\alpha}^{(i)}(\textbf{r}) =  
   \frac{\overline{\rho}_{i}}{N_{i}Q_{i}}
   \int\limits_{\alpha(s) = \alpha} 
   q(\textbf{r},s) q^{\dagger}(\textbf{r},s)

in which :math:`\overline{\rho}_{i}` is the average overall volume
fraction of molecular species :math:`i` within the mixture. 

**Solvent Species**

Each solvent species :math:`i` is associated with a specific monomer type 
:math:`\alpha` and a volume :math:`i`. A "monomer" type that is assigned 
to a solvent species may or may not also be contained within one or more 
of the polymeric species. In the single molecule problem, each solvent of
monomer type :math:`\alpha` is subjected to an effective energy penalty
:math:`k_{B}T N_{i}\omega_{\alpha}(\textbf{r})` at location 
:math:`\textbf{r}`, giving a concentration proportional to 
:math:`\exp(-N_{i}\omega_{\alpha}(\textbf{r})`. The contribution of 
solvent species 


