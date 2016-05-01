
.. _theory-page:

***************************************
Appendix: Self-Consistent Field Theory
***************************************

PSCF solves self-consistent field theory (SCFT) equations for an 
incompressible mixture of any number of linear block copolymer species and 
point-like solvent molecular species. 

SCFT for a liquid of flexible polymers is based on a mean-field 
approximation that allows us to predict properties of an interacting 
liquid by considering considering the behavior of a corresponding gas of 
noninteracting molecules in a spatially inhomogeneous chemical potential 
landscape. In what follows, let :math:`k_{B}T\omega_{\alpha}(\textbf{r})` 
denote the chemical potential field for monomers of type :math:`\alpha`, 
which gives the free energy cost of placing such a monomer at location
:math:`\textbf{r}`. Let :math:`\rho_{\alpha}(\textbf{r})` denote the
corresponding average volume fraction of monomers of type :math:`\alpha`
at position :math:`\textbf{r}`.

PSCF implements a version of SCFT for incompressible liquids in which
each species of molecule occupies a well-defined volume, independent 
of composition and (modest) changes in pressure. In what follows, 
:math:`N_{i}` is dimensionless measure of the volume per molecule (or size)
of species :math:`i`, which is defined for both polymeric and point-like 
molecules as the ratio of the volume occupied by one molecule of that 
species to an monomer reference volume. The length of each block 
within a block copolymer is specified similarly, as a ratio of the
block volume to a monomer reference volume. 

In what follows, let :math:`C` be the number of distinct monomer or solvent 
types in the system. Let :math:`P` denote the number of polymer species and 
:math:`S` be the number of solvent species.  Here, we use a convention in 
which integer indices :math:`\alpha, \beta = 1, \ldots, C` indicate monomer 
types, and indices :math:`i, j` denote molecular species.  Species indices 
are ordered with all polymeric species listed first, so that species index 
values in the range :math:`i, j = 1, \ldots, P` denote polymeric species, 
and values in the range :math:`P+1,\ldots, P+S` denote solvent species.

**Polymer Species**

Each polymeric species with species index :math:`i` and overall 
chain length :math:`N_{i}` is treated as a random walk characterized by 
a contour :math:`\textbf{R}(s)`, where :math:`s` is a contour variable
with a range :math:`0 \leq s \leq N_{i}`. For each such species, we
define a pair of functions :math:`q_{i}(\textbf{r}, s)` and 
:math:`q^{\dagger}_{i}(\textbf{r}, s)`. The functions :math:`q_{i}` and
:math:`q^{\dagger}` are normalized partition functions for chain segments 
corresponding to contour variable domains :math:`[0,s]` and :math:`[s,N_{i}]`, 
respectively, when the monomer at contour position :math:`s` is constrained 
to position :math:`\textbf{R}(s) = \textbf{r}`. These functions obey a
pair of modified diffusion equations

.. math::

  \frac{\partial q_{i}}{\partial s} =  -H_{\alpha(s)}q_{i} 

  \frac{\partial q_{i}^{\dagger}}{\partial s} = +H_{\alpha(s)}q_{i}^{\dagger}

in which :math:`\alpha(s)` is the monomer type of the block containing 
monomer :math:`s` of polymer species :math:`i`, and in which :math:`H_{\alpha}` 
is a linear diferential operator

.. math::

  H_{\alpha} = -\frac{b_{\alpha}^{2}}{6}\nabla^{2} 
             + \omega_{\alpha}(\textbf{r})

in which :math:`b_{\alpha}` is the statistical segment length for monomers of
type :math:`\alpha`, and :math:`\omega_{\alpha}` is the corresponding chemical
potential field. These equations must be solved for :math:`0 < s < N_{i}` 
subject to an initial condition

.. math::

   q_{i}(\textbf{r},s=0) = q^{\dagger}_{i}(\textbf{r},s=N) = 1

for all :math:`\textbf{r}`, and boundary condition requiring that :math:`q` 
and :math:`q^{\dagger}` be periodic functions of :math:`\textbf{r}` with 
the periodicity of some specified Bravais lattice. 

The quantity :math:`Q_{i}` is a normalized overall partition function
for chains of species :math:`i`, given by an integral

.. math::

   Q_{i} = \frac{1}{V}\int \! d\textbf{r} \; q(\textbf{r},s=N)

in which the integral is taken over one unit cell of a periodic structure,
and :math:`V` denotes the generalized volume per unit cell of a structure
that periodic in :math:`D` dimensions (i.e., the volume per unit cell for 
a 3D crystal, area per 2D unit cell for a 2D structure such as the hexagonal 
cylinder phase, and the length per unit cell in a 1D lamellar phase).

The probability of finding a specific monomer :math:`s` of species 
:math:`i` at position :math:`\textbf{r}` is proportional to the product 
:math:`q_{i}(\textbf{r},s) q^{\dagger}_{i}(\textbf{r},s)` of the constrained
partition functions for the two chain segments that meet at monomer s.
Let :math:`\rho_{\alpha}^{(i)}` denote the contribution to the local
volume fraction of :math:`\alpha` monomers from monomers of a polymer 
species :math:`i` that contains at least one block of monomer type
:math:`\alpha`. This quantity is given by a product

.. math::

   \rho_{\alpha}^{(i)}(\textbf{r}) =  
   \frac{\overline{\phi}_{i}}{N_{i}Q_{i}}
   \int\limits_{\alpha(s)=\alpha} \! ds \;
   q(\textbf{r},s) q^{\dagger}(\textbf{r},s)

in which :math:`\overline{\phi}_{i}` is the average overall volume
fraction of molecule species :math:`i` within the mixture, and the
integral with respect to :math:`s` is taken only over blocks of 
monomer type :math:`\alpha`.

**Solvent Species**

Each solvent species :math:`i` is associated with a specific monomer type 
:math:`\alpha` and a volume :math:`i`. A "monomer" type that is assigned 
to a solvent species may or may not also be contained within one or more 
of the polymeric species. 

In the single molecule problem for solvent species, the free energy penalty 
for a solvent molecule of monomer type :math:`\alpha` to be located at position 
:math:`\textbf{r}` is given by :math:`k_{B}T N_{i}\omega_{\alpha}(\textbf{r})`.
This yields a solvent concentration 
:math:`\rho_{i}(\textbf{r}) \propto \exp(-N_{i}\omega_{\alpha}(\textbf{r}))`. 

The normalized overall partition for such a point-like species is given by 
an integral 

.. math::

   Q_{i} = \frac{1}{V}\int \! d\textbf{r} \; \exp(-N_{i}\omega_{\alpha}(\textbf{r}))

The contribution of solvent species :math:`i` of type :math:`\alpha` to 
the local volume fraction of :math:`\alpha` is given by a ratio

.. math::

   \rho_{\alpha}^{(i)}(\textbf{r}) = 
   \frac{\overline{\phi}_{i}}{Q_{i}} 
   \exp(-N_{i}\omega_{\alpha}(\textbf{r}))

in which :math:`\overline{\phi}_{i}` is the overall volume fraction of
species :math:`i` within the mixture. 

The total volume fraction :math:`\rho_{\alpha}(\textbf{r})` for each
monomer type :math:`\alpha` is simply given by the sum of contributions
from all polymeric species that contain a block or blocks of type 
:math:`\alpha` plus the contribution of any solvent of type :math:`\alpha`.

**Self-Consistent Field Equations**

The monomer chemical potential fields are given, within the standard 
approximation for excess free energies in terms of binary Flory-Huggins 
interaction parameters, as functions

.. math::

   \omega_{\alpha}(\textbf{r}) = \sum_{\beta = 1}^{C}
   \chi_{\alpha\beta} \rho_{\beta}(\textbf{r}) + \xi(\textbf{r})

in which :math:`\chi_{\alpha\beta}` is a binary interaction parameter for
interactions between monomers of types :math:`\alpha` and :math:`\beta`,
and :math:`\xi(\textbf{r})` is a Lagrange multiplier pressure 
field.  The interaction parameters in PSCF satisfy
:math:`\chi_{\alpha\alpha}=0` and (obviously)
:math:`\chi_{\alpha\beta} = \chi_{\beta\alpha}`.

The field :math:`\xi(\textbf{r})` must be chosen such that the monomer 
concentrations satisfy the incompressibility constraint

.. math::

   1 = \sum_{\alpha=1}^{C} \rho_{\alpha}(\textbf{r})

**Thermodynamic Properties**

The Helmholtz free energy :math:`f` per monomer reference volume, as given
in the output file, is given by a sum

.. math::

    \frac{f}{k_{B}T} 
    & = 
    \sum_{i=1}^{P+S} \frac{\overline{\phi}_{i}}{N_{i}} 
    \left [ \ln ( \overline{\phi}_{i} / Q_{i}) - 1 \right ] \\
    & -  \frac{1}{V}
          \sum_{\alpha=1}^{C} 
          \int \! d\textbf{r} \; 
          \omega_{\alpha}(\textbf{r})
          \rho_{\alpha}(\textbf{r}) \\
    & +  \frac{1}{2V} 
          \sum_{\alpha, \beta =1}^{C}  \chi_{\alpha\beta}
          \int \! d\textbf{r} \; 
          \rho_{\alpha}(\textbf{r})
          \rho_{\beta}(\textbf{r})

Note that the sum over species in the first line is a sum over all species,
including polymeric and solvent species, with different ways of defining 
:math:`Q_{i}` for different types of molecule.

The corresponding chemical potential :math:`\mu_{i}` for species :math:`i` 
is given by

.. math::

    \frac{\mu_{i}}{k_{B}T} = \ln(\overline{\phi}_{i}/Q_{i})

The value given in the output file is :math:`\mu_{i}/k_{B}T`.

The macroscopic physical pressure :math:`P` is computed from the identity

.. math::

    P = - \frac{f}{v} + \sum_{i=1}\frac{\mu_{i}\overline{\phi}_{i}}{N_{i}v} 
      
in which :math:`v` is the monomer reference volume and :math:`f` is 
the Helmholtz free energy per reference volume. Note that :math:`f/v`
is the Helmholtz free energy per volume and 
:math:`\overline{\phi}_{i}/(N_{i}v)` is the average number of 
molecules of species :math:`i` per unit volume. The value given in the 
output file is the dimensionless value :math:`Pv/k_{B}T`.

**Ensembles**

PSCF can be carry out calculations using either canonical ensemble 
or grand-canonical ensemble. 

In canonical ensemble a value of the overall volume fraction 
:math:`\overline{\phi}_{i}` must be given for each species in 
the input parameter file, and values of chemical potential are 
computed from the solution.

In grand canonical ensemble, a value of the normalized chemical 
potential :math:`\mu_{i}/k_{B}T` must be given for each species in
the input parameter file, and average volume fractions for each
species are computed.

In grand-canonical ensemble, values for the Lagrange multplier
field :math:`\xi(\textbf{r})` and the macroscopic pressure :math:`P`
are uniquely determined by the values for the chemical potentials. 

In canonical ensemble, the value of the Lagrange multplier field 
:math:`\xi(\textbf{r})` is defined only to within a arbitrary
spatially homogeneous constant. As a result, the chemical potentials
and the macroscopic pressure :math:`P` are also undefined in this
ensemble, unless an additional constraint is imposed. PSCF resolves 
this ambiguity by requiring, as a matter of convention, that the 
spatial average of :math:`\xi(\textbf{r})` vanish. In this ensemble,
PSCF also outputs values for the pressure, chemical potentials, and 
:math:`\omega` fields that are all consistent with this convention
for the average value of :math:`\xi`. Values for the Hemholtz free
energy density of an incompressible liquid can, however, be shown 
to be independent of changes in the value of :math:`\xi` by a 
homogeneous constant, and are thus independent of this choice of
convention.

