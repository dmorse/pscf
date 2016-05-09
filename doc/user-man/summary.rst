
.. _summary-page:

*******************
Output Summary File
*******************

When the self-consistent field problem is solved for a specific set 
of parameters, in response to an ITERATE command in the parameter file, 
PSCF outputs an output summary file named "out", in addition to files
containing the final omega and rho fields. The format of this summary
file is similar to that of the parameter file, except that it 
contains some additional information about thermodynamic properties 
and computational details near the end of the file.

**Example**

Here is an example of an output summary file produced by a simulation of
the BCC of a diblock copolymer melt. The majority of the file, up until
the FINISH line, has the same format as a parameter file, so that this
part of the file can be used as the basis of a parameter file for a new
simulation. The THERMO and STATISTICS sections appear after the FINISH 
line, and contain additional information about, respectively, thermodynamic 
properties of the converged solution and about computational aspeces of
the simulation.

::

    format  1  0
   
   MONOMERS            
   N_monomer           
                      2
   kuhn                
       1.0000000000E+00    1.0000000000E+00
   
   CHAINS              
   N_chain             
                      1
   N_block             
                      2
   block_monomer       
                      1                   2
   block_length        
       2.0000000000E-01    8.0000000000E-01
   
   COMPOSITION         
   ensemble            
                      0
   phi_chain           
       1.0000000000E+00
   
   INTERACTION         
   interaction_type    
                  'chi'
   chi                 
       2.4000000000E+01
   
   UNIT_CELL           
   dim                 
                      3
   crystal_system      
                'cubic'
   N_cell_param        
                      1
   cell_param          
       1.8447588337E+00
   
   DISCRETIZATION      
   ngrid               
                     24                  24                  24
   chain_step          
       1.0000000000E-02
   
   BASIS               
   group_name          
             'I m -3 m'
   
   ITERATE             
   input_filename      
             'in.omega'
   output_prefix       
                 'out/'
   max_itr             
                     20
   error_max           
       1.0000000000E-08
   domain              
                      T
   itr_algo            
                   'NR'
   N_cut               
                     95
   
   FINISH              
   
   THERMO              
   f_Helmholtz         
       2.7451703490E+00
   f_homo              
       2.8400000000E+00
   pressure            
       3.7365224085E+00
   mu_chain            
       6.4816927575E+00
   stress              
      -6.0758291562E-11
   
   STATISTICS          
   N_star              
                    231
   final_error         
       6.0758291562E-09
   iterations          
                     10
   basis_time          
       5.6800000000E+00
   scf_time            
       1.3674400000E+02


.. _summary-thermodynamics-sec:

**THERMO Section**

The thermodynamics section contains final values for the following
thermodynamic properties of the converged solution.  All quantities
that involve energy (i.e., f_Helmholtz and mu) are output in units 
in which thermal energy is set to kT=1. The pressure is output in 
units in which kT=1 and in which the monomer reference volume is =1. 

 
   =============== ====================================================
   Variable        Description
   =============== ====================================================
   f_Helmholtz     Helmholtz free energy per monomer / kT
   f_homo          f_Helhmoltz of a hypothetical homogeneous mixture
   pressure        Macroscopic pressure x monomer volume / kT
   mu_chain        chemical potential / kT, for each polymer species
   mu_solvent      chemical potential / kT, for each solvent species
   phi_chain       total volume fraction for each polymer species 
   phi_solvent     total volume fraction for each solvent species
   stress          derivatives of free energy per monomer / kT
   =============== ====================================================

Values of mu_chain and mu_solvent appear in this section if and only if 
the calculation was carried out in canonical ensemble (ensemble == 0), 
in which case the corresponding species volume fractions are given as 
input parameters in the COMPOSITION section. Conversely, values of
phi_chain and phi_solvent vectors appear in this output section only 
if the calculation was carried out in grand-canonical ensemble 
(ensemble == 1), in which case the corresponding chemical potential
values are given as inputs in the COMPOSITION section.

**Units**: The Helmholtz free energy f_Helmholtz is given in this
section is a dimensionless free energy per monomer, normalized by kT. 
In the simple case of a single component block copolymer melt, the 
free energy per chain is then given by the product of f_Helmholtz 
and the overall chain length (sum of the block lengths given in the 
parameter file).  Similarly, the reported pressure is a dimensionless 
value obtained by multiplying the pressure by the monomer reference
volume and then dividing by kT. Values of mu_chain and (if present) 
mu_solvent are instead free energies per molecule, normalized by kT.

In a system with more than one chain and/or solvent component, a value 
would be given for the chemical potential of each species, with one 
value per line. The mu_solvent (in canonical ensemble) or phi_solvent
(in grand-canonical ensemble) array appears only if the parameter file
has a SOLVENTS input section with one or more solvent species. 

Each element of the array of values of "stress" is the derivative of
the dimensionless free energy per monomer f_Helmholtz with respect to
one of the unit cell parameters. The number of elements is thus equal
to N_cell_param, the number of parameters required to describe the
unit cell. All components of this array should be very close to zero 
at the end of a computation with a flexible unit cell (domain == T). 

Please refer to the :ref:`theory-page` for documentaiton of the 
mathematical expressions used to compute free energies and other
physical properties.

.. _summary-statistics-sec:

**STATISTICS Section**

The statistics section contains information about the size N_star of 
the basis used to approximate the rho and omega fields, the number of 
iterations required, the final error, and the amount of time taken 
for different parts of the computation. Times are given in seconds.
