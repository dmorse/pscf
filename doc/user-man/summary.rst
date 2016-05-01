
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
       2.5000000000E-01    7.5000000000E-01
   
   COMPOSITION         
   ensemble            
                      0
   phi_chain           
       1.0000000000E+00
   
   INTERACTION         
   interaction_type    
                  'chi'
   chi                 
       2.0000000000E+01
   
   UNIT_CELL           
   dim                 
                      3
   crystal_system      
                'cubic'
   N_cell_param        
                      1
   cell_param          
       1.9254998725E+00
   
   DISCRETIZATION      
   ngrid               
                     10                  10                  10
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
                     28
   
   FINISH              
   
   THERMO              
   f_Helmholtz         
       2.6137508664E+00
   f_homo              
       2.7500000000E+00
   pressure            
       3.5507319562E+00
   mu_chain            
       6.1644828226E+00
   stress              
      -3.4708983737E-12
   
   STATISTICS          
   N_star              
                     28
   Final Error         
       3.4708983737E-10
   Iterations          
                      3
   Basis Time          
       2.0689300000E-01
   SCF Time            
       4.1129500000E+00


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
   stress          derivatives of free energy per monomer / kT
   =============== ====================================================

Note that the Helmholtz free energy is output per monomer, normalized
by kT. In the simple case of a single component block copolymer melt, 
the free energy per chain is then given by the product of f_Helmholtz 
and the overall chain length (sum of the block lengths given in the 
parameter file).  Similarly, the reported pressure is a dimensionless 
value obtained by multiplying the pressure by the monomer volume and
then dividing by kT. Values of mu_chain and (if present) mu_solvent
are instead free energies per molecule, normalized by kT.

In a system with more than one component, a value would be given for
the chemical potential of each species, with one value per line. The
mu_solvent array appears only if there is one or more solvent species.

Please refer to the :ref:`theory-page` for the precise mathematical 
expressions used to obtain these quantities.

.. _summary-statistics-sec:

**STATISTICS Section**

The statistics section contains information about the size N_star of 
the basis used to approximate the rho and omega fields, the number of 
iterations required, the final error, and the amount of time taken 
for different parts of the computation. Times are given in seconds.
