#!/usr/bin/env python2.4

from io import *
from version import *
import string 

class OutFile(object):
    '''
    An OutFile object contains the data in a *.out output
    file produced by F90 program pscf_pd. 
    '''
   
    def __init__(self,filename):
        '''
        PURPOSE
          Read and parse out file filename. 
          Create an OutFile object.
        ARGUMENT
          filename - string 
        COMMENT
          File is opened and closed within body of method
        '''
        self.file = open(filename,'r')
	self.io   = IO()
	file      = self.file
	self.N_chain   = 0
	self.N_solvent = 0

	# Dictionary flags indicates which sections are present
	# Keys are flag string values ('CHEMISTRY' etc.), values are all 1
	self.flags = {}

        # Read version line
	self.version = Version(self.file)

        # Read input script
        next = 1
        while next :
            flag = self.file.readline().strip()
	    if flag == '':
                continue

	    # Set key in self.flags dictionary
	    if flag == 'GRID':
                self.flags['DISCRETIZATION'] = 1
	    elif flag == 'SPECTRAL_BASIS':
                self.flags['BASIS'] = 1
            else:
                self.flags[flag] = 1

            if flag == 'CHEMISTRY':
                self.input_chemistry()
                self.flags['MONOMERS'] = 1
                self.flags['CHAINS'] = 1
                self.flags['SOLVENTS'] = 1
                self.flags['COMPOSITION'] = 1
                self.flags['INTERACTION'] = 1
            if flag == 'MONOMERS':
                self.input_monomers()
            if flag == 'CHAINS':
                self.input_chains()
            if flag == 'SOLVENTS':
                self.input_solvents()
            if flag == 'COMPOSITION':
                self.input_composition()
            if flag == 'INTERACTION':
                self.input_interaction()
            elif flag == 'UNIT_CELL':
                self.input_unit_cell()
            elif flag == 'GRID':
                self.N_grid = self.input_vec('int')
            elif flag == 'DISCRETIZATION':
                self.N_grid = self.input_vec('int')
                self.chain_step = self.input_var('real')
            elif flag == 'FILE_PREFIXES':
                self.input_prefix  = self.input_var('char')
                self.output_prefix = self.input_var('char')
            elif flag == 'BASIS':
                self.group_name = self.input_var('char')
            elif flag == 'SPECTRAL_BASIS':
                self.group_name = self.input_var('char')
                self.Gabs_max = self.input_var('real',comment='Gabs_max',f='N')
                self.output_prefix = self.input_var('char',comment='output_prefix',f='N')
		self.flags['FILE_PREFIXES'] = 1
            elif flag == 'READ_OMEGA':
		if not self.__dict__.has_key('input_prefix'):
                    self.input_prefix = self.input_var('char')
		    self.flags['FILE_PREFIXES'] = 1
                self.reorder_omega = self.input_var('logic')
            elif flag == 'ITERATE':
                self.max_itr   = self.input_var('int')
		comment = self.file.readline().strip()
		if comment == 'ds' or comment == 'chain_step':
                    self.chain_step = self.input_var('real',f='N')
		    comment = self.file.readline().strip()
		if comment == 'extrapolation_order' or comment == 'extr_order':
                    self.extrapolation_order = self.input_var('int',f='N')
		    comment = self.file.readline().strip()
		if comment == 'N_cut':
                    self.N_cut = self.input_var('int',f='N')
		    comment = self.file.readline().strip()
                self.error_max = self.input_var('real',f='N')
                self.domain    = self.input_var('logic')
            elif flag == 'SWEEP':
                self.s_max     = self.input_var('real')
                self.input_increments()
            elif flag == 'FINISH':
                next = 0

        # Read output post-script, older format (version 0 9)
	version = self.version
	if version.eq(0,9):

	    # Loop over postscript lines
            next = 1
            while next:
                line = self.file.readline()
                if line:
                    if line.strip() == '':
                       continue
                    line = line.split()
                    name = line[0]
    		    if name == 'phi': 
                        name = 'phi_chain'
    		    if name == 'mu': 
                        name = 'mu_chain'
                    if name == 'mu_chain' and self.ensemble == 0 :
                        self.mu_chain = \
                           self.input_vec('real',self.N_chain,s='C',f='N')
                    if name == 'mu_solvent' and self.ensemble == 0 :
                        self.mu_solvent = \
                           self.input_vec('real',self.N_solvent,s='C',f='N')
                    if name == 'phi_chain' and self.ensemble == 1 :
                        self.phi_chain = \
                           self.input_vec('real',self.N_chain,s='C',f='N')
                    if name == 'phi_solvent' and self.ensemble == 1 :
                        self.phi_solvent = \
                           self.input_vec('real',self.N_solvent,s='C',f='N')
                    if name == 'N_star':
                        self.N_star = int(line[2])
                    if name == 'f_Helmholtz':
	                self.flags['THERMO'] = 1
                        self.f_Helmholtz = float(line[2])
                    if name == 'f_homo':
                        self.f_homo = float(line[2])
                    if name == 'pressure' :
                        self.pressure = float(line[2])
                    if name == 'stress' :
                        self.stress = [ float(x) for x in line[2:] ]
                    if name == 'Final':
                        self.final_error = float(line[3])
                    if name == 'Iterations':
	                self.flags['STATISTICS'] = 1
                        self.iteration = int(line[2])
                    if name == 'Basis':
                        self.basis_time = float(line[3])
                    if name == 'SCF':
                        self.scf_time = float(line[3])
                else:
                    next = 0

        # Read output post-script, newer format (version 0 9)
        else:

            # Read input script
            next = 1
            while next :
                line = self.file.readline()
		if not line:
                    break
                flag = self.file.readline().strip()
    	        if flag == '':
                    continue
                self.flags[flag] = 1

                if flag == 'THERMO':
                    self.f_Helmholtz = self.input_var('real')
                    self.f_homo      = self.input_var('real')
                    self.pressure    = self.input_var('real')
                    if self.ensemble == 0 :
                        if self.N_chain > 0 :
                            self.mu_chain = \
                               self.input_vec('real',self.N_chain,s='C')
                        if self.N_solvent > 0 :
                            self.mu_solvent = \
                               self.input_vec('real',self.N_solvent,s='C')
                    if self.ensemble == 1 :
                        if self.N_chain > 0 :
                            self.phi_chain = \
                               self.input_vec('real',self.N_chain,s='C')
                        if self.N_solvent > 0 :
                            self.phi_solvent = \
                               self.input_vec('real',self.N_solvent,s='C')
                    self.stress = self.input_vec('real',self.N_cell_param)

                if flag == 'DECOMPOSE':
                    self.overlap_AB  = self.input_var('real')
                    if self.N_monomer > 2 :
                        self.overlap_BC  = self.input_var('real')
                        self.overlap_CA  = self.input_var('real')
                    self.f_enthalpy  = self.input_var('real')
                    self.f_head      = self.input_var('real')
                    self.f_tail      = self.input_var('real')
                    self.f_excess    = self.input_var('real')

		elif flag == 'STATISTICS':
                    self.N_star      = self.input_var('int')
                    self.final_error = self.input_var('real')
                    self.iteration   = self.input_var('int')
                    self.basis_time  = self.input_var('real')
                    self.scf_time    = self.input_var('real')

        self.file.close()
	self.file = None

        # Make sorted list of attribute names
        self.att_names = self.__dict__.keys()
        self.att_names.sort()

    def write(self, file, major=1, minor=0):
	self.file = file
	self.version.major = major
	self.version.minor = minor
	self.version.write(file)
        if self.version.eq(0,9):
            if self.flags.has_key('CHEMISTRY'):
                file.write("\n%-15s\n" % 'CHEMISTRY')
                self.output_chemistry()
        else:
            if self.flags.has_key('MONOMERS'):
                file.write("\n%-15s\n" % 'MONOMERS')
                self.output_monomers()
            if self.flags.has_key('CHAINS'):
                file.write("\n%-15s\n" % 'CHAINS')
                self.output_chains()
            if self.flags.has_key('SOLVENTS'):
                file.write("\n%-15s\n" % 'SOLVENTS')
                self.output_solvents()
            if self.flags.has_key('COMPOSITION'):
                file.write("\n%-15s\n" % 'COMPOSITION')
                self.output_composition()
            if self.flags.has_key('INTERACTION'):
                file.write("\n%-15s\n" % 'INTERACTION')
                self.output_interaction()
        if self.flags.has_key('UNIT_CELL'):
            file.write("\n%-15s\n" % 'UNIT_CELL')
            self.output_unit_cell()
        if self.flags.has_key('DISCRETIZATION'):
            file.write("\n%-15s\n" % 'DISCRETIZATION')
            self.output_vec( 'int', 'N_grid', self.dim)
            self.output_var( 'real', 'chain_step')
            #self.output_var( 'int', 'extrapolation_order')
        if self.flags.has_key('FILE_PREFIXES'):
            file.write("\n%-15s\n" % 'FILE_PREFIXES')
            self.output_var('char', 'input_prefix')
            self.output_var('char', 'output_prefix')
        if self.flags.has_key('BASIS'):
            file.write("\n%-15s\n" % 'BASIS')
            self.output_var('char', 'group_name')
        #if self.flags.has_key('READ_OMEGA'):
        #    file.write("\n%-15s\n" % 'READ_OMEGA')
        #    output_var(file,'logic', self.reorder_omega, 'reorder_omega')
	if self.flags.has_key('ITERATE'):
            file.write("\n%-15s\n" % 'ITERATE')
            self.output_var( 'int', 'max_itr')
            self.output_var( 'int', 'N_cut')
            self.output_var( 'real', 'error_max')
            self.output_var( 'logic', 'domain')
        if self.flags.has_key('SWEEP'):
            file.write("\n%-15s\n" % 'SWEEP')
            self.output_var( 'real', 's_max')
            self.output_increments()
        file.write("\n%-15s\n" % 'FINISH')

        if self.version.eq(0,9):

           # Write phi or mu
           file.write("\n\n")
           if self.__dict__.has_key('mu_chain') and self.ensemble == 0 :
              self.output_vec('real','mu_chain',self.N_chain,s='C',f='A')
           if self.__dict__.has_key('mu_solvent') and self.ensemble == 0 :
              self.output_vec('real','mu_solvent',self.N_solvent,s='C',f='A')
           if self.__dict__.has_key('phi_chain') and self.ensemble == 1 :
               self.output_vec('real','phi_chain',self.N_chain,s='C',f='A')
           if self.__dict__.has_key('phi_solvent') and self.ensemble == 1 :
               self.output_vec('real','phi_solvent',self.N_solvent,s='C',f='A')
           file.write("\n" + '************************************' +"\n\n")
   
           # Write poscript
           if self.__dict__.has_key('N_star'):
               file.write('N_star      =  %20d' % self.N_star + "\n")
           if self.__dict__.has_key('f_Helmholtz'):
               file.write('f_Helmholtz =  %20.7E' % self.f_Helmholtz + "\n")
           if self.__dict__.has_key('f_homo'):
               file.write('f_homo      =  %20.7E' % self.f_homo + "\n")
           if self.__dict__.has_key('pressure'):
               file.write('pressure    =  %20.7E' % self.pressure + "\n")
           if self.__dict__.has_key('cell_param'):
               file.write('cell_param  =  ')
               for i in range(self.N_cell_param):
                   file.write('%20.7E' % self.cell_param[i])
               file.write("\n")
           if self.__dict__.has_key('stress'):
               file.write('stress      =  ')
               for i in range(self.N_cell_param):
                   file.write('%20.7E' % self.stress[i])
               file.write("\n")
           if self.__dict__.has_key('final_error'):
               file.write('Final Error =  %20.7E' % self.final_error + "\n")
           if self.__dict__.has_key('iteration'):
               file.write('Iterations  =  %20d' % self.iteration + "\n")
           if self.__dict__.has_key('basis_time'):
               file.write('Basis time  =  %20.7E' % self.basis_time + "\n")
           if self.__dict__.has_key('scf_time'):
               file.write('SCF time    =  %20.7E' % self.scf_time + "\n")

        else:

            if self.flags.has_key('THERMO'):
                file.write("\n%-15s\n" % 'THERMO')
                self.output_var('real','f_Helmholtz')
                self.output_var('real','f_homo')
                self.output_var('real','pressure')
                if self.ensemble == 0 :
                    if self.N_chain > 0 :
                        self.output_vec('real', \
			                'mu_chain',self.N_chain,s='C')
                    if self.N_solvent > 0 :
                        self.output_vec('real', \
			                'mu_solvent',self.N_solvent,s='C')
                elif self.ensemble == 1 :
                    if self.N_chain > 0 :
                        self.output_vec('real', \
			                'phi_chain',self.N_chain,s='C')
                    if self.N_solvent > 0 :
                        self.output_vec('real', \
			                'phi_solvent',self.N_solvent,s='C')
                self.output_vec('real','stress',self.N_cell_param)

            if self.flags.has_key('DECOMPOSE'):
                file.write("\n%-15s\n" % 'DECOMPOSE')
                self.output_var('real','overlap_AB')
                if self.N_monomer > 2 :
                    self.output_var('real','overlap_BC')
                    self.output_var('real','overlap_CA')
                self.output_var('real','f_enthalpy')
                self.output_var('real','f_head')
                self.output_var('real','f_tail')
                self.output_var('real','f_excess')

            if self.flags.has_key('STATISTICS') :
                file.write("\n%-15s\n" % 'STATISTICS')
                self.output_var('int','N_star')
                self.output_var('real','final_error')
                self.output_var('int','iteration')
                self.output_var('real','basis_time')
                self.output_var('real','scf_time')

	file.close()
	self.file = None


    def input_chemistry(self):
        ''' Analog of subroutine input_chemistry in chemistry_mod.f '''
	# Monomers 
        self.N_monomer = self.input_var('int')
        N_monomer      = self.N_monomer
        self.names     = self.input_vec('char', comment='names') # optional
        self.kuhn      = self.input_vec('real')

	#Interaction
        self.interaction_type  = self.input_var('char')
	print 'Reading interaction type =' + self.interaction_type
        if  self.interaction_type == 'B' or self.interaction_type == 'b':
            self.interaction_type = 'chi'
	    print 'Changing interaction type to chi'
        elif self.interaction_type == 'T' or self.interaction_type == 't':
            self.interaction_type = 'chi_T'
	    print 'Changing interaction type to chi_T'
        if  self.interaction_type == 'chi':
            self.chi = self.input_mat('real',N_monomer,N_monomer,s='L')
        elif self.interaction_type == 'chi_T':
            self.chi_A = self.input_mat('real',N_monomer,N_monomer,s='L')
            self.chi_B = self.input_mat('real',N_monomer,N_monomer,s='L')
            self.Temperature = self.input_var('real')

	# Chains 
        self.N_chain = self.input_var('int',f='A')
	if self.N_chain:
            self.N_block = self.input_vec('int',n=self.N_chain,s='C')
            self.file.readline()
            self.block_monomer = []
            for j in range(self.N_chain):
                self.block_monomer.append( self.input_vec('int',f='N') )
            self.file.readline()
            self.block_length = []
            for j in range(self.N_chain):
                self.block_length.append( self.input_vec('real',f='N') )
        else:
            self.N_chain = 0

	# Solvents (if any)
        self.N_solvent = self.input_var('int',comment='N_solvent', f='A')
	if self.N_solvent:
            self.solvent_monomer = self.input_vec('int', self.N_solvent, s='C')
            self.solvent_size    = self.input_vec('real',self.N_solvent, s='C')
        else:
            self.N_solvent = 0

	# Ensemble and composition
        self.ensemble = self.input_var('int',f='A')
        if self.ensemble == 0:
            if self.N_chain > 0:
                self.phi_chain = self.input_vec('real',n=self.N_chain,s='C',f='A') 
	    if self.N_solvent > 0:
                self.phi_solvent = self.input_vec('real',n=self.N_solvent,s='C',f='A') 
        elif self.ensemble == 1:
            if self.N_chain > 0:
                self.mu_chain = self.input_vec('real', n=self.N_chain,s='C',f='A') 
            if self.N_solvent > 0:
                self.mu_solvent = self.input_vec('real',n=self.N_solvent,s='C',f='A') 

    def input_monomers(self):
        ''' Analog of subroutine input_monomers in chemistry_mod.f '''
	# Monomers 
        self.N_monomer = self.input_var('int')
        N_monomer      = self.N_monomer
        self.kuhn      = self.input_vec('real')

    def input_interaction(self):
        ''' Analog of subroutine input_interaction in chemistry_mod.f '''
        self.interaction_type = self.input_var('char')
	N_monomer = self.N_monomer
        if self.interaction_type == 'chi':
            self.chi = self.input_mat('real',N_monomer,N_monomer,s='L')
        elif self.interaction_type == 'chi_T':
            self.chi_A = self.input_mat('real',N_monomer,N_monomer,s='L')
            self.chi_B = self.input_mat('real',N_monomer,N_monomer,s='L')
            self.Temperature = self.input_var('real')

    def input_chains(self):
        self.N_chain = self.input_var('int',f='A')
	if self.N_chain:
            self.N_block = self.input_vec('int',n=self.N_chain,s='C')
            self.file.readline()
            self.block_monomer = []
            for j in range(self.N_chain):
                self.block_monomer.append( self.input_vec('int',f='N') )
            self.file.readline()
            self.block_length = []
            for j in range(self.N_chain):
                self.block_length.append( self.input_vec('real',f='N') )
        else:
            self.N_chain = 0

    def input_solvents(self):
        self.N_solvent = self.input_var('int',f='A')
	if self.N_solvent:
            self.solvent_monomer = self.input_vec('int', self.N_solvent, s='C')
            self.solvent_size    = self.input_vec('real',self.N_solvent, s='C')
        else:
            self.N_solvent = 0

    def input_composition(self):
        self.ensemble = self.input_var('int',f='A')
	N_chain   = self.N_chain
	N_solvent = self.N_solvent
        if self.ensemble == 0:
            if self.N_chain > 0:
                self.phi_chain = self.input_vec('real',n=N_chain,s='C',f='A') 
	    if self.N_solvent > 0:
                self.phi_solvent = self.input_vec('real',n=N_solvent,s='C',f='A') 
        elif self.ensemble == 1:
            if self.N_chain > 0:
                self.mu_chain = self.input_vec('real', n=N_chain,s='C',f='A') 
            if self.N_solvent > 0:
                self.mu_solvent = self.input_vec('real',n=N_solvent,s='C',f='A') 

    def output_chemistry(self):
        ''' Analog of subroutine input_chemistry in chemistry_mod.f '''
        self.output_var( 'int', 'N_monomer' )
	N_monomer = self.N_monomer
        #self.output_vec('char', self.names, N_monomer, 'names')
        self.output_vec('real', 'kuhn', N_monomer )
        self.output_var('char', 'chi_flag' )
        if  self.chi_flag == 'B' or self.chi_flag == 'b':
            self.output_mat('real','chi',N_monomer,N_monomer,s='L')
        elif self.chi_flag == 'T' or self.chi_flag == 't':
            self.output_mat('real','chiA',N_monomer,N_monomer,s='L')
            self.output_mat('real','chiB',N_monomer,N_monomer,s='L')
            self.output_var('real', 'Temperature')
        self.output_var( 'int', 'N_chain')
	N_chain = self.N_chain
	if N_chain > 0:
            self.output_vec( 'int', 'N_block', N_chain, s='C')
            self.file.write('block_monomer'+"\n")
            for j in range(self.N_chain):
               self.io.output_vec(self.file,'int',self.block_monomer[j],self.N_block[j],f='N') 
            self.file.write('block_length'+"\n")
            for j in range(self.N_chain):
                self.io.output_vec(self.file,'real',self.block_length[j],self.N_block[j],f='N') 
        self.output_var('int','N_solvent')
	N_solvent = self.N_solvent
        if self.N_solvent > 0:
            self.output_vec('int','solvent_monomer',N_solvent,s='C')
            self.output_vec('real','solvent_size',N_solvent,s='C')
        self.output_var('int', 'ensemble')
        if self.ensemble == 0:
            if N_chain > 0:
                self.output_vec('real','phi_chain',N_chain,s='C',f='A') 
            if N_solvent > 0:
                self.output_vec('real','phi_solvent',N_solvent,s='C',f='A') 
        elif self.ensemble == 1:
            if N_chain > 0:
                self.output_vec('real','mu_chain',N_chain,s='C',f='A') 
            if N_solvent > 0:
                self.output_vec('real','mu_solvent',N_solvent,s='C',f='A') 

    def output_monomers(self):
        ''' Analog of subroutine output_monomers in chemistry_mod.f '''
        self.output_var( 'int', 'N_monomer' )
        self.output_vec('real', 'kuhn', self.N_monomer )

    def output_interaction(self):
        ''' Analog of subroutine output_interaction in chemistry_mod.f '''
	N_monomer = self.N_monomer
        self.output_var('char', 'interaction_type' )
        if  self.interaction_type == 'chi':
            self.output_mat('real','chi',N_monomer,N_monomer,s='L')
        if  self.interaction_type == 'chi_T':
            self.output_mat('real','chiA',N_monomer,N_monomer,s='L')
            self.output_mat('real','chiB',N_monomer,N_monomer,s='L')
            self.output_var('real', 'Temperature')

    def output_chains(self):
        self.output_var( 'int', 'N_chain')
	N_chain = self.N_chain
	if N_chain > 0:
            self.output_vec( 'int', 'N_block', N_chain, s='C')
            self.file.write('block_monomer'+"\n")
            for j in range(self.N_chain):
               self.io.output_vec(self.file,'int',self.block_monomer[j],self.N_block[j],f='N') 
            self.file.write('block_length'+"\n")
            for j in range(self.N_chain):
                self.io.output_vec(self.file,'real',self.block_length[j],self.N_block[j],f='N') 

    def output_solvents(self):
        self.output_var('int','N_solvent')
	N_solvent = self.N_solvent
        if self.N_solvent > 0:
            self.output_vec('int','solvent_monomer',N_solvent,s='C')
            self.output_vec('real','solvent_size',N_solvent,s='C')

    def output_composition(self):
        self.output_var('int', 'ensemble')
        N_chain   = self.N_chain
        N_solvent = self.N_solvent
        if self.ensemble == 0:
            if N_chain > 0:
                self.output_vec('real','phi_chain',N_chain,s='C',f='A') 
            if N_solvent > 0:
                self.output_vec('real','phi_solvent',N_solvent,s='C',f='A') 
        elif self.ensemble == 1:
            if N_chain > 0:
                self.output_vec('real','mu_chain',N_chain,s='C',f='A') 
            if N_solvent > 0:
                self.output_vec('real','mu_solvent',N_solvent,s='C',f='A') 

    def input_unit_cell(self):
        ''' Analog of subroutine input_unit_cell in unit_cell_mod.f '''
        self.dim = self.input_var('int')
        self.crystal_system = self.input_var('char')
        self.N_cell_param = self.input_var('int')
        self.cell_param = self.input_vec('real',self.N_cell_param)

    def output_unit_cell(self):
        ''' Analog of subroutine output_unit_cell in unit_cell_mod.f '''
        self.output_var('int','dim')
        self.output_var('char','crystal_system')
        self.output_var('int','N_cell_param')
        self.output_vec('real','cell_param',self.N_cell_param)

    def input_increments(self):
        ''' Analog of subroutine input_increments in sweep_mod.f '''
	self.increments = {}
        next = 1
        while next:
            comment = self.file.readline().strip()
            self.increments[comment] = 1
            if comment == 'd_kuhn':
                 self.d_kuhn = self.input_vec('real',f='N')
            elif comment == 'd_chi':
                 self.d_chi = \
                 self.input_mat('real',self.N_monomer,self.N_monomer,f='N',s='L')
            elif comment == 'd_temperature':
                 self.d_temperature = self.input_var('real',f='N')
            elif comment == 'd_block_length':
                 self.d_block_length = []
                 for i in range(self.N_chain):
                     self.d_block_length.append(self.input_vec('real',f='N'))
            elif comment == 'd_phi' or comment == 'd_phi_chain':
                 self.increments['d_phi_chain'] = 1
                 self.d_phi_chain = []
                 for i in range(self.N_chain):
                     self.d_phi_chain.append(self.input_var('real',f='N'))
            elif comment == 'd_mu' or comment == 'd_mu_chain':
                 self.increments['d_mu_chain'] = 1
                 self.d_mu_chain = []
                 for i in range(self.N_chain):
                     self.d_mu_chain.append(self.input_var('real',f='N'))
            elif comment == 'd_cell_param':
                 self.d_cell_param(self.input_vec('real',f='N'))
            elif comment == 'end_increments':
                 next = 0

    def output_increments(self):
        ''' Analog of subroutine output_increments in sweep_mod.f '''
	N_mon = self.N_monomer
        if self.increments.has_key('d_kuhn'):
             self.output_vec('real','d_kuhn',f='A')
        if self.increments.has_key('d_chi'):
             self.output_mat('real','d_chi',N_mon,N_mon,f='A',s='L')
        if self.increments.has_key('d_temperature'):
             self.output_var('real','temperature',f='A')
        if self.increments.has_key('d_block_length'):
             self.file.write('d_block_length' + "\n")
             for i in range(self.N_chain):
		 d_block_length = self.d_block_length[i]
		 N_block         = self.N_block[i]
                 self.output_vec('real','d_block_length',N_block,f='N')
        if self.increments.has_key('d_phi_chain'):
             self.file.write('d_phi_chain' + "\n")
             for i in range(self.N_chain):
                self.output_var('real',self.d_phi_chain[i],'d_phi_chain',f='N')
        if self.increments.has_key('d_mu_chain'):
             for i in range(self.N_chain):
                self.output_var('real',self.d_mu_chain[i],'d_mu_chain',f='A')
        if self.increments.has_key('d_cell_param'):
                self.output_var('real',self.d_cell_param,self.N_cell_param,'d_cell_param',f='A')
        self.file.write('end_increments' + "\n")

    # Input methods (wrapper for self.io.input_... methods of IO)
    def input_var(self, type, comment = None, f='A'):
        return self.io.input_var(self.file, type, comment, f)

    def input_vec(self, type, n=None, comment=None, s='R',f='A'):
        return self.io.input_vec(self.file, type, n, comment, s, f)

    def input_mat(self, type, m, n=None, comment=None, s='L', f='A'):
        return self.io.input_mat(self.file, type, m, n, comment, s, f)

    # Output methods (output by name)
    def output_var(self, type, name, f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_var(self.file, type, data, name, f)

    def output_vec(self, type, name, n=None, s='R', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_vec(self.file, type, data, n, name, s, f)

    def output_mat(self, type, name, m, n=None, s='L', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_mat(self.file, type, data, m, n, name, s, f)


    def __getitem__(self,key):
        return self.__dict__[key]

    def __str__(self):
        s = []
        for key in self.att_names:
            s.append( key +  ' : ' + str( self[key] ) )
        return string.join(s,'\n')

    def eval(self,expr1):
        '''
	Returns the value of a python expression calculated
	by using the key names of attributes of an Outfile
	as variable names. 
	'''
        for key in self.__dict__.keys():
            exec( key + '= self.' + key )
        return eval(expr1)
         
if __name__ == '__main__':

    '''
    USAGE
      scf_out.py filename
    PURPOSE
      Reads and parses scf output file with path filename.
      Prints the dictionary of attributes
    '''
    import sys 

    # Read file "filename"
    filename = sys.argv[1]
    x = OutFile(filename)
    print x

    # Print dictionary of attributes
    # for key in x.__dict__.keys():
    #    print key, ':', x.__dict__[key]
