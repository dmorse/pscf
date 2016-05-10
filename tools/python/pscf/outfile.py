from io import *
from version import *
import string 
import sys

class OutFile(object):
    '''
    An OutFile object contains the data in a output summary
    file produced by F90 program pscf.
    '''
   
    def __init__(self,filename):
        '''
        PURPOSE
          Read and parse out file filename. 
          Create an OutFile object.
        ARGUMENT
          filename - name of PSCF output summary file (string)
        COMMENT
          The file is opened and closed within body of method
        '''
        self.file = open(filename, 'r')
        self.io   = IO()
        file = self.file
        self.N_chain   = 0
        self.N_solvent = 0

        # Dictionary flags indicates which sections are present
        # Keys are flag string values ('MONOMERS' etc.), values are all 1
        self.flags = {}

        # List sections lists section names (flags) in order read
        self.sections = []

        # Read version line
        self.version = Version(self.file)

        # Read input parameter file sections
        next = 1
        while next:
            next = self.read_param_section(file)

        # Read additional output sections
        next = 1
        while next :
            next = self.read_output_section()

        self.file.close()
        self.file = None

        # Make sorted list of attribute names
        self.att_names = self.__dict__.keys()
        self.att_names.sort()

    def write(self, file, major=1, minor=0):
        # If file argument is a string, open a file of that name
        if type(file) == type('thing'):
            temp = open(file,'w')
            file = temp
        self.file = file
        self.version.major = major
        self.version.minor = minor
        self.version.write(file)
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
            self._output_vec( 'int', 'N_grid', self.dim)
            self._output_var( 'real', 'chain_step')
            #self._output_var( 'int', 'extrapolation_order')
        if self.flags.has_key('BASIS'):
            file.write("\n%-15s\n" % 'BASIS')
            self._output_var('char', 'group_name')
        if self.flags.has_key('ITERATE'):
            file.write("\n%-15s\n" % 'ITERATE')
            self.output_iterate()
        if self.flags.has_key('SWEEP'):
            file.write("\n%-15s\n" % 'SWEEP')
            self._output_var( 'real', 's_max')
            self.output_increments()
        file.write("\n%-15s\n" % 'FINISH')

        if self.flags.has_key('THERMO'):
            file.write("\n%-15s\n" % 'THERMO')
            self.output_thermo()

        if self.flags.has_key('DECOMPOSE'):
            file.write("\n%-15s\n" % 'DECOMPOSE')
            self.output_decompose()

        if self.flags.has_key('STATISTICS') :
            file.write("\n%-15s\n" % 'STATISTICS')
            self.output_statistics()

        file.close()
        self.file = None

    def read_param_section(self, file):

        # Read the next non-empty line
        hasFlag = False
        while not hasFlag:
            line = self.file.readline()
            if not line:
                next = 0
                return next
            flag = line.strip()
            if flag != '':
                hasFlag = True

        # Set key in self.flags dictionary
        self.flags[flag] = 1
        self.sections.append(flag)

        next = 1
        if flag == 'MONOMERS':
            self.input_monomers()
        elif flag == 'CHAINS':
            self.input_chains()
        elif flag == 'SOLVENTS':
            self.input_solvents()
        elif flag == 'COMPOSITION':
            self.input_composition()
        elif flag == 'INTERACTION':
            self.input_interaction()
        elif flag == 'UNIT_CELL':
            self.input_unit_cell()
        elif flag == 'DISCRETIZATION':
            self.N_grid = self._input_vec('int')
            self.chain_step = self._input_var('real')
        elif flag == 'BASIS':
            self.group_name = self._input_var('char')
        elif flag == 'ITERATE':
            self.input_iterate()
        elif flag == 'FINISH':
            next = 0
        else:
            next = 0

        return next

    def read_output_section(self):
        next = 1

        # Read next non-empty line
        hasFlag = False
        while not hasFlag:
            line = self.file.readline()
            if not line:
                next = 0
                return next
            flag = line.strip()
            if flag != '':
                hasFlag = True

        # Read output section
        self.flags[flag] = 1
        if flag == 'THERMO':
            self.input_thermo()
        if flag == 'DECOMPOSE':
            self.input_decompose()
        elif flag == 'STATISTICS':
            self.input_statistics()
        else:
            next = 0

        return next


    def input_monomers(self):
        ''' Analog of subroutine input_monomers in chemistry_mod.f '''
        # Monomers 
        self.N_monomer = self._input_var('int')
        N_monomer      = self.N_monomer
        self.kuhn      = self._input_vec('real')

    def input_chains(self):
        self.N_chain = self._input_var('int',f='A')
        if self.N_chain:
            self.N_block = self._input_vec('int',n=self.N_chain,s='C')
            self.file.readline()
            self.block_monomer = []
            for j in range(self.N_chain):
                self.block_monomer.append( self._input_vec('int',f='N') )
            self.file.readline()
            self.block_length = []
            for j in range(self.N_chain):
                self.block_length.append( self._input_vec('real',f='N') )
        else:
            self.N_chain = 0

    def input_solvents(self):
        self.N_solvent = self._input_var('int',f='A')
        if self.N_solvent:
            self.solvent_monomer = self._input_vec('int', self.N_solvent, s='C')
            self.solvent_size    = self._input_vec('real',self.N_solvent, s='C')
        else:
            self.N_solvent = 0

    def input_composition(self):
        self.ensemble = self._input_var('int',f='A')
        N_chain   = self.N_chain
        N_solvent = self.N_solvent
        if self.ensemble == 0:
            if self.N_chain > 0:
                self.phi_chain = self._input_vec('real',n=N_chain,s='C',f='A') 
            if self.N_solvent > 0:
                self.phi_solvent = self._input_vec('real',n=N_solvent,s='C',f='A') 
        elif self.ensemble == 1:
            if self.N_chain > 0:
                self.mu_chain = self._input_vec('real', n=N_chain,s='C',f='A') 
            if self.N_solvent > 0:
                self.mu_solvent = self._input_vec('real',n=N_solvent,s='C',f='A') 

    def input_interaction(self):
        ''' Analog of subroutine input_interaction in chemistry_mod.f '''
        self.interaction_type = self._input_var('char')
        N_monomer = self.N_monomer
        if self.interaction_type == 'chi':
            self.chi = self._input_mat('real',N_monomer,N_monomer,s='L')
        elif self.interaction_type == 'chi_T':
            self.chi_A = self._input_mat('real',N_monomer,N_monomer,s='L')
            self.chi_B = self._input_mat('real',N_monomer,N_monomer,s='L')
            self.Temperature = self._input_var('real')

    def output_monomers(self):
        ''' Analog of subroutine output_monomers in chemistry_mod.f '''
        self._output_var( 'int', 'N_monomer' )
        self._output_vec('real', 'kuhn', self.N_monomer )

    def output_chains(self):
        self._output_var( 'int', 'N_chain')
        N_chain = self.N_chain
        if N_chain > 0:
            self._output_vec( 'int', 'N_block', N_chain, s='C')
            self.file.write('block_monomer'+"\n")
            for j in range(self.N_chain):
               self.io.output_vec(self.file, 'int',self.block_monomer[j],self.N_block[j],f='N') 
            self.file.write('block_length'+"\n")
            for j in range(self.N_chain):
                self.io.output_vec(self.file, 'real',self.block_length[j],self.N_block[j],f='N') 

    def output_solvents(self):
        self._output_var('int', 'N_solvent')
        N_solvent = self.N_solvent
        if self.N_solvent > 0:
            self._output_vec('int', 'solvent_monomer',N_solvent,s='C')
            self._output_vec('real', 'solvent_size',N_solvent,s='C')

    def output_composition(self):
        self._output_var('int', 'ensemble')
        N_chain   = self.N_chain
        N_solvent = self.N_solvent
        if self.ensemble == 0:
            if N_chain > 0:
                self._output_vec('real', 'phi_chain',N_chain,s='C',f='A') 
            if N_solvent > 0:
                self._output_vec('real', 'phi_solvent',N_solvent,s='C',f='A') 
        elif self.ensemble == 1:
            if N_chain > 0:
                self._output_vec('real', 'mu_chain',N_chain,s='C',f='A') 
            if N_solvent > 0:
                self._output_vec('real', 'mu_solvent',N_solvent,s='C',f='A') 

    def output_interaction(self):
        ''' Analog of subroutine output_interaction in chemistry_mod.f '''
        N_monomer = self.N_monomer
        self._output_var('char', 'interaction_type' )
        if  self.interaction_type == 'chi':
            self._output_mat('real', 'chi',N_monomer,N_monomer,s='L')
        if  self.interaction_type == 'chi_T':
            self._output_mat('real', 'chiA',N_monomer,N_monomer,s='L')
            self._output_mat('real', 'chiB',N_monomer,N_monomer,s='L')
            self._output_var('real', 'Temperature')

    def input_unit_cell(self):
        ''' Analog of subroutine input_unit_cell in unit_cell_mod.f '''
        self.dim = self._input_var('int')
        self.crystal_system = self._input_var('char')
        self.N_cell_param = self._input_var('int')
        self.cell_param = self._input_vec('real',self.N_cell_param)

    def output_unit_cell(self):
        ''' Analog of subroutine output_unit_cell in unit_cell_mod.f '''
        self._output_var('int', 'dim')
        self._output_var('char', 'crystal_system')
        self._output_var('int', 'N_cell_param')
        self._output_vec('real', 'cell_param', self.N_cell_param)

    def input_iterate(self):
        self.input_filename = self._input_var('char')
        self.output_prefix = self._input_var('char')
        self.max_itr = self._input_var('int')
        self.error_max = self._input_var('real')
        self.domain = self._input_var('logic')
        self.itr_algo = self._input_var('char')
        if self.itr_algo == 'NR':
            self.N_cut = self._input_var('int')
        if self.itr_algo == 'AM':
            self.N_history = self._input_var('int')

    def output_iterate(self):
        self._output_var('char', 'input_filename')
        self._output_var('char', 'output_prefix')
        self._output_var('int', 'max_itr')
        self._output_var('real', 'error_max')
        self._output_var('logic', 'domain')
        self._output_var('char', 'itr_algo')
        if self.itr_algo == 'NR':
            self._output_var('int', 'N_cut')
        if self.itr_algo == 'AM':
            self._output_var('int', 'N_history')

    def input_increments(self):
        ''' Analog of subroutine input_increments in sweep_mod.f '''
        self.increments = {}
        next = 1
        while next:
            comment = self.file.readline().strip()
            self.increments[comment] = 1
            if comment == 'd_kuhn':
                 self.d_kuhn = self._input_vec('real',f='N')
            elif comment == 'd_chi':
                 self.d_chi = \
                 self._input_mat('real',self.N_monomer,self.N_monomer,f='N',s='L')
            elif comment == 'd_temperature':
                 self.d_temperature = self._input_var('real',f='N')
            elif comment == 'd_block_length':
                 self.d_block_length = []
                 for i in range(self.N_chain):
                     self.d_block_length.append(self._input_vec('real',f='N'))
            elif comment == 'd_phi' or comment == 'd_phi_chain':
                 self.increments['d_phi_chain'] = 1
                 self.d_phi_chain = []
                 for i in range(self.N_chain):
                     self.d_phi_chain.append(self._input_var('real',f='N'))
            elif comment == 'd_mu' or comment == 'd_mu_chain':
                 self.increments['d_mu_chain'] = 1
                 self.d_mu_chain = []
                 for i in range(self.N_chain):
                     self.d_mu_chain.append(self._input_var('real',f='N'))
            elif comment == 'd_cell_param':
                 self.d_cell_param(self._input_vec('real',f='N'))
            elif comment == 'end_increments':
                 next = 0

    def output_increments(self):
        ''' Analog of subroutine output_increments in sweep_mod.f '''
        N_mon = self.N_monomer
        if self.increments.has_key('d_kuhn'):
            self._output_vec('real', 'd_kuhn',f='A')
        if self.increments.has_key('d_chi'):
            self._output_mat('real', 'd_chi',N_mon,N_mon,f='A',s='L')
        if self.increments.has_key('d_temperature'):
            self._output_var('real', 'temperature',f='A')
        if self.increments.has_key('d_block_length'):
            self.file.write('d_block_length' + "\n")
            for i in range(self.N_chain):
                d_block_length = self.d_block_length[i]
                N_block = self.N_block[i]
                self._output_vec('real', 'd_block_length',N_block,f='N')
        if self.increments.has_key('d_phi_chain'):
            self.file.write('d_phi_chain' + "\n")
            for i in range(self.N_chain):
                self._output_var('real',self.d_phi_chain[i], 'd_phi_chain',f='N')
        if self.increments.has_key('d_mu_chain'):
             for i in range(self.N_chain):
                self._output_var('real',self.d_mu_chain[i], 'd_mu_chain',f='A')
        if self.increments.has_key('d_cell_param'):
            self._output_var('real',self.d_cell_param,self.N_cell_param, 'd_cell_param',f='A')
        self.file.write('end_increments' + "\n")

    def input_thermo(self):
        self.f_Helmholtz = self._input_var('real')
        self.f_homo = self._input_var('real')
        self.pressure = self._input_var('real')
        if self.ensemble == 0 :
            if self.N_chain > 0 :
                self.mu_chain = \
                   self._input_vec('real',self.N_chain,s='C')
            if self.N_solvent > 0 :
                self.mu_solvent = \
                   self._input_vec('real',self.N_solvent,s='C')
        if self.ensemble == 1 :
            if self.N_chain > 0 :
                self.phi_chain = \
                   self._input_vec('real',self.N_chain,s='C')
            if self.N_solvent > 0 :
                self.phi_solvent = \
                   self._input_vec('real',self.N_solvent,s='C')
        self.stress = self._input_vec('real',self.N_cell_param)

    def output_thermo(self):
        self._output_var('real', 'f_Helmholtz')
        self._output_var('real', 'f_homo')
        self._output_var('real', 'pressure')
        if self.ensemble == 0 :
            if self.N_chain > 0 :
                self._output_vec('real', \
                                 'mu_chain',self.N_chain,s='C')
            if self.N_solvent > 0 :
                self._output_vec('real', \
                                 'mu_solvent',self.N_solvent,s='C')
        elif self.ensemble == 1 :
            if self.N_chain > 0 :
                self._output_vec('real', \
                                 'phi_chain',self.N_chain,s='C')
            if self.N_solvent > 0 :
                self._output_vec('real', \
                                 'phi_solvent',self.N_solvent,s='C')
        self._output_vec('real', 'stress',self.N_cell_param)

    def input_statistics(self):
        self.N_star      = self._input_var('int')
        self.final_error = self._input_var('real')
        self.iterations  = self._input_var('int')
        self.basis_time  = self._input_var('real')
        self.scf_time    = self._input_var('real')

    def output_statistics(self):
        self._output_var('int', 'N_star')
        self._output_var('real', 'final_error')
        self._output_var('int', 'iterations')
        self._output_var('real', 'basis_time')
        self._output_var('real', 'scf_time')

    def input_decompose(self):
        self.overlap_AB  = self._input_var('real')
        if self.N_monomer > 2 :
            self.overlap_BC  = self._input_var('real')
            self.overlap_CA  = self._input_var('real')
        self.f_enthalpy  = self._input_var('real')
        self.f_head      = self._input_var('real')
        self.f_tail      = self._input_var('real')
        self.f_excess    = self._input_var('real')

    def output_decompose(self):
        self._output_var('real', 'overlap_AB')
        if self.N_monomer > 2 :
            self._output_var('real', 'overlap_BC')
            self._output_var('real', 'overlap_CA')
        self._output_var('real', 'f_enthalpy')
        self._output_var('real', 'f_head')
        self._output_var('real', 'f_tail')
        self._output_var('real', 'f_excess')

    # "Private" methods

    # Input methods (wrapper for self.io.input_... methods of IO)
    def _input_var(self, type, comment = None, f='A'):
        return self.io.input_var(self.file, type, comment, f)

    def _input_vec(self, type, n=None, comment=None, s='R',f='A'):
        return self.io.input_vec(self.file, type, n, comment, s, f)

    def _input_mat(self, type, m, n=None, comment=None, s='L', f='A'):
        return self.io.input_mat(self.file, type, m, n, comment, s, f)

    # Output methods (output by name)
    def _output_var(self, type, name, f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
            self.io.output_var(self.file, type, data, name, f)

    def _output_vec(self, type, name, n=None, s='R', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
            self.io.output_vec(self.file, type, data, n, name, s, f)

    def _output_mat(self, type, name, m, n=None, s='L', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
            self.io.output_mat(self.file, type, data, m, n, name, s, f)


    def __getitem__(self,key):
        return self.__dict__[key]

    def __str__(self):
        s = []
        for key in self.att_names:
            s.append( key +  ' : ' + str( self[key] ) )
        return string.join(s, '\n')

    def eval(self, expr1):
        '''
        Returns the value of a python expression calculated
        by using the key names of attributes of an Outfile
        as variable names. 
        '''
        for key in self.__dict__.keys():
            exec( key + '= self.' + key )
        return eval(expr1)
         
