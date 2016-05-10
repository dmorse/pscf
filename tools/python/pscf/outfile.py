from io import *
from version import *
from paramfile import *
import string 
import sys

class OutFile(ParamFile):
    '''
    An OutFile object contains the data in a output summary
    file produced by F90 program pscf. 

    Class Outfile is dervied from ParamFile because the output 
    file format begins with a parameter file section, to which 
    a few output sections are added. 
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
        elif flag == 'DECOMPOSE':
            self.input_decompose()
        elif flag == 'STATISTICS':
            self.input_statistics()
        else:
            next = 0

        return next

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

#    def eval(self, expr1):
#        '''
#        Returns the value of a python expression calculated
#        by using the key names of attributes of an Outfile
#        as variable names. 
#        '''
#        for key in self.__dict__.keys():
#            exec( key + '= self.' + key )
#        return eval(expr1)
#         
