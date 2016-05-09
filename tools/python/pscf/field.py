from io import *
from version import *
import string 
import sys

class Field(object):
    '''
    An Field object contains the data in a field file
    produced by F90 program pscf 
    '''
   
    def __init__(self,filename):
        '''
        PURPOSE
          Read and parse field file.
          Create a Field object.
        ARGUMENT
          filename - string 
        COMMENT
          File is opened and closed within body of method
        '''
        self.file = open(filename, 'r')
	self.io   = IO()
	file = self.file

        # Read version line
	self.version = Version(self.file)

        self.__input_unit_cell()
        self.group_name = self.__input_var('char')
        self.N_monomer = self.__input_var('int')
        self.N_star = self.__input_var('int')

        # Define empty lists
        self.fields = []
        self.waves = []
        self.counts = []

        for i in range(self.N_star):

            data = file.readline().split()
            if len(data) != self.N_monomer + self.dim + 1:
                raise IoException('Incorrect number of elements in field line')
            j = 0

            # Read field coefficients
            self.fields.append([])
            for k in range(self.N_monomer):
                value = float(data[j])
                self.fields[i].append(value)
                j += 1

            # Read field coefficients
            self.waves.append([])
            for k in range(self.dim):
                value = int(data[j])
                self.waves[i].append(value)
                j += 1

            # Read star_count
            self.counts.append(int(data[j]))

        self.file.close()
	self.file = None

    def write(self, file, major=1, minor=0):
	self.file = file
	self.version.major = major
	self.version.minor = minor
	self.version.write(file)

        self.__output_unit_cell()
        self.__output_var('char', 'group_name')
        self.__output_var( 'int', 'N_monomer')
        self.__output_var( 'int', 'N_star')

        for i in range(self.N_star):
            for k in range(self.N_monomer):
                file.write('%20.12E' % self.fields[i][k])
            file.write('    ')
            for k in range(self.dim):
                file.write('%4d' % self.waves[i][k])
            file.write('%6d' % self.counts[i])
            file.write("\n")

	file.close()
	self.file = None

    # "Private" methods

    # Input methods (wrapper for self.io.input_... methods of IO)
    def __input_var(self, type, comment = None, f='A'):
        return self.io.input_var(self.file, type, comment, f)

    def __input_vec(self, type, n=None, comment=None, s='R',f='A'):
        return self.io.input_vec(self.file, type, n, comment, s, f)

    # Output methods (output by name)
    def __output_var(self, type, name, f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_var(self.file, type, data, name, f)

    def __output_vec(self, type, name, n=None, s='R', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_vec(self.file, type, data, n, name, s, f)

    def __input_unit_cell(self):
        ''' Analog of subroutine __input_unit_cell in unit_cell_mod.f '''
        self.dim = self.__input_var('int')
        self.crystal_system = self.__input_var('char')
        self.N_cell_param = self.__input_var('int')
        self.cell_param = self.__input_vec('real',self.N_cell_param)

    def __output_unit_cell(self):
        ''' Analog of subroutine __output_unit_cell in unit_cell_mod.f '''
        self.__output_var('int', 'dim')
        self.__output_var('char', 'crystal_system')
        self.__output_var('int', 'N_cell_param')
        self.__output_vec('real', 'cell_param', self.N_cell_param)

