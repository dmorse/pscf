from io import *
from version import *
import string
import sys

class Field(object):
    '''
    A Field object contains the data in a field file in 
    the PSCF symmetry-adapated Fourier expansion format.
    '''

    # "Public" methods

    def __init__(self,filename):
        '''
        PURPOSE
          Constructor: Read a PSCF symmetry-adapted field 
          file, and create a new Field object to hold data.
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

        self._input_unit_cell()
        self.group_name = self._input_var('char')
        self.N_monomer = self._input_var('int')
        self.N_star = self._input_var('int')

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
        '''
        PURPOSE
           Write field to file in PSCF symmetry-adapted format.
        ARGUMENTS
           file  - file object or file name string
           major - major file format version number
           minor - minor file format version number
        COMMENT
           if file is a field object, it must be open for writing
        '''

        # If file argument is a string, open a file of that name
        if type(file) == type('thing'):
            temp = open(file,'w')
            file = temp
        self.file = file
           
	self.version.major = major
	self.version.minor = minor
	self.version.write(file)

        self._output_unit_cell()
        self._output_var('char', 'group_name')
        self._output_var( 'int', 'N_monomer')
        self._output_var( 'int', 'N_star')

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

    def addMonomer(self, value = 0.0):
        ''' 
        PURPOSE
           Add a field with a coefficients set to common value
        ARGUMENTS
           value - value of all coefficients for new momoner
        COMMENT
            N_momomer is incremented by 1, and new field is last
        '''
        self.N_monomer += 1
        for k in range(self.N_star):
            self.fields[k].append(value)

    def duplicateMonomer(self, i):
        ''' 
        PURPOSE
           Add a field by duplicating field i
        ARGUMENTS
           i - index in range [0, N_monomer-1]
        COMMENT
            N_momomer is incremented, and duplicate field is last
        '''
        self.N_monomer += 1
        for k in range(self.N_star):
            self.fields[k].append(self.fields[k][i])

    def switchMonomers(self, i, j):
        '''
        PURPOSE
           Switch coefficients of fields i and j
        ARGUMENTS
           i - index in range [0, N_monomer-1]
           j - index in range [0, N_monomer-1]
        '''
        for k in range(self.N_star):
            temp = self.fields[k][i]
            self.fields[k][i] = self.fields[k][j]
            self.fields[k][j] = temp

    # "Private" methods

    # Wrappers for input_... output_.. methods of IO)

    def _input_var(self, type, comment = None, f='A'):
        return self.io.input_var(self.file, type, comment, f)

    def _input_vec(self, type, n=None, comment=None, s='R',f='A'):
        return self.io.input_vec(self.file, type, n, comment, s, f)

    # Output methods (output by name)
    def _output_var(self, type, name, f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_var(self.file, type, data, name, f)

    def _output_vec(self, type, name, n=None, s='R', f='A'):
        if self.__dict__.has_key(name):
            data = self.__dict__[name]
	    self.io.output_vec(self.file, type, data, n, name, s, f)

    def _input_unit_cell(self):
        ''' Analog of subroutine _input_unit_cell in unit_cell_mod.f '''
        self.dim = self._input_var('int')
        self.crystal_system = self._input_var('char')
        self.N_cell_param = self._input_var('int')
        self.cell_param = self._input_vec('real',self.N_cell_param)

    def _output_unit_cell(self):
        ''' Analog of subroutine _output_unit_cell in unit_cell_mod.f '''
        self._output_var('int', 'dim')
        self._output_var('char', 'crystal_system')
        self._output_var('int', 'N_cell_param')
        self._output_vec('real', 'cell_param', self.N_cell_param)

