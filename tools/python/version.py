#!/usr/bin/env python2.4

from io import *
import string 

class Version(object):
    '''
    Abstract representation of a file format version tag
    '''
   
    def __init__(self,file):
        '''
        PURPOSE
        ARGUMENT
          file - file handle, opened for reading
        COMMENT
          File is opened and closed within body of method
        '''

        # Read first line to determine file format
	line = file.readline().strip()
        if line == '':
            self.major = 0
            self.minor = 9
        else:
  	    line = line.split()
	    if line[0] == 'format':
                self.major = int(line[1])
                self.minor = int(line[2])
            else: 
                raise 'Invalid file format line' 

    def write(self, file, major=None, minor=None):
        if major is None:
            major = self.major
        if minor is None:
            minor = self.minor
        line = 'format %2d %-3d' % (major, minor)
        file.write(line + "\n")

    def eq(self, major, minor):
        if (self.major == major and self.minor == minor):
            return 1
        else:
            return 0

    def ge(self, major, minor):
        if (self.major >= major and self.minor >= minor):
            return true
        else:
            return false

