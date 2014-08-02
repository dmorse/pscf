#!/usr/bin/env python2.4

from   outfile import OutFile
import sys, os, string

class Sweep(object):
    '''
    PURPOSE
       A Sweep object contains the data in all the *.out output files 
       in a directory produced by a scf sweep. 

       The class emulates a list of OutFile objects by overloading the
       __getitem__ method: If s is a Sweep object constructed from 
       directory 'dir', then s[i] is the OutFile object constructed 
       by parsing the file 'dir/i.out' , where i is an integer.
    INSTANCE ATTRIBUTES
       self.outfiles -- a list of OutFile objects. Item self.outfile[i] 
                        is produced from directory/i.out
    '''
   
    def __init__(self,directory):
        '''
        PURPOSE
           Create a Sweep object by reading all of the integer.out
           files in a specified directory created by an scf sweep.
        ARGUMENT
           directory - string, path name for directory
        '''
        old_cwd = os.getcwd()
        os.chdir(directory)

        #Find filenames of *.out files
        filenames = {}
        for file in os.listdir('.'):
            if file.endswith(r".out"):
                filenames[file[:-4]] = file
        n = len(filenames)

        # Make list self.outfiles of scf_out objects
        self.outfiles = []
        i = 0
        while i < n:
            key = str(i)
            if key in filenames.keys():
                filename = filenames[key]
                self.outfiles.append(OutFile(filename))
            else:
                self.outfiles.append(None)
                n = n + 1
            i = i + 1

        # Return to original current working directory
        os.chdir(old_cwd)

    def __getitem__(self,i):
        return self.outfiles[i]

    def write(self,expr1,expr2):
        s = []
        for i in range(len(self.outfiles)):
            if self[i]:
                # Make local variables: For each attribute 
                # x[i].key of OutFile object self[i], make a
                # corresponding local variable named "key", 
                # by executing: key = self[i].key
                for key in self[i].__dict__.keys():
                    exec( key + '= self[i].' + key )
                # Print results of expressions expr1 and expr2
                line = str(eval(expr1)) +  '   ' + str(eval(expr2)) 
                s.append(line)
        return string.join(s,"\n")


if __name__ == '__main__':

    '''
    SCRIPT
       sweep.py directory [ expr1 expr2 ]
    PURPOSE
       Reads outfiles from the directory with the specified path, 
       creates a corresponding Sweep object, and prints a two-column
       summary of results. The format of the output depends upon 
       whether optional command line arguments expr1 and expr2 are 
       present:

       1) If expr1 and expr2 are absent, each entry in column 1 of 
       the output is an integer i corresponding to the integer in
       a filename directory/i.out, and column 2 is the Helmholtz 
       free energy reported in that output file.

       2) If present, arguments expr1 and expr2 are Python
       expressions that may refer to local variables whose names 
       and values correspond to the names and values of attributes 
       of the current OutFile object. The names of the allowed 
       variables also correspond to names of corresponding 
       variables in the scf code. 
    EXAMPLES
       1) To produce a two-column output in which the first column is
       the [0][1] element of the chi matrix (i.e., the value of chi
       for a system with two monomer types) and the second is the
       Helmholtz free energy, enter:

           > sweep.py directory chi[0][1] f_Helmholtz 

       2) To produce a file in which the second column is the difference
       between the Helmholtz free energy and that of a homogeneous
       system, enter:

           > sweep.py directory chi[0][1] 'f_Helmholtz - f_homo'

       Quotes are required around expressions that contain spaces, to 
       force the unix shell to pass the expression to the script as a 
       single argument.  
    '''

    if len(sys.argv) == 1:
        print "Error: Required directory argument not supplied"
    # Read directory and create Sweep object
    x = Sweep(sys.argv[1])
    # Read two column summary
    if len(sys.argv) == 4:
        print x.write(sys.argv[2],sys.argv[3])
    else:
        print x.write('i','f_Helmholtz')
