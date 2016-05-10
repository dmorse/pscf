from outfile import OutFile
import sys, os, string

class Sweep(object):
    """
    Contains data in all summary output files produced by a PSCF sweep.

    The constructor reads in, parses and stores all of the files with
    names of the form i.out in a directory, where i is an integer in
    a sequence range [0, n]. 

    The class emulates a list of Outfile objects: If s is a Sweep object 
    constructed from directory 'dir', then s[i] returns the OutFile object 
    constructed by parsing the file 'dir/i.out' , where i is an integer.

    Implementation: Use of the [] subscript operator is implemented by 
    overloading the __getitem__ method.
    """
   
    def __init__(self, directory):
        """
        Create a Sweep object by reading all output files in a directory.

        Arguments:
           directory -- path name for directory [string]
        """
        old_cwd = os.getcwd()
        os.chdir(directory)

        # Find filenames of *.out files
        filenames = {}
        for file in os.listdir('.'):
            if file.endswith(r".out"):
                filenames[file[:-4]] = file
        n = len(filenames)

        # Make list self._outfiles of OutFile objects
        self._outfiles = []
        i = 0
        while i < n:
            key = str(i)
            if key in filenames.keys():
                filename = filenames[key]
                self._outfiles.append(OutFile(filename))
            else:
                self._outfiles.append(None)
                n = n + 1
            i = i + 1

        # Return to original current working directory
        os.chdir(old_cwd)

    def write(self, expr1, expr2):
        """
        Return a two column string containing values of python expressions.

        The arguments expr1 and expr2 are python expressions with 
        numerical (integer or real) values that are constructed using 
        the names of attributes of individual OutFile objects as variables. 
        The value of each expression is evaluated using the eval method of
        the ParamFile class, which is inherited by OutFile. 

        This class is intended to make it easy to plot output pairs of
        values in a form suitable for plotting or further analysis. The
        simplest application is to give expr1 (the first column) as the
        value of some control parameter that is varied in the sweep, and
        expr2 (the second column) as f_Helmholtz, the Helmholtz free 
        energy per monomer.
        """
        s = []
        for i in range(len(self._outfiles)):
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

    def __getitem__(self,i):
        """
        Allow a Sweep to emulate a list of OutFile objects.
        """
        return self._outfiles[i]
    
