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

