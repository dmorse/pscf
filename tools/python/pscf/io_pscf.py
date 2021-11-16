class IO:
    """
    Provide interfaces for formatted file IO of variables.

    The class IO is analogous to the io_mod module of scf. Class
    methods provides generic interfaces for reading in scalars, 
    vectors, and matrices for several data types.  Because variables 
    in python are not typed, only one input and one output function 
    is required for io of scalar variables, one for vectors, and one 
    for matrices. In both input and output functions, the variable 
    type is passed passed as a string argument 'type', which can have 
    values 'int', 'real', 'char' and (for scalars) 'logical'.
    """

    def __init__(self):
        """
        Default construct an initially empty object.
        """
        self.comment = None

    def input_comment(self, file, comment=None):
        """
        Read the comment associated with a variable.
        """
        if not self.comment:
            self.comment = file.readline().strip()
        if comment:
            if comment == self.comment:
                self.comment = None
                return 1
            else:
                return 0
        else:
           self.comment = None
           return 1
            
    def input_var(self,file,type,comment=None,f='A'):
        """ 
        Read and return a scalar variable from file. 
     
        Arguments:
          file -- file object (must be opened for reading)
          type -- string, 'int', 'real', 'char', or 'logic'
          f    -- 'A' -> comment string on line above data
               -- 'N' -> no comment string
        """
        if f == 'A':
            if not self.input_comment(file, comment):
                return None
        data = file.readline().strip()
        if type == 'int':
            return int(data)
        elif type == 'real':
            return float(data)
        elif type == 'char':
            return strip_quotes(data)
        elif type == 'logic':
            if data in ['T','true', 'TRUE', 'True']:
                return 1
            elif data in ['F','false','FALSE','FALSE']:
                return 0
            else:
                raise("Invalid logical variable:", data)
        else:
            raise 'Illegal type in input_var'
    
    def input_vec(self,file,type,n=None, comment=None, s='R', f='A'):
        """
        Read and return a vector/array of n variables from file.

        Arguments:
          file -- file object (must be opened for reading)
          type -- string, = 'int', 'real', 'char', or 'logic'
          n    -- number of elements in vector
          f    -- 'A' -> comment on line above data
               -- 'N' -> no comment 
          s    -- 'R' -> row vector (one line)
               -- 'C' -> column vector (n lines)
        """
        if f == 'A':
            if not self.input_comment(file,comment):
                return None
        if s == 'R':   # Row Vector
            data = file.readline().split()
            if n:
                data = data[:n]
        elif s == 'C': # Column Vector
            if not n: raise 'No n value in input_int_vec for column vector'
            data = []
            for i in range(n):
                data.append( file.readline().split()[0] )
        if type == 'int':
            return [ int(x) for x in data ]
        elif type == 'real':
            return [ float(x) for x in data ]
        elif type == 'char':
            return [ strip_quotes(x) for x in data ]
        else:
            raise 'Illegal type in input_vec'
    
    # Matrices
        
    def input_mat(self, file, type, m, n=None, comment=None, s=None,f='A'):
        """
        Read and return an m x n matrix from file.

        Arguments:
        file -- file object (must be opened for reading)
        type -- string, = 'int' or 'real'
        m    -- number of rows in matrix
        n    -- number of columns (set to m by default)
        f    -- 'A' -> comment on line above data
             -- 'N' -> no comment 
        s    -- 'N' -> full matrix
             -- 'S' -> symmetric
             -- 'L' -> symmetric lower triangular (0 diagonals)
        """
        if f == 'A':
            if not self.input_comment(file,comment):
                return None
        # Emulate initialization of m x n array
        if not n:
           n = m
        data = []
        for i in range(m):
            data.append([])
            for j in range(n):
                if type == 'int':
                    data[i].append(0) 
                elif type == 'real':
                    data[i].append(0.0) 
        # Read matrix
        if s == 'N' or s == 'S':
            min = 0
        elif s == 'L':
            min = 1
        else:
            raise 'Invalid style in input_mat: s=' + str(s)
        for i in range(min,m) :
            line = file.readline().split()
            if s == 'N':
                line = line[:n]
            elif s == 'S':
                line = line[:i+1]
            elif s == 'L':
                line = line[:i+1]
            for j in range(len(line)):
                if type == 'int':
                    datum = int(line[j])
                elif type == 'real':
                    datum = float(line[j])
                else:
                    raise 'Illegal type in input_mat'
                data[i][j] = datum
                if s == 'S' or s == 'L':
                    data[j][i] = datum
        return data
    
    
    def format_var(self, type, data):
        """
        Returns a formatted string representation of data.

        Arguments:
        type -- string identifying type (int/real/char/logic)
        data -- variable of specified type
        """
        if type == 'int':
            return '%20d' % data
        elif type == 'real':
            return '%20.10E' % data
        elif type == 'char':
            data = "'" + data + "'" 
            return "%20s" % data 
        elif type == 'logic':
            if data:
                data = 'T'
            else:
                data = 'F'
            return "%20s" % data 
        else:
            raise 'Illegal type in format_var'
      
    def output_comment(self, file,  comment, n=20):
        """
        Write comment (variable name) to file.

        Arguments:
        file     -- output file object
        comment  -- comment string
        n        -- field size (default = 20)
        """
        comment = comment.strip()
        comment = comment.ljust(n)
        comment = file.write(comment + "\n")
    
    def output_var(self, file, type, data, comment, f='A'):
        """
        Output a scalar variable to a file

        Arguments:
        file -- output file object (must be opened for writing)
        type -- type identifer string: 'int', 'real', 'char', or 'logic'
        f    -- 'A' -> comment string on line above data
             -- 'N' -> no comment string
        """
        if f == 'A':
            self.output_comment(file,comment)
        file.write(self.format_var(type,data) + "\n")
    
    def output_vec(self, file, type, data, n=None, comment=None, s='R', f='A'):
        """
        Output a vector/array to a file.
       
        Arguments:
        file -- file object (must be opened for reading)
        type -- string, = 'int', 'real', 'char', or 'logic'
        n    -- number of elements in vector
        f    -- 'A' -> comment on line above data
             -- 'N' -> no comment 
        s    -- 'R' -> row vector (one line)
             -- 'C' -> column vector (n lines)
        """
        if not n:
            n = len(data)
        if f == 'A':
            self.output_comment(file, comment)
        for i in range(n):
            file.write(self.format_var(type,data[i]))
            if s == 'C':
                file.write("\n")
        if s == 'R':
            file.write("\n")
    
    def output_mat(self, file, type, data, m=None, n=None, 
                   comment=None, s='L',f='A'):
        """
        Output an m x n matrix to file.

        Arguments:
        file -- file object (must be opened for reading)
        type -- string, = 'int' or 'real'
        m    -- number of rows in matrix
        n    -- number of columns (set to m by default)
        f    -- 'A' -> comment on line above data
             -- 'N' -> no comment 
        s    -- 'N' -> full matrix
             -- 'S' -> symmetric
             -- 'L' -> symmetric lower triangular (0 diagonals)
        """
        if not m:
           m = len(data)
        if not n:
           n = m
        if f == 'A':
            self.output_comment(file, comment)
        # Write matrix element values
        if s == 'N' or s == 'S':
            min = 0
        elif s == 'L':
            min =1
        for i in range(min,m) :
            if s == 'N':
                max = n
            elif s == 'S':
                max = i+1
            elif s == 'L':
                max = i
            vec = data[i][:max]
            self.output_vec(file, type, vec, max, None, s='R', f='N')


def strip_quotes(q_string):
    """
    Strip quote marks from a string, if any, and return stripped string.

    Argument:
    q_string -- upstripped string, which may contain final and initial quotes.
    """
    q_string = q_string.strip()
    if q_string[0] == "'" and q_string[-1] == "'":
        return q_string[1:-1]
    elif q_string[0] == '"' and q_string[-1] == '"':
        return q_string[1:-1]
    else:
        return q_string

class IoException(Exception):
    """
    Exception for a file syntax error.
    """

    def __init__(self, message):
        """
        Constructor, which stores a comment.
        """
        super(IoException, self).__init__(message)
