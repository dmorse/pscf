
# Copyright (c) 2002 Trent Mick
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""Test preprocessing of inputs/... with preprocess.py."""

import sys
import os
import unittest
import difflib
import pprint

import testsupport



#----- test cases

class PreprocessInputsTestCase(unittest.TestCase):
    def setUp(self):
        if not os.path.exists('tmp'):
            os.mkdir('tmp')

    def tearDown(self):
        testsupport.rmtree('tmp')


def _testOneInputFile(self, fname):
    import preprocess

    infile = os.path.join('inputs', fname) # input
    reffile = os.path.join('outputs', fname) # expected output
    outfile = os.path.join('tmp', fname) # actual output
    errfile = os.path.join('outputs', fname+'.err')  # expected error
    optsfile = os.path.join('inputs', fname+'.opts') # input options

    # Determine input options to use, if any.
    opts = {}
    if os.path.exists(optsfile):
        for line in open(optsfile, 'r').readlines():
            if line[-1] == "\n": line = line[:-1]
            name, value = line.split('=', 1)
            try:
                value = eval(value)
            except NameError:
                pass
            opts[name] = value
        #print "options from '%s': %s" % (optsfile, pprint.pformat(opts))

    # If there is no reference output file this means that processing
    # this file is expected to fail.
    if os.path.exists(reffile):
        ref = open(reffile, 'r').readlines()
        if not sys.platform.startswith("win"):
            ref = [line.replace('\\','/') for line in ref] # use Un*x paths
        preprocess.preprocess(infile, outfile, **opts)
        out = open(outfile, 'r').readlines()
        if ref != out:
            diff = list(difflib.ndiff(ref, out))
            self.fail("%r != %r:\n%s"\
                      % (reffile, outfile, pprint.pformat(diff)))
    elif os.path.exists(errfile):
        err = open(errfile, 'r').read()
        if not sys.platform.startswith("win"):
            err = err.replace('\\','/') # use Un*x paths
        try:
            preprocess.preprocess(infile, outfile, **opts)
        except preprocess.PreprocessError, ex:
            #print "XXX ex: %s" % str(ex).strip()
            self.failUnlessEqual(err.strip(), str(ex).strip())
        else:
            self.fail("No PreprocessError when expected one.")
    else:
        self.fail("No reference or error file for '%s'." % infile)

    # Ensure next test file gets a clean preprocess.
    del sys.modules['preprocess']
        

for fname in os.listdir('inputs'):
    if fname.endswith(".opts"): continue # skip input option files
    testFunction = lambda self, fname=fname: _testOneInputFile(self, fname)
    name = 'test_'+fname
    setattr(PreprocessInputsTestCase, name, testFunction)


#---- mainline

def suite():
    """Return a unittest.TestSuite to be used by test.py."""
    return unittest.makeSuite(PreprocessInputsTestCase)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(sys.stdout, verbosity=2)
    result = runner.run(suite())

