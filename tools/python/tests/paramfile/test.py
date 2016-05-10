#!/usr/bin/env python2.7

from pscf.paramfile import ParamFile

print "Reading file 'in'"
thing = ParamFile('in')
print "Writing 'out1'"
thing.write('out1')
print "Reading 'out1'"
thing2 = ParamFile('out1')
print "Writing 'out2'"
thing2.write('out2')

print "\nTesting eval method (expression evaluation):\n"
expr = 'block_length[0][0]+block_length[0][1]'
print expr + ' = ' + str(thing.eval(expr))
