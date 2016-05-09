#!/usr/bin/env python2.7

from pscf.outfile import OutFile

print "Reading file 'in'"
thing = OutFile('in')
print "Writing 'out1'"
thing.write('out1')
print "Reading 'out1'"
thing2 = OutFile('out1')
print "Writing 'out2'"
thing2.write('out2')

print "\nTesting eval method (expression evaluation):\n"
expr = 'block_length[0][0]+block_length[0][1]'
print expr + ' = ' + str(thing.eval(expr))
