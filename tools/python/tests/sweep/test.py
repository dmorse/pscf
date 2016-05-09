#!/usr/bin/env python2.7

from pscf.sweep import *
sweep = Sweep('out')
print
print sweep[5].chi[0][1]
print
print sweep.write('block_length[0][0]','f_Helmholtz')
