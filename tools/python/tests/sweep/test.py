#!/usr/bin/env python2.4

from sweep import *
sweep = Sweep('ABbcp_lam_1')
print
print sweep[5].chi[0][1]
print
print sweep.write('block_length[0][0]','f_Helmholtz')
