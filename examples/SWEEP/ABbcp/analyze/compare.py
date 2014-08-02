#!/usr/bin/env python2.4

from sweep import *
import sys, os

# Prefix (i.e., directory path) for output files
prefix = 'ref/'

# Read lamellar phase
lam = Sweep('../lam/1/ref')
string = lam.write('block_length[0][0]','f_Helmholtz')
file = open(prefix + 'lam','w')
file.write(string)
file.close()

# Read hex phase
hex = Sweep('../hex/1/ref')
string = hex.write('block_length[0][0]','f_Helmholtz')
file = open(prefix + 'hex','w')
file.write(string)
file.close()

# Read gyr phase
gyr = Sweep('../gyr/1/ref')
string = gyr.write('block_length[0][0]','f_Helmholtz')
file = open(prefix + 'gyr','w')
file.write(string)
file.close()

# Read bcc phase
bcc = Sweep('../bcc/1/ref')
string = bcc.write('block_length[0][0]','f_Helmholtz')
file = open(prefix + 'bcc','w')
file.write(string)
file.close()
