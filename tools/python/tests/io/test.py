#!/usr/bin/env python2.4

from io import *
infile  = open('in','r')
outfile = open('out','w')
io = IO()

int_var = io.input_var(infile,'int')
io.output_var(outfile,'int',int_var,'int_var')

real_var = io.input_var(infile,'real')
io.output_var(outfile,'real',real_var,'real_var')

char_var = io.input_var(infile,'char')
io.output_var(outfile,'char',char_var,'char_var')

int_vec_R = io.input_vec(infile,'int')
io.output_vec(outfile,'int',int_vec_R, None, 'int_vec_R')

real_vec_R = io.input_vec(infile,'real')
io.output_vec(outfile,'real',real_vec_R, None, 'real_vec_R')

int_vec_C = io.input_vec(infile,'int', 2, s='C')
io.output_vec(outfile,'int',int_vec_C, 2, 'int_vec_C',s='C')

real_vec_C = io.input_vec(infile,'real', 2, s='C')
io.output_vec(outfile,'real',real_vec_C, 2, 'real_vec_C',s='C')
io.output_vec(outfile,'real',real_vec_C, None, 'real_vec_C',s='C')

real_mat_L = io.input_mat(infile,'real', 2, 2, s='L')
io.output_mat(outfile,'real',real_mat_L, 2, 2, 'real_mat_L',s='L')

real_mat_L = io.input_mat(infile,'real', 3, 3, s='L')
io.output_mat(outfile,'real',real_mat_L, 3, 3, 'real_mat_L',s='L')

