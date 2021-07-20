#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 16:56:56 2021

@author: mjh1n20
"""

from sympy import *
import numpy as np

# Read in the algebra in sympy form from file
svs_str = [["" for i in range(4)] for j in range(5)]
fvs_str = [[["" for i in range(4)] for j in range(5)] for k in range(3)]
infile = open('ISFullOutput.txt','r')
print(infile.readline())

for i in range(5):
    for j in range(4):
        #print(infile.readline())
        #print(type(infile.readline()))
        svs_str[i][j] = str(infile.readline())
        
print(infile.readline())
for i in range(5):
    for j in range(4):
        for k in range(3):
            fvs_str[k][i][j] = str(infile.readline())
infile.close()

#print(fvs_str)

# Use string replacement to get C++ style expressions for the state & flux vectors
func_of = '(t, x, y, z)'
components = ', i, j, k)]'
prim_vars = ['vx', 'vy', 'vz', 'p', 'n', 'rho']
for i in range(5):
    for j in range(4):
        #print(svs_str[i][j])
        for prim_var in prim_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
        print(svs_str[i][j])