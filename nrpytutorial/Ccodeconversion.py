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
        svs_str[i][j] = str(infile.readline())
        
print(infile.readline())
for i in range(5):
    for j in range(4):
        for k in range(3):
            fvs_str[k][i][j] = str(infile.readline())

infile.close()

# Use string replacement to get C++ style expressions for the state & flux vectors
func_of = '(t, x, y, z)'
components = ', i, j, k)]'
prim_vars = ['vx', 'vy', 'vz', 'p', 'n', 'rho']
diss_vars = ['qx', 'qy', 'qz', 'Pi', 'pixx', 'pixy', 'pixz', 'piyx', 'piyy', 'piyz', 'pizx', 'pizy', 'pizz']
diss1s = ['qx1', 'qy1', 'qz1', 'Pi1', 'pixx1', 'pixy1', 'pixz1', 'piyx1', 'piyy1', 'piyz1', 'pizx1', 'pizy1', 'pizz1']
aux_vars = ['T', 'W', 'qv', 'pitt', 'pitx', 'pity', 'pitz']
dissNSs = ['qxNS', 'qyNS', 'qzNS', 'PiNS', 'pixxNS', 'pixyNS', 'pixzNS', 'piyxNS', 'piyyNS', 'piyzNS', 'pizxNS', 'pizyNS', 'pizzNS']
for i in range(5):
    for j in range(4):
        #print(svs_str[i][j])
        for prim_var in prim_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
            for k in range(3):
                fvs_str[k][i][j] = svs_str[i][j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)

        for diss_var in diss_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)
            for k in range(3):
                fvs_str[k][i][j] = svs_str[i][j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)

        for aux_var in aux_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(aux_var)+func_of,'prims[ID(Prims::'+str(aux_var)+components)
            for k in range(3):
                fvs_str[k][i][j] = svs_str[i][j].replace(str(aux_var)+func_of,'prims[ID(Prims::'+str(aux_var)+components)

state_vec = ['D', 'Sx', 'Sy', 'Sz', 'E']

outfile = open('CcodeOutput.txt','w')

outfile.write('STATE VECTOR STARTS\n')
for i in range(5):
    outfile.write('\n'+state_vec[i]+'\n')
    for j in range(4):
        outfile.write(svs_str[i][j])

outfile.write('\nFLUX VECTOR STARTS \n')
for i in range(3):
    for j in range(5):
        for k in range(4):
            outfile.write(fvs_str[i][j][k])
        outfile.write('\n')
    outfile.write('\n')

outfile.close()