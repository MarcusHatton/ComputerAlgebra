#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 17:17:47 2021

@author: mjh1n20
"""

import sympy as sp
import numpy as np

# Reads in the tau_X series expansion for the state & flux vectors, 
# as well as the LO forms of the diss variables and does the temporal
# and spatial derivatives on the stae & flux components, respectively
svs = [["" for i in range(4)] for j in range(5)]
fvs = [[["" for i in range(4)] for j in range(5)] for k in range(3)]
LO = ["" for i in range(13)]
infile = open('ISFullOutput.txt','r')

t, x, y, z = sp.symbols('t x y z')
X = [x, y, z]

# Read first line which is text
print(infile.readline())
# Read in state vector components
for i in range(5):
    for j in range(4):
        svs[i][j] = sp.sympify(infile.readline().replace('LO',''))
        #svs[i][j] = sp.simplify(sp.expand(svs[i][j].diff(t)))
        svs[i][j] = sp.simplify(sp.expand(svs[i][j].diff(t)))
    
# Another line of flex before the flux vector        
print(infile.readline())
# Read in flux vector components
for i in range(3):
    for j in range(5):
        for k in range(4):
            fvs[i][j][k] = sp.sympify(infile.readline())
            #fvs[i][j][k] = sp.simplify(sp.expand(fvs[i][j][k].diff(X[i])))
            fvs[i][j][k] = fvs[i][j][k].diff(X[i])

# Another line of text before LO CE Correction expressions
print(infile.readline())
for i in range(13):
    LO[i] = sp.sympify(infile.readline())

infile.close()


# for i in range(5):
#     for j in range(1,4):
#         for k in range(len(jac_vars)):
#             if svs[i][j] == 0:
#                 continue
#             #print(svs[i][j])
#             #print(jac_vars[k].diff(t))
#             svs[i][j] = svs[i][j].subs(jac_vars[k].diff(t),dt_jac_vars[k])

"""
Outstanding issue: For Jacobian solving, expressions for state and flux vectors
are read in from state_flux.txt ... where its q1, q2, q3 etc. but then for the
substitution into the time-derived and timescale-expanded forms (svs, fvs) its
q1LO, q2LO etc. ... probably a string replacement is easiest
"""




