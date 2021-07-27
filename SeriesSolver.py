#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 17:17:47 2021

@author: mjh1n20
"""

import sympy as sp

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
        svs[i][j] = sp.sympify(infile.readline())
        svs[i][j] = sp.simplify(sp.expand(svs[i][j].diff(t)))
    
# Another line of flex before the flux vector        
print(infile.readline())
# Read in flux vector components
for i in range(3):
    for j in range(5):
        for k in range(4):
            fvs[i][j][k] = sp.sympify(infile.readline())
            fvs[i][j][k] = sp.simplify(sp.expand(fvs[i][j][k].diff(X[i])))

# Another line of text before LO CE Correction expressions
print(infile.readline())
for i in range(13):
    LO[i] = sp.sympify(infile.readline())

infile.close()

print(svs[1][1])
print(fvs[0][1][2])