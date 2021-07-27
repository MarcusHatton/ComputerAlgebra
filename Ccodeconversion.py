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
LO_str = ["" for i in range(13)]
infile = open('ISFullOutput.txt','r')

# Read first line which is text
print(infile.readline())
# Read in state vector components
for i in range(5):
    for j in range(4):
        svs_str[i][j] = str(infile.readline())
    
# Another line of flex before the flux vector        
print(infile.readline())
# Read in flux vector components
for i in range(3):
    for j in range(5):
        for k in range(4):
            fvs_str[i][j][k] = str(infile.readline())

# Another line of text before LO CE Correction expressions
print(infile.readline())
for i in range(13):
    LO_str[i] = str(infile.readline())

infile.close()

# Use string replacement to get C++ style expressions for the state & flux vectors
func_of = '(t, x, y, z)'
components = ', i, j, k)]'
prim_vars = ['v1', 'v2', 'v3', 'p', 'n', 'rho']
diss_vars = ['q1', 'q2', 'q3', 'Pi', 'pi11', 'pi12', 'pi13', 'pi21', 'pi22', 'pi23', 'pi31', 'pi32', 'pi33']
diss1s = ['q1LO', 'q2LO', 'q3LO', 'PiLO', 'pi11LO', 'pi12LO', 'pi13LO', 'pi21LO', 'pi22LO', 'pi23LO', 'pi31LO', 'pi32LO', 'pi33LO']
aux_vars = ['T', 'W', 'qv', 'pi00', 'pi01', 'pi02', 'pi03']
dissNSs = ['q1NS', 'q2NS', 'q3NS', 'PiNS', 'pi11NS', 'pi12NS', 'pi13NS', 'pi21NS', 'pi22NS', 'pi23NS', 'pi31NS', 'pi32NS', 'pi33NS']
for i in range(5):
    for j in range(4):
        #print(svs_str[i][j])
        for prim_var in prim_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
            for k in range(3):
                fvs_str[k][i][j] = fvs_str[k][i][j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
            
            
        for diss_var in diss_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)
            for k in range(3):
                fvs_str[k][i][j] = fvs_str[k][i][j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)

        for aux_var in aux_vars:
            svs_str[i][j] = svs_str[i][j].replace(str(aux_var)+func_of,'aux[ID(Aux::'+str(aux_var)+components)
            for k in range(3):
                fvs_str[k][i][j] = fvs_str[k][i][j].replace(str(aux_var)+func_of,'aux[ID(Aux::'+str(aux_var)+components)

        for dissNS in dissNSs:
            svs_str[i][j] = svs_str[i][j].replace(str(dissNS)+func_of,'aux[ID(Aux::'+str(dissNS)+components)
            for k in range(3):
                fvs_str[k][i][j] = fvs_str[k][i][j].replace(str(dissNS)+func_of,'aux[ID(Aux::'+str(dissNS)+components)

        for diss1 in diss1s:
            svs_str[i][j] = svs_str[i][j].replace(str(diss1)+func_of,'aux[ID(Aux::'+str(diss1)+components)
            for k in range(3):
                fvs_str[k][i][j] = fvs_str[k][i][j].replace(str(diss1)+func_of,'aux[ID(Aux::'+str(diss1)+components)

for i in range(13):
    for prim_var in prim_vars:
        LO_str[i] = LO_str[i].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', x)', \
                                      '((prims[ID(Prims::'+str(prim_var)+', i+1, j, k)] - prims[ID(Prims::'+str(prim_var)+', i-1, j, k)])/(d->dx))')
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', y)', \
                                      '((prims[ID(Prims::'+str(prim_var)+', i, j+1, k)] - prims[ID(Prims::'+str(prim_var)+', i, j-1, k)])/(d->dy))')
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', z)', \
                                      '((prims[ID(Prims::'+str(prim_var)+', i, j, k+1)] - prims[ID(Prims::'+str(prim_var)+', i, j, k-1)])/(d->dz))')
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(prim_var)+', i, j, k)]')


    for diss_var in diss_vars:
        LO_str[i] = LO_str[i].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(diss_var)+components+', x)', \
                                      '((prims[ID(Prims::'+str(diss_var)+', i+1, j, k)] - prims[ID(Prims::'+str(diss_var)+', i-1, j, k)])/(d->dx))')
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(diss_var)+components+', y)', \
                                      '((prims[ID(Prims::'+str(diss_var)+', i, j+1, k)] - prims[ID(Prims::'+str(diss_var)+', i, j-1, k)])/(d->dy))')
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(diss_var)+components+', z)', \
                                      '((prims[ID(Prims::'+str(diss_var)+', i, j, k+1)] - prims[ID(Prims::'+str(diss_var)+', i, j, k-1)])/(d->dz))')
        LO_str[i] = LO_str[i].replace('Derivative(prims[ID(Prims::'+str(diss_var)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(diss_var)+', i, j, k)]')

    for aux_var in aux_vars:
        LO_str[i] = LO_str[i].replace(str(aux_var)+func_of,'aux[ID(Aux::'+str(aux_var)+components)
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', x)', \
                                      '((aux[ID(Aux::'+str(aux_var)+', i+1, j, k)] - aux[ID(Aux::'+str(aux_var)+', i-1, j, k)])/(d->dx))')
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', y)', \
                                      '((aux[ID(Aux::'+str(aux_var)+', i, j+1, k)] - aux[ID(Aux::'+str(aux_var)+', i, j-1, k)])/(d->dy))')
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', z)', \
                                      '((aux[ID(Aux::'+str(aux_var)+', i, j, k+1)] - aux[ID(Aux::'+str(aux_var)+', i, j, k-1)])/(d->dz))')
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(aux_var)+', i, j, k)]')
    for dissNS in dissNSs:
        LO_str[i] = LO_str[i].replace(str(dissNS)+func_of,'aux[ID(Aux::'+str(dissNS)+components)
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(dissNS)+components+', x)', \
                                      '((aux[ID(Aux::'+str(dissNS)+', i+1, j, k)] - aux[ID(Aux::'+str(dissNS)+', i-1, j, k)])/(d->dx))')
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(dissNS)+components+', y)', \
                                      '((aux[ID(Aux::'+str(dissNS)+', i, j+1, k)] - aux[ID(Aux::'+str(dissNS)+', i, j-1, k)])/(d->dy))')
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(dissNS)+components+', z)', \
                                      '((aux[ID(Aux::'+str(dissNS)+', i, j, k+1)] - aux[ID(Aux::'+str(dissNS)+', i, j, k-1)])/(d->dz))')
        LO_str[i] = LO_str[i].replace('Derivative(aux[ID(Aux::'+str(dissNS)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(dissNS)+', i, j, k)]')


state_vec = ['D', 'Sx', 'Sy', 'Sz', 'E']

outfile = open('CcodeOutput.txt','w')

outfile.write('STATE VECTOR STARTS\n')
for i in range(5):
    outfile.write('\n'+state_vec[i]+'\n')
    for j in range(4):
        outfile.write(svs_str[i][j].replace('aux[ID(Aux::W, i, j, k)]**2','sqr(aux[ID(Aux::W, i, j, k)])')\
                      .replace('prims[ID(Prims::vx, i, j, k)]**2','sqr(prims[ID(Prims::vx, i, j, k)])')\
                      .replace('prims[ID(Prims::vy, i, j, k)]**2','sqr(prims[ID(Prims::vy, i, j, k)])')\
                      .replace('prims[ID(Prims::vz, i, j, k)]**2','sqr(prims[ID(Prims::vz, i, j, k)])'))
            

dirs = ['x', 'y', 'z']
outfile.write('\nFLUX VECTOR STARTS \n')
for i in range(3):
    outfile.write(dirs[i]+'\n')
    for j in range(5):
        outfile.write(state_vec[j]+'\n')
        for k in range(4):
            outfile.write(fvs_str[i][j][k].replace('aux[ID(Aux::W, i, j, k)]**2','sqr(aux[ID(Aux::W, i, j, k)])')\
                      .replace('prims[ID(Prims::vx, i, j, k)]**2','sqr(prims[ID(Prims::vx, i, j, k)])')\
                      .replace('prims[ID(Prims::vy, i, j, k)]**2','sqr(prims[ID(Prims::vy, i, j, k)])')\
                      .replace('prims[ID(Prims::vz, i, j, k)]**2','sqr(prims[ID(Prims::vz, i, j, k)])'))
        outfile.write('\n')
    outfile.write('\n')

outfile.write('\nDiss LO Corrections START'+'\n')
for i in range(13):
    outfile.write('aux[ID(Aux::'+diss1s[i]+components+' = '+LO_str[i].replace('aux[ID(Aux::W, i, j, k)]**2','sqr(aux[ID(Aux::W, i, j, k)])')\
                      .replace('prims[ID(Prims::vx, i, j, k)]**2','sqr(prims[ID(Prims::vx, i, j, k)])')\
                      .replace('prims[ID(Prims::vy, i, j, k)]**2','sqr(prims[ID(Prims::vy, i, j, k)])')\
                      .replace('prims[ID(Prims::vz, i, j, k)]**2','sqr(prims[ID(Prims::vz, i, j, k)])'))
    outfile.write('\n')

outfile.close()





