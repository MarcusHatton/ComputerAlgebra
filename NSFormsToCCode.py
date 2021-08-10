# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 00:56:38 2021

@author: marcu
"""

from sympy import *
import numpy as np

# Read in the algebra in sympy form from file
svs_str = ["" for j in range(5)]
fvs_str = [["" for j in range(5)] for k in range(3)]

infile = open('state_flux_NS.txt','r')
#infile = open('state_NS_dt.txt','r')

# Read in state vector components
for i in range(5):
    svs_str[i] = str(infile.readline())
    
# Read in flux vector components
for i in range(3):
    for j in range(5):
        fvs_str[i][j] = str(infile.readline())

infile.close()

# Use string replacement to get C++ style expressions for the state & flux vectors
func_of = '(t, x, y, z)'
components = ', i, j, k)]'
prim_vars = ['v1', 'v2', 'v3', 'p', 'n', 'rho']
diss_vars = ['q1', 'q2', 'q3', 'Pi', 'pi11', 'pi12', 'pi13', 'pi21', 'pi22', 'pi23', 'pi31', 'pi32', 'pi33']
diss1s = ['q1LO', 'q2LO', 'q3LO', 'PiLO', 'pi11LO', 'pi12LO', 'pi13LO', 'pi21LO', 'pi22LO', 'pi23LO', 'pi31LO', 'pi32LO', 'pi33LO']
aux_vars = ['T', 'W', 'qv', 'pi00', 'pi01', 'pi02', 'pi03']
dissNSs = ['q1NS', 'q2NS', 'q3NS', 'PiNS', 'pi11NS', 'pi12NS', 'pi13NS', 'pi21NS', 'pi22NS', 'pi23NS', 'pi31NS', 'pi32NS', 'pi33NS']
for i in range(3):
    for j in range(5):
        #print(svs_str[i])
        for prim_var in prim_vars:
            svs_str[j] = svs_str[j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
            svs_str[j] = svs_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(prim_var)+', i, j, k)]')
            fvs_str[i][j] = fvs_str[i][j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
            NS_str = ["" for j in range(13)]

            
        for diss_var in diss_vars:
            svs_str[j] = svs_str[j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)
            fvs_str[i][j] = fvs_str[i][j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)

        for aux_var in aux_vars:
            svs_str[j] = svs_str[j].replace(str(aux_var)+func_of,'aux[ID(Aux::'+str(aux_var)+components)
            svs_str[j] = svs_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(aux_var)+', i, j, k)]')
            fvs_str[i][j] = fvs_str[i][j].replace(str(aux_var)+func_of,'aux[ID(Aux::'+str(aux_var)+components)

        for dissNS in dissNSs:
            svs_str[j] = svs_str[j].replace(str(dissNS)+func_of,'aux[ID(Aux::'+str(dissNS)+components)
            svs_str[j] = svs_str[j].replace('Derivative(aux[ID(Aux::'+str(dissNS)+components+', t)', \
                                      'tderivs[ID(TDerivs::dt'+str(dissNS)+', i, j, k)]')
            fvs_str[i][j] = fvs_str[i][j].replace(str(dissNS)+func_of,'aux[ID(Aux::'+str(dissNS)+components)

        for diss1 in diss1s:
            svs_str[j] = svs_str[j].replace(str(diss1)+func_of,'aux[ID(Aux::'+str(diss1)+components)
            fvs_str[i][j] = fvs_str[i][j].replace(str(diss1)+func_of,'aux[ID(Aux::'+str(diss1)+components)

state_vec = ['D', 'Sx', 'Sy', 'Sz', 'E']

outfile = open('NSinCCode.txt','w')

outfile.write('STATE VECTOR STARTS\n')
for i in range(5):
    #outfile.write('\n'+state_vec[i]+'\n')
    outfile.write(svs_str[i].replace('aux[ID(Aux::W, i, j, k)]**2','sqr(aux[ID(Aux::W, i, j, k)])')\
                      .replace('prims[ID(Prims::vx, i, j, k)]**2','sqr(prims[ID(Prims::vx, i, j, k)])')\
                      .replace('prims[ID(Prims::vy, i, j, k)]**2','sqr(prims[ID(Prims::vy, i, j, k)])')\
                      .replace('prims[ID(Prims::vz, i, j, k)]**2','sqr(prims[ID(Prims::vz, i, j, k)])')\
                      .replace('21','12').replace('31','13').replace('32','23'))
            

dirs = ['x', 'y', 'z']
outfile.write('\nFLUX VECTOR STARTS \n')
for i in range(3):
    outfile.write(dirs[i]+'\n')
    for j in range(5):
        #outfile.write(state_vec[j]+'\n')
        outfile.write(fvs_str[i][j].replace('aux[ID(Aux::W, i, j, k)]**2','sqr(aux[ID(Aux::W, i, j, k)])')\
                      .replace('prims[ID(Prims::vx, i, j, k)]**2','sqr(prims[ID(Prims::vx, i, j, k)])')\
                      .replace('prims[ID(Prims::vy, i, j, k)]**2','sqr(prims[ID(Prims::vy, i, j, k)])')\
                      .replace('prims[ID(Prims::vz, i, j, k)]**2','sqr(prims[ID(Prims::vz, i, j, k)])')\
                      .replace('21','12').replace('31','13').replace('32','23'))
    outfile.write('\n')

outfile.close()





