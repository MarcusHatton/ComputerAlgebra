# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 15:55:43 2021

@author: mjh1n20
"""

from sympy import *
import numpy as np

num_lines = 5

# Read in the algebra in sympy form from file
NS_str = ["" for j in range(num_lines)]

infile = open('BDNK_simple.txt','r')
#infile = open('jac_sol_simple.txt','r')

# Read in state vector components
for i in range(num_lines):
    NS_str[i] = str(infile.readline())

infile.close()

# Use string replacement to get C++ style expressions for the state & flux vectors
func_of = '(t, x, y, z)'
components = ', i, j, k)]'
prim_vars = ['v1', 'v2', 'v3', 'p', 'n', 'rho']
diss_vars = ['q1', 'q2', 'q3', 'Pi', 'pi11', 'pi12', 'pi13', 'pi21', 'pi22', 'pi23', 'pi31', 'pi32', 'pi33']
diss1s = ['q1LO', 'q2LO', 'q3LO', 'PiLO', 'pi11LO', 'pi12LO', 'pi13LO', 'pi21LO', 'pi22LO', 'pi23LO', 'pi31LO', 'pi32LO', 'pi33LO']
aux_vars = ['T', 'W', 'qv', 'pi00', 'pi01', 'pi02', 'pi03']
dissNSs = ['q1NS', 'q2NS', 'q3NS', 'PiNS', 'pi11NS', 'pi12NS', 'pi13NS', 'pi21NS', 'pi22NS', 'pi23NS', 'pi31NS', 'pi32NS', 'pi33NS']
cons_vars = ['D', 'S1', 'S2', 'S3', 'E']
all_prim_vars = prim_vars + diss_vars + cons_vars # pretend
all_aux_vars = aux_vars + dissNSs + diss1s
for j in range(num_lines):
    for prim_var in all_prim_vars:
        NS_str[j] = NS_str[j].replace(str(prim_var)+func_of,'prims[ID(Prims::'+str(prim_var)+components)
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', t)', \
                                      'aux[ID(Aux::dt'+str(prim_var)+', i, j, k)]')
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', x)', \
                                      '((prims[ID(Prims::'+str(prim_var)+', i+1, j, k)] - prims[ID(Prims::'+str(prim_var)+', i-1, j, k)])/(2*d->dx))')
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', y)', \
                                      '((prims[ID(Prims::'+str(prim_var)+', i, j+1, k)] - prims[ID(Prims::'+str(prim_var)+', i, j-1, k)])/(2*d->dy))')
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', z)', \
                                      '((prims[ID(Prims::'+str(prim_var)+', i, j, k+1)] - prims[ID(Prims::'+str(prim_var)+', i, j, k-1)])/(2*d->dz))')

        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', t, x)', \
                                      '((aux[ID(Aux::dt'+str(prim_var)+', i+1, j, k)] - aux[ID(Aux::dt'+str(prim_var)+', i-1, j, k)])/(2*d->dx))')            
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', t, y)', \
                                      '((aux[ID(Aux::dt'+str(prim_var)+', i, j+1, k)] - aux[ID(Aux::dt'+str(prim_var)+', i, j-1, k)])/(2*d->dy))')            
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', t, z)', \
                                      '((aux[ID(Aux::dt'+str(prim_var)+', i, j, k+1)] - aux[ID(Aux::dt'+str(prim_var)+', i, j, k-1)])/(2*d->dz))')            
        NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', (t, 2))', str(0)) 
                    # not great simply setting second time derivs to zero...
            
    # for diss_var in diss_vars:
    #     NS_str[j] = NS_str[j].replace(str(diss_var)+func_of,'prims[ID(Prims::'+str(diss_var)+components)
    #     NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', x)', \
    #                                   '((prims[ID(Prims::'+str(prim_var)+', i+1, j, k)] - prims[ID(Prims::'+str(prim_var)+', i-1, j, k)])/(d->dx))')
    #     NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', y)', \
    #                                   '((prims[ID(Prims::'+str(prim_var)+', i, j+1, k)] - prims[ID(Prims::'+str(prim_var)+', i, j-1, k)])/(d->dy))')
    #     NS_str[j] = NS_str[j].replace('Derivative(prims[ID(Prims::'+str(prim_var)+components+', z)', \
                                      
    for aux_var in all_aux_vars:
        NS_str[j] = NS_str[j].replace(str(aux_var)+func_of,'aux[ID(Aux::'+str(aux_var)+components)
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', t)', \
                                      'aux[ID(Aux::dt'+str(aux_var)+', i, j, k)]')
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', x)', \
                                      '((aux[ID(Aux::'+str(aux_var)+', i+1, j, k)] - aux[ID(Aux::'+str(aux_var)+', i-1, j, k)])/(d->dx))')
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', y)', \
                                      '((aux[ID(Aux::'+str(aux_var)+', i, j+1, k)] - aux[ID(Aux::'+str(aux_var)+', i, j-1, k)])/(d->dy))')
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', z)', \
                                      '((aux[ID(Aux::'+str(aux_var)+', i, j, k+1)] - aux[ID(Aux::'+str(aux_var)+', i, j, k-1)])/(d->dz))')

        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', t, x)', \
                                      '((aux[ID(Aux::dt'+str(aux_var)+', i+1, j, k)] - aux[ID(Aux::dt'+str(aux_var)+', i-1, j, k)])/(d->dx))')            
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', t, y)', \
                                      '((aux[ID(Aux::dt'+str(aux_var)+', i, j+1, k)] - aux[ID(Aux::dt'+str(aux_var)+', i, j-1, k)])/(d->dy))')            
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', t, z)', \
                                      '((aux[ID(Aux::dt'+str(aux_var)+', i, j, k+1)] - aux[ID(Aux::dt'+str(aux_var)+', i, j, k-1)])/(d->dz))')            
        NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(aux_var)+components+', (t, 2))', str(0))



    # for dissNS in dissNSs:
    #     NS_str[j] = NS_str[j].replace(str(dissNS)+func_of,'aux[ID(Aux::'+str(dissNS)+components)
    #     NS_str[j] = NS_str[j].replace('Derivative(aux[ID(Aux::'+str(dissNS)+components+', t)', \
    #                               'tderivs[ID(TDerivs::dt'+str(dissNS)+', i, j, k)]')

    # for diss1 in diss1s:
    #     NS_str[j] = NS_str[j].replace(str(diss1)+func_of,'aux[ID(Aux::'+str(diss1)+components)


outfile = open('BDNKInCCode.txt','w')

for i in range(num_lines):
    outfile.write(NS_str[i].replace('aux[ID(Aux::W, i, j, k)]**2','sqr(aux[ID(Aux::W, i, j, k)])')\
                  .replace('aux[ID(Aux::W, i, j, k)]**3','aux[ID(Aux::W, i, j, k)]*sqr(aux[ID(Aux::W, i, j, k)])')
                  .replace('prims[ID(Prims::v1, i, j, k)]**2','sqr(prims[ID(Prims::v1, i, j, k)])')\
                  .replace('prims[ID(Prims::v2, i, j, k)]**2','sqr(prims[ID(Prims::v2, i, j, k)])')\
                  .replace('prims[ID(Prims::v3, i, j, k)]**2','sqr(prims[ID(Prims::v3, i, j, k)])')\
                  .replace('21','12').replace('31','13').replace('32','23')\
                  .replace('Gamma','d->gamma'))
outfile.close()