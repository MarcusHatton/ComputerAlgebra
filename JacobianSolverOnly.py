# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 17:39:57 2021

@author: mjh1n20
"""

import sympy as sp
import numpy as np

t, x, y, z = sp.symbols('t x y z')

D, S1, S2, S3, E = sp.Function('D')(t, x, y, z), sp.Function('S1')(t, x, y, z), \
    sp.Function('S2')(t, x, y, z), sp.Function('S3')(t, x, y, z), sp.Function('E')(t, x, y, z)
cons = [D, S1, S2, S3, E]
dtcons = [D.diff(t), S1.diff(t), S2.diff(t), S3.diff(t), E.diff(t)]
# Prims
v1, v2, v3 = sp.Function('v1')(t, x, y, z), sp.Function('v2')(t, x, y, z), sp.Function('v3')(t, x, y, z)
vs = [v1, v2, v3]
W = sp.Function('W')(t, x, y, z)

n = sp.Function('n')(t, x, y, z)
rho = sp.Function('rho')(t, x, y, z)
p = sp.Function('p')(t, x, y, z)
#p = sp.Function('p')(rho, n)
Gamma = sp.symbols('Gamma')
prim_vars = [v1, v2, v3, p, n, rho]

# state vector simplest
# D = n
# S1 = (rho + p)*v1
# S2 = (rho + p)*v2
# S3 = (rho + p)*v3
# E = rho

# better
D = (rho - p/(Gamma-1))*W # n is substituted here
S1 = (rho + p)*v1*W**2
S2 = (rho + p)*v2*W**2
S3 = (rho + p)*v3*W**2
E = (rho + p)*W**2 - p

sv = [D, S1, S2, S3, E]

# Declare the vars that we need to calculate the Jacobian of the state vector wrt
# (Essentially the variables that have time derivatives in the CE expansion)
jac_vars = [p, rho, v1, v2, v3]
#jac_vars = [W, v1, v2, v3]

# Convert the state and flux vectors into Matrices so that the Jacobian
# function can be used 
jac_vars_Mat = sp.Matrix(jac_vars)

# Convert state vector into sympy Matrix
sv_Mat = sp.Matrix(sv)
# Calculate Jacobian of state vector wrt jac_vars
sv_Jac = sv_Mat.jacobian(jac_vars_Mat)
sv_Jac_inv = sv_Jac.inv()

#jac_sol_simple_outfile = open('jac_sol_simple.txt','w')
jac_sol_full_outfile = open('jac_sol_full.txt','w')

for i in range(len(jac_vars)):
    print(str(jac_vars[i].diff(t))+': \n')
    print((sv_Jac_inv*sp.Matrix(dtcons))[i])
    print('\n')
    jac_sol_full_outfile.write(str((sv_Jac_inv*sp.Matrix(dtcons))[i])+'\n')

#jac_sol_simple_outfile.close()
jac_sol_full_outfile.close()

