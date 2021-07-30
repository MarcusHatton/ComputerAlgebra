#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 17:17:47 2021

@author: mjh1n20
"""

import sympy as sp
import numpy as np

sv = ["" for i in range(5)]
fv = [["" for i in range(5)] for j in range(3)]

# Read in original state and flux vector forms (undifferentiated)
sf_infile = open('state_flux.txt','r')
for i in range(5):
    sv[i] = sp.sympify(sf_infile.readline())
for i in range(3):
    for j in range(5):
        fv[i][j] = sp.sympify(sf_infile.readline())
sf_infile.close()

# Setup indep vars
t, x, y, z = sp.symbols('t x y z')
X = [x, y, z]

# Conserveds just for clarity
D, S1, S2, S3, E = sp.Function('D')(t, x, y, z), sp.Function('S1')(t, x, y, z), \
    sp.Function('S2')(t, x, y, z), sp.Function('S3')(t, x, y, z), sp.Function('E')(t, x, y, z)
cons = [D, S1, S2, S3, E]

# Dependent vars
# Diss vars
q1, q2, q3 = sp.Function('q1')(t, x, y, z), sp.Function('q2')(t, x, y, z), sp.Function('q3')(t, x, y, z)
qs = [q1, q2, q3]
Pi =  sp.Function('Pi')(t, x, y, z)
pi11, pi12, pi13, pi21, pi22, pi23, pi31, pi32, pi33 = \
    sp.Function('pi11')(t, x, y, z), sp.Function('pi12')(t, x, y, z), sp.Function('pi13')(t, x, y, z), \
    sp.Function('pi21')(t, x, y, z), sp.Function('pi22')(t, x, y, z), sp.Function('pi23')(t, x, y, z), \
    sp.Function('pi31')(t, x, y, z), sp.Function('pi32')(t, x, y, z), sp.Function('pi33')(t, x, y, z)
diss_vars = [q1, q2, q3, Pi, pi11, pi12, pi13, pi21, pi22, pi23, pi31, pi32, pi33]

# Prims
v1, v2, v3 = sp.Function('v1')(t, x, y, z), sp.Function('v2')(t, x, y, z), sp.Function('v3')(t, x, y, z)
vs = [v1, v2, v3]
p, n, rho  = sp.Function('p')(t, x, y, z), sp.Function('n')(t, x, y, z),sp.Function('rho')(t, x, y, z)
prim_vars = [v1, v2, v3, p, n, rho]
                  
# Aux
T, W = sp.Function('T')(t, x, y, z), sp.Function('W')(t, x, y, z)
qv, pi00, pi01, pi02, pi03 = sp.Function('qv')(t, x, y, z), sp.Function('pi00')(t, x, y, z), \
    sp.Function('pi01')(t, x, y, z), sp.Function('pi02')(t, x, y, z), sp.Function('pi03')(t, x, y, z)
aux_vars = [T, W, qv, pi00, pi01, pi02, pi03]

# Declare the vars that we need to calculate the Jacobian of the state vector wrt
# (Essentially the variables that have time derivatives in the CE expansion)
jac_vars = [n, W]
for diss_var in diss_vars:
    jac_vars.append(diss_var)

# Convert the state and flux vectors into Matrices so that the Jacobian
# function can be used 
jac_vars_Mat = sp.Matrix(jac_vars)

# Convert state vector into sympy Matrix
sv_Mat = sp.Matrix(sv)
# Calculate Jacobian of state vector wrt jac_vars
sv_Jac = sv_Mat.jacobian(jac_vars_Mat)
# Now do the same for each of the flux vector components (x, y, z) in turn
# Not actually sure we need these yet...
fvx_Mat = sp.Matrix(fv[0])
fvy_Mat = sp.Matrix(fv[1])
fvz_Mat = sp.Matrix(fv[2])
fvx_Jac = fvx_Mat.jacobian(jac_vars_Mat)
fvy_Jac = fvy_Mat.jacobian(jac_vars_Mat)
fvz_Jac = fvz_Mat.jacobian(jac_vars_Mat)

# Form a list of the time derivatives of the required variables
dt_jac_vars = np.zeros_like(jac_vars,dtype=type(jac_vars[0]))
for i in range(len(jac_vars)):
    #dt_jac_vars[i] = jac_vars[i].diff(t)
    for j in range(len(cons)):
        # If the derivative of the conserved wrt the jac_var (n, W, etc.)
        # is zero then its reciprocal should be set to zero (not inf)
        if sv_Jac[i+j*len(jac_vars)] == 0:
            continue
        # Little numbering hack picks out all the required partial derivs
        # for each of the conserveds e.g. dn/dt = dD/dt*dn/dD + dS1/dt*dn/dS1 + ...
        dt_jac_vars[i] += cons[j].diff(t)*(1/sv_Jac[i+j*len(jac_vars)])




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
        svs[i][j] = sp.sympify(infile.readline())
        #svs[i][j] = sp.simplify(sp.expand(svs[i][j].diff(t)))
        svs[i][j] = svs[i][j].diff(t)
    
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
