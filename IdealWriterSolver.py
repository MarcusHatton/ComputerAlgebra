#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 16:15:27 2021

@author: mjh1n20
"""

import sympy as sp
import numpy as np
#import pickle

# Setup symbols
kappa, tau_q, zeta, tau_Pi, eta, tau_pi = sp.symbols('kappa tau_q zeta tau_Pi eta tau_pi', real=True, positive=True)
strengths = [kappa,zeta,eta]
# This looks v silly but it makes some loops nicer
timescales = [tau_q, tau_q, tau_q, tau_Pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi]

# Using epsilon?
epsilon = sp.symbols('epsilon',positive=True,real=True)
epsilon=1
for i in range(len(strengths)):
    strengths[i] *= epsilon

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

# DissNSs - will need to swap Diss vars for these for Jacobian calcs
q1NS, q2NS, q3NS = sp.Function('q1NS')(t, x, y, z), sp.Function('q2NS')(t, x, y, z), sp.Function('q3NS')(t, x, y, z)
PiNS =  sp.Function('PiNS')(t, x, y, z)
pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS = \
    sp.Function('pi11NS')(t, x, y, z), sp.Function('pi12NS')(t, x, y, z), sp.Function('pi13NS')(t, x, y, z), \
    sp.Function('pi21NS')(t, x, y, z), sp.Function('pi22NS')(t, x, y, z), sp.Function('pi23NS')(t, x, y, z), \
    sp.Function('pi31NS')(t, x, y, z), sp.Function('pi32NS')(t, x, y, z), sp.Function('pi33NS')(t, x, y, z)
dissNSs = [q1NS, q2NS, q3NS, PiNS, pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS]

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

# For writing these out to be read in for Jacobian calcs in SeSo.py
state_flux_outfile = open('state_flux_ideal.txt','w')
# state vector
D = n*W
Sx = (rho + p)*W**2*v1
Sy = (rho + p)*W**2*v2
Sz = (rho + p)*W**2*v3
E = (rho + p)*W**2 - p
sv = [D, Sx, Sy, Sz, E]
for i in range(len(sv)):
    state_flux_outfile.write(str(sv[i])+'\n')

# flux vector
fv = np.zeros((3,len(sv)),dtype=type(sv[0]))
for i in range(3):
    fv[i][0] = sv[0]*vs[i]
    fv[i][1] = sv[1]*vs[i]
    fv[i][2] = sv[2]*vs[i]
    fv[i][3] = sv[3]*vs[i]
    fv[i][i+1] += p
    fv[i][4] = (sv[4] + p)*vs[i]
    for j in range(len(sv)):
        fv[i][j] = sp.simplify(sp.expand(fv[i][j]))
        state_flux_outfile.write(str(fv[i][j])+'\n')

state_flux_outfile.close()

# Declare the vars that we need to calculate the Jacobian of the state vector wrt
# (Essentially the variables that have time derivatives in the CE expansion)
jac_vars = [p, n, v1, v2, v3]

# Convert the state and flux vectors into Matrices so that the Jacobian
# function can be used 
jac_vars_Mat = sp.Matrix(jac_vars)

# Convert state vector into sympy Matrix
sv_Mat = sp.Matrix(sv)
# Calculate Jacobian of state vector wrt jac_vars
sv_Jac = sv_Mat.jacobian(jac_vars_Mat).inv()
# Now do the same for each of the flux vector components (x, y, z) in turn
# Not actually sure we need these yet...
fvx_Mat = sp.Matrix(fv[0])
fvy_Mat = sp.Matrix(fv[1])
fvz_Mat = sp.Matrix(fv[2])
fvx_Jac = fvx_Mat.jacobian(jac_vars_Mat)
fvy_Jac = fvy_Mat.jacobian(jac_vars_Mat)
fvz_Jac = fvz_Mat.jacobian(jac_vars_Mat)

TURN INTO A MATRIX CALC NOT SUM!

# Form a list of the time derivatives of the required variables
dt_jac_vars = np.zeros_like(jac_vars,dtype=type(jac_vars[0]))
for i in range(len(jac_vars)):
    #dt_jac_vars[i] = jac_vars[i].diff(t)
    for j in range(len(cons)):
        # If the derivative of the conserved wrt the jac_var (n, W, etc.)
        # is zero then its reciprocal should be set to zero (not inf)
        #if sv_Jac[i+j*len(jac_vars)] == 0:
        #    continue
        # Little numbering hack picks out all the required partial derivs
        # for each of the conserveds e.g. dn/dt = dD/dt*(dD/dn)^-1 + dS1/dt*dn/dS1 + ... + dTau/dt*dn/dTau
        dt_jac_vars[i] += cons[j].diff(t)*(sv_Jac[i+j*len(jac_vars)])

# Define NS vector forms (H & F)
svNS = np.zeros_like(sv,dtype=type(sv[0]))
svNS[0] = 0
svNS[1] = Pi*W**2*v1 + (q1 + qv*v1)*W + pi01
svNS[2] = Pi*W**2*v2 + (q2 + qv*v2)*W + pi02
svNS[3] = Pi*W**2*v3 + (q3 + qv*v3)*W + pi03
svNS[4] = Pi*(W**2 - 1) + 2*qv*W + pi00

fvNS = np.zeros_like(fv,dtype=type(sv[0]))
for i in range(3):
    fvNS[i][0] = svNS[0]*vs[i]
    fvNS[i][1] = svNS[1]*vs[i] + W*(qs[i]*v1 - qv*vs[i]*v1) 
    fvNS[i][2] = svNS[2]*vs[i] + W*(qs[i]*v2 - qv*vs[i]*v2) 
    fvNS[i][3] = svNS[3]*vs[i] + W*(qs[i]*v3 - qv*vs[i]*v3) 
    fvNS[i][i+1] += Pi
    fvNS[i][4] = svNS[4]*vs[i] + W*(qs[i] - qv*vs[i]) 
    for j in range(len(sv)):
        fvNS[i][j] = sp.simplify(sp.expand(fv[i][j]))
















