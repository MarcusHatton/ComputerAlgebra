# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 01:09:42 2021

@author: marcu
"""

from sympy import *
import numpy as np

from outputC import lhrh         # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import finite_difference as fin  # NRPy+: Finite difference C code generation moduleimport finite_difference as fd

# Setup symbols
kappa, tauq, zeta, tauPi, eta, taupi = symbols('kappa tauq zeta tauPi eta taupi', real=True, positive=True)
strengths = [kappa,zeta,eta]
# This looks v silly but it makes some loops nicer
timescales = [tauq, tauq, tauq, tauPi, taupi, taupi, taupi, taupi, taupi, taupi, taupi, taupi, taupi]

# Using epsilon?
epsilon = symbols('epsilon',positive=True,real=True)
epsilon=1
for i in range(len(strengths)):
    strengths[i] *= epsilon

# Setup indep vars
t, x, y, z = symbols('t x y z')
X = [x, y, z]

# Dependent vars
# Diss vars
qx, qy, qz = Function('qx')(t, x, y, z), Function('qy')(t, x, y, z), Function('qz')(t, x, y, z)
qs = [qx, qy, qz]
Pi =  Function('Pi')(t, x, y, z)
pixx, pixy, pixz, piyx, piyy, piyz, pizx, pizy, pizz = \
    Function('pixx')(t, x, y, z), Function('pixy')(t, x, y, z), Function('pixz')(t, x, y, z), \
    Function('piyx')(t, x, y, z), Function('piyy')(t, x, y, z), Function('piyz')(t, x, y, z), \
    Function('pizx')(t, x, y, z), Function('pizy')(t, x, y, z), Function('pizz')(t, x, y, z)
diss_vars = [qx, qy, qz, Pi, pixx, pixy, pixz, piyx, piyy, piyz, pizx, pizy, pizz]

# Prims
vx, vy, vz = Function('vx')(t, x, y, z), Function('vy')(t, x, y, z), Function('vz')(t, x, y, z)
vs = [vx, vy, vz]
p, n, rho  = Function('p')(t, x, y, z), Function('n')(t, x, y, z),Function('rho')(t, x, y, z)
prim_vars = [vx, vy, vz, p, n, rho]
                  
# Aux
T, W = Function('T')(t, x, y, z), Function('W')(t, x, y, z)
qv, pitt, pitx, pity, pitz = Function('qv')(t, x, y, z), Function('pitt')(t, x, y, z), \
    Function('pitx')(t, x, y, z), Function('pity')(t, x, y, z), Function('pitz')(t, x, y, z)
aux_vars = [T, W, qv, pitt, pitx, pity, pitz]

# May choose not to define explicitly for readability
# vsqrd = vx**2 + vy**2 + vz**2
# W = 1/sqrt(1-vsqrd)
# qv = qx*vx + qy*vy + qz*vz
# pitt = pixx + piyy + pizz
# pitx = vx*pixx + vy*pixy + vz*pixz
# pity = vx*piyz + vy*piyy + vz*piyz
# pitz = vx*pizx + vy*pizy + vz*pizz

# CE LO Corrections
qx1, qy1, qz1 = Function('qx1')(t, x, y, z), Function('qy1')(t, x, y, z), Function('qz1')(t, x, y, z)
Pi1 =  Function('Pi1')(t, x, y, z)
pixx1, pixy1, pixz1, piyx1, piyy1, piyz1, pizx1, pizy1, pizz1 = \
    Function('pixx1')(t, x, y, z), Function('pixy1')(t, x, y, z), Function('pixz1')(t, x, y, z), \
    Function('piyx1')(t, x, y, z), Function('piyy1')(t, x, y, z), Function('piyz1')(t, x, y, z), \
    Function('pizx1')(t, x, y, z), Function('pizy1')(t, x, y, z), Function('pizz1')(t, x, y, z)
diss1s = [qx1, qy1, qz1, Pi1, pixx1, pixy1, pixz1, piyx1, piyy1, piyz1, pizx1, pizy1, pizz1]

# state vector
D = n*W
Sx = (rho + p + Pi)*W**2*vx + (qx + qv*vx)*W + pitx
Sy = (rho + p + Pi)*W**2*vy + (qy + qv*vy)*W + pity
Sz = (rho + p + Pi)*W**2*vz + (qz + qv*vz)*W + pitz
E = (rho + p + Pi)*W**2 - (p + Pi) + 2*qv*W + pitt
state_vec = [D, Sx, Sy, Sz, E]

# flux vector
flux_vec = np.zeros((3,len(state_vec)),dtype=type(state_vec[0]))
for i in range(3):
    flux_vec[i][0] = state_vec[0]*vs[i]
    flux_vec[i][1] = state_vec[1]*vs[i] + W*(qs[i]*vx - qv*vs[i]*vx) 
    flux_vec[i][2] = state_vec[2]*vs[i] + W*(qs[i]*vy - qv*vs[i]*vy) 
    flux_vec[i][3] = state_vec[3]*vs[i] + W*(qs[i]*vz - qv*vs[i]*vz) 
    flux_vec[i][i+1] += (p + Pi)
    flux_vec[i][4] = (state_vec[4] + p)*vs[i] + W*(qs[i] - qv*vs[i]) 
    for j in range(len(state_vec)):
        flux_vec[i][j] = simplify(expand(flux_vec[i][j]))

# source vector
# Navier-Stokes equ forms
qxNS, qyNS, qzNS = Function('qxNS')(t, x, y, z), Function('qyNS')(t, x, y, z), Function('qzNS')(t, x, y, z)
PiNS =  Function('PiNS')(t, x, y, z)
pixxNS, pixyNS, pixzNS, piyxNS, piyyNS, piyzNS, pizxNS, pizyNS, pizzNS = \
    Function('pixxNS')(t, x, y, z), Function('pixyNS')(t, x, y, z), Function('pixzNS')(t, x, y, z), \
    Function('piyxNS')(t, x, y, z), Function('piyyNS')(t, x, y, z), Function('piyzNS')(t, x, y, z), \
    Function('pizxNS')(t, x, y, z), Function('pizyNS')(t, x, y, z), Function('pizzNS')(t, x, y, z)
dissNSs = [qxNS, qyNS, qzNS, PiNS, pixxNS, pixyNS, pixzNS, piyxNS, piyyNS, piyzNS, pizxNS, pizyNS, pizzNS]

# Will need to define these at some point...
#dissNSs[0] = -strengths[0]*T.diff(x)
#dissNSs[1] = -strengths[1]*v.diff(x)
#dissNSs[2] = -2*strengths[2]*(4/3)*v.diff(x)

source_vec = np.zeros_like(state_vec)
diss_sources = np.zeros_like(diss_vars)
for i in range(len(diss1s)):
    diss_sources[i] = (n/timescales[i])*(dissNSs[i] - diss_vars[i])
    #print(diss_sources[i])
    #print('\n')

# Assemble full system
IS_sys = np.zeros_like(state_vec)
IS_sys_LO = np.zeros_like(state_vec)
source_vec_LO = source_vec
diss_eqsLO = np.zeros_like(diss_sources)
diss1sLO = np.zeros_like(diss1s)
diss_sourcesLO = np.zeros_like(diss_sources)
for i in range(len(diss_sources)):
    # Calc diss1s to LO
    diss_sourcesLO[i] = diss_sources[i].subs(diss_vars[i],dissNSs[i] + timescales[i]*diss1s[i])
    diss_eqsLO[i] = simplify(Eq((D*dissNSs[i]).diff(t) + (D*dissNSs[i]*vx).diff(x) + \
                    (D*dissNSs[i]*vy).diff(y) + (D*dissNSs[i]*vz).diff(z), diss_sourcesLO[i]))
    #diss_eqsLO[i] = diss_eqs[i].subs(diss_vars[i],dissNSs[i] + timescales[i]*diss1s[i])
    #diss_eqsLO[i] = diss_eqs[i].subs(timescales[i],0)
    diss1sLO[i] = simplify(solve(diss_eqsLO[i],diss1s[i])[0])
    #print(diss1sLO[i])
    #print('\n')     

# Copy these so the substitution loop works
state_vec_LO = state_vec
flux_vec_LO = flux_vec
for i in range(len(IS_sys)):
    # Untouched IS system
    IS_sys[i] = simplify(Eq(state_vec[i].diff(t) + flux_vec[0][i].diff(x) \
                    + flux_vec[1][i].diff(y) + flux_vec[2][i].diff(z), source_vec[i]))
    # Substitute CE expansion
    for j in range(len(diss_vars)):
        state_vec_LO[i] = simplify(state_vec_LO[i].subs(diss_vars[j],dissNSs[j] + timescales[j]*diss1s[j]))# + timescales[j]*diss1sLO[j])
        #print(state_vec_LO[i])
        #print('\n')
        for k in range(len(X)):
            # substitute uncalculated 1s here for viewing clarity
            flux_vec_LO[k][i] = simplify(flux_vec_LO[k][i].subs(diss_vars[j],dissNSs[j] + timescales[j]*diss1s[j]))#+ timescales[j]*diss1sLO[j])
    # Form the equation with the expansion
    IS_sys_LO[i] = simplify(Eq(state_vec_LO[i].diff(t) + flux_vec_LO[0][i].diff(x) \
                    + flux_vec_LO[1][i].diff(y) + flux_vec_LO[2][i].diff(z), 0))

# Separate out the zeroth order parts of the system and each timescale order
# 4 for the 3 timescales plus one for not one of them (zeroth order)
IS_sys_series = np.zeros((len(IS_sys),4),dtype=type(IS_sys_LO[0]))
state_vec_series = np.zeros((len(state_vec_LO),4),dtype=type(state_vec_LO[0]))
flux_vec_series = np.zeros(shape=(flux_vec_LO.shape[0],flux_vec_LO.shape[1],4),dtype=type(IS_sys_LO[0]))
for i in range(len(IS_sys)):
    IS_sys_series[i][0] = simplify(expand(IS_sys_LO[i].lhs).as_independent(tauq)[0].as_independent(tauPi)[0].as_independent(taupi)[0])
    state_vec_series[i][0] = simplify(expand(state_vec_LO[i]).as_independent(tauq)[0].as_independent(tauPi)[0].as_independent(taupi)[0])
    for k in range(len(X)):
        flux_vec_series[k][i][0] = simplify(expand(flux_vec_LO[k][i]).as_independent(tauq)[0].as_independent(tauPi)[0].as_independent(taupi)[0])
    #print(IS_sys_series[i][0])
    for j in range(3):
    #    IS_sys_LO[i] = IS_sys_LO[i].subs(timescales[j],0)
        # +2 hack ensures one of each timescale
        IS_sys_series[i][j+1] = simplify(expand(IS_sys_LO[i].lhs).as_independent(timescales[j+2])[1])
        state_vec_series[i][j+1] = simplify(expand(state_vec_LO[i]).as_independent(timescales[j+2])[1])
        # a catch for the 1 returned when the timescale is not present in that component (should be zero therefore)
        if (state_vec_series[i][j+1] == 1):
            IS_sys_series[i][j+1] = 0
        if (state_vec_series[i][j+1] == 1):
            state_vec_series[i][j+1] = 0
        for k in range(len(X)):
            flux_vec_series[k][i][j+1] = simplify(expand(flux_vec_LO[k][i]).as_independent(timescales[j+2])[1])
            if(flux_vec_series[k][i][j+1] == 1):
                flux_vec_series[k][i][j+1] = 0

# for i in range(5): # D, Sx, Sy, Sz, E
#     for j in range(4): # tau^0, tau^q, tau^Pi, tau^pi
#         print(i, j)
#         #print(IS_sys_series[i][j])
#         print("state vector")
#         print(state_vec_series[i][j])
#         print("flux vectors")
#         for k in range(len(X)):
#             print(k)
#             print(flux_vec_series[k][i][j])
#         print('\n')
#     print('\n')


# Write the un-differentiated state and flux vectors to file
outfile = open('ISFullOutput.txt','w')

outfile.write('State Vector (5x4): 5 components, 4 timescale separations \n')
for i in range(5):
    for j in range(4):
        outfile.write(str(state_vec_series[i][j])+'\n')
        
outfile.write('Flux Vector (3x5x4): 3 Directions, 5 components, 4 timescale separations \n')
for i in range(3):
    for j in range(5):
        for k in range(4):
            outfile.write(str(flux_vec_series[i][j][k])+'\n')

outfile.write('First order CE corrections \n')
for i in range(len(diss1sLO)):
    outfile.write(str(diss1sLO[i])+'\n')

outfile.close()


