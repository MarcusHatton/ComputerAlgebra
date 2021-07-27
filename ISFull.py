# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 01:09:42 2021

@author: marcu
"""

from sympy import *
import numpy as np
#import pickle

# Setup symbols
kappa, tau_q, zeta, tau_Pi, eta, tau_pi = symbols('kappa tau_q zeta tau_Pi eta tau_pi', real=True, positive=True)
strengths = [kappa,zeta,eta]
# This looks v silly but it makes some loops nicer
timescales = [tau_q, tau_q, tau_q, tau_Pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi]

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
q1, q2, q3 = Function('q1')(t, x, y, z), Function('q2')(t, x, y, z), Function('q3')(t, x, y, z)
qs = [q1, q2, q3]
Pi =  Function('Pi')(t, x, y, z)
pi11, pi12, pi13, pi21, pi22, pi23, pi31, pi32, pi33 = \
    Function('pi11')(t, x, y, z), Function('pi12')(t, x, y, z), Function('pi13')(t, x, y, z), \
    Function('pi21')(t, x, y, z), Function('pi22')(t, x, y, z), Function('pi23')(t, x, y, z), \
    Function('pi31')(t, x, y, z), Function('pi32')(t, x, y, z), Function('pi33')(t, x, y, z)
diss_vars = [q1, q2, q3, Pi, pi11, pi12, pi13, pi21, pi22, pi23, pi31, pi32, pi33]

# Prims
v1, v2, v3 = Function('v1')(t, x, y, z), Function('v2')(t, x, y, z), Function('v3')(t, x, y, z)
vs = [v1, v2, v3]
p, n, rho  = Function('p')(t, x, y, z), Function('n')(t, x, y, z),Function('rho')(t, x, y, z)
prim_vars = [v1, v2, v3, p, n, rho]
                  
# Aux
T, W = Function('T')(t, x, y, z), Function('W')(t, x, y, z)
qv, pi00, pi01, pi02, pi03 = Function('qv')(t, x, y, z), Function('pi00')(t, x, y, z), \
    Function('pi01')(t, x, y, z), Function('pi02')(t, x, y, z), Function('pi03')(t, x, y, z)
aux_vars = [T, W, qv, pi00, pi01, pi02, pi03]

# May choose not to define explicitly for readability
# vsqrd = v1**2 + v2**2 + v3**2
# W = 1/sqrt(1-vsqrd)
qv = q1*v1 + q2*v2 + q3*v3
pi00 = pi11 + pi22 + pi33
pi01 = v1*pi11 + v2*pi12 + v3*pi13
pi02 = v1*pi21 + v2*pi22 + v3*pi23
pi03 = v1*pi31 + v2*pi32 + v3*pi33

# CE LO Corrections
q1LO, q2LO, q3LO = Function('q1LO')(t, x, y, z), Function('q2LO')(t, x, y, z), Function('q3LO')(t, x, y, z)
PiLO =  Function('PiLO')(t, x, y, z)
pi11LO, pi12LO, pi13LO, pi21LO, pi22LO, pi23LO, pi31LO, pi32LO, pi33LO = \
    Function('pi11LO')(t, x, y, z), Function('pi12LO')(t, x, y, z), Function('pi13LO')(t, x, y, z), \
    Function('pi21LO')(t, x, y, z), Function('pi22LO')(t, x, y, z), Function('pi23LO')(t, x, y, z), \
    Function('pi31LO')(t, x, y, z), Function('pi32LO')(t, x, y, z), Function('pi33LO')(t, x, y, z)
diss1s = [q1LO, q2LO, q3LO, PiLO, pi11LO, pi12LO, pi13LO, pi21LO, pi22LO, pi23LO, pi31LO, pi32LO, pi33LO]

# state vector
D = n*W
Sx = (rho + p + Pi)*W**2*v1 + (q1 + qv*v1)*W + pi01
Sy = (rho + p + Pi)*W**2*v2 + (q2 + qv*v2)*W + pi02
Sz = (rho + p + Pi)*W**2*v3 + (q3 + qv*v3)*W + pi03
E = (rho + p + Pi)*W**2 - (p + Pi) + 2*qv*W + pi00
state_vec = [D, Sx, Sy, Sz, E]

# flux vector
flux_vec = np.zeros((3,len(state_vec)),dtype=type(state_vec[0]))
for i in range(3):
    flux_vec[i][0] = state_vec[0]*vs[i]
    flux_vec[i][1] = state_vec[1]*vs[i] + W*(qs[i]*v1 - qv*vs[i]*v1) 
    flux_vec[i][2] = state_vec[2]*vs[i] + W*(qs[i]*v2 - qv*vs[i]*v2) 
    flux_vec[i][3] = state_vec[3]*vs[i] + W*(qs[i]*v3 - qv*vs[i]*v3) 
    flux_vec[i][i+1] += (p + Pi)
    flux_vec[i][4] = (state_vec[4] + p)*vs[i] + W*(qs[i] - qv*vs[i]) 
    for j in range(len(state_vec)):
        flux_vec[i][j] = simplify(expand(flux_vec[i][j]))

# source vector
# Navier-Stokes equ forms
q1NS, q2NS, q3NS = Function('q1NS')(t, x, y, z), Function('q2NS')(t, x, y, z), Function('q3NS')(t, x, y, z)
PiNS =  Function('PiNS')(t, x, y, z)
pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS = \
    Function('pi11NS')(t, x, y, z), Function('pi12NS')(t, x, y, z), Function('pi13NS')(t, x, y, z), \
    Function('pi21NS')(t, x, y, z), Function('pi22NS')(t, x, y, z), Function('pi23NS')(t, x, y, z), \
    Function('pi31NS')(t, x, y, z), Function('pi32NS')(t, x, y, z), Function('pi33NS')(t, x, y, z)
dissNSs = [q1NS, q2NS, q3NS, PiNS, pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS]

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
    diss_eqsLO[i] = simplify(Eq((D*dissNSs[i]).diff(t) + (D*dissNSs[i]*v1).diff(x) + \
                    (D*dissNSs[i]*v2).diff(y) + (D*dissNSs[i]*v3).diff(z), diss_sourcesLO[i]))
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
    IS_sys_series[i][0] = simplify(expand(IS_sys_LO[i].lhs).as_independent(tau_q)[0].as_independent(tau_Pi)[0].as_independent(tau_pi)[0])
    state_vec_series[i][0] = simplify(expand(state_vec_LO[i]).as_independent(tau_q)[0].as_independent(tau_Pi)[0].as_independent(tau_pi)[0])
    for k in range(len(X)):
        flux_vec_series[k][i][0] = simplify(expand(flux_vec_LO[k][i]).as_independent(tau_q)[0].as_independent(tau_Pi)[0].as_independent(tau_pi)[0])
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
#pickle_out = open('ISFullPickleOutput.txt','w')

outfile.write('State Vector (5x4): 5 components, 4 timescale separations \n')
for i in range(5):
    for j in range(4):
        outfile.write(str(state_vec_series[i][j])+'\n')
        #pickle.dump(state_vec_series[i][j],pickle_out)
        
outfile.write('Flux Vector (3x5x4): 3 Directions, 5 components, 4 timescale separations \n')
for i in range(3):
    for j in range(5):
        for k in range(4):
            outfile.write(str(flux_vec_series[i][j][k])+'\n')
            #pickle.dump(flux_vec_series[i][j][k],pickle_out)


outfile.write('First order CE corrections \n')
for i in range(len(diss1sLO)):
    outfile.write(str(diss1sLO[i])+'\n')
    #pickle_out.dump(diss1sLO[i],pickle_out)

outfile.close()






















