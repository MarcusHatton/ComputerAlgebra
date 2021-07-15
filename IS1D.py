# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 16:27:53 2021

@author: marcu
"""

from sympy import *
import numpy as np

# Setup symbols
kappa, tauq, zeta, tauPi, eta, taupi = symbols('kappa tauq zeta tauPi eta taupi', real=True)
strengths = [kappa,zeta,eta]
timescales = [tauq, tauPi, taupi]

# Using epsilon?
epsilon = symbols('epsilon',positive=True,real=True)
epsilon=1
for i in range(len(strengths)):
    strengths[i] *= epsilon

# Setup indep vars
t, x, y, z = symbols('t x y z')

# Dependent vars
# Diss vars
q, Pi, pi = Function('q')(x, t), Function('Pi')(x, t), Function('pi')(x, t)
diss_vars = [q, Pi, pi]
# Prims
v, p, n, rho, T, W = Function('v')(x, t), Function('p')(x, t), Function('n')(x, t), \
                  Function('rho')(x, t), Function('T')(x, t), Function('W')(x, t), 

# Aux
qv, pi00, pi01 = Function('qv')(x, t), Function('pi00')(x, t), Function('pi01')(x, t)
vsqrd = v**2
# May choose not to define explicitly for readability
#W = 1/sqrt(1-v**2)
#qv = q*v
#pi00 = pi
#pi01 = v*pi

# CE LO Corrections
q1, Pi1, pi1 = Function('q1')(x, t), Function('Pi1')(x, t), Function('pi1')(x, t)
diss1s = [q1, Pi1, pi1]

# state vector
D = n*W
S = (rho + p + Pi)*W**2*v + (q + qv*v)*W + pi01
E = (rho + p + Pi)*W**2 - (p + Pi) + 2*qv*W + pi00
state_vec = [D, S, E]

# flux vector
flux_vec = np.zeros_like(state_vec)
flux_vec[0] = state_vec[0]*v
flux_vec[1] = state_vec[1]*v + (p + Pi) + W*(q*v - qv*v*v) + pi
flux_vec[2] = (state_vec[2] + p)*v + W*(q - qv*v) + pi01

# source vector
# Navier-Stokes equ forms
qNS, PiNS, piNS = Function('qNS')(x, t), Function('PiNS')(x, t), Function('piNS')(x, t)
diss_NSs = np.zeros_like(diss_vars)
diss_NSs = [qNS, PiNS, piNS]
#diss_NSs[0] = -strengths[0]*T.diff(x)
#diss_NSs[1] = -strengths[1]*v.diff(x)
#diss_NSs[2] = -2*strengths[2]*(4/3)*v.diff(x)
source_vec = np.zeros_like(state_vec)
diss_sources = np.zeros_like(diss_vars)
for i in range(len(diss1s)):
    diss_sources[i] = (n/timescales[i])*(diss_NSs[i] - diss_vars[i])
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
    diss_sourcesLO[i] = diss_sources[i].subs(diss_vars[i],diss_NSs[i] + timescales[i]*diss1s[i])
    diss_eqsLO[i] = simplify(Eq((D*diss_NSs[i]).diff(t) + (D*diss_NSs[i]*v).diff(x), diss_sourcesLO[i]))
    #diss_eqsLO[i] = diss_eqs[i].subs(diss_vars[i],diss_NSs[i] + timescales[i]*diss1s[i])
    #diss_eqsLO[i] = diss_eqs[i].subs(timescales[i],0)
    diss1sLO[i] = simplify(solve(diss_eqsLO[i],diss1s[i])[0])
    #print(diss1sLO[i])
    #print('\n')     

# Copy these so the substitution loop works
state_vec_LO = state_vec
flux_vec_LO = flux_vec
for i in range(len(IS_sys)):
    # Untouched IS system
    IS_sys[i] = simplify(Eq(state_vec[i].diff(t) + flux_vec[i].diff(x),source_vec[i]))
    # Substitute CE expansion
    for j in range(len(diss_vars)):
        #print(state_vec_LO[i])
        #print('\n')
        # substitute uncalculated 1s here for viewing clarity
        state_vec_LO[i] = simplify(state_vec_LO[i].subs(diss_vars[j],diss_NSs[j] + timescales[j]*diss1s[j]))# + timescales[j]*diss1sLO[j])
        flux_vec_LO[i] = simplify(flux_vec_LO[i].subs(diss_vars[j],diss_NSs[j] + timescales[j]*diss1s[j]))#+ timescales[j]*diss1sLO[j])
    # Form the equation with the expansion
    IS_sys_LO[i] = simplify(Eq(state_vec_LO[i].diff(t) + flux_vec_LO[i].diff(x),0))

# Separate out the zeroth order parts of the system and each timescale order
IS_sys_series = np.zeros((len(IS_sys),len(timescales)+1),dtype=type(IS_sys_LO[0]))
for i in range(len(IS_sys)):
    IS_sys_series[i][0] = simplify(IS_sys_LO[i].lhs.as_independent(timescales[0])[0].as_independent(timescales[1])[0].as_independent(timescales[2])[0])
    #print(IS_sys_series[i][0])
    for j in range(len(diss_vars)):
    #    IS_sys_LO[i] = IS_sys_LO[i].subs(timescales[j],0)
        IS_sys_series[i][j+1] = simplify(expand(IS_sys_LO[i].lhs).as_independent(timescales[j])[1])

for i in range(3):
    for j in range(4):
        print(IS_sys_series[i][j])
        print('\n')
    print('\n')










