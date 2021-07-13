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
qNS, PiNS, piNS = Function('qNS')(x, t), Function('PiNS')(x, t), Function('piNS')(x, t)
diss_vars = [q, Pi, pi]
# Prims
v, p, n, rho, T, W = Function('v')(x, t), Function('p')(x, t), Function('n')(x, t), \
                  Function('rho')(x, t), Function('T')(x, t), Function('W')(x, t), 
# Aux
vsqrd = v**2
#W = 1/sqrt(1-v**2)
qv = q*v
pi00 = pi
pi01 = v*pi
# CE LO Corrections
q1, Pi1, pi1 = Function('q1')(x, t), Function('Pi1')(x, t), Function('pi1')(x, t)
diss1s = [q1, Pi1, pi1]

# state vector
D = n*W
S = (rho + p + Pi)*W**2*v + (q + qv*v)*W + pi01
E = (rho + p + Pi)*W**2 - (p + Pi) + 2*qv*W + pi00
#Y = D*q
#U = D*Pi
#Z = D*pi
state_vec = [D, S, E]

# flux vector
flux_vec = np.zeros_like(state_vec)
flux_vec[0] = state_vec[0]*v
flux_vec[1] = state_vec[1]*v + (p + Pi) + W*(q*v - qv*v*v) + pi
flux_vec[2] = (state_vec[2] + p)*v + W*(q - qv*v) + pi01
#flux_vec[3] = state_vec[3]*v 
#flux_vec[4] = state_vec[4]*v
#flux_vec[5] = state_vec[5]*v

# source vector
# Navier-Stokes equ forms
diss_NSs = np.zeros_like(diss_vars)
diss_NSs = [qNS, PiNS, piNS]
#diss_NSs[0] = -strengths[0]*T.diff(x)
#diss_NSs[1] = -strengths[1]*v.diff(x)
#diss_NSs[2] = -2*strengths[2]*(4/3)*v.diff(x)
source_vec = np.zeros_like(state_vec)
diss_sources = np.zeros_like(diss_vars)
for i in range(len(diss1s)):
    diss_sources[i] = (n/timescales[i])*(diss_NSs[i] - diss_vars[i])

# Assemble full system
IS_sys = np.zeros_like(state_vec)
IS_sys_LO = np.zeros_like(state_vec)
state_vec_LO = state_vec
flux_vec_LO = flux_vec
source_vec_LO = source_vec
diss_eqs = np.zeros_like(diss_sources)
diss1sLO = np.zeros_like(diss1s)
for i in range(len(diss_sources)):
    diss_eqs[i] = Eq((D*diss_vars[i]).diff(t) + (D*diss_vars[i]*v).diff(x), diss_sources[i])
    diss_eqs[i] = diss_eqs[i].subs(diss_vars[i],diss_NSs[i] + timescales[i]*diss1s[i])
    diss_eqs[i] = diss_eqs[i].subs(timescales[i],0)
    diss1sLO[i] = solve(diss_eqs[i],diss1s[i])[0]
    #print(diss1sLO[i])
    #print('\n')     

for i in range(len(IS_sys)):
    IS_sys[i] = Eq(state_vec[i].diff(t) + flux_vec[i].diff(x),source_vec[i])
    # Substitute CE expansion
    for j in range(len(diss_vars)):
        state_vec_LO[i] = state_vec_LO[i].subs(diss_vars[j],diss_NSs[j] + timescales[j]*diss1sLO[j])
        flux_vec_LO[i] = flux_vec_LO[i].subs(diss_vars[j],diss_NSs[j] + timescales[j]*diss1sLO[j])
    # Form the equation with the expansion
    IS_sys_LO[i] = Eq(state_vec_LO[i].diff(t) + flux_vec_LO[i].diff(x),0)
    print(IS_sys_LO[i])
    print('\n')

for i in range(len(IS_sys)):
    for j in range(len(diss_vars)):
        IS_sys_LO[i] = IS_sys_LO[i].subs(timescales[j],0)












