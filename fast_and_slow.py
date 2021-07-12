# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 00:54:44 2021

@author: marcu
"""

from sympy import *
import numpy as np

kappa, tauq, zeta, tauPi = symbols('kappa tauq zeta tauPi', real=True)
epsilon = symbols('epsilon',positive=True,real=True)
epsilon=1
strengths = [kappa,zeta]
for i in range(len(strengths)):
    strengths[i] *= epsilon
timescales = [tauq, tauPi]
t, x = symbols('t x')
q, Pi = Function('q')(x, t), Function('Pi')(x, t)
T, v = Function('T')(x, t), Function('v')(x, t)
q1, Pi1 = Function('q1')(x, t), Function('Pi1')(x, t)

fasts = [q, Pi]
slows = [T, v]
fast1s = [q1, Pi1]

dxslows = np.zeros_like(slows)
fastNSs = np.zeros_like(fasts)

for i in range(len(fasts)):
    dxslows[i] = slows[i].diff(x)
    fastNSs[i] = -strengths[i]*dxslows[i]
    fasts[i] = fastNSs[i] + timescales[i]*fast1s[i]


dtfasts = np.zeros_like(fasts)
dxfasts = np.zeros_like(fasts)

for i in range(len(fasts)):
    dtfasts[i] = fasts[i].diff(t)
    dxfasts[i] = fasts[i].diff(x)

beqs = np.zeros_like(fasts)
beqsLO = np.zeros_like(fasts)
bfluxes = [0, fasts[1]*dxslows[1] + slows[1]*dxfasts[1]]
dtslows = np.zeros_like(slows)
for i in range(len(fasts)):
    beqs[i] = Eq(dtfasts[i] + bfluxes[i], (1/timescales[i])*(fastNSs[i] - fasts[i]))
    beqsLO[i] = beqs[i].subs(timescales[i],0)
    #beqsLO[i] = beqs[i].subs(epsilon,0)
    fast1s[i] = solve(beqsLO[i],fast1s[i])[0]
    fasts[i] = fastNSs[i] + timescales[i]*fast1s[i]
    dxfasts[i] = fasts[i].diff(x)
    dtslows[i] = slows[i].diff(t)

#print(fast1s)

aeqs = np.zeros_like(fasts)
aexprs = np.zeros_like(fasts)
aeqsLO = np.zeros_like(fasts)
afluxes = [dxfasts[0], 2*slows[1]*dxslows[1] + dxfasts[1]]
dtslowsLO = np.zeros_like(fasts)
dtslowsNLO = np.zeros_like(fasts)
dtslowsols = np.zeros_like(fasts)

for i in range(len(fasts)):
    aeqs[i] = Eq(dtslows[i] + afluxes[i], 0)
    aexprs[i] = dtslows[i] + afluxes[i]
    aeqsLO[i] = aeqs[i].subs(timescales[i], 0)
    dtslowsLO[i] = solve(aeqsLO[i], dtslows[i])[0]
    #print(aeqs[i])
    aeqs[i] = aeqs[i].subs(dtslows[i],dtslowsLO[i])
    dtslowsNLO[i] = simplify(aexprs[i].subs(dtslows[i],dtslowsLO[i]))
    dtslowsols[i] = dtslowsLO[i] + dtslowsNLO[i]

print(dtslowsols)
#print(expand(dtslowsols[1]).as_ordered_terms())


