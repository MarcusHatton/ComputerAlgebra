# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 00:54:44 2021

@author: marcu
"""

from sympy import *

kappa, tauq, zeta, tauPi = symbols('kappa tauq zeta tauPi')
strengths = [kappa,zeta]
timescales = [tauq, tauPi]
t, x = symbols('t x')
q, Pi = Function('q')(x, t), Function('Pi')(x, t)
T, v = Function('T')(x, t), Function('v')(x, t)
q1, Pi1 = Function('q1')(x, t), Function('Pi1')(x, t)

fasts = [q, Pi]
slows = [T, v]
fast1s = [q1, Pi1]

dxT = T.diff(x)
dxv = v.diff(x)
qNS = -kappa*dxT
PiNS = -zeta*dxv
fastNSs = [qNS, PiNS]
for i in len(fasts):
    fast[i] = fastNSs[i] + timescales[i]*fast1s[i]
    
q = qNS + tauq*q1
dtq = q.diff(t)

eq10b =  Eq(-dtq, (1/tauq)*(-qNS + q))
eq10bLO = eq10b.subs(tauq,0)
#print(eq10b)

q1 = solve(eq10bLO,q1)[0]
#print(q1)

q = qNS + tauq*q1
dxq = q.diff(x)
dtT = T.diff(t)
eq10a = Eq(dtT + dxq, 0)
expr10a = dtT + dxq
eq10aLO = eq10a.subs(tauq,0)
#print(expr10a)

dtTLO = solve(eq10aLO,dtT)[0]
#print(dtTLO)

eq10a = eq10a.subs(dtT,dtTLO)
#print(simplify(eq10a))

dtTNLO = simplify(expr10a.subs(dtT,dtTLO))
#print(simplify(dtTNLO))

dtTsol = dtTLO + dtTNLO
print(dtTsol)

#dtTsol = solve(eq10a,dtT)[0]
#print(dtTsol)
