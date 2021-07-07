# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 16:10:02 2021

@author: marcu
"""

from sympy import *

kappa, tauq, t, x = symbols('kappa tauq t x')
q = Function('q')(x, t)
T = Function('T')(x, t)
q1 = Function('q1')(x, t)
#q, T = symbols('q T', cls=Function)(x, t)
#qNS = Function('qNS')


#epsilon = symbols('epsilon', positive=True)
#kappa = Function('kappa')(T)
#tauq = Function('tau_q')(T)

dxT = T.diff(x)
qNS = -kappa*dxT
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
