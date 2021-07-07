# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 16:10:02 2021

@author: marcu
"""

from sympy import *

kappa, tauq, t, x = symbols('kappa tauq t x')
q = Function('q')(x, t)
T = Function('T')(x, t)
#q, T = symbols('q T', cls=Function)(x, t)
#qNS = Function('qNS')

#kappa = 1
#tauq = 0

dxT = T.diff(x)

dtq = q.diff(t)

qNS = -kappa*dxT

eq10b =  Eq(dtq, 0)

if (tauq != 0):
    eq10b = Eq(dtq, (1/tauq)*(qNS - q))
    qsol = solve(eq10b,q)
    q = qsol[0]
else:
    q = qNS # limit case
    
dtT = T.diff(t)
dxq = q.diff(x)

eq10a = Eq(dtT + dxq, 0)

dtTsol = solve(eq10a,dtT)
print(dtTsol)

