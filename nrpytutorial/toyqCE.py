# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 16:10:02 2021

@author: marcu
"""

from sympy import *

from outputC import lhrh         # NRPy+: Core C code output module
import NRPy_param_funcs as par   # NRPy+: parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import finite_difference as fin  # NRPy+: Finite difference C code generation moduleimport finite_difference as fd


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

#dtTsol = solve(eq10a,dtT)[0]
print(dtTsol)
#print(str(dtTsol).replace('T','U'))

#fdcoeffs, fdstencil = fd.compute_fdcoeffs_fdstencl("dDD01",FDORDER=4)
#print(fdcoeffs, fdstencil)

### NrPy stuff starts ###

# Set the spatial dimension to 1
par.set_paramsvals_value("grid::DIM = 1")

# Register the input gridfunction "phi" and the gridfunction to which data are output, "output":
T, dtT_NrPy = gri.register_gridfunctions("AUX",["T","dtT_NrPy"])

# Declare phi_dDD as a rank-2 indexed expression: phi_dDD[i][j] = \partial_i \partial_j phi
T_dDD = ixp.declarerank2("T_dDD","nosym")
#U_dDD = ixp.declarerank2("U_dDD","nosym")

# Set output to \partial_0^2 phi
dtT_NrPy = T_dDD[0][0] #+ U_dDD[0][0]

# Output to the screen the core C code for evaluating the finite difference derivative
fin.FD_outputC("stdout",lhrh(lhs=gri.gfaccess("out_gf","dtT_NrPy"),rhs=dtT_NrPy))