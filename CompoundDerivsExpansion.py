# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 13:27:51 2024

@author: mjh1n20
"""

# dvxdt = -gamma*vx*dEdt/(gamma*W**2*p + gamma*W**2*e - gamma*p - gamma*e + p + e) + (gamma*vx - vx)*dDdt/(gamma*W**3*p + gamma*W**3*e - gamma*W*p - gamma*W*e + W*p + W*e) + dSxdt/(W**2*p + W**2*e);


from sympy import *
import numpy as np

# Setup indep vars
t, x, y, z = symbols('t x y z')
X = [x, y, z]

# Primitives
vx, vy, vz = Function('vx')(t, x, y, z), Function('vy')(t, x, y, z), Function('vz')(t, x, y, z)
vs = [vx, vy, vz]
p, n, e  = Function('p')(t, x, y, z), Function('n')(t, x, y, z),Function('e')(t, x, y, z)
prim_vars = [vx, vy, vz, p, n, e]

# Auxiliary
T, W = Function('T')(t, x, y, z), Function('W')(t, x, y, z)
W = 1 # simplification

dDdt, dSxdt, dSydt, dSzdt, dEdt = Function('dDdt')(t, x, y, z), Function('dSxdt')(t, x, y, z),\
    Function('dSydt')(t, x, y, z), Function('dSzdt')(t, x, y, z), Function('dEdt')(t, x, y, z)

dedt, dpdt, dndt, dvxdt, dvydt, dvzdt  = Function('dedt')(t, x, y, z), Function('dpdt')(t, x, y, z),\
    Function('dndt')(t, x, y, z), Function('dvxdt')(t, x, y, z), Function('dvydt')(t, x, y, z), Function('dvzdt')(t, x, y, z)

# EoS
gamma = symbols('gamma',positive=True,real=True)

# Use EoS to write dtn
dndt = dedt - (1/(gamma-1))*dpdt;
# Now use chain rule to write dtW(dtv1,dtv2,dtv3...)
dWdt = W*W*W*(vx*dvxdt + vy*dvydt + vz*dvzdt);
dTdt = (1/n)*dpdt - (p/(n**2))*dndt;

dpdt = (gamma - 1)*dEdt/(gamma*W**2 - gamma + 1) + (-gamma*W + W)*dDdt/(gamma*W**2 - gamma + 1)
dedt = (gamma*W**2 - gamma - W**2 + 1)*dDdt/(gamma*W**3 - gamma*W + W) + dEdt/(gamma*W**2 - gamma + 1)
dvxdt = -gamma*vx*dEdt/(gamma*W**2*p + gamma*W**2*e - gamma*p - gamma*e + p + e) + (gamma*vx - vx)*dDdt/(gamma*W**3*p + gamma*W**3*e - gamma*W*p - gamma*W*e + W*p + W*e) + dSxdt/(W**2*p + W**2*e)
dvydt = -gamma*vy*dEdt/(gamma*W**2*p + gamma*W**2*e - gamma*p - gamma*e + p + e) + (gamma*vy - vy)*dDdt/(gamma*W**3*p + gamma*W**3*e - gamma*W*p - gamma*W*e + W*p + W*e) + dSydt/(W**2*p + W**2*e)
dvzdt = -gamma*vz*dEdt/(gamma*W**2*p + gamma*W**2*e - gamma*p - gamma*e + p + e) + (gamma*vz - vz)*dDdt/(gamma*W**3*p + gamma*W**3*e - gamma*W*p - gamma*W*e + W*p + W*e) + dSzdt/(W**2*p + W**2*e)

dXdts = {"dpdt" : dpdt, "dedt": dedt, "dvxdt": dvxdt, "dvydt": dvydt, "dvzdt": dvzdt, "dndt": dndt, "dWdt": dWdt, "dTdt": dTdt}
dXs = {"dx": x, "dy": y, "dz": z}

# print(diff(dvxdt, x))

outfile = open('CompoundDerivsExpansion.txt','w')

outfile.write('8 dXdts, 3 deriv directions for each \n')

# for dX in dXs.keys():
#     print(dpdt.diff(dXs[dX]))
#print(dXs.keys())

deriv_replacement_strs = ["vx", "vy", "vz", "p", "n", "e", "dDdt", "dEdt", "dSxdt", "dSydt", "dSzdt", "dndt", "dWdt", "dTdt",
                          "dvxdt", "dvydt", "dvzdt", "dpdt", "dedt"]


for dXdt in dXdts.keys():
    for dX in dXs.keys():
        outfile.write(dX+dXdt+' = '+str(dXdts[dXdt].diff(dXs[dX]))+'\n')

outfile.write('Now write in C++-friendly code\n')

for dXdt in dXdts.keys():
    for dX in dXs.keys():
        original_str = dX+dXdt+' = '+str(dXdts[dXdt].diff(dXs[dX]))
        mod_str = original_str.replace('(t, x, y, z)', "")
        for deriv_str in deriv_replacement_strs:
            mod_str = mod_str.replace('Derivative('+deriv_str+', '+dX[-1]+')', dX+deriv_str)
            mod_str = mod_str.replace(deriv_str+'**2', 'sqr('+deriv_str+')')
        outfile.write(mod_str+";\n")

outfile.close()

















