# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 01:11:51 2021

@author: marcu
"""

import sympy as sp
import numpy as np
#import pickle

# Setup symbols
kappa, tau_q, zeta, tau_Pi, eta, tau_pi = sp.symbols('kappa tau_q zeta tau_Pi eta tau_pi', real=True, positive=True)
strengths = [kappa,zeta,eta]
# This looks v silly but it makes some loops nicer
timescales = [tau_q, tau_q, tau_q, tau_Pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi, tau_pi]

# Using epsilon?
epsilon = sp.symbols('epsilon',positive=True,real=True)
epsilon=1
for i in range(len(strengths)):
    strengths[i] *= epsilon

# Setup indep vars
t, x, y, z = sp.symbols('t x y z')
X = [x, y, z]

# Conserveds just for clarity
D, S1, S2, S3, E = sp.Function('D')(t, x, y, z), sp.Function('S1')(t, x, y, z), \
    sp.Function('S2')(t, x, y, z), sp.Function('S3')(t, x, y, z), sp.Function('E')(t, x, y, z)
cons = [D, S1, S2, S3, E]
dtcons = [D.diff(t), S1.diff(t), S2.diff(t), S3.diff(t), E.diff(t)]

# Dependent vars
# Diss vars
q1, q2, q3 = sp.Function('q1')(t, x, y, z), sp.Function('q2')(t, x, y, z), sp.Function('q3')(t, x, y, z)
qs = [q1, q2, q3]
Pi =  sp.Function('Pi')(t, x, y, z)
pi11, pi12, pi13, pi21, pi22, pi23, pi31, pi32, pi33 = \
    sp.Function('pi11')(t, x, y, z), sp.Function('pi12')(t, x, y, z), sp.Function('pi13')(t, x, y, z), \
    sp.Function('pi21')(t, x, y, z), sp.Function('pi22')(t, x, y, z), sp.Function('pi23')(t, x, y, z), \
    sp.Function('pi31')(t, x, y, z), sp.Function('pi32')(t, x, y, z), sp.Function('pi33')(t, x, y, z)
diss_vars = [q1, q2, q3, Pi, pi11, pi12, pi13, pi21, pi22, pi23, pi31, pi32, pi33]

# DissNSs - will need to swap Diss vars for these for Jacobian calcs
q1NS, q2NS, q3NS = sp.Function('q1NS')(t, x, y, z), sp.Function('q2NS')(t, x, y, z), sp.Function('q3NS')(t, x, y, z)
PiNS =  sp.Function('PiNS')(t, x, y, z)
pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS = \
    sp.Function('pi11NS')(t, x, y, z), sp.Function('pi12NS')(t, x, y, z), sp.Function('pi13NS')(t, x, y, z), \
    sp.Function('pi21NS')(t, x, y, z), sp.Function('pi22NS')(t, x, y, z), sp.Function('pi23NS')(t, x, y, z), \
    sp.Function('pi31NS')(t, x, y, z), sp.Function('pi32NS')(t, x, y, z), sp.Function('pi33NS')(t, x, y, z)
dissNSs = [q1NS, q2NS, q3NS, PiNS, pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS]

# Prims
v1, v2, v3 = sp.Function('v1')(t, x, y, z), sp.Function('v2')(t, x, y, z), sp.Function('v3')(t, x, y, z)
vs = [v1, v2, v3]
p, n, rho  = sp.Function('p')(t, x, y, z), sp.Function('n')(t, x, y, z),sp.Function('rho')(t, x, y, z)
prim_vars = [v1, v2, v3, p, n, rho]
                  
# Aux
T, W, Theta = sp.Function('T')(t, x, y, z), sp.Function('W')(t, x, y, z), sp.Function('Theta')(t, x, y, z)
qv, pi00, pi01, pi02, pi03 = sp.Function('qv')(t, x, y, z), sp.Function('pi00')(t, x, y, z), \
    sp.Function('pi01')(t, x, y, z), sp.Function('pi02')(t, x, y, z), sp.Function('pi03')(t, x, y, z)
aux_vars = [T, W, qv, pi00, pi01, pi02, pi03]

""" WRITE THE IDEAL PARTS OF THE STATE & FLUX VECTORS TO FILE """

# For writing these out to be read in for Jacobian calcs in SeSo.py
state_flux_outfile = open('state_flux_ideal.txt','w')
# state vector
D = n
S1 = (rho + p)*v1
S2 = (rho + p)*v2
S3 = (rho + p)*v3
E = rho
sv = [D, S1, S2, S3, E]
for i in range(len(sv)):
    state_flux_outfile.write(str(sv[i])+'\n')

# flux vector
fv = np.zeros((3,len(sv)),dtype=type(sv[0]))
for i in range(3):
    fv[i][0] = sv[0]*vs[i]
    fv[i][1] = sv[1]*vs[i]
    fv[i][2] = sv[2]*vs[i]
    fv[i][3] = sv[3]*vs[i]
    fv[i][i+1] += p
    fv[i][4] = (sv[4] + p)*vs[i]
    for j in range(len(sv)):
        fv[i][j] = sp.simplify(sp.expand(fv[i][j]))
        state_flux_outfile.write(str(fv[i][j])+'\n')

state_flux_outfile.close()

""" WRITE THE NS PARTS OF THE STATE & FLUX VECTORS TO FILE """

# May choose not to define explicitly for readability
#vsqrd = v1**2 + v2**2 + v3**2
#W = 1/sqrt(1-vsqrd)
qvNS = q1NS*v1 + q2NS*v2 + q3NS*v3
qsNS = [q1NS, q2NS, q3NS]
pi00NS = pi11NS + pi22NS + pi33NS
pi01NS = v1*pi11NS + v2*pi12NS + v3*pi13NS
pi02NS = v1*pi12NS + v2*pi22NS + v3*pi23NS
pi03NS = v1*pi13NS + v2*pi23NS + v3*pi33NS

svNS_out = open('state_flux_NS.txt','w')

# Define NS vector forms (H & F)
# NOTE ALL TERMS SHOULD BE NS
svNS = np.zeros_like(sv,dtype=type(sv[0]))
svNS[0] = 0
svNS[1] = PiNS*v1 + q1NS + pi01NS
svNS[2] = PiNS*v2 + q2NS + pi02NS
svNS[3] = PiNS*v3 + q3NS + pi03NS
svNS[4] = 2*qvNS + pi00NS
for i in range(len(svNS)):
    if svNS[i] == 0:
        svNS_out.write(str(0)+'\n')
        continue
    svNS_out.write(str(svNS[i].diff(t))+'\n')

fvNS = np.zeros_like(fv,dtype=type(sv[0]))
for i in range(3):
    fvNS[i][0] = svNS[0]*vs[i]
    fvNS[i][1] = svNS[1]*vs[i] + W*(qsNS[i]*v1 - qvNS*vs[i]*v1) 
    fvNS[i][2] = svNS[2]*vs[i] + W*(qsNS[i]*v2 - qvNS*vs[i]*v2) 
    fvNS[i][3] = svNS[3]*vs[i] + W*(qsNS[i]*v3 - qvNS*vs[i]*v3) 
    fvNS[i][i+1] += PiNS
    fvNS[i][4] = svNS[4]*vs[i] + W*(qsNS[i] - qvNS*vs[i]) 
    for j in range(len(svNS)):
        fvNS[i][j] = sp.simplify(sp.expand(fvNS[i][j]))
        svNS_out.write(str(fvNS[i][j])+'\n')

svNS_out.close()

C = sp.symbols('C')

# Define NS forms
q1NS = -kappa * T.diff(x) 
q2NS = -kappa * T.diff(y) 
q3NS = -kappa * T.diff(z)

Theta = v1.diff(x) + v2.diff(y) + v3.diff(z)
PiNS = -zeta * ( Theta )

pi11NS = -2*eta*( 2*(v1).diff(x) - C*Theta )
pi12NS = -2*eta*( (v2).diff(x) + (v1).diff(y) )
pi13NS = -2*eta*( (v3).diff(x) + (v1).diff(z) )
pi22NS = -2*eta*( 2*(v2).diff(y) - C*Theta )
pi23NS = -2*eta*( (v3).diff(y) + (v2).diff(z) )
pi33NS = -2*eta*( 2*(v3).diff(x) - C*Theta )
pi21NS = pi12NS
pi31NS = pi13NS
pi32NS = pi23NS

# Subsitute it all into the time-differentiated expressions
""" CALCULATE THE TIME DERIVATIVE OF THE NS STATE VECTOR """

dtsvNS = np.zeros_like(sv,dtype=type(sv[0]))
dtsvNS[0] = 0
dtsvNS[1] = (PiNS*v1 + q1NS + pi01NS).diff(t)
dtsvNS[2] = (PiNS*v2 + q2NS + pi02NS).diff(t)
dtsvNS[3] = (PiNS*v3 + q3NS + pi03NS).diff(t)
dtsvNS[4] = (2*qvNS + pi00NS).diff(t)

dtsvNS_out = open('state_NS_dt.txt','w')

for i in range(len(svNS)):
    #dtsvNS_out.write(str(sp.simplify(sp.expand(svNS[i].diff(t))))+'\n')
    dtsvNS_out.write(str(dtsvNS[i])+'\n')

dtsvNS_out.close()

dissNSs = [q1NS, q2NS, q3NS, PiNS, pi11NS, pi12NS, pi13NS, pi21NS, pi22NS, pi23NS, pi31NS, pi32NS, pi33NS]
dtsvNS_out = open('dt_NSs.txt','w')

for i in range(len(dissNSs)):
    #dtsvNS_out.write(str(sp.simplify(sp.expand(svNS[i].diff(t))))+'\n')
    dtsvNS_out.write(str(sp.simplify(sp.expand(dissNSs[i].diff(t))))+'\n')

dtsvNS_out.close()




# Declare the vars that we need to calculate the Jacobian of the state vector wrt
# (Essentially the variables that have time derivatives in the CE expansion)
jac_vars = [rho, n, v1, v2, v3]
#jac_vars = [W, v1, v2, v3]

# Convert the state and flux vectors into Matrices so that the Jacobian
# function can be used 
jac_vars_Mat = sp.Matrix(jac_vars)

# Convert state vector into sympy Matrix
sv_Mat = sp.Matrix(sv)
# Calculate Jacobian of state vector wrt jac_vars
sv_Jac = sv_Mat.jacobian(jac_vars_Mat)
sv_Jac_inv = sv_Jac.inv()
#sv_Jac_psrinv = sv_Jac.transpose()*(sv_Jac*sv_Jac.transpose()).inv()

# Form a list of the time derivatives of the required variables
dt_jac_vars = np.zeros_like(jac_vars,dtype=type(jac_vars[0]))
for i in range(len(jac_vars)):
    #dt_jac_vars[i] = jac_vars[i].diff(t)
    for j in range(len(cons)):
        # If the derivative of the conserved wrt the jac_var (n, W, etc.)
        # is zero then its reciprocal should be set to zero (not inf)
        if sv_Jac[i+j*len(jac_vars)] == 0:
            continue
        # Little numbering hack picks out all the required partial derivs
        # for each of the conserveds e.g. dn/dt = dD/dt*(dD/dn)^-1 + dS1/dt*dn/dS1 + ... + dTau/dt*dn/dTau
        #dt_jac_vars[i] += cons[j].diff(t)*(sv_Jac_inv[i*len(jac_vars)+j])
        dt_jac_vars[i] += cons[j].diff(t)*(1/sv_Jac[i+j*len(jac_vars)])

dt_jac_vars = sv_Jac_inv*sp.Matrix(dtcons)

# Write individual time differentials out so they can be pasted for use in
# the NS expression in the flux term in the model correction...
dts_out = open('dt_prims.txt','w')
for i in range(len(jac_vars)):
    dts_out.write(str(jac_vars[i].diff(t))+' = '+str(dt_jac_vars[i])+'\n')
dts_out.close()

dtsvNS_out = open('state_NS_dt.txt','a')
dtsvNS_out.write('\n'+'SUBSTITUTED EXPRESSION BEGINS'+'\n')
# Now make all the necessary substitutions into the time-diff'd state vector
for i in range(len(dtsvNS)):
    for j in range(len(jac_vars)):
        if dtsvNS[i] == 0:
            continue
        dtsvNS[i] = dtsvNS[i].subs(jac_vars[j].diff(t),dt_jac_vars[j])
    dtsvNS_out.write(str(dtsvNS[i])+'\n')
        
dtsvNS_out.close()


