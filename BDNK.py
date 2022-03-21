# -*- coding: utf-8 -*-
"""
Created on Wed March 16 2022

@author: mjh1n20
"""

import sympy as sp
import numpy as np

t, x, y, z = sp.symbols('t x y z')

D, S1, S2, S3, E = sp.Function('D')(t, x, y, z), sp.Function('S1')(t, x, y, z), \
    sp.Function('S2')(t, x, y, z), sp.Function('S3')(t, x, y, z), sp.Function('E')(t, x, y, z)
cons = [D, S1, S2, S3, E]
f_D, f_S1, f_S2, f_S3, f_E = sp.Function('f_D')(t, x, y, z), sp.Function('f_S1')(t, x, y, z), \
    sp.Function('f_S2')(t, x, y, z), sp.Function('f_S3')(t, x, y, z), sp.Function('f_E')(t, x, y, z)
fluxes = [3*[f_D, f_S1, f_S2, f_S3, f_E]]

# Prims
v1, v2, v3 = sp.Function('v1')(t, x, y, z), sp.Function('v2')(t, x, y, z), sp.Function('v3')(t, x, y, z)
# v2, v3 = 0, 0 # big big simplification
vs = [v1, v2, v3]
W = sp.Function('W')(t, x, y, z)
# W = 1/(1-v1**2 + v2**2 + v3**2)**(1/2) # optional setting

n = sp.Function('n')(t, x, y, z)
rho = sp.Function('rho')(t, x, y, z)
p = sp.Function('p')(t, x, y, z)
#p = sp.Function('p')(rho, n)
Gamma = sp.symbols('Gamma')
prim_vars = [v1, v2, v3, p, n, rho]

# state vector simplest
# D = (rho - p/(1-Gamma))
# S1 = (rho + p)*v1
# S2 = (rho + p)*v2
# S3 = (rho + p)*v3
# E = rho

# better
D = (rho - p/(1-Gamma))*W # n is substituted here with eos
S1 = (rho + p)*v1*W**2
S2 = (rho + p)*v2*W**2
S3 = (rho + p)*v3*W**2
E = (rho + p)*W**2 - p

sv = [D, S1, S2, S3, E]
#sv = [D, S1, E] # 1D version

# make flux vector
fv = np.zeros((3,5),dtype=type(sv[0]))
for i in range(3):
    fv[i,0] = D*vs[i]
    for j in range(3):
        fv[i,j+1] = sv[j+1]*vs[i]
        if (i == j):
            fv[i,j+1] += p
    fv[i,4] = (sv[4] + p)*vs[i]

# split into components
fv_x, fv_y, fv_z = fv[0], fv[1], fv[2]

# pure 1D version
# fv_x = [D*v1, S1*v1 + p, (E + p)*v1]
# fv_y, fv_z = fv_x, fv_x # silly hack, unused vectors

# Declare the vars that we need to calculate the Jacobian of the state vector wrt
# (Essentially the variables that have time derivatives in the CE expansion)
jac_vars = [p, rho, v1, v2, v3]
jac_vars_diff_x = [p.diff(x), rho.diff(x), v1.diff(x), v2.diff(x), v3.diff(x)]
jac_vars_diff_y = [p.diff(y), rho.diff(y), v1.diff(y), v2.diff(y), v3.diff(y)]
jac_vars_diff_z = [p.diff(z), rho.diff(z), v1.diff(z), v2.diff(z), v3.diff(z)]
#jac_vars = [W, v1, v2, v3]

# Convert the state and flux vectors into Matrices so that the Jacobian
# function can be used 
jac_vars_Mat = sp.Matrix(jac_vars)

# Convert state vector into sympy Matrix
sv_Mat = sp.Matrix(sv)
# Calculate Jacobian of state vector wrt jac_vars
sv_Jac = sv_Mat.jacobian(jac_vars_Mat)
sv_Jac_inv = sv_Jac.inv()

# do the same for flux vector without inversion
fvx_Mat, fvy_Mat, fvz_Mat = sp.Matrix(fv_x), sp.Matrix(fv_y), sp.Matrix(fv_z)
fvx_Jac, fvy_Jac, fvz_Jac = fvx_Mat.jacobian(jac_vars_Mat), fvy_Mat.jacobian(jac_vars_Mat), fvz_Mat.jacobian(jac_vars_Mat)

#print(sv_Jac)

outfile = open('BDNK_3D_full.txt','w')

result = []
for i in range(len(jac_vars)):
    result.append((-sv_Jac_inv*( fvx_Jac*sp.Matrix(jac_vars_diff_x) + fvy_Jac*sp.Matrix(jac_vars_diff_y) + \
                                fvz_Jac*sp.Matrix(jac_vars_diff_z) ))[i])

for i in range(len(jac_vars)):
    print(str(jac_vars[i].diff(t))+'  =  ')
    print(result[i])
    print('\n')
    outfile.write(str(jac_vars[i].diff(t))+'  =  ')
    outfile.write(str(sp.simplify(sp.expand(result[i])))+'\n')


# for i in range(len(jac_vars)):
#     print(str(jac_vars[i].diff(t))+'  =  ')
#     print((-sv_Jac_inv*(fvx_Jac*sp.Matrix(jac_vars_diff_x)))[i])
#     print('\n')
#     outfile.write(str(jac_vars[i].diff(t))+'  =  ')
#     outfile.write(str(sp.simplify(sp.expand((-sv_Jac_inv*fvx_Jac*sp.Matrix(jac_vars_diff_x))[i])))+'\n')

outfile.close()

