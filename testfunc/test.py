#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 22:56:11 2024

@author: mr
"""
import sympy as sp
sp.init_printing(use_unicode=True)

# Define the symbols
x, y, z, r, theta, phi = sp.symbols('x y z r theta phi')

# Define the Cartesian to spherical coordinate transformation
r_expr = sp.sqrt(x**2 + y**2 + z**2)
theta_expr = sp.atan(sp.sqrt(y**2 + z**2) / x)
phi_expr = sp.atan(z/ y)

# Define the function U in terms of x, y, and z
U = sp.Function('U')(r, theta, phi)

# Calculate the derivatives of r, theta, and phi with respect to x, y, and z
dr_dx = sp.diff(r_expr, x)
dtheta_dx = sp.diff(theta_expr, x)
dphi_dx = sp.diff(phi_expr, x)

dr_dy = sp.diff(r_expr, y)
dtheta_dy = sp.diff(theta_expr, y)
dphi_dy = sp.diff(phi_expr, y)

dr_dz = sp.diff(r_expr, z)
dtheta_dz = sp.diff(theta_expr, z)
dphi_dz = sp.diff(phi_expr, z)

# Calculate the derivatives of U with respect to x, y, and z
dU_dr = sp.diff(U, r)
dU_dtheta = sp.diff(U, theta)
dU_dphi = sp.diff(U, phi)


# Apply the chain rule to calculate dU/dx
dU_dx = dU_dr * dr_dx + dU_dtheta * dtheta_dx + dU_dphi * dphi_dx
dU_dy = dU_dr * dr_dy + dU_dtheta * dtheta_dy + dU_dphi * dphi_dy
dU_dz = dU_dr * dr_dz + dU_dtheta * dtheta_dz + dU_dphi * dphi_dz

#chain rule for second derivative
dU_dx = sp.simplify(dU_dx.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dU_dy = sp.simplify(dU_dy.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dU_dz = sp.simplify(dU_dz.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))

dr_dx = sp.simplify(dr_dx.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dtheta_dx = sp.simplify(dtheta_dx.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dphi_dx = sp.simplify(dphi_dx.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))

dr_dy = sp.simplify(dr_dy.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dtheta_dy = sp.simplify(dtheta_dy.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dphi_dy = sp.simplify(dphi_dy.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))

dr_dz = sp.simplify(dr_dz.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dtheta_dz = sp.simplify(dtheta_dz.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))
dphi_dz = sp.simplify(dphi_dz.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)}))


#dU_dxx = sp.diff(dU_dx,r) *dr_dx + sp.diff(dU_dx,theta)  * dtheta_dx + sp.diff(dU_dx, phi) * dphi_dx
#dU_dyy = sp.diff(dU_dy,r) *dr_dy + sp.diff(dU_dy,theta)  * dtheta_dy + sp.diff(dU_dy, phi) * dphi_dy
#dU_dzz = sp.diff(dU_dz,r) *dr_dz + sp.diff(dU_dz,theta)  * dtheta_dz + sp.diff(dU_dz, phi) * dphi_dz

#dU_dxy = sp.diff(dU_dx,r) *dr_dy + sp.diff(dU_dx,theta)  * dtheta_dy + sp.diff(dU_dx, phi) * dphi_dy
#dU_dxz = sp.diff(dU_dx,r) *dr_dz + sp.diff(dU_dx,theta)  * dtheta_dz + sp.diff(dU_dx, phi) * dphi_dz
dU_dyz = sp.diff(dU_dy,r) *dr_dz + sp.diff(dU_dy,theta)  * dtheta_dz + sp.diff(dU_dy, phi) * dphi_dz


AAA = sp.simplify(dU_dyz)
#BBB = sp.trigsimp(dU_dyy)
#CCC = sp.trigsimp(dU_dzz)

print(sp.latex(sp.expand(sp.simplify(sp.expand_trig(AAA)))))

#print(sp.latex(sp.simplify(AAA+BBB+CCC)))

#print(sp.latex(sp.simplify(dU_dxx)))

#print("The derivative of U with respect to x is:", sp.latex(dU_dx))
#print("The derivative of U with respect to y is:", sp.latex(dU_dy))
#print("The derivative of U with respect to z is:", dU_dz)

#AAA = dU_dz.subs({x:r * sp.cos(theta),y:r* sp.sin(theta)* sp.cos(phi),z:r* sp.sin(theta)* sp.sin(phi)})
#aaa = sp.trigsimp(AAA)
#print(sp.latex((aaa)))