#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 16:08:01 2024

@author: mr
"""
import sympy as sp
sp.init_printing(use_unicode=True)

A, theta, phi, omega = sp.symbols('A theta phi omega')

u10 = 2/5 * (-2*A**5 + 6*A**4-5*A**3 +1)
u12 = 4/5 * (1-A)**3 *A**2

u1 = u10 + u12 * 1/2 * (sp.cos(theta)**2 -1)
u = omega**2 * u1

Ua = sp.diff(u, A)
Uaa = sp.diff(Ua, A)
Ut = sp.diff(u, theta)
Utt = sp.diff(Ut, theta) 


lhs = Uaa + 2/A * Ua + 1/(A**2 * (1-A)**2) * (Utt + Ut * sp.cot(theta)) 
rhs = -36 * omega**2 * A*(1-A)**2 *sp.sin(theta) * 1/(1+A*u)**7
print(sp.diff(u , A,A))