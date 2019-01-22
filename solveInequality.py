# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 17:03:15 2019

solveInequality.py

@author: lvanhulle
"""

from sympy.abc import x, y
from sympy import Poly
from sympy.solvers.inequalities import solve_rational_inequalities
from sympy import solveset, S

#res = solve_rational_inequalities([[
#        ((Poly(x+y-2), Poly(1,x)), '>'),
##        ((Poly(x-10), Poly(1,x)), '<')
#        ]])

res = solveset([x + 7,  y - 10], dict=True)

print(res)