# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 21:10:53 2019

tensileStrength.py

@author: Myself
"""

import  csv
from pint import UnitRegistry

ureg = UnitRegistry()

tensileStrengths = dict()

with open('Wire Tensile Strengths.txt', 'r') as f:
    next(f)
    fReader = csv.reader(f)
    headers = next(fReader)
    headers[0] = 'WireDia'
    columns = [[] for _ in headers]
    for line in fReader:
        for header, l, c in zip(headers, line, columns):
            if header == headers[0]:
                l = float(l)*ureg.inch
                l.ito(ureg.mm)
            else:
                l = float(l) * 1e3 * ureg.psi
                l.ito(ureg.Pa)
            c.append(l)
    for header, column in zip(headers, columns):
        tensileStrengths[header] = column

"""
#Below is code used to create the Wire Tensile Strengths.txt file

import csv

fname = 'Min Tensile Strength MW HD OT.txt'

plainSteelTens = []

with open(fname, 'r') as f:
    prev = ['0']*4
    for line in reversed(list(csv.reader(f))):
        line.extend([None]*(4-len(line)))
        curr = []
        for p, l in zip(prev, line):
            if l:
                curr.append(l)
            else:
                curr.append(p)
        plainSteelTens.append(curr)
        prev = curr
        
ssfname = 'Min Tensile Strength 302 SST 17-7 ST.txt'

SSteelTens = []

with open(ssfname, 'r') as f:
    prev = ['0']*3
    for line in reversed(list(csv.reader(f))):
        line.extend([None]*(3-len(line)))
        curr = []
        for p, l in zip(prev, line):
            if l:
                curr.append(l)
            else:
                curr.append(p)
        SSteelTens.append(curr)
        prev = curr
        
plainSteelTens.reverse()
SSteelTens.reverse()

with open('Wire Tensile Strengths.txt', 'w', newline='') as f:
    fWriter = csv.writer(f)
    plainIter = iter(plainSteelTens)
    ssIter = iter(SSteelTens)
    fWriter.writerow(next(plainIter)[:1])
    next(ssIter)
    fWriter.writerow(next(plainIter) + next(ssIter)[1:])
    plain = next(plainIter)
    ss = next(ssIter)
    while(True):
        try:
            if plain[0] == ss[0]:
                fWriter.writerow(plain + ss[1:])
                plain = next(plainIter)
                ss = next(ssIter)
            elif plain[0] <= ss[0]:
                fWriter.writerow(plain + ss[1:])
                plain = next(plainIter)
            else:
                fWriter.writerow(ss[:1] + plain[1:] + ss[1:])            
                ss = next(ssIter)
        except:
            break
"""

























