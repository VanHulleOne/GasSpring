# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 20:29:50 2019

springs.py

@author: Myself
"""

from collections import namedtuple as nt
from math import pi

COMP_FILE = 'CS Comp cat.txt'
SPRING_FILE = 'CSCompSprings.txt'
IMP = 0
METRIC = 1
MIN = 0
MAX = 1
INF = float('inf')
#                 Min Tens Mod of Elasticity
#                   (psi)
matProperties = {'MW': (280000, 11500000), # MW
                 'SST':(210000, 10000000)} # SST

offsets = {'OD':0, 'length':3, 'ID':5, 'rate':7, 'deflection':9, 'load':11,
           'solidLength':13, 'wireDia':15} # offsets to switch between units
options = ['material', 'ends', 'finish']

Spring = nt('Spring',
            'OD_Imp \
            OD_Metric \
            num \
            length_Imp \
            length_Metric \
            ID_Imp \
            ID_Metric \
            rate_Imp \
            rate_Metric \
            deflection_Imp \
            deflection_Metric \
            load_Imp \
            load_Metric \
            solidLength_Imp \
            solidLength_Metric \
            wireDia_Imp \
            wireDia_Metric \
            totalCoils \
            material \
            ends \
            finish')

springs = []

with open(SPRING_FILE, 'r') as f:
    for line in f:
        temp = line.split(sep=' ')
        if len(temp) == 21:            
            springs.append(Spring(*([float(x) for x in temp[:2]] + 
                                     [temp[2]] +
                                     [float(x) for x in temp[3:18]] + 
                                     temp[18:])))
        else:
            print(temp)

def filterer(spring, *, units=IMP,
             OD=None,
             length=None,
             ID=None,
             rate=None,
             deflection=None,
             load=None,
             solidLength=None,
             wireDia=None,
             totalCoils=None,
             material=None,
             ends=None,
             finish=None):
    
    for key, value in locals().items():
        if value is not None:
            if key in options:
                if getattr(spring, key) not in value:
                    return False
            elif key in offsets:
                springValue = spring[offsets[key]+units]
                if value[MIN] > springValue < value[MAX]:
                    return False
            elif key == 'totalCoils':
                if value[MIN] > spring.totalCoils < value[MAX]:
                    return False
    return True

#print(filterer(springs[0], OD=(0, 0.04)))

def getSprings(*,
               units=IMP,
               OD=None,
               ID=None,
               length1=None,
               length2=None,
               force1=None,
               force2=None,
               tol=0.2,
               material=None,
               ends=None):
    rate = (force2-force1)/(length1-length2)
    filteredSprings = [spring for spring in springs if filterer(spring,
                                                                units=units,
                                                                OD=OD,
                                                                ID=ID,
                                                                solidLength=[0, length2],
                                                                length=[length1, INF],
                                                                rate=[(1-tol)*rate, (1+tol)*rate],
                                                                material=material
                                                                )]
    doubleFilter = []
    for spring in filteredSprings:
        rateAct = spring[offsets['rate']+units]
        freeLength = spring[offsets['length']+units]
        f1Act = (freeLength-length1)*rateAct
        if (1-tol)*force1 < f1Act < (1+tol)*force1:
            doubleFilter.append(spring)
            
    fatigueFiltered = [s for s in doubleFilter if fatigueLimit(s, length1, length2)[0]]
            
    return fatigueFiltered

def fatigueLimit(spring, l1, l2, units=IMP):
    _id = spring[offsets['ID']+units]
    wireDia = spring[offsets['wireDia']+units]
    tensileStrength = matProperties[spring.material][0]
    freeLength = spring[offsets['length']+units]
    rate = spring[offsets['rate']+units]

    K = (4.0*_id/wireDia - 1)/(4.0*_id/wireDia - 4) + 0.615*wireDia/_id
    
    f1 = (freeLength - l1)*rate
    f2 = (freeLength - l2)*rate
    
    sig1 = 8*f1*_id*K/(pi*wireDia**3)
    sig2 = 8*f2*_id*K/(pi*wireDia**3)
    
    R1 = sig1/tensileStrength
    R2 = sig2/tensileStrength
    
    if R1 > 0.34 and R2 < 0.48:
        return True, R1, R2
    if R2 < (0.4*R1 + 0.34):
        return True, R1, R2
    return False, R1, R2

fs = getSprings(ID=[0.4,0.6], length1=2, length2=1, force1=2, force2=10, material='SST', ends='CG')
    
    










