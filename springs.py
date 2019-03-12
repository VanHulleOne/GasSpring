# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 20:29:50 2019

springs.py

@author: Myself
"""

from collections import namedtuple as nt
import bisect
from math import pi

from pint import UnitRegistry

from tensileStrength import tensileStrengths

ureg = UnitRegistry()

spring_rate = ureg('N/mm')
imp_spring_rate = ureg('lbf/in')

COMP_FILE = 'CS Comp cat.txt'
SPRING_FILE = 'CSCompSprings.txt'
IMP = 0
METRIC = 1
MIN = 0
MAX = 1
INF = float('inf')

NOT_POSSIBLE = -1
NOT_RECOMMENDED = 0
L_10e6 = 1
L_MILLION = 2
L_INF = 3


#                 Min Tens Mod of Elasticity
#                   (psi)
#matProperties = {'MW': (280_000, 11_500_000), # MW
#                 'SST':(210_000, 10_000_000)} # SST
#
#offsets = {'OD':0, 'length':3, 'ID':5, 'rate':7, 'deflection':9, 'load':11,
#           'solidLength':13, 'wireDia':15} # offsets to switch between units
#options = ['material', 'ends', 'finish']

#Spring = nt('Spring',
#        0    'OD_Imp \
#        1    OD_Metric \
#        2    num \
#        3    length_Imp \
#        4    length_Metric \
#        5    ID_Imp \
#        6    ID_Metric \
#        7    rate_Imp \
#        8    rate_Metric \
#        9    deflection_Imp \
#       10    deflection_Metric \
#       11    load_Imp \
#       12    load_Metric \
#       13    solidLength_Imp \
#       14    solidLength_Metric \
#       15    wireDia_Imp \
#       16    wireDia_Metric \
#       17    totalCoils \
#       18    material \
#       19    ends \
#       20    finish')

springs = []

strenghtReductionFactors = {'MW':0.45, 'HD':0.40, 'OT':0.45, 'SST':0.30, '17-7':0.45}

class Spring():
    def __init__(self, catRow):
        self.name = catRow[2]
        self.OD = float(catRow[1]) * ureg.mm
        self.freeLength = float(catRow[4]) * ureg.mm
        self.ID = float(catRow[6]) * ureg.mm
        self.wireDia = float(catRow[16]) * ureg.mm
        self.rate = float(catRow[8]) * spring_rate
        self.maxDeflection = float(catRow[10]) * ureg.mm
        self.maxLoad = float(catRow[12]) * ureg.N
        self.solidLength = float(catRow[14]) * ureg.mm
        self.numCoils = float(catRow[17])
        self.material = catRow[18]
        self.ends = catRow[19]
        
        self.checks()
        
    def checks(self):
        if abs(self.OD - self.ID - 2*self.wireDia) > .1 * self.wireDia:
            print(self.name, self.OD, self.ID, self.wireDia)
        
    def getForce(self, length):
        return (self.freeLength - length) * self.rate
    
    def getStress(self, deflection):
        D = self.OD - self.wireDia
        C = D/self.wireDia
        K = (4*C-1)/(4*C-4) + 0.615/C
        stress = 8*self.rate*D*K*deflection/(pi*self.wireDia**3)
        return stress
    
    def expectedLife(self, deflection):
        if self.freeLength - self.solidLength > deflection:
            return NOT_POSSIBLE
        if deflection > self.maxDeflection:
            return NOT_RECOMMENDED
        
        offset = bisect.bisect_right(tensileStrengths['WireDia'], self.wireDia)
        minTens = tensileStrengths[self.material][offset]
        reductionFactor = strenghtReductionFactors[self.material]
        
        stress = self.getStress(deflection)
        
        if stress > minTens * reductionFactor:
            return L_10e6
        if stress > minTens * (reductionFactor + 0.1):
            return L_MILLION
        return L_INF
        
        

with open(SPRING_FILE, 'r') as f:
    for line in f:
        temp = line.split(sep=' ')
        if len(temp) == 21:            
            springs.append(Spring(temp))
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
                if getattr(spring, key) != value:
                    return False
            elif key in offsets:
                springValue = spring[offsets[key]+units]
                if springValue < value[MIN] or springValue > value[MAX]:#value[MIN] > springValue > value[MAX]:
                    return False
            elif key == 'totalCoils':
                raise Exception('Total Coils not implimented')
#                if value[MIN] > spring.totalCoils > value[MAX]:
#                    return False
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
                                                                material=material,
                                                                ends=ends
                                                                )]
    print('Num filtered:', len(filteredSprings))
    doubleFilter = []
    for spring in filteredSprings:
        rateAct = spring[offsets['rate']+units]
        freeLength = spring[offsets['length']+units]
        f1Act = (freeLength-length1)*rateAct
        if (1-tol)*force1 < f1Act < (1+tol)*force1:
            doubleFilter.append(spring)
            
    print('Num double:', len(doubleFilter))
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

#fs = getSprings(ID=[0.4,0.6], length1=2, length2=1, force1=2, force2=10, material='SST', ends='CG')
    
#stripSprings = getSprings(ID=[0.323, INF], OD=[0, 0.49], length1=2.046, length2=1.862,
#                          force1=4, force2=8, material='SST', ends='CG', tol=0.3)
#print('NumSprings', len(stripSprings))

#retSprings = getSprings(ID=[0.403, INF], OD=[0, 0.58], length1=0.531, length2=0.374,
#                       force1=0.7, force2=1.0, material='SST', ends='CG', tol=0.25)
#print('NumSprings', len(retSprings))









