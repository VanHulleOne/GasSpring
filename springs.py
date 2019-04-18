# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 20:29:50 2019

springs.py

@author: Myself
"""

from collections import namedtuple as nt
import bisect
from math import pi

#from pint import UnitRegistry

from tensileStrength import tensileStrengths, ureg

#ureg = UnitRegistry()

spring_rate = ureg('N/mm')
imp_spring_rate = ureg('lbf/in')

COMP_FILE = 'CS Comp cat.txt'
#SPRING_FILE = 'CSCompSprings.txt'
SPRING_FILE = 'correctedSpring.txt'

IMP = 0
METRIC = 1
MIN = 0
MAX = 1
INF = float('inf')

NOT_POSSIBLE = -1
NOT_RECOMMENDED = 0
L_10e5 = 1
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

strenghtReductionFactors = {'MW':0.45, 'SPR': 0.45, 'HD':0.40, 'OT':0.45,
                            'SST':0.30, '17-7':0.45, 'BC':0.45, 'PB':0.40}
elasticModuli = {'MW':11.5e6, 'SPR': 11.5e6, 'HD':11.5e6,
                 'OT':11.5e6, 'SST':10e6, '17-7':11e6,
                 'BC':18.5e6, 'PB':15e6,
                 }
activeCoilsDiff = {'O':0, 'OG':1, 'C':2, 'CG':2}

class Spring():
    def __init__(self, catRow):
        self.name = catRow[2]
        self.OD = float(catRow[0]) 
        self.freeLength = float(catRow[3]) 
        self.ID = float(catRow[5]) 
        self.wireDia = float(catRow[15]) 
        self.rate = float(catRow[7])
        self.maxDeflection = float(catRow[9]) 
        self.maxLoad = float(catRow[11]) 
        self.solidLength = float(catRow[13]) 
        self.numCoils = float(catRow[17])
        self.material = catRow[18]
        self.ends = catRow[19]
        self.safeSolidLength = bool(self.freeLength - self.maxDeflection < 1.01*self.solidLength)
        
        if self.material == 'SPR':
            self.material = 'HD'

        self.activeCoils = self.numCoils - activeCoilsDiff[self.ends]
        
        self.possibleDataErrors = False        
        self.checks()
        
    def checks(self):
        if abs(self.OD - self.ID - 2*self.wireDia) > .1 * self.wireDia:
#            print(self.name, self.OD, self.ID, self.wireDia)
            self.possibleDataErrors = True
        if abs(self.maxDeflection*self.rate - self.maxLoad) > .1 * self.maxLoad:
            self.possibleDataErrors = True
#            error = (self.maxDeflection*self.rate-self.maxLoad)/self.maxLoad
#            print(self.name, self.maxDeflection, self.rate, self.maxLoad, error, sep=': ')
            
        calcRate = elasticModuli[self.material]*self.wireDia**4/(8*self.activeCoils*(self.OD-self.wireDia)**3)

        if abs(calcRate-self.rate)/self.rate > 0.25:
#            print(self.name, self.rate, calcRate, sep='\t')
            self.possibleDataErrors = True
        
    def getForce(self, length):
        return (self.freeLength - length) * self.rate
    
    def getStress(self, length):
        deflection = self.freeLength - length
        D = self.OD - self.wireDia
        C = D/self.wireDia
        K = (4*C-1)/(4*C-4) + 0.615/C # Wahl stress correction factor: Table 5.2 Pg 190
        stress = 8*self.rate*D*K*deflection/(pi*self.wireDia**3)
        return stress
    
    def getMinTensileStrength(self):
        offset = bisect.bisect_right(tensileStrengths['WireDia'], self.wireDia-.0001)
        minTens = tensileStrengths[self.material][offset]
        return minTens
    
    def _getMillionDefl(self):
        minTens = self.getMinTensileStrength()
        reductionFactor = strenghtReductionFactors[self.material]
        
        stress = minTens * reductionFactor
        
        D = self.OD - self.wireDia
        C = D/self.wireDia
        K = (4*C-1)/(4*C-4) + 0.615/C # Wahl stress correction factor: Table 5.2 Pg 190
        
        return pi*self.wireDia**3*stress/(8*self.rate*D*K)
        
    def _getnCycles(self):
        deflection = self._getMillionDefl()
        minTens = self.getMinTensileStrength()
        reductionFactor = strenghtReductionFactors[self.material]
        
        K_s1 = 0
        K_s2 = self.getStress(self.freeLength-deflection)/minTens
        
        K_U = 0.56
        C_S = 0.5546
        M = -0.009
        C_E = 0.662
        Y = -0.0622
        
        K_E = K_U*(K_s2 - K_s1)/(2*K_U - (K_s2+K_s1))
                
        n = (C_E/K_E)**(1/-Y)
        return n
    
    def expectedLife(self, length):
        deflection = self.freeLength - length
        if deflection > self.freeLength - self.solidLength:
            return NOT_POSSIBLE
        if deflection > self.maxDeflection:
            return NOT_RECOMMENDED
        
        offset = bisect.bisect_right(tensileStrengths['WireDia'], self.wireDia)
        minTens = tensileStrengths[self.material][offset]
        reductionFactor = strenghtReductionFactors[self.material]
        
        stress = self.getStress(length)
        
        if stress > minTens * reductionFactor:
            return L_10e5
        if stress > minTens * (reductionFactor - 0.1):
            return L_MILLION
        return L_INF
    
    def __repr__(self):
        return self.name
        

with open(SPRING_FILE, 'r') as f:
    for line in f:
        temp = line.split(sep=' ')
        if len(temp) == 21:            
            springs.append(Spring(temp))
        else:
            print(temp)

def minMaxFilter(inSprings, minMaxDict):
    for key, values in minMaxDict.items():
        if values[0] == 0 and values[1] == INF:
            continue
        inSprings = [spring for spring in inSprings if values[0] <= getattr(spring, key) <= values[1]]
    print('After minmaxFilter:', len(inSprings))
    return inSprings

def optionsFilter(inSprings, options):
    for key, value in options.items():
        if value is None:
            continue
        inSprings = [spring for spring in inSprings if getattr(spring, key) == value]
    print('After optionsFilter:', len(inSprings))
    return inSprings

def lifeFilter(inSprings, length1, length2, minLife):
    values = sorted(l for l in [length1, length2] if l is not None)
    if values:
        inSprings = [s for s in inSprings if s.expectedLife(values[0]) >= minLife]
    print('After lifeFilter:', len(inSprings))
    return inSprings
        
def forceFilter(inSprings, lengths, forces, tolerances):
    for length, force, tol in zip(lengths, forces, tolerances):
        if length is None:
            continue
        length=length
        force = force
        inSprings = [s for s in inSprings if (1-tol)*force <= s.getForce(length) <= (1+tol)*force]        
    print('After forceFilter:', len(inSprings))
    return inSprings
    

        
            
def getSprings2(thisSprings,*,
                minOD=0,
                maxOD=INF,
                minID=0,
                maxID=INF,
                minFreeLength=0,
                maxFreeLength=INF,
                length1=None,
                length2=None,
                force1=None,
                force2=None,
                safeSolidLength=None, # If the max deflection is at the solid length the spring can't be wrecked by over compression
                material=None,
                ends='CG',
                finish=None,
                tolerance1=0.2,
                tolerance2=0.2,
                numResults=INF,
                minLife=L_INF,
                forceFactor=2,
                ):
    minMaxDict = {'OD':[minOD, maxOD],
                  'ID':[minID, maxID],
                  'freeLength':[minFreeLength, maxFreeLength],
                  }
    optionsDict = {'material':material,
                   'ends':ends,
                   'finish': finish,
                   'safeSolidLength':safeSolidLength,
                   }
    
    thisSprings = optionsFilter(thisSprings, optionsDict)
    thisSprings = minMaxFilter(thisSprings, minMaxDict)
    thisSprings = lifeFilter(thisSprings, length1, length2, minLife)
    
    lessForceSprings = forceFilter(thisSprings, [length1, length2], [force1/forceFactor, force2/forceFactor], [tolerance1, tolerance2])
    targetForceSprings = forceFilter(thisSprings, [length1, length2], [force1, force2], [tolerance1, tolerance2])
    moreForceSprings = forceFilter(thisSprings, [length1, length2], [force1*forceFactor, force2*forceFactor], [tolerance1, tolerance2])

#    results = getInformativeResults(lessForceSprings, targetForceSprings, moreForceSprings)
    return lessForceSprings, targetForceSprings, moreForceSprings

s = getSprings2(springs, maxOD=0.55, minID=0.4, length1=0.7, length2=0.38,
                force1=0.49, force2=1.7, material='SST', tolerance1=0.8, safeSolidLength=None)

lengths = [len(x) for x in s]
print('Lengths:', lengths)





#def fatigueLimit(spring, l1, l2, units=IMP):
#    _id = spring[offsets['ID']+units]
#    wireDia = spring[offsets['wireDia']+units]
#    tensileStrength = matProperties[spring.material][0]
#    freeLength = spring[offsets['length']+units]
#    rate = spring[offsets['rate']+units]
#
#    K = (4.0*_id/wireDia - 1)/(4.0*_id/wireDia - 4) + 0.615*wireDia/_id
#    
#    f1 = (freeLength - l1)*rate
#    f2 = (freeLength - l2)*rate
#    
#    sig1 = 8*f1*_id*K/(pi*wireDia**3)
#    sig2 = 8*f2*_id*K/(pi*wireDia**3)
#    
#    R1 = sig1/tensileStrength
#    R2 = sig2/tensileStrength
#    
#    if R1 > 0.34 and R2 < 0.48:
#        return True, R1, R2
#    if R2 < (0.4*R1 + 0.34):
#        return True, R1, R2
#    return False, R1, R2

