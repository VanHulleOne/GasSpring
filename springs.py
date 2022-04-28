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

from tensileStrength import tensileStrengths#, #ureg

#ureg = UnitRegistry()

# spring_rate = ureg('N/mm')
# imp_spring_rate = ureg('lbf/in')

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

lifeDict = {'10^5': L_10e5,
            '10^6': L_MILLION,
            'Inf': L_INF,
            }

GoodmanLine = nt('GoodmanLine', 'life y_intercept max_y x_at_max_y')


'''
Goodman Graph data is estimated from SAE Spring Design Manual AE-21 1996, Fig. 5.1 (page 192)

"Fig 5.1-Fatigue strength diagram for round wire helical compression springs which are not preset.
All stresses are Wahl corrected with diagram representing a B-10 fatigue life.
Diagram applicable to springs which are not preset and of the following materials:
    Music Steel Spring Wire and Springs-SAE J178
    Hard Drawn Carbon Vale Spring Quality Wiare and Springs-SAE J172
    Oil Tempered Cabon Walve Spring Quality Wire and Springs-SAE J351
    Oil Tempered Chromium-Vanadium Vale Spring Quality Wire and Springs-SAE J132"
    
'''
goodmanLines = [GoodmanLine('10^5', .414, .500, .318),
                GoodmanLine('10^6', .378, .491, .342),
                GoodmanLine('10^7', .345, .480, .351),
                ]

springs = []

# From Century Spring catalogue Appendix A, Material Properties of Common Spring Materials
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
        self.minLength = self.freeLength - self.maxDeflection
        
        self.lastLength = self.freeLength
        
        if self.material == 'SPR':
            self.material = 'HD'

        self.activeCoils = self.numCoils - activeCoilsDiff[self.ends]
        
        self.possibleDataErrors = False        
        self.checks()
        
        self.minTens = self.getMinTensileStrength()
        self.reductionFactor = strenghtReductionFactors[self.material]
        
    def summary(self):
        s = [self.name]
        s.append('OD: ' + str(self.OD))
        s.append('ID: ' + str(self.ID))
        s.append('Free Length: ' + str(self.freeLength))
        s.append('Rate: ' + str(self.rate))
        
        return '\t'.join(s)
        
    def checks(self):
        if abs(self.OD - self.ID - 2*self.wireDia) > .1 * self.wireDia:
            self.possibleDataErrors = True
        if abs(self.maxDeflection*self.rate - self.maxLoad) > .1 * self.maxLoad:
            self.possibleDataErrors = True
            
        calcRate = elasticModuli[self.material]*self.wireDia**4/(8*self.activeCoils*(self.OD-self.wireDia)**3)

        if abs(calcRate-self.rate)/self.rate > 0.25:
            self.possibleDataErrors = True
        
    def getForce(self, length):
        return (self.freeLength - length) * self.rate
    
    def getDeflection(self, force):
        deflection = force/self.rate
        if deflection > self.maxDeflection:
            return None 
        return deflection
    
    def getMaxGoodman_Ks2(self, Ks1, goodmanLine):
        x1 = 0
        x2 = goodmanLine.x_at_max_y    
        y1 = goodmanLine.y_intercept
        y2 = goodmanLine.max_y
        
        y = (y2-y1)/(x2-x1)*(Ks1-x1)+y1
        
        return min(y, goodmanLine.max_y)
        
    
    def getGoodmanLife(self, length1, length2=None):
        if length2 is None:
            length1, length2 = self.freeLength, length1
        
        minTens = self.getMinTensileStrength()
        initStress = self.getStress(length1)
        maxStress = self.getStress(length2)
        
        Ks1 = initStress/minTens
        Ks2 = maxStress/minTens
        
        for goodmanLine in goodmanLines[::-1]:
            maxKs2 = self.getMaxGoodman_Ks2(Ks1, goodmanLine)
            if Ks2 < maxKs2:
                print(self, goodmanLine.life)
                return
        print('Not much')
        
    
    def getStress(self, length=None, deflection = None):
        if deflection is None:
            deflection = self.freeLength - length
        D = self.OD - self.wireDia
        C = D/self.wireDia
        K = (4*C-1)/(4*C-4) + 0.615/C # Wahl stress correction factor: Table 5.2 Pg 190
        stress = 8*self.rate*D*K*deflection/(pi*self.wireDia**3)
        return stress
    
    def getMinTensileStrength(self):
        offset = bisect.bisect_right(tensileStrengths['WireDia'], self.wireDia,
                                     hi=len(tensileStrengths['WireDia'])-1)
        minTens = tensileStrengths[self.material][offset]
        return minTens
    
    def getMillionDefl(self):
        stress = self.minTens * self.reductionFactor
        
        D = self.OD - self.wireDia
        C = D/self.wireDia
        K = (4*C-1)/(4*C-4) + 0.615/C # Wahl stress correction factor: Table 5.2 Pg 190
        
        return pi*self.wireDia**3*stress/(8*self.rate*D*K)
        
    # def getnCycles(self, deflection = None):
    #     '''
    #     This spring life calculation should be accurate for the life of plain carbon
    #     springs when the number of cycles is between 5e5 and 1e7 [Stone]. Above 
    #     2e6 these springs should experience near infinite life.
    #     For stainless springs these values become less accurate.
        
    #     '''
    #     if deflection is None:
    #         deflection = self._getMillionDefl()

    #     stress = self.getStress(deflection = deflection)
    #     minTens = self.getMinTensileStrength()
    #     reductionFactor = strenghtReductionFactors[self.material]
        
    #     K_s1 = 0
    #     K_s2 = stress/minTens
        
    #     K_U = 0.56
    #     C_S = 0.5546
    #     M = -0.009
    #     C_E = 0.662
    #     Y = -0.0622
        
    #     K_E = K_U*(K_s2 - K_s1)/(2*K_U - (K_s2+K_s1))
                
    #     n = (C_E/K_E)**(1/-Y)
    #     return n

    def getLifeFOS(self, minLife):
        if self.lastLength >= self.freeLength:
            return '-'
        
        life = lifeDict[minLife]
        
        stress = self.getStress(self.lastLength)
        factor = 10*(-stress/self.minTens + self.reductionFactor)
        
        if life == L_10e5:
            if factor > 0:
                return self.expectedLife(self.lastLength) - L_10e5 + factor
            else:
                millDeflection = self.getMillionDefl()
                currDeflection = self.freeLength - self.lastLength
                return millDeflection/currDeflection
        
        if life == L_MILLION:
            return factor
        return factor - 1

    def expectedLife(self, length):
        self.lastLength = length
        deflection = self.freeLength - length
        if deflection > self.freeLength - self.solidLength:
            return NOT_POSSIBLE
        if deflection > self.maxDeflection:
            return NOT_RECOMMENDED
        
        stress = self.getStress(length)
        
        if stress < self.minTens * (self.reductionFactor - 0.1):
            return L_INF
        if stress < self.minTens * self.reductionFactor:
            return L_MILLION
        return L_10e5
    
    def __repr__(self):
        return 'Spring(' + self.name + ')'

        

with open(SPRING_FILE, 'r') as f:
    for line in f:
        temp = line.split(sep=' ')
        if len(temp) == 21:            
            springs.append(Spring(temp))
        else:
            print(temp)

def getFourDeflections(spring, stroke, minF1, maxF1, minF2, maxF2):
    d1 = spring.getDeflection(minF1)    
    d2 = min(x for x in (spring.getDeflection(maxF1), spring.maxDeflection) if x is not None)     
    d3 = spring.getDeflection(minF2)
    d4 = min(x for x in (spring.getDeflection(maxF2), spring.maxDeflection) if x is not None)
    
    return d1, d2, d3, d4

def strokeFilter(inSprings, stroke, minF1, maxF1, minF2, maxF2):
    outSprings = []
    for s in inSprings:
        d1, d2, d3, d4 = getFourDeflections(s, stroke, minF1, maxF1, minF2, maxF2)
        if d1 is None or d3 is None:
            continue
        
        minDeltaD = d3-d2
        maxDeltaD = d4-d1
        if minDeltaD < stroke and maxDeltaD > stroke:
            outSprings.append(s)
    return outSprings

springDict = {s.name:s for s in springs}

def rangeAttrFilter(lsprings, attr, minN, maxN):
    return [s for s in lsprings if minN <= s.__getattribute__(attr) <= maxN]

def lengthFilter(lsprings, length):
    return [s for s in lsprings if s.minLength <= length <= s.freeLength]

def forceFilter(lsprings, length, minF, maxF):
    return [s for s in lsprings if minF <= s.getForce(length) <= maxF]

def selectionFilter(lsprings, attr, value):
    return [s for s in lsprings if s.__getattribute__(attr) == value]

def lifeFilter(inSprings, minLife, *, length1=None, length2=None, stroke=None, forces=None):
    if stroke and forces:
        outSprings = []
        for spring in inSprings:
            d1, d2, d3, d4 = getFourDeflections(spring, stroke, *forces)
            if d1 + stroke > d3:
                length1 = spring.freeLength - d1
                length2 = spring.freeLength - (d1 + stroke)
            else:
                length2 = spring.freeLength - d3
                length1 = length2 + stroke
            if spring.expectedLife(length2) >= lifeDict[minLife]:
                outSprings.append(spring)
        return outSprings
    try:
        minLength = min(l for l in [length1, length2] if l is not None)
        return [s for s in inSprings if s.expectedLife(minLength) >= lifeDict[minLife]]
    except Exception:
        return inSprings
