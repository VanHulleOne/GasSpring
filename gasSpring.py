# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 08:32:03 2019

gasSpring.py

@author: Myself
"""

from pint import UnitRegistry
import math
import numpy as np

import matplotlib.pyplot as plt

ureg = UnitRegistry()

lengthMax = 20.6 * ureg.inch
lengthMin = 12.76 * ureg.inch

forceMax = 157 * ureg.lbs
forceMin = 100 * ureg.lbs

springRate = (forceMax - forceMin) / (lengthMax - lengthMin)

numSprings = 2

assemblyMass = 141 * ureg.lbs
assemblyCOM_X = 6.0 * ureg.inch
assemblyCOM_Y = 1.5 * ureg.inch
assemblyLenth = 14 * ureg.inch

assemblyCOM_radius = (assemblyCOM_X**2 + assemblyCOM_Y**2)**0.5
assemblyStartAngle = math.atan2(assemblyCOM_Y.magnitude, assemblyCOM_X.magnitude) * ureg.rad

movingMnt_X = 5.875 * ureg.inch
movingMnt_Y = -1.1 * ureg.inch
movingMnt_radius = (movingMnt_X**2 + movingMnt_Y**2)**0.5
movingMntStartAngle = math.atan2(movingMnt_Y.magnitude, movingMnt_X.magnitude) * ureg.rad

fixedMnt_X = -2 * ureg.inch # -2
fixedMnt_Y = -18.9 * ureg.inch

deg = np.linspace(0, -90, 19) * ureg.deg
radians = deg.to('rad')

aCOM_X = assemblyCOM_radius * np.cos(assemblyStartAngle + radians)
aCOM_Y = assemblyCOM_radius * np.sin(assemblyStartAngle + radians)

mX = movingMnt_radius * np.cos(movingMntStartAngle + radians)
mY = movingMnt_radius * np.sin(movingMntStartAngle + radians)

springLength = np.sqrt((mX-fixedMnt_X)**2 + (mY-fixedMnt_Y)**2)
springForce = (lengthMax - springLength)*springRate + forceMin

assemblyTorque = aCOM_X * assemblyMass

springMomentArm = (mX*fixedMnt_Y - mY*fixedMnt_X)/springLength
springTorque = springMomentArm * springForce * numSprings

gearboxTorque = assemblyTorque + springTorque


fig, ax = plt.subplots()

#ax.plot(deg, assemblyTorque / assemblyLenth, label = 'Assembly')
#ax.plot(deg, springTorque / assemblyLenth, label = 'Spring')
#ax.plot(deg, gearboxTorque / assemblyLenth, label = 'Gearbox')
#
#ax.set(xlabel='Angle (deg)', ylabel='Force (lb)',
#       title='Table Edge Force Values')
#ax.grid()
#
#plt.gca().invert_xaxis()
#plt.legend()
#plt.show()

def validFixed(step=0.2, box=6, tol = 0.125*ureg.inch):
    goodVals = []
    xVals = np.linspace(fixedMnt_X.magnitude-box, fixedMnt_X.magnitude+box, 2*box/step+1)*ureg.inch
    yVals = np.linspace(fixedMnt_Y.magnitude-box, fixedMnt_Y.magnitude+box, 2*box/step+1)*ureg.inch
    for xVal in xVals:
        for yVal in yVals:
            
            springLength = np.sqrt((mX-xVal)**2 + (mY-yVal)**2)
            _min = np.amin(springLength)
            _max = np.amax(springLength)
            
            if tol+lengthMin < _min and lengthMax-tol > _max:
                goodVals.append((xVal, yVal))
                
    return goodVals

f = validFixed()
fX = [x[0].magnitude for x in f]
fY = [y[1].magnitude for y in f]
ax.scatter(fX, fY, c='blue')
plt.show()

