# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 08:32:03 2019

gasSpring.py

@author: Myself
"""

from pint import UnitRegistry
import math
import numpy as np

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

assemblyCOM_radius = (assemblyCOM_X**2 + assemblyCOM_Y**2)**0.5
assemblyStartAngle = math.atan2(assemblyCOM_Y.magnitude, assemblyCOM_X.magnitude) * ureg.rad

movingMnt_X = 5.875 * ureg.inch
movingMnt_Y = -1.1 * ureg.inch
movingMnt_radius = (movingMnt_X**2 + movingMnt_Y**2)**0.5
movingMntStartAngle = math.atan2(movingMnt_Y.magnitude, movingMnt_X.magnitude) * ureg.rad

fixedMnt_X = -2 * ureg.inch
fixedMnt_Y = -18.9 * ureg.inch

deg = np.linspace(0, -90, 19) * ureg.deg
radians = deg.to('rad')

aCOM_X = assemblyCOM_radius * np.cos(assemblyStartAngle + radians)
aCOM_Y = assemblyCOM_radius * np.sin(assemblyStartAngle + radians)

mX = movingMnt_radius * np.cos(movingMntStartAngle + radians)
mY = movingMnt_radius * np.sin(movingMntStartAngle + radians)

springLength = np.sqrt((mX-fixedMnt_X)**2 + (mY-fixedMnt_Y)**2)

