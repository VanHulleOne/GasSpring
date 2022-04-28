# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 08:08:09 2021

testSprings.py

@author: lvanhulle
"""

import unittest
import springs as sp

class Tests(unittest.TestCase):
    
    def testRangeAttrFilter(self):
        l = [2,3,4,5,
             1,0,-1,
             6,7,8]
        res = sp.rangeAttrFilter(l, 'numerator', 2,5)
        self.assertEqual(l[:4], res)


if __name__ == '__main__':
    unittest.main()