#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 15:28:13 2021

@author: christian
"""

import bisect

test_array = [(1,2),(3,4),(5,9),(5,6),(5,7000),(7,8),(9,10)]

class KeyList(object):
    # bisect doesn't accept a key function, so we build the key into our sequence.
    def __init__(self, l, key):
        self.l = l
        self.key = key
    def __len__(self):
        return len(self.l)
    def __getitem__(self, index):
        return self.key(self.l[index])
    
    
left = bisect.bisect_left(KeyList(test_array, key=lambda x: x[0]), 4.9)
right = bisect.bisect_left(KeyList(test_array, key=lambda x: x[0]), 5.1)

print(test_array[left:right])

