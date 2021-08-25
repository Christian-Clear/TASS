#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:54:13 2021

@author: christian
"""

import pandas as pd

def split_tol(test_list, tol):
    res = []
    last = test_list[0]
    for ele in test_list:
        if ele-last > tol:
            yield res
            res = []
        res.append(ele)
        last = ele
    yield res
    

strans_levs = list(pd.read_csv('ni2_input_levhams.lev', dtype={'parity':float}).transpose().to_dict().values()) 
energies = [x['energy'] for x in strans_levs[:36]]

lines = list(pd.read_csv('ni2_input_cutdown.lin', dtype={'parity':float}).transpose().to_dict().values()) 
lines = [x['wavenumber'] for x in lines]


pred_lines = []

for energy in energies:
    for line in lines:
        pred_lines.append(abs(energy - line))
        
        

pred_lines = list(dict.fromkeys(sorted(pred_lines))) 

res = list(split_tol(pred_lines, 0.01))

# print(res)

for lev in res:
    if len(lev) > 2:
        print(lev)