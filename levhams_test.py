#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:54:13 2021

@author: christian
"""

import pandas as pd

def levhams_match_tol(line_list, tol):
    """Generator function to split predicted lines groups with each element sperarated from its neighbour 
    by <= a given tolerance.
    """
    matches = []
    last = line_list[0]['wavenumber']
    
    for element in line_list:
        if element['wavenumber'] - last > tol:
            yield matches
            matches = []
            
        matches.append(element)
        last = element['wavenumber']
    yield matches
    

strans_levs = list(pd.read_csv('ni2_input_levhams.lev', dtype={'parity':float}).transpose().to_dict().values()) 
levels = [x for x in strans_levs[:36]]

lines = list(pd.read_csv('ni2_input_cutdown.lin', dtype={'parity':float}).transpose().to_dict().values()) 
lines = [x['wavenumber'] for x in lines]

pred_lines = []



for level in levels:
    for line in lines:
        pred_lines.append({'wavenumber': level['energy'] - line, 'level':level['label']})
        pred_lines.append({'wavenumber': level['energy'] + line, 'level':level['label']})


pred_lines = sorted(pred_lines, key=lambda k: k['wavenumber'])

pred_levels = list(levhams_match_tol(pred_lines, 0.001))

for level in pred_levels:
    if len(level) > 5:
        sum_lines = 0.
        print('-------------------------------------')
        for line in level:
            print(line)
            #print('\n')
            sum_lines += line['wavenumber']
        print('                                     ' + str(sum_lines/len(level)) + '   ' + str(abs(level[0]['wavenumber'] - level[-1]['wavenumber'])))
    
