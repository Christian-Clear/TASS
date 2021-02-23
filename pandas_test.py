#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 15:49:01 2021

@author: christian
"""

import numpy as np
import pandas as pd
import bisect

class KeyList(object):
    # bisect doesn't accept a key function, so we build the key into our sequence.
    def __init__(self, l, key):
        self.l = l
        self.key = key
    def __len__(self):
        return len(self.l)
    def __getitem__(self, index):
        return self.key(self.l[index])





# df = pd.read_csv('ni2_input.lin')
# df.to_pickle('panda_test.pkl')
strans_lev_file = 'ni2_input.lev'
df = pd.read_pickle('panda_test.pkl')
#df.set_index('wavenumber')


strans_levs = np.loadtxt(strans_lev_file, dtype=str, delimiter=',', skiprows=1)
strans_levs_even = [x for x in strans_levs if x[-1]=='1']
strans_levs_even = sorted(strans_levs_even, key=lambda x: x[1])
strans_levs_odd = [x for x in strans_levs if x[-1]=='0']
strans_levs_odd = sorted(strans_levs_odd, key=lambda x: x[1])

strans_wn_discrim = .05

df['main_desig'] = np.empty((len(df), 0)).tolist()

desig_list = df[['wavenumber', 'main_desig']].values.tolist()

#print(strans_levs_even)
  
for even_lev in strans_levs_even:    
    j_even = float(even_lev[1])
    
    left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: float(x[1])), j_even - 1)
    right_j = bisect.bisect_right(KeyList(strans_levs_odd, key=lambda x: float(x[1])), j_even + 1)
                    
    for odd_lev in strans_levs_odd[left_j:right_j]:
        j_odd = float(odd_lev[1])                    
        label_even = even_lev[0]
        label_odd = odd_lev[0]
        energy_even = float(even_lev[2])
        energy_odd = float(odd_lev[2])                    
        match_wn = abs(energy_even - energy_odd)
            
        left = bisect.bisect_left(KeyList(desig_list, key=lambda x: x[0]), match_wn - strans_wn_discrim)
        right = bisect.bisect_left(KeyList(desig_list, key=lambda x: x[0]), match_wn + strans_wn_discrim)
            
        for matched_line in desig_list[left:right]:
            matched_line[1].append({'even_level': label_even, 'odd_level':label_odd})
    
df.update(desig_list)

print(df.tail())
#print(desig_list)
            
            
            
            
            
            
            

















            
            
            # matching_lines.append((match_wn, label_even, label_odd))
            
            # left_i = df['wavenumber'].searchsorted(match_wn - strans_wn_discrim, side = 'left')
            # right_i = df['wavenumber'].searchsorted(match_wn - strans_wn_discrim, side = 'right')
            
            
            # # #print(left_i, right_i)
            # for x in df.iloc[left_i : right_i]:
            #     print(x['main_desig'])
            
            
            
            #matched_df = df[(df.wavenumber >= match_wn - strans_wn_discrim) & (df.wavenumber <= match_wn + strans_wn_discrim)].copy()  # returns dataframe of matching lines
            # matched_df.main_desig = matched_df.main_desig.apply(lambda x: x+[{'level_even':label_even, 'level_odd':label_odd}])
            # df.update(matched_df)
 








           
# matching_lines = [x for x in matching_lines if x[0] >= df['wavenumber'].min() - strans_wn_discrim and x[0] <= df['wavenumber'].max() + strans_wn_discrim]

# # for matching_line in matching_lines:
# #       matched_df = df.loc[(df.wavenumber >= match_wn - strans_wn_discrim) & (df.wavenumber <= match_wn + strans_wn_discrim)].copy()  # returns dataframe of matching lines
# #       #matched_df.main_desig = matched_df.main_desig.apply(lambda x: x+[{'level_even':label_even, 'level_odd':label_odd}])
# #       #df.update(matched_df)
      
    
           
# #print(matched_df)
# print(df)
#print(df.loc[df.main_desig.str.len() > 0 ])

# lines = df.loc[df.main_desig.str.len() > 0 ].values.tolist()

# for line in lines:
#     level_string = ''
#     for lev_pair in line[-1]:
#         if not level_string:
#             sep = ''
#         else:
#             sep = ',     '
#         level_string += sep + lev_pair['level_even'] + ' - ' + lev_pair['level_odd']
        
#     line[-1] = level_string
    



























