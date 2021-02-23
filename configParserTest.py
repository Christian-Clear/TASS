#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 13:28:56 2021

@author: christian
"""

import configparser

config = configparser.ConfigParser()

config.read('test.ini')

# print(config.getfloat('strans', 'wn_discrim'))

# config.add_section('project')
# config['project']['project_config'] = 'test.ini'

# config.add_section('strans')
# config['strans']['wn_discrim'] = '0.05'


# with open('tame.ini', 'w') as configfile:
#     config.write(configfile)
    
    
config['files']['other_lev_files'] = 'ni1,ni1_input.lev\nhe1,he1_input.lev'


with open('test.ini', 'w') as configfile:
    config.write(configfile)