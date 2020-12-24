#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 20:16:59 2020

@author: keegan
"""

import random
import string
from joblib import load

Events = load('Events.joblib')

seed = ''.join(random.choices(string.digits, k=8)) 

f = open('Seed_History.txt','a+')

f.write(seed + ' \n')

f.close()

r = open('r2.in','w+') 
    
r.writelines(['/run/initialize \n \n', \
              '/SetSeed/set {} \n \n'.format(seed), \
              '/gun/particle mu- \n \n', \
              '#/vis/scene/endOfEventAction accumulate \n \n', \
              '/run/beamOn {} \n \n'.format(Events)])

r.close()

