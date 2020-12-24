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

r = open('r2.in','w+') #beware!!! using the same seed for real and sky
    
r.writelines(['/run/initialize \n \n', \
              '/SetSeed/set {} \n \n'.format(''.join(random.choices(string.digits, k=10))), \
              '/gun/particle mu- \n \n', \
              '#/vis/scene/endOfEventAction accumulate \n \n', \
              '/run/beamOn {} \n \n'.format(Events)])

r.close()

