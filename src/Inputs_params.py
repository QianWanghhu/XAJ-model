import numpy as np
import pandas as pd
import json
from runoff_generation import *
from routing import *


fpath = '../data/'
fn = 'data.csv'
weights = np.loadtxt(f'{fpath}raingauge_weights.txt')    
P, E, Q_obs = read_inputs(fpath + fn)

# The initial states of the catchment and runoff generation parameters
f_parms = 'runoff_params.json'
f_initials = 'init_states.json'
S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD = \
    read_parms(f_parms, f_initials)

# routing parameters
fn_routing = 'routing_params.json'
CS, CI, CG, L, n, K, xx, dt, C0, C1, C2, periods, AA, NT = routing_params(fn_routing)
