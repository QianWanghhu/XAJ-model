import numpy as np
import pandas as pd
import json

# prepare input data
def read_inputs(fn, weights=None):
    """
    Read the input of P, E
    """
    print('================Load the required input data===============')
    data_inputs = pd.read_csv(fn)
    rain_col = [col for col in data_inputs.columns if 'P' in col]
    if len(rain_col) > 1:
        assert weights.shape[0] == len(rain_col), "Each rain gauge should have a weight"
    
    P, E, Q_obs = data_inputs.loc[:, rain_col].values, \
        data_inputs.loc[:, 'E'].values, data_inputs.loc[:, 'Q'].values

    P_weighted = (P * weights).sum(axis=1)

    return P_weighted, E, Q_obs

fn = 'data.csv'
weights = np.loadtxt('raingauge_weights.txt')    
P, E, Q_obs = read_inputs(fn, weights)

from runoff_generation import *
f_parms = 'runoff_params.json'
f_initials = 'init_states.json'
S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD = \
    read_parms(f_parms, f_initials)


WM, WMM, MS, EP, KI, KG = \
    calculate_capacities(WDM, WUM, WLM, SM, b, EX, kc, KI,  KG, E)

WL, WU, WD, S0 = soil_variable(WLM, WUM, WDM, WU, WL, WD, S0, SM)


RSS, RII, RGG, e_list, WU1, WL1, WD1, R, Rb = \
    runoff_production(P, E, S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD)