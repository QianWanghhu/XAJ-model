import numpy as np
import pandas as pd
import json

user_parms = {
    'im': im,
    'WDM': WDM,
    'WUM': WUM,
    'WLM': WLM,
    'SM': SM,
    'b' : b,
    'EX': EX,
    'kc': kc,
    'C': C,
    'KI': KI,
    'KG': KG,
    'WU': WU,
    'WL': WL,
    'WD': WD
}

inital_state = {
    'S0': S0,
    'FR0': FR0,
    'WU': WU,
    'WL': WL,
    'WD': WD
}

json.dump(user_parms, open(f'../data/input_parms.json', 'wb'),
        indent = 2)

json.dump(inital_state, open(f'../data/initial_states.json', 'wb'),
        indent = 2)
