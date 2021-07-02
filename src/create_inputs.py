import numpy as np
import pandas as pd
import json

pr_parms = {
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

routing_parms = {
    'CS': CS,
    'CI': CI,
    'CG': CG,
    'L': L,
    'x1': x1,
    'K': K,
    'xx': xx,
    'periods': periods,
    'AA': AA
}

json.dump(pr_parms, open(f'../data/pr_parms.json', 'wb'),
        indent = 2)

json.dump(inital_state, open(f'../data/initial_states.json', 'wb'),
        indent = 2)

json.dump(routing_parms, open(f'../data/routing_parms.json', 'wb'),
        indent = 2)
