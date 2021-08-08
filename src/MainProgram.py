import numpy as np
import pandas as pd
import json
from runoff_generation import *
from routing import *
from Inputs_params import *

def run_XAJ():
    
    NP = P.shape[1] # 雨量站数目
    FloodTime = np.zeros(shape=(NT, NP))
    QCE = np.zeros(NT)
    for i in range(NP):
        print(f'===================Subcatchment {i}===================')
        RSS, RII, RGG = \
        runoff_production(P[:, i], E, S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD)
        
        QC1 = routing(Q_obs, NT, CS, CI, CG, L, n[i], C0, C1, C2, periods, RII, RGG, RSS, AA)
        QCE = QCE + QC1 * weights[i]
    np.savetxt('Q.txt', QCE)
    # End for

if __name__ == "__main__":
    run_XAJ()