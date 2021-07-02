# import packages 
import json
import numpy as np

def routing_parms(f_parms):
    routing_parms = json.load(open(f_parms, 'rb'))
    CS = routing_parms['CS']
    CI = routing_parms['CI']
    CG = routing_parms['CG']
    L = np.floor(routing_parms['L'])
    n = np.floor(routing_parms['x1'])
    K = routing_parms['K']
    xx = routing_parms['xx']
    periods = routing_parms['periods']
    dt = K
    AA = routing_parms['AA']

    C0 = (0.5 * dt - K * xx)/(0.5 * dt + K - K * xx)
    C1 = (0.5 * dt + K * xx)/(0.5 * dt + K - K * xx)
    C2 = (-0.5 * dt + K - K * xx)/(0.5 * dt + K - K * xx)
    CI = CI ** (1.0 / 24.0)
    CG = CG ** (1.0 / 24.0)
    return CS, CI, CG, L, n, K, xx, dt, C0, C1, C2, periods, AA

def slope_runoff(QS, QI, QG, RS, RI, RG, CI, CG, AA, periods, NT):
    scale_factor = AA / 3.6 / periods
    for tt in range(1, NT):
        QS[tt] = RS[tt] * scale_factor
        QI[tt] = CI * QI[tt - 1] + (1.0 - CI) * RI[tt] * scale_factor
        QG[tt] = CG * QG[tt] + (1.0 - CG) *RG[tt] * scale_factor
    
    return QS, QI, QG


def network_runoff(NT, L, QQ, QT, CS):
    for tt  in range(1, NT - 1):
        if ((tt - L) > 1):
            QQ[tt + 1] = CS * QQ[tt] + (1.0 - CS) * QT[tt - L]
    return QQ


def muskingum():
    pass

def routing(Q_obs, NT, CS, CI, CG, L, n, K, xx, dt, C0, C1, C2, periods, RI, RG, RS, AA):
    QS, QI, QG, QQ, QT = np.zeros(NT), np.zeros(NT), \
        np.zeros(NT), np.zeros(NT), np.zeros(NT)
    QG[0] = Q_obs[0]
    QQ[0] = Q_obs[0]
    
    QS, QI, QG = slope_runoff(QS, QI, QG, RS, RI, RG, CI, CG, AA, periods, NT)

    QT = QS +QI + QG

    QQ = network_runoff(NT, L, QQ, QT, CS)

    # Muskingum


