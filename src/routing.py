# import packages 
import json
import numpy as np

def routing_params(f_parms):
    routing_parameters = json.load(open(f_parms, 'rb'))

    CS = routing_parameters['CS'] # 地面径流消退系数
    CI = routing_parameters['CI'] # 壤中流消退系数
    CG = routing_parameters['CG'] # 地下径流消退系数
    L = int(np.floor(routing_parameters['L']))
    n = np.floor(routing_parameters['x1'])
    n = [int(i) for i in n]
    K = routing_parameters['K'] # 槽蓄系数
    xx = routing_parameters['xx'] # 流量比重系数
    periods = routing_parameters['periods']
    dt = K
    AA = routing_parameters['AA'] # 流域面积
    NT = routing_parameters['NT']

    C0 = (0.5 * dt - K * xx)/(0.5 * dt + K - K * xx)
    C1 = (0.5 * dt + K * xx)/(0.5 * dt + K - K * xx)
    C2 = (-0.5 * dt + K - K * xx)/(0.5 * dt + K - K * xx)
    CI = CI ** (1.0 / 24.0)
    CG = CG ** (1.0 / 24.0)
    return CS, CI, CG, L, n, K, xx, dt, C0, C1, C2, periods, AA, NT

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
            # breakpoint()
            QQ[tt + 1] = CS * QQ[tt] + (1.0 - CS) * QT[int(tt - L)]
    return QQ


def muskingum(Q_1, QQ, n, C0, C1, C2):
    """
    Parameters:
    ===========
    Q_1: the first observed Q
    QQ: simulated Q

    Return:
    =======
    QC: 
    """
    if n > 0:
        QC2 = np.zeros(QQ.shape[0])
        QC3 = np.zeros(QQ.shape[0])
        QC3[0] = Q_1
        QC2 = QQ
        for i in range(n): # TODO: check the value of n whether starting from 0 or 1
            for t in range(QQ.shape[0] - 1):
                # breakpoint()
                QC3[t+1]= C0 * QC2[t+1] + C1 * QC2[t] + C2 * QC3[t]
            QC2 = QC3
        QC = QC2        
    else:
        QC = QQ
    return QC
# End muskingum()

def routing(Q_obs, NT, CS, CI, CG, L, n, C0, C1, C2, periods, RI, RG, RS, AA):
    QS, QI, QG, QQ, QT = np.zeros(NT), np.zeros(NT), \
        np.zeros(NT), np.zeros(NT), np.zeros(NT)
    QG[0] = Q_obs[0]
    QQ[0] = Q_obs[0]
    print('------------Run slope_runoff--------------')
    QS, QI, QG = slope_runoff(QS, QI, QG, RS, RI, RG, CI, CG, AA, periods, NT)

    QT = QS + QI + QG
    print('------------Run network_runoff--------------')
    QQ = network_runoff(NT, L, QQ, QT, CS)
    print('------------Run muskingum--------------')
    QC = muskingum(Q_obs[0], QQ, n, C0, C1, C2)
    
    return QC

# End routing()


