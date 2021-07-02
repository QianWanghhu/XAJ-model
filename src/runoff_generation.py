# import packages
import numpy as np
import json
import pandas as pd

def read_parms(f_parms, f_initials):
    print('================Load the required parameters===============')
    user_parms = json.load(open(f_parms, 'rb')) # f'../data/input_parms.json'
    print('================Load the required initial state variables===============')
    initial_states = json.load(open(f_initials, 'rb'))
    # read parameter values
    im = user_parms['im']
    WDM = user_parms['WDM']
    WUM = user_parms['WUM']
    WLM = user_parms['WLM']
    SM = user_parms['SM']
    b = user_parms['b']
    EX = user_parms['EX']
    kc = user_parms['kc']
    C = user_parms['C']
    KI = user_parms['KI']
    KG = user_parms['KG']

    # read initial states
    S0 = initial_states['S0']
    FR0 = initial_states['FR0']
    WU = initial_states['WU']
    WL = initial_states['WL']
    WD = initial_states['WD']
    return S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD


def calculate_capacities(WDM, WUM, WLM, SM, b, EX, kc, KI,  KG, E):
    WM = WLM + WUM + WDM # 流域平均蓄水容量
    WMM = WM * (b + 1.0) #流域单点最大蓄水量
    MS = SM * (EX + 1.0) #流域最大自由水蓄量
    EP = kc * E #计算蒸散发能力

    KIt = (1.0 - (1.0 - KI - KG) ** (1.0 / (24.0 * 1.0)))/(1.0 + KG/KI) # 划分小段后的参数值
    KGt = KG * KIt / KI
    KI = KIt
    KG = KGt

    return WM, WMM, MS, EP, KI, KG

def soil_variable(WLM, WUM, WDM, WU, WL, WD, S0, SM):
    """
    Analyze the inputs for soil water to make sure that the input variables are of consistence.
    """
    # user_parms = json.load(open(f'../data/input_parms.json', 'rb'))
    # WM, WMM, MS, EP, KIt, KGt, KI, KG  = calculate_capacities(user_parms)

    if WU > WUM:
        WL = WL + WU - WUM
        WU = WUM

    if WL > WLM:
        WD = WD + WL - WLM
        WL = WLM
    
    if WD > WDM: WD = WDM

    if S0 > SM: S0 = SM

    return WL, WU, WD, S0
# END soil_variable()

def update_soil_state(PEt, WUt, WLt, WDt, Rt, WUM, WLM, WDM, W0):

    """
    Update the soil water states.
    Parameters:
    ------------

    Return:
    -------
    WUt, WLt, WDt: float, t时刻对应的上层，下层和深层蓄水容量
    """
    if WUt + PEt - Rt > WUM:
        if (WUt + WLt + PEt - Rt - WUM > WLM):
            WUt = WUM
            WLt = WLM
            WDt = W0 + PEt - Rt - WUt - WLt
            if WDt > WDM: WDt = WDM
            # END if
        else:
            WLt = WUt + WLt + PEt - Rt - WUM
            WUt = WUM
        # END if-else
    else:
        WUt = WUt  + PEt - Rt

    return WUt, WLt, WDt

def runoff_single_time(PEt, im, WMM, W0, WM, b):
    """
    t时刻流域蓄满产流。
    Parameters:
    ------------
    PEt: float, t时刻的净雨
    im: float, 不透水面积比例
    WMM: float, 流域单点最大蓄水量
    W0: float, 上层土壤含水量
    WM: float, 流域平均蓄水容量
    b: float, 反映流域包气带蓄水容量分布的不均匀性

    Return:
    --------
    Rbt: runoff at the impervious area.
    Rt: runoff generated in the pervious region.
    """
    if PEt > 0:
        Rbt = im * PEt   #不透水面积产流
        a_t = WMM * (1.0 - (1.0 - W0 / WM) ** (1.0 / (1.0 + b))) #W0对应的纵坐标
        if ((a_t + PEt) < WMM):
            Rt = PEt - WM + W0 + WM * ((1 - (PEt + a_t) / WMM) ** (b + 1.0))
        else:
            Rt = PEt - (WM - W0)
        # END if-else
    else:
        Rt = 0.0
        Rbt = 0.0
    # END if-else

    return Rbt, Rt

def separate_source(Rt, PEt, S0t, FR0t, KI, KG, MS, EX):
    RS, RI, RG = 0.0, 0.0, 0.0
    if Rt > 0:
        FRt = Rt / PEt
        S = S0t * FR0t / FRt
        QD = Rt / FRt
        N = np.floor(QD / 5) + 1 # 将每个计算时段的入流R，按5mm分成N段
        KIt = (1.0 - (1.0 - KI - KG)**(1.0 / (N * 1.0))) / (1.0 + KG / KI)
        KGt = KG * KIt / KI
        QD = QD / N
        SMM1 = MS
        SM1=SMM1 / (1.0 + EX)
        SS=S

        for i in range(N):
            if (S > SM1):
                AU = SMM1
            else:
                AU = SMM1 * (1.0 - (1.0 - S/SM1)**(1.0 / (1.0 + EX))) #S0对应的纵坐标
            # END if-else

            if (QD + AU < SMM1):
                RS = FRt * (QD + S - SM1 + SM1 * ((1.0 - (QD + AU) / SMM1)**(EX + 1.0))) + RS
            else:
                RS = FRt * (QD + S - SM1) + RS
            # END if-else

            S = i * QD - RS / FRt + S
            RI = KIt * FRt * S + RI
            RG = KGt * FRt * S + RG
            S = SS + i * QD - (RS + RI + RG) / FRt

            if S > SM1: S = SM1

            if S < 0: S = 0.0
        # END for
        S0t = S, FR0t = FRt
    else:
        RG = S0t * KG * FR0t
        RI = S0t * KI * FR0t
        S0t = S0t - (RG + RI) / FR0t
    # END if-else

    return RS, RI, RG, S0t, FR0t

def soil_evpa(Pt, EPt, WUt, WLt, WDt, C, WLM):
    """
    三层蒸散发模型。
    Parameters: 
    -----------
    Pt: float, t时刻的降水量 
    EPt: float, t时刻对应的蒸发能力
    WUt, WLt, WDt: float, t时刻对应的地表，浅层和深层土壤含水量
    C: float, 浅层蒸散发折算系数
    WLM: float, 下层土壤水含水容量

    Return:
    --------
    EUt, EDt, ELt: float, t时刻对应的地表，浅层和深层蒸发水量
    """
    if Pt - EPt < 0.0:
        if ((WUt + Pt) > EPt):
            EUt = EPt
            ELt = 0.0
            EDt = 0.0
            WUt = WUt + Pt -EPt
        else:
            EUt = WUt + Pt
            WUt = 0.0
            if (WLt > C * WLM):
                ELt = (EPt - EUt) * WLt / WLM
                EDt = 0.0
                WLt = WLt - ELt
            else:
                if (WLt > C * (EPt - EUt)):
                    ELt = C * (EPt - EUt)
                    EDt = 0.0
                    WLt = WLt - ELt
                else:
                    ELt = WLt
                    WLt = 0.0
                    EDt = C * (EPt - EUt) - WLt
                    WDt = WDt - EDt
                    if WDt <= 0:
                        EDt = EDt + WDt
                        WDt = 0.0
                    # END if-else
                # END if-else    
            # END if-else
        # END if-else
    # END if-else
    
    return EUt, ELt, EDt

def runoff_production(P, E, S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD):
    # T = 1
    # read inputs and calculate related parameters
    WM, WMM, MS, EP, KI, KG = calculate_capacities(WDM, WUM, WLM, SM, b, EX, kc, KI,  KG, E)
    WL, WU, WD, S0 = soil_variable(WLM, WUM, WDM, WU, WL, WD, S0, SM)

    # set the intial state values
    WL1, WD1, WU1 = np.zeros_like(P), np.zeros_like(P), np.zeros_like(P)
    RSS, RII, RGG = np.zeros_like(P), np.zeros_like(P), np.zeros_like(P) 
    WL1[0], WD1[0], WU1[0] = WL, WU, WD
    RSS[0]=0.0; RII[0]=0.0; RGG[0]=0.0

    for t in range(1, P.shape[0]):
        W0 = WL + WU + WD  # t时段的初始土壤含水量
        if (W0 > WM): W0 = WM

        # 流域蒸散发，三层蒸发模型
        if (P[t] - EP[t]) < 0.0:
            EU, EL, ED = soil_evpa(P[t], EP[t], WU, WL, WD, C, WLM)
            UE = EU + ED + EL
            PE = P[t] - UE #扣除蒸发的净雨
        else:
            PE = P[t] -EP[t]

        if PE < 0: PE = 0.0

        # 流域产流蓄满产流
        Rb, R = runoff_single_time(PE, im, WMM, W0, WM, b)
        if (np.abs(R) < 1E-10): R = 0.0

        # UPDATE SOIL WATER
        if PE > 0:
            WU, WL, WD = update_soil_state(PE, WU, WL, WD, R, WUM, WLM, WDM, W0)

        # 分水源
        # TO　ADD　separate_source()
        RS, RI, RG, S0, FR0 = separate_source(R, PE, S0, FR0, KI, KG, MS, EX)

        RS = Rb + RS * (1.0 - im)
        RG = RG * (1.0 - im)
        RI = RI * (1.0 - im)

        WU1[t] = WU
        WL1[t] = WL
        WD1[t] = WD
        RSS[t] = RS
        RII[t] = RI
        RGG[t] = RG

    return RSS, RII, RGG
    # END runoff_production()


def read_inputs(fn):
    """
    Read the input of P, E
    """
    print('================Load the required input data===============')
    PE_inputs = np.loadtxt(fn, skiprows=1, usecols=[1, 2])
    P, E = PE_inputs[:, 0], PE_inputs[:, 1]

    return P, E

def save_runoff(fn, f_parms, f_initials, runoff_file):
    # read data
    P, E = read_inputs(fn)
    S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD = \
        read_parms(f_parms, f_initials)

    RSS, RII, RGG = \
        runoff_production(P, E, S0, FR0, im, WDM, WUM, WLM, SM, b, EX, kc, C, KI,  KG, WU, WL, WD)
    df = pd.DataFrame(data=np.array([RSS, RII, RGG]).T, columns=['RSS', 'RII', 'RGG'] )
    df.to_csv(runoff_file)
    return df