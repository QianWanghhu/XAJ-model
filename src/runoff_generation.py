# import packages
import numpy as np
import json

def calculate_capacities(user_parms):
    WM = user_parms['WLM'] + user_parms['WUM'] + user_parms['WDM'] # 流域平均蓄水容量
    WMM = WM * (user_parms['b'] + 1.0) #流域单点最大蓄水量
    MS = user_parms['SM'] * (user_parms['EX'] + 1.0) #流域最大自由水蓄量
    EP = user_parms['kc'] * user_parms['E'] #计算蒸散发能力

    KIt = (1.0 - (1.0 - user_parms['KI'] - user_parms['KG']) ** (1.0 / (24.0 * 1.0)))/(1.0 + user_parms['KG']/user_parms['KI']) # 划分小段后的参数值
    KGt = user_parms['KG'] * KIt / user_parms['KI']
    KI = KIt; KG = KGt

    return WM, WMM, MS, EP, KIt, KGt, KI, KG

def soil_variable(user_parms, WLM, WUM, WDM, WU, WL, WD, S0, SM):
    """
    Analyze the inputs for soil water to make sure that the input variables are of consistence.
    """
    # user_parms = json.load(open(f'../data/input_parms.json', 'rb'))
    WM, WMM, MS, EP, KIt, KGt, KI, KG  = calculate_capacities(user_parms)

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

def initial_state():
    """
    Assign the initial state of variables from inputs.
    """
    pass

def runoff_single_time():
    pass

def separate_source():
    pass

def soil_evpa():
    pass

def runoff_production():
    pass

