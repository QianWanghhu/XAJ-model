# import packages
import numpy as np


def soil_variable(WLM, WUM, WDM, WU, WL, WD, S0, SM):
    """
    Analyze the inputs to make sure that the input variables are of consistence.
    """
    if WU>WUM:
        WL=WL+WU-WUM
        WU=WUM

    if WU > WLM:
        WD=WD+WL-WLM
        WL=WLM
    
    if WD>WDM: WD=WDM

    if S0 > SM: S0 = SM

    return WL, WU, WD, S0
# END soil_variable()

def runoff_single_time():
    pass

def separate_source():
    pass

def soil_evpa():
    pass

def runoff_production():
    pass

