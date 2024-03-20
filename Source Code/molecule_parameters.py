import sympy as sy
import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j

c = 29979.2458 #cm/us

params_general = {
'mu_B': 1.399624494, #MHz/Gauss
'g_S': 2.0023,
'g_L': 1,
'2_e0c': 2*3.76730314*10**(10),# (V/cm)^2/(W/um^2)
'mu_N': 7.62259323*10**(-4), #MHz/Gauss
}


### YbOH Parameters ###

# params_174X000 = {
# 'Be': 7348.40053,
# 'Gamma_SR': -81.150,
# 'bF': 4.80,
# 'c': 2.46,
# 'b': (4.80-2.46/3),
# 'D': 0.006084,
# 'Gamma_D': 0.00476,
# 'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
# }

# params_174X000 = { #YbF
# 'Be': 7233.8271,
# 'Gamma_SR': -13.41679,
# 'bF': 170.26374,
# 'c': 85.4028,
# 'muE': 3.91*0.503412 #Debye in MHz/(V/cm)
# }

params_174X000 = { #RaF
'Be': 5755.56,
'Gamma_SR': 175.38,
'bF': 96.3,
'c': 19,
'muE': 3.91*0.503412, #Debye in MHz/(V/cm)
'D': 1.4e-7*c
}

params_174X010 = {
'Be': 7328.48,
'Gamma_SR': -87.69,
'Gamma_Prime': 0,
'bF': 4.80,
'c': 2.46,
'b': (4.80-2.46/3),
'q_lD': 12.20, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
'p_lD': 12.09,
'muE': 2.15*0.503412,
#Origin 319.909053
}

# #YbOH With gamma prime
# params_174X010 = {
# 'Be': 7328.644,
# 'Gamma_SR': -88.660,
# 'Gamma_Prime': 17.38, #optical 15.61,detuned raman 17.38
# 'bF': 4.08, #4.07 splitting fit, 4.08 combined fit
# 'c': 3.41, #3.49 splitting fit, 3.41 combined fit
# 'q_lD': -12.03, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
# 'p_lD': -11.30, #optical -10.73, detuned raman -11.30
# 'muE': 2.15*0.503412,
# 'Origin': 319.90901,
# 'g_l': 0.0,
# 'g_S':2.07,
# }

#BaOH
# params_174X010 = {
# 'Be': 6485.2640,
# 'Gamma_SR': 68.65,
# 'Gamma_Prime': 0, #optical 15.61,detuned raman 17.38
# 'bF': 4.08, #YbOH Value
# 'c': 3.41, #YbOH Value
# 'q_lD': -9.4932, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
# 'p_lD': 2.33, #optical -10.73, detuned raman -11.30
# 'muE': 1.43*0.503412,
# 'Origin': 	341.6,
# 'g_l': 0.0,
# 'g_S':2.0023,
# }

#226RaF
params_174A000 = {
'Be': 5726.48,
'ASO': 0*2067.6*c, #Fixed from 1350 cm^-1
'a': -1.6/2,         # YbF value, using h1/2 = a- (bf+2c/3) = A||/2
'bF': 0,     #
'c': 0,     #
'd': -4.6, #YbF value
'p+2q': -0.41071*c,
'q':0,
'D': 1.4e-7*c,
'p2q_D': 1.9e-7*c,
'g_lp': -0.724,#-0.865,
'muE': 0.43*0.503412,
'g_S': 2.0023,
'Origin': 13284.427+0*2067.6/2, #Pi1/2 origin + ASO(=1350 cm-1)
}

# 174YbOH
# params_174A000 = {
# 'Be': 7586.3+2*0.006952*0,
# 'ASO': 4.04719818*10**7, #Fixed from 1350 cm^-1
# 'a': 0.01,         # extrapolated from YbF
# 'bF': 0.06985,     # extrapolated from YbF
# 'c': 0.1799,     # extrapolated from YbF
# 'p+2q': -13133-0.1139*0,
# 'q':0,
# 'D': 0.006952,
# 'p2q_D': 0.1139,
# 'g_lp': -0.724,#-0.865,
# 'muE': 0.43*0.503412,
# 'g_S': 1.86,
# }

params_173X000 = { # all units MHz except for muE
'Be': 7351.24,
'Gamma_SR': -81.06,
'bFYb': -1883.21,
'cYb': -81.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': -3318.70,
'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
}

params_171X000 = {
'Be': 7359.81,
'Gamma_SR': -80.85,
'bFYb': 6823.58,
'cYb': 233.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': 0,
'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
}

params_171X010 = {
'Be': 7359.81,
'Gamma_SR': -93.593,
'Gamma_Prime':0,
'bFYb': 6823.58,
'cYb': 233.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': 0,
'q_lD': -14.752,
'muE': 1.09*0.503412 #Debye in MHz/(V/cm)
}

params_173X010 = { # all units MHz except for muE
'Be': 7351.2,
'Gamma_SR': -93.593,
'Gamma_Prime':0,
'bFYb': -1883.21,
'cYb': -81.84,
'bFH': 4.80,
'cH': 2.46,
'e2Qq0': -3318.7,
'q_lD': -14.752,
'muE': 2.12*0.503412 #Debye in MHz/(V/cm)
}

params_173A000 = {
'Be': 7590.30,
'ASO': 0, #Actually 4.047*10**7,
'h1/2Yb': -126.51,
'dYb': -261.72,
'bFH': 0.07,     # extrapolated from YbF
'cH': -0.18,     # extrapolated from YbF
'e2Qq0': -1924.66,
'p+2q': -13141.63,
'muE': 0.43*0.503412
}

params_171A000 = {
'Be': 7597.79,
'ASO': 0, #Actually 4.047*10**7,
'h1/2Yb': 443.69,
'dYb': 959.04,
'bFH': 48.8, #0.07,     # extrapolated from YbF
'cH': 24.46, #-0.18,     # extrapolated from YbF
'e2Qq0': 0,
'p+2q': -13150.91,
'muE': 0.43*0.503412
}


### CaOH Parameters ###
#X(000) Taken from Louis Baum thesis and Steimle papers from 90s
#Vibrational states from Fletcher et al

params_40X000 = {
'Be': 10023.0841,
'D': 1.154*10**-2,
'Gamma_SR': 34.7593,
'bF': 2.602,
'c': 2.053,
'b': (2.602-2.053/3),
'muE': 1.465*0.503412, #Debye in MHz/(V/cm)
# 'g_N': 5.253736,
}

#Whenever possible, constants are taken from Fletcher et all, Milimeter Wave Hydroxide paper
# params_40X010 = {
# 'Be': 9996.7518,
# 'D': 0.0117696,
# 'Gamma_SR': 35.051,
# 'Gamma_Prime': 0,
# 'bF': 2.244, #2.602, #2.29 fit?
# 'c': 2.607, #2.053, #2.52 fit?
# # 'b': (2.29-2.52/3),
# 'p_lD': -0.05,
# 'q_lD': -21.6492,
# 'q_lD_D': 6.4*10**-5,
# 'muE': 1.465*0.503412,
# 'azz': 3.33441*10**(-8)*5.525**2/25.875, #Debye^2/MHz in units of MHz/(V/cm)^2. Using Lan values
# 'axxyy': 9.46049*10**(-8)*6.165**2/29.34
# }

#Coxon parameters for comparison:
params_40X010 = {
'Be': 9996.82,
'D': 0.008823,
'Gamma_SR': 35.5,
'Gamma_Prime': 0,
'bF': 2.45,#2.247, #2.293,#2.2445, #2.602
'c': 2.6,#2.601,#2.522,##2.6074, #2.053
# 'b': (2.29-2.52/3),
'p_lD': -0.00,
'q_lD': -21.53,
# 'q_lD_D': 6.4*10**-5,
'muE': 1.465*0.503412,
# 'g_N': 5.253736,
'azz': 3.5555*10**(-8), #From Lan for X(000)
'axxyy': 1.1718*10**(-7) #From Lan for X(000)
# 'azz': 3.33441*10**(-8)*5.525**2/25.875, #Calc from Lan's static #Debye^2/MHz in units of MHz/(V/cm)^2. Using Lan values
# 'axxyy': 9.46049*10**(-8)*6.165**2/29.34 #Calc from Lan's static
}


params_40A000 = {
'Be': 10229.52,
'ASO': 2.00316*10**6,
'a': 0,         # extrapolated from YbF
'bF': 0.07,     # extrapolated from YbF
'c': -0.18,     # extrapolated from YbF
'p+2q': -1305 ,
'q': -9.764,
'g_lp': -0.865, #Unknown
'muE': 0.836*0.503412
}

params_40B000 = {
'Be': 10175.2,
'Gamma_SR': -1307.54,
'bF': 2.602*0.04, #Hyperfine extrapolated from CaF A state to X state ratio
'c': 2.053*0.04,
'b': (2.602-2.053/3)*0.04,
'muE': 0.744*0.503412 #Debye in MHz/(V/cm)
}

YbOH_params = {
'174X000':{**params_general,**params_174X000},
'174X010':{**params_general,**params_174X010},
'173X000':{**params_general,**params_173X000},
'173X010':{**params_general,**params_173X010},
'174A000':{**params_general,**params_174A000},
'173A000':{**params_general,**params_173A000},
'171A000':{**params_general,**params_171A000},
'171X000':{**params_general,**params_171X000},
'171X010':{**params_general,**params_171X010},
}

CaOH_params = {
'40X000':{**params_general,**params_40X000},
'40X010':{**params_general,**params_40X010},
'40A000':{**params_general,**params_40A000},
'40B000':{**params_general,**params_40B000},
}

all_params = {
'YbOH': YbOH_params,
'CaOH': CaOH_params}
