import sympy as sy
import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j
from formalism import convert_formalism


def bF_2_b(bF_value,c_value):
    b = bF_value - c_value/3
    return b

def b_2_bF(b_value, c_value):
    bF = b_value + c_value/3
    return bF

def abinitio_2_effective_hyperfine(A_par,A_perp,state='Sigma'):
    if 'Sigma' in state:
        b = A_perp
        c = A_par - A_perp
        bF = b_2_bF(b,c)
    return bF,c



params_general = {
'mu_B': 1.399624494, #MHz/Gauss
'g_S': 2.0023,
'g_L': 1,
'2_e0c': 2*3.76730314*10**(10),# (V/cm)^2/(W/um^2)
'mu_N': 7.62259323*10**(-4), #MHz/Gauss,
'c': 29979.2458 #cm/us
}

c=params_general['c']

# Dictionary to hold all molecule parameters

species_names = ['YbOH','CaOH','BaOH','BaF','RaF','RaOH','YbF','DyO']

# Build nested structure with fermion/boson branches
molecules = {molecule_name: {'boson': {}, 'fermion': {}} for molecule_name in species_names}


molecules['RaF']['boson']['X0'] = {
    'Be': 5755.56,
    'Gamma_SR': 175.38,
    'bF': 96.3,
    'c': 19,
    'muE': 3.91*0.503412, #Debye in MHz/(V/cm)
    'D': 1.4e-7*c
    }
molecules['RaF']['boson']['A0'] = {
    'formalism': 'N2',   # raw paper (N²) values; converter → engine R² (spec)
    'Lambda': 1,         # A²Π, declared case-(a) basis Λ
    'Be': 0.191015*c,    # B′₀(N²) = 0.191015(5)[15] cm⁻¹. Udrescu et al.,
                         # Nat. Phys. 2024, DOI 10.1038/s41567-023-02296-w,
                         # Table I, ²²⁶Ra¹⁹F A²Π₁/₂ v=0. Was 5743.96 MHz —
                         # wrong source, ~17.5 MHz high vs this paper (Task 6).
    'ASO': 1350*c,       # Fixed from 1350 cm^-1
    'h1/2': 0,           # Calc relayed from Silviu, h1/2 = a-(bf+2c/3)=A||/2
    'a':19/2,
    'bF': 0,
    'c': 0,
    'd': -9,             # Silviu; PGopher convention, B+C needs minus sign
    'p+2q': -0.41071*c,  # X(N²) paper value
    'q':0,
    'D': 1.4e-7*c,
    'p2q_D': 1.9e-7*c,   # X_D for p+2q
    'g_lp': -0.41071/(2*0.191015),#-0.865,  # 0.191015 = B′₀(N²) cm⁻¹ (Udrescu 2024)
    'g_l': -6000/(2*0.191015*c),
    'muE': 1*0.503412,
    'g_S': 2.0023,
    'Origin': 13284.427, # PHYSICAL band origin (N²); converter adds R² G-row.
                         # NO hand +B/c term (that was the +11.6 MHz bug).
    }


molecules['BaF']['boson']['X0'] = {
    'formalism': 'N2',      # convention declaration (Sigma: convention-independent; converter no-ops at Lambda=0)
    'Lambda': 0,            # X^2Sigma+: |Lambda|=0
    'Be': 6473.9588,        # Ryzlewicz1980 [22]; 138Ba19F; 0.21594802 cm^-1
    'Gamma_SR': 80.984,     # Ryzlewicz1980 [22]; 0.0027013 cm^-1
    'bF': b_2_bF(63.509, 8.224),  # =66.2503; Frosch-Foley b,c Ernst1986 [23]; bF=b+c/3 (B&C Eq.8.511 p.605/p.845; Steimle PRA84 012508 TblIII bF(F)=66.25 direct)
    'c': 8.224,             # Ernst1986 [23]; 0.0002743 cm^-1
    'D': 0.0055250,         # Effantin1990 [18]; 1.84294e-7 cm^-1
    'muE': 3.170*0.503412,  # 3.170(3) D, Steimle PRA84 012508 Table V; Debye->MHz/(V/cm)
    }
molecules['BaF']['boson']['A0'] = {
    'Be': 6347.847,         # Steimle2011 [24], 0.2117414 cm^-1, held fixed in fit (paper N^2 value). Converter applies the N^2->R^2 B-row (-2*Lam^2*D = -0.012 MHz, immaterial) since this dict is now formalism:'N2'.
    'ASO': 18955512.5,      # This work arXiv:2511.06986; 632.287838 cm^-1
    'h1/2': 0,
    'a': 26.55 - 0.5*(-5.54),  # =29.32; arXiv:2511.06986 Eq.(4)/B&C p.845 h1/2=a-(b+c)/2; Denis2022 [39]. Gated on Table II 21.8 MHz (plan Task 2)
    'bF': 0,
    'c': 0,
    'd': -3.58,             # Denis2022 [39]; pgopher->B+C needs minus (RaF A0 precedent)
    'p+2q': -7713.96 + 2*(-2.52),  # =-7719.00 (spec/§4 wrote -7718.99: hand-add slip, -7713.96-5.04=-7719.00; 0.01 MHz, immaterial); Effantin1990 [18] combined p+2q (op-grounded hamiltonian_builders.py:147; Steimle (p+2q)=-7721 cross-check)
    'q': -2.52,             # Effantin1990 [18]; -8.40e-5 cm^-1
    'p2q_D': -0.00699,      # Effantin1990 [18] p_D; -2.332e-7 cm^-1
    'D': 0.006007,          # Effantin1990 [18]; 2.0036e-7 cm^-1
    'A_D': 0.93,            # Steimle2011 [24]; consumed iff A2Pi builder reads it (plan Task 2 confirms)
    'muE': 1.50*0.503412,   # 1.50(2) D A2Pi1/2, Steimle PRA84 012508 Table V; Debye->MHz/(V/cm)
    'g_S': 2.0023,
    'g_lp': -0.536,         # MEASURED parity-dependent anisotropic g-factor g_l'(A2Pi,v=0), 138BaF: -0.536(23). Steimle, Frey, Le, DeMille, Rahmlow, Linton, Phys. Rev. A 84 012508 (2011) [24], Table V (read citation-grade; p.22: Curl-analog estimate -0.2575/(2*0.21954)=-0.587 "in very good agreement with the experimental value (=-0.536(23))"). Same p/(2B) convention as the validated RaF A0 g_lp (below) & H_even_A ZeemanParityZ. Dimensionless; Zeeman-only -> ZERO effect on B=0 spectra. (g_S fixed 2.002 and perp g_l constrained 0 in Steimle's A2Pi fit, consistent with this dict.)
    'Origin': 11946.109609,  # PHYSICAL N^2 band origin T0,0 [cm^-1] (arXiv:2511.06986 Tbl III, "This work"; pgopher-DEFAULT N^2 fit). The +Be_A/c R^2 G-row is now applied by the converter (formalism:'N2' below), NOT by hand. Verified deltas vs the old hand path: Origin -0.006, Be -0.012, p+2q -0.007, ASO +0.93 MHz (on ~1.9e7) -- all << 1-2 MHz validation tols. Spec docs/superpowers/specs/2026-05-16-xa-spectra-toolkit-design.md §2; thesis App A.3.1.
    'formalism': 'N2',        # paper fit is pgopher-default N^2; converter -> engine R^2 (adds Be*Lam^2/c G-row + Table 7.2 terms)
    'Lambda': 1,              # A^2Pi: |Lambda|=1 (REQUIRED whenever 'formalism' set)
    }

bF_225, c_225 = abinitio_2_effective_hyperfine(A_par=-0.5692*c, A_perp=-0.5445*c, state='Sigma')

molecules['RaF']['fermion']['X0'] = {
    'Be': 5755.56,
    'Gamma_SR': 175.38*1.004, # scaled from 226
    'bFYb': bF_225,
    'cYb': c_225,
    'bFH': 96.3,
    'cH': 19,
    'e2Qq0': 0,
    'muE': 3.91*0.503412, #Debye in MHz/(V/cm)
    'D': 1.4e-7*c
    }

molecules['RaF']['fermion']['A0'] = {
    'formalism': 'N2',
    'Lambda': 1,
    'Be': 5729.03,       # B(N²) paper value
    'ASO': 1350*c,
    'h1/2Yb': -2852/2,   # From Skripnikov theory
    'dYb': -1*-0.076*c,  # Wilkins expt; pgopher convention → add - sign
    'h1/2H': 19/2,
    'aH':19/2,
    'bFH': 0,
    'cH': 0,
    'dH': -9,
    'p+2q': -0.4109*c,   # X(N²) paper value
    'q':0,
    'e2Qq0': 0,
    'D': 1.4e-7*c,
    'p2q_D': 1.9e-7*c,
    'g_lp': -0.4109/(2*0.191100),#-0.865,
    'muE': 1.0*0.503412,
    'g_S': 2.0023,
    'Origin': 13284.427, # PHYSICAL band origin (N²); converter adds R² G-row
    }

molecules['YbOH']['fermion']['171A000'] = {
    'Be': 7597.79,
    'ASO': 0, #Actually 4.047*10**7,
    'h1/2Yb': 443.69,
    'dYb': 959.04,
    'bFH': 48.8, #0.07,     # extrapolated from YbF
    'cH': 24.46, #-0.18,     # extrapolated from YbF
    'e2Qq0': 0,
    'p+2q': -13150.91,
    'muE': 0.43*0.503412,
    'g_S': 2.0023
    }

molecules['RaOH']['fermion']['X010'] = {
    'Be': 5814.3, # X000 value
    'Gamma_SR': 151, # X000 value
    'Gamma_Prime':0,
    'bFYb': bF_225,
    'cYb': c_225,
    'bFH': 4.08, # YbOH Value
    'cH': 3.41, # YbOH Value
    'e2Qq0': 0,
    'q_lD': -15,
    'muE': 2.15*0.503412 # Debye in MHz/(V/cm)
    }



molecules['YbOH']['boson']['X000'] = {
    'Be': 7348.40053,
    'Gamma_SR': -81.150,
    'bF': 4.80,
    'c': 2.46,
    'b': bF_2_b(4.80,2.46),
    'D': 0.006084,
    'Gamma_D': 0.00476,
    'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
    }
molecules['YbOH']['boson']['X010'] = {
    'Be': 7328.644,
    'Gamma_SR': -88.660,
    'Gamma_Prime': 17.38, #optical 15.61,detuned raman 17.38
    'bF': 4.08, #4.07 splitting fit, 4.08 combined fit
    'c': 3.41, #3.49 splitting fit, 3.41 combined fit
    'b': bF_2_b(4.08,3.41),
    'q_lD': -12.03, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
    'p_lD': -11.30, #optical -10.73, detuned raman -11.30
    'D': 0.006952,
    'Gamma_D': 0.00476,
    'muE': 2.15*0.503412,
    'Origin': 319.90901,
    'g_l': 0.0,
    'g_S':2.07,
    }
molecules['YbOH']['boson']['A000'] = {
    'Be': 7590.30,
    'ASO': 4.04719818*10**7, #Actually 4.047*10**7,
    'h1/2': -126.51,
    'd': -261.72,
    'a': 0.01,         # extrapolated from YbF
    'bF': 0.06985,     # extrapolated from YbF
    'c': 0.1799,     # extrapolated from YbF
    'b': bF_2_b(0.06985,0.1799),
    'p+2q': -13133,
    'q': -9.4932,
    'D': 0.006952,
    'p2q_D': 0.1139,
    'g_lp': -0.724,#-0.865,
    'muE': 0.43*0.503412,
    'g_S': 1.86,
    }


molecules['YbOH']['fermion']['173X000'] = {
    'Be': 7351.24,
    'Gamma_SR': -81.06,
    'bFYb': -1883.21,
    'cYb': -81.84,
    'bFH': 4.80,
    'cH': 2.46,
    'e2Qq0': -3318.70,
    'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
    }
molecules['YbOH']['fermion']['173X010'] = {
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
molecules['YbOH']['fermion']['173A000'] = {
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

molecules['YbOH']['fermion']['171X000'] = {
    'Be': 7359.81,
    'Gamma_SR': -80.85,
    'bFYb': 6823.58,
    'cYb': 233.84,
    'bFH': 4.80,
    'cH': 2.46,
    'e2Qq0': 0,
    'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
    }
molecules['YbOH']['fermion']['171X010'] = {
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
molecules['YbOH']['fermion']['171A000'] = {
    'Be': 7597.79,
    'ASO': 0, #Actually 4.047*10**7,
    'h1/2Yb': 443.69,
    'dYb': 959.04,
    'bFH': 48.8, #0.07,     # extrapolated from YbF
    'cH': 24.46, #-0.18,     # extrapolated from YbF
    'e2Qq0': 0,
    'p+2q': -13150.91,
    'muE': 0.43*0.503412,
    'g_S': 2.0023
    }


molecules['YbF']['boson']['X0']={
    'Be': 7233.8271,
    'Gamma_SR': -13.41679,
    'bF': 170.26374,
    'c': 85.4028,
    'muE': 3.91*0.503412 #Debye in MHz/(V/cm)
    }


molecules['CaOH']['boson']['40X000'] = {
    'Be': 10023.0841,
    'D': 1.154*10**-2,
    'Gamma_SR': 34.7593,
    'bF': 2.602,
    'c': 2.053,
    'b': bF_2_b(2.602,2.053),
    'muE': 1.465*0.503412, #Debye in MHz/(V/cm)
    # 'g_N': 5.253736,
    }
molecules['CaOH']['boson']['40X010'] = {
    'Be': 9996.82,
    'D': 0.008823,
    'Gamma_SR': 35.5,
    'Gamma_Prime': 0,
    'bF': 2.45,#2.247, #2.293,#2.2445, #2.602
    'c': 2.6,#2.601,#2.522,##2.6074, #2.053
    'b': bF_2_b(2.45,2.6),
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
molecules['CaOH']['boson']['40A000'] = {
    'Be': 10229.52,
    'ASO': 2.00316*10**6,
    'a': 0,         # extrapolated from YbF
    'bF': 0.07,     # extrapolated from YbF
    'c': -0.18,     # extrapolated from YbF
    'b': bF_2_b(0.07,-0.18),
    'p+2q': -1305 ,
    'q': -9.764,
    'g_lp': -0.865, #Unknown
    'muE': 0.836*0.503412
    }
molecules['CaOH']['boson']['40B000'] = {
    'Be': 10175.2,
    'Gamma_SR': -1307.54,
    'bF': 2.602*0.04, #Hyperfine extrapolated from CaF A state to X state ratio
    'c': 2.053*0.04,
    'b': bF_2_b(2.602*0.04,2.053*0.04),
    'muE': 0.744*0.503412 #Debye in MHz/(V/cm)
    }


molecules['BaOH']['boson']['X000'] = {
    'Be': 6485.2640,
    'Gamma_SR': 68.65,
    'bF': 4.08, #YbOH Value
    'c': 3.41, #YbOH Value
    'b': bF_2_b(4.08,3.41),
    'D': 0.006952,
    'muE': 1.43*0.503412,
    'Origin': 	341.6,
    'g_l': 0.0,
    'g_S':2.0023,
    }
molecules['BaOH']['boson']['X010'] = {
    'Be': 6485.2640,
    'Gamma_SR': 68.65,
    'bF': 4.09, #YbOH Value
    'c': 3.41, #YbOH Value
    'b': bF_2_b(4.08,3.41),
    'D': 0.006952,
    'muE': 1.43*0.503412,
    'Origin': 	341.6,
    'g_l': 0.0,
    'g_S':2.0023,
    }


molecules['DyO']['boson']['X000'] = {
    'Be': 0.33*0,
    'ASO': 0, #Actually 4.047*10**7,
    'h1/2Yb': 1000,
    'dYb': -0,
    'bFH': 0,     # extrapolated from YbF
    'cH': -0.0,     # extrapolated from YbF
    'e2Qq0': 0,
    'p+2q': 0,
    'muE': 2*0.503412
    }
molecules['DyO']['boson']['A000'] = {
    'Be': 0.33*0,
    'ASO': 0, #Actually 4.047*10**7,
    'h1/2Yb': 1000,
    'dYb': -0,
    'bFH': 0,     # extrapolated from YbF
    'cH': -0.0,     # extrapolated from YbF
    'e2Qq0': 0,
    'p+2q': 0,
    'muE': 2*0.503412
    }




def get_molecule_params(molecule_name, elec_state,vib_state, fermion_or_boson=None,
                        include_general=True, overrides=None):

    vibronic_label = f'{elec_state}{vib_state}'

    if molecule_name not in molecules:
        raise KeyError(f"Unknown molecule '{molecule_name}'. Available: {list(molecules)}")

    molecule = molecules[molecule_name]

    if fermion_or_boson is None:
        fermion_or_boson='boson'
    
    vibronic_states = molecule[fermion_or_boson]

    if vibronic_label not in vibronic_states:
        raise KeyError(f"Unknown state '{vibronic_label}' for '{molecule_name}' ({fermion_or_boson}). "
                       f"Available: {list(vibronic_states)}")

    state_params = vibronic_states[vibronic_label]

    if include_general:
        merged = dict(params_general)
        merged.update(state_params)
    else:
        merged = state_params

    if overrides:
        merged.update(overrides)
    
    # Engine wants R²: convert only an explicitly N²-tagged set; R²/untagged
    # passes through unchanged (copied). convert_formalism flips the tag.
    if merged.get('formalism') == 'N2':
        return convert_formalism(merged, c_cm=c)
    return dict(merged)


### Legacy Code ###


# ### YbOH Parameters ###


# params_174X000 = { #YbF
# 'Be': 7233.8271,
# 'Gamma_SR': -13.41679,
# 'bF': 170.26374,
# 'c': 85.4028,
# 'muE': 3.91*0.503412 #Debye in MHz/(V/cm)
# }

# # params_174X000 = { #RaF
# # 'Be': 5755.56,
# # 'Gamma_SR': 175.38,
# # 'bF': 96.3,
# # 'c': 19,
# # 'muE': 3.91*0.503412, #Debye in MHz/(V/cm)
# # 'D': 1.4e-7*c
# # }

# # params_174X010 = {
# # 'Be': 7328.48,
# # 'Gamma_SR': -87.69,
# # 'Gamma_Prime': 0,
# # 'bF': 4.80,
# # 'c': 2.46,
# # 'b': (4.80-2.46/3),
# # 'q_lD': 12.20, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
# # 'p_lD': 12.09,
# # 'muE': 2.15*0.503412,
# # #Origin 319.909053
# # }

# # #YbOH With gamma prime
# # params_174X010 = {
# # 'Be': 7328.644,
# # 'Gamma_SR': -88.660,
# # 'Gamma_Prime': 17.38, #optical 15.61,detuned raman 17.38
# # 'bF': 4.08, #4.07 splitting fit, 4.08 combined fit
# # 'c': 3.41, #3.49 splitting fit, 3.41 combined fit
# # 'q_lD': -12.03, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
# # 'p_lD': -11.30, #optical -10.73, detuned raman -11.30
# # 'muE': 2.15*0.503412,
# # 'Origin': 319.90901,
# # 'g_l': 0.0,
# # 'g_S':2.07,
# # }

# #BaOH
# # params_174X010 = {
# # 'Be': 6485.2640,
# # 'Gamma_SR': 68.65,
# # 'Gamma_Prime': 0, #optical 15.61,detuned raman 17.38
# # 'bF': 4.08, #YbOH Value
# # 'c': 3.41, #YbOH Value
# # 'q_lD': -9.4932, #Should be minus if neg parity lower and parity is (-1)^(J-l-S), but for modeling I treat this as positive for now....
# # 'p_lD': 2.33, #optical -10.73, detuned raman -11.30
# # 'muE': 1.43*0.503412,
# # 'Origin': 	341.6,
# # 'g_l': 0.0,
# # 'g_S':2.0023,
# # }

# #226RaF
# # params_174A000 = {
# # 'Be': 5726.48 - 2*1.4e-7*c, #R^2 form
# # 'ASO': 2067.6*c, #Fixed from 1350 cm^-1
# # 'h1/2': 0,         # Calc relayed from Silviu, using h1/2 = a- (bf+2c/3) = A||/2
# # 'a':19/2,
# # 'bF': 0,     #
# # 'c': 0,     #
# # 'd': -9, #Calc relayed from Silviu. They use PGopher convention, we use B+C, so need minus sign
# # 'p+2q': -0.41071*c+1.9e-7*c, #R^2 form
# # 'q':0,
# # 'D': 1.4e-7*c,
# # 'p2q_D': 1.9e-7*c,
# # 'g_lp': -0.724,#-0.865,
# # 'muE': 0.43*0.503412,
# # 'g_S': 2.0023,
# # 'Origin': 13284.427+0*2067.6/2 + 5726.48/c, #Pi1/2 origin + ASO(=1350 cm-1), also with B offset due to code being R^2
# # }

# # 174YbOH
# # params_174A000 = {
# # 'Be': 7586.3+2*0.006952*0,
# # 'ASO': 4.04719818*10**7, #Fixed from 1350 cm^-1
# # 'a': 0.01,         # extrapolated from YbF
# # 'bF': 0.06985,     # extrapolated from YbF
# # 'c': 0.1799,     # extrapolated from YbF
# # 'p+2q': -13133-0.1139*0,
# # 'q':0,
# # 'D': 0.006952,
# # 'p2q_D': 0.1139,
# # 'g_lp': -0.724,#-0.865,
# # 'muE': 0.43*0.503412,
# # 'g_S': 1.86,
# # }

# params_173X000 = { # all units MHz except for muE
# 'Be': 7351.24,
# 'Gamma_SR': -81.06,
# 'bFYb': -1883.21,
# 'cYb': -81.84,
# 'bFH': 4.80,
# 'cH': 2.46,
# 'e2Qq0': -3318.70,
# 'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
# }

# params_171X000 = {
# 'Be': 7359.81,
# 'Gamma_SR': -80.85,
# 'bFYb': 6823.58,
# 'cYb': 233.84,
# 'bFH': 4.80,
# 'cH': 2.46,
# 'e2Qq0': 0,
# 'muE': 1.9*0.503412 #Debye in MHz/(V/cm)
# }

# params_171X010 = {
# 'Be': 7359.81,
# 'Gamma_SR': -93.593,
# 'Gamma_Prime':0,
# 'bFYb': 6823.58,
# 'cYb': 233.84,
# 'bFH': 4.80,
# 'cH': 2.46,
# 'e2Qq0': 0,
# 'q_lD': -14.752,
# 'muE': 1.09*0.503412 #Debye in MHz/(V/cm)
# }

# params_173X010 = { # all units MHz except for muE
# 'Be': 7351.2,
# 'Gamma_SR': -93.593,
# 'Gamma_Prime':0,
# 'bFYb': -1883.21,
# 'cYb': -81.84,
# 'bFH': 4.80,
# 'cH': 2.46,
# 'e2Qq0': -3318.7,
# 'q_lD': -14.752,
# 'muE': 2.12*0.503412 #Debye in MHz/(V/cm)
# }

# params_173A000 = { #DyO
# 'Be': 0.33*c,
# 'ASO': 0, #Actually 4.047*10**7,
# 'h1/2Yb': 1000,
# 'dYb': -0,
# 'bFH': 0,     # extrapolated from YbF
# 'cH': -0.0,     # extrapolated from YbF
# 'e2Qq0': 0,
# 'p+2q': 0,
# 'muE': 2*0.503412
# }

# # params_173A000 = {
# # 'Be': 7590.30,
# # 'ASO': 0, #Actually 4.047*10**7,
# # 'h1/2Yb': -126.51,
# # 'dYb': -261.72,
# # 'bFH': 0.07,     # extrapolated from YbF
# # 'cH': -0.18,     # extrapolated from YbF
# # 'e2Qq0': -1924.66,
# # 'p+2q': -13141.63,
# # 'muE': 0.43*0.503412
# # }

# params_171A000 = {
# 'Be': 7597.79,
# 'ASO': 0, #Actually 4.047*10**7,
# 'h1/2Yb': 443.69,
# 'dYb': 959.04,
# 'bFH': 48.8, #0.07,     # extrapolated from YbF
# 'cH': 24.46, #-0.18,     # extrapolated from YbF
# 'e2Qq0': 0,
# 'p+2q': -13150.91,
# 'muE': 0.43*0.503412
# }


# ### CaOH Parameters ###
# #X(000) Taken from Louis Baum thesis and Steimle papers from 90s
# #Vibrational states from Fletcher et al

# # params_40X000 = {
# # 'Be': 10023.0841,
# # 'D': 1.154*10**-2,
# # 'Gamma_SR': 34.7593,
# # 'bF': 2.602,
# # 'c': 2.053,
# # 'b': (2.602-2.053/3),
# # 'muE': 1.465*0.503412, #Debye in MHz/(V/cm)
# # # 'g_N': 5.253736,
# # }

# #Whenever possible, constants are taken from Fletcher et all, Milimeter Wave Hydroxide paper
# # params_40X010 = {
# # 'Be': 9996.7518,
# # 'D': 0.0117696,
# # 'Gamma_SR': 35.051,
# # 'Gamma_Prime': 0,
# # 'bF': 2.244, #2.602, #2.29 fit?
# # 'c': 2.607, #2.053, #2.52 fit?
# # # 'b': (2.29-2.52/3),
# # 'p_lD': -0.05,
# # 'q_lD': -21.6492,
# # 'q_lD_D': 6.4*10**-5,
# # 'muE': 1.465*0.503412,
# # 'azz': 3.33441*10**(-8)*5.525**2/25.875, #Debye^2/MHz in units of MHz/(V/cm)^2. Using Lan values
# # 'axxyy': 9.46049*10**(-8)*6.165**2/29.34
# # }

# # #Coxon parameters for comparison:
# # params_40X010 = {
# # 'Be': 9996.82,
# # 'D': 0.008823,
# # 'Gamma_SR': 35.5,
# # 'Gamma_Prime': 0,
# # 'bF': 2.45,#2.247, #2.293,#2.2445, #2.602
# # 'c': 2.6,#2.601,#2.522,##2.6074, #2.053
# # # 'b': (2.29-2.52/3),
# # 'p_lD': -0.00,
# # 'q_lD': -21.53,
# # # 'q_lD_D': 6.4*10**-5,
# # 'muE': 1.465*0.503412,
# # # 'g_N': 5.253736,
# # 'azz': 3.5555*10**(-8), #From Lan for X(000)
# # 'axxyy': 1.1718*10**(-7) #From Lan for X(000)
# # # 'azz': 3.33441*10**(-8)*5.525**2/25.875, #Calc from Lan's static #Debye^2/MHz in units of MHz/(V/cm)^2. Using Lan values
# # # 'axxyy': 9.46049*10**(-8)*6.165**2/29.34 #Calc from Lan's static
# # }


# # params_40A000 = {
# # 'Be': 10229.52,
# # 'ASO': 2.00316*10**6,
# # 'a': 0,         # extrapolated from YbF
# # 'bF': 0.07,     # extrapolated from YbF
# # 'c': -0.18,     # extrapolated from YbF
# # 'p+2q': -1305 ,
# # 'q': -9.764,
# # 'g_lp': -0.865, #Unknown
# # 'muE': 0.836*0.503412
# # }

# # params_40B000 = {
# # 'Be': 10175.2,
# # 'Gamma_SR': -1307.54,
# # 'bF': 2.602*0.04, #Hyperfine extrapolated from CaF A state to X state ratio
# # 'c': 2.053*0.04,
# # 'b': (2.602-2.053/3)*0.04,
# # 'muE': 0.744*0.503412 #Debye in MHz/(V/cm)
# # }

