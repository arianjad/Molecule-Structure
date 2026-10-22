"""BaF X-A simulated spectrum vs the 08/27/2025 fluorescence scan.

Simulation is generated from molecular structure via gen_spectra
(no pasted arrays). Freq offset is 348.66 THz.
"""
from config_path import add_to_sys_path
add_to_sys_path()                       # walk up to "Source Code"

import numpy as np
import matplotlib.pyplot as plt
from Energy_Levels import MoleculeLevels
import gen_spectra as xs

# --- measured data (unchanged: [freq(MHz)-offset, counts, stdev]) ---
# Data format is [freq (MHz) - offset, fluorescence counts (arb), and count stdev]
# Freq offset is 348.66 THz
# Data is from 08272025
fluor_data = [np.array([1205., 1210., 1215., 1220., 1225., 1230., 1235., 1240., 1245.,
        1250., 1255., 1260., 1265., 1270., 1275., 1280., 1285., 1290.,
        1295., 1300., 1305., 1310., 1105., 1110., 1115., 1120., 1125.,
        1130., 1135., 1140., 1145., 1150., 1155., 1160., 1165., 1170.,
        1175., 1180., 1185., 1190., 1195., 1200.]),
 np.array([  954.43592634,  1087.61957589,   953.54723214,   829.73112723,
         1185.73243304,  1662.5906529 ,  1761.31720982,  2429.26184152,
         3533.5843378 ,  3789.80350446,  5360.14004464,  7503.88844866,
         7876.46177455,  8462.56964286,  7924.8356808 ,  7821.04029018,
         9434.2538058 ,  7251.5007952 ,  3839.25515625,  4594.57449777,
         2883.02945312,  2267.90547991,  3057.64712054,  4432.69589286,
         6285.13082589,  6199.64004464, 11302.64309152, 16451.95758929,
        14443.29041295, 17100.06632813, 19282.46429688, 16190.796875  ,
        13476.84639509,  7493.23444196,  7256.00859375,  3795.53261161,
         3969.40377232,  2422.1769308 ,  2310.4871317 ,  1852.63741071,
         1451.77878348,  1077.76672991]),
 np.array([ 235.00589891,  128.99784432,   81.46605303,  243.85854923,
         196.81471109,  473.50919145,  429.74031751,  474.4173157 ,
         498.28412148,  879.07844073, 1567.73056767, 1999.87169295,
        1864.54178033, 4916.22051105, 1947.09522348, 1789.61127126,
        2178.55245282, 2076.69114048, 1049.67959235,  464.17120557,
         864.92703449,  932.33154183,  941.21718282,  946.48199856,
        1895.57958143,  971.23372669, 3450.89241376, 5705.77957989,
        1244.18776877, 4230.17815931, 4184.97357791, 3499.85325979,
        2948.30779245, 2537.25328151, 3409.19722747, 1110.35611458,
        1115.43950672,  349.21045873,  800.05779797,  349.06214641,
         407.58383444,  144.66573808])]

offset = 348660000                      # 348.66 THz, the x-axis reference
tweak = -19                             # laser-frequency calibration shift
yscale = 1.3

def build(elec, N_list):
    return MoleculeLevels.initialize_state(
        molecule_name='BaF', elec_state=elec, vib_state=0,
        N_list=np.array(N_list), fermion_or_boson='boson',
        M_sublevels='none', I_nuclei=[0, 1/2], isotope=138,
        round=8, params=None, P_values=[1/2])

g = build('X', [1])                      # the N=1 feature in the 08/27 scan
e = build('A', [1])                      # A 2Pi1/2 J=1/2
g.eigensystem(0, 0); e.eigensystem(0, 0)

gidx = g.select_q({'N': 1})
eidx = e.select_q({'J': 0.5}, parity='+')
lines = xs.line_list(g, e, gidx, eidx, origin=e.parameters['Origin'])
lines = lines.assign(freq=lines['freq'] - offset)        # to the plot axis

fig, ax = plt.subplots()
xs.plot_spectrum(
    lines=lines, ax=ax, sticks=True,
    broaden_kw=dict(shape='gaussian', fwhm=40.0, cluster=True),
    normalize=True,
    experimental=dict(freq=fluor_data[0], signal=fluor_data[1],
                      err=fluor_data[2], tweak=tweak, yscale=yscale),
    label='Simulation')
ax.set_xlabel('Laser Frequency (MHz) - 348.66 THz')
ax.set_ylabel('Fluorescence (arb)')
fig.savefig('BaF_X_A_N1_sim_vs_data.pdf')
print("wrote BaF_X_A_N1_sim_vs_data.pdf;", len(lines), "lines")
