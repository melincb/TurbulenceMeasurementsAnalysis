'''With this file we take PSD from shear1 and shear2 from priprava_podatkov_za_izracun_epsilon_shear.py, compute epsilon
from both PSDs using various methods and prepare this epsilon for plotting with
Risanje_epsilon_shear_na_vseh_tockah_ob_isti_ruti.py.'''

from fitting_procedures_for_Nasmyth_curve import fit_Nasmyth_MLE, fit_Nasmyth_Lueck
import numpy as np
import pickle
import scipy.stats

save_outputs_as_pickle = False      # If you want to save the output, set this to true
datumska_mapa = 'all_2009_04_20'    # date of measurements ('all_2008_08_21', 'all_2008_11_18', 'all_2009_04_20)

# Load the data
outputs_shear = pickle.load(open(f'shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p', 'rb'))

# THE BASIC MANUAL SETTINGS END HERE


prvi_slovar = {}    # the first dictionary
for ioutput in range(len(outputs_shear['tempcor'])): # Go through all data
    prg = outputs_shear['postaja_ruta_globina'][ioutput]    # station-route-depth
    postaja_ruta = prg[:12] # station-route
    globina = float(prg[13:])   # depth in dbar

    # Compute decimal logarithm of epsilon from shear sensor 1 using MLE with 12 degrees of freedom
    log_epsilon1_MLE = fit_Nasmyth_MLE(psd=outputs_shear['useful_PSD1'][ioutput], wavenumbers=outputs_shear['wavenumbers'][ioutput],min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][ioutput], DoF=12)
    # Compute decimal logarithm of epsilon from shear sensor 2 using MLE with 12 degrees of freedom
    log_epsilon2_MLE = fit_Nasmyth_MLE(psd=outputs_shear['useful_PSD2'][ioutput], wavenumbers=outputs_shear['wavenumbers'][ioutput],min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][ioutput], DoF=12)
    # Compute decimal logarithm of mean of MLE epsilons from shear sencors 1 and 2
    log_epsilon_MLE_mean = np.log10((10**log_epsilon1_MLE + 10**log_epsilon2_MLE) / 2)

    # Compute decimal logarithm of epsilon from shear sensor 1 using Lueck's integral method IM
    log_epsilon1_IM = fit_Nasmyth_Lueck(psd=outputs_shear['useful_PSD1'][ioutput], wavenumbers=outputs_shear['wavenumbers'][ioutput], min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][ioutput])
    # Compute decimal logarithm of epsilon from shear sensor 2 using IM
    log_epsilon2_IM = fit_Nasmyth_Lueck(psd=outputs_shear['useful_PSD2'][ioutput], wavenumbers=outputs_shear['wavenumbers'][ioutput], min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][ioutput])
    # Compute decimal logarithm of mean of IM epsilons from shear sencors 1 and 2
    log_epsilon_IM_mean= np.log10((10 ** log_epsilon1_IM + 10 ** log_epsilon2_IM) / 2)

    nu = 1.702747 - 0.05126103 * outputs_shear['tempcor'][ioutput] + 0.0005918645 * outputs_shear['tempcor'][ioutput] ** 2  # MSSpro user manual
    nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual, nu is viscosity
    Re_b = 10**log_epsilon_MLE_mean / (nu * outputs_shear['N2'][ioutput])   # buoyancy Reynolds number

    # Append to dictionaries by station-route
    if postaja_ruta in prvi_slovar.keys():
        prvi_slovar[postaja_ruta][0].append(globina)
        prvi_slovar[postaja_ruta][1].append(log_epsilon_MLE_mean)
        prvi_slovar[postaja_ruta][2].append(Re_b)
        prvi_slovar[postaja_ruta][3].append(outputs_shear['peps'][ioutput])
        prvi_slovar[postaja_ruta][4].append(log_epsilon1_MLE)
        prvi_slovar[postaja_ruta][5].append(log_epsilon2_MLE)
        prvi_slovar[postaja_ruta][6].append(log_epsilon_IM_mean)
        prvi_slovar[postaja_ruta][7].append(log_epsilon1_IM)
        prvi_slovar[postaja_ruta][8].append(log_epsilon2_IM)
    else:
        print(postaja_ruta)
        prvi_slovar[postaja_ruta] = [[globina], [log_epsilon_MLE_mean], [Re_b], [outputs_shear['peps'][ioutput]], [log_epsilon1_MLE], [log_epsilon2_MLE], [log_epsilon_IM_mean], [log_epsilon1_IM], [log_epsilon2_IM]]


drugi_slovar = {}   # the second dictionary
for postaja_ruta in prvi_slovar.keys(): # Look at each combination station-depth
    seznam_seznamov = np.array(prvi_slovar[postaja_ruta])   # the list of lists
    razvrscen_seznam_seznamov = seznam_seznamov[:, seznam_seznamov[0].argsort()]    # list of lists, sorted by depth
    drugi_slovar[postaja_ruta] = {'Press':razvrscen_seznam_seznamov[0], \
                                  'log_epsilon_MLE_mean':razvrscen_seznam_seznamov[1], \
                                  'Re_b':razvrscen_seznam_seznamov[2], \
                                  'log_peps':razvrscen_seznam_seznamov[3], \
                                  'log_epsilon_MLE_1':razvrscen_seznam_seznamov[4], \
                                  'log_epsilon_MLE_2':razvrscen_seznam_seznamov[5], \
                                  'log_epsilon_IM_mean':razvrscen_seznam_seznamov[6], \
                                  'log_epsilon_IM_1':razvrscen_seznam_seznamov[7], \
                                  'log_epsilon_IM_2':razvrscen_seznam_seznamov[8]}

if save_outputs_as_pickle:
    pickle.dump(drugi_slovar, open(f'slovar_epsilon_za_risat_NA_1m_{datumska_mapa}.p', 'wb'))
# The keys of the dumped dictionary are 'LK01_RUTA_01', 'LK02_RUTA_02', etc.
# The values are also dictionaries.