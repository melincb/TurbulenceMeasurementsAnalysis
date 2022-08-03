'''By this script we calculate epsilon from temperature profiles. To save the results uncomment the last line'''
import numpy as np
import pickle
from fitting_procedures_for_Batchelor_curve import iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral, fit_Batchelorjeve_krivulje_na_razliko_logaritmov

# Define which outputs to use
datumska_mapa = 'all_2008_11_18'    # 'all_2008_08_21' or 'all_2008_11_18' or 'all_2009_04_20'
outputs = pickle.load(open('grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m' + datumska_mapa + '.p', 'rb'))
# also works with outputs = pickle.load(open('grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji' + datumska_mapa + '.p', 'rb'))


# Let's apply the MLE fit on outputs!

new_outputs_na_1m = {'epsilon_shear':[], 'useful_grad_T_PSD':[], 'peps':[], 'tempcor':[], 'wavenumbers':[], 'thermdiss_MSSpro':[], 'st_zdruzenih':[], 'postaja_ruta_globina':[], 'N2':[], 'epsilon_grad_T_DoF6':[], 'epsilon_grad_T_DoF12':[], 'thermdiss_DoF12':[], 'thermdiss_DoF6':[], 'min_fitted_wn':[]}


for ioutput in range(len(outputs['tempcor'])):
    print(ioutput, len(outputs['tempcor']))

    success = True
    try:
        # The initial guess for epsilon comes from LSnTD (i.e. fit_Batchelorjeve_krivulje_na_razliko_logaritmov)
        grad_T_fit_start = fit_Batchelorjeve_krivulje_na_razliko_logaritmov(psd=outputs['useful_grad_T_PSD'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput])[0]
        # Let's try MLE with 6 degrees of freedom
        grad_T_fit_final_DoF6 = iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(psd=outputs['useful_grad_T_PSD'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput], starting_epsilon=10**grad_T_fit_start, DoF=6)
        # Let's try MLE with 12 degrees of freedom
        grad_T_fit_final_DoF12 = iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(psd=outputs['useful_grad_T_PSD'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput], starting_epsilon=10**grad_T_fit_start, DoF=12)
        epsilon_DoF6 = np.log10(grad_T_fit_final_DoF6[0])
        thermdiss_DoF6 = np.log10(grad_T_fit_final_DoF6[1])
        epsilon_DoF12 = np.log10(grad_T_fit_final_DoF12[0])
        thermdiss_DoF12 = np.log10(grad_T_fit_final_DoF12[1])
        min_wn = grad_T_fit_final_DoF12[2]
        print(grad_T_fit_start, epsilon_DoF12, epsilon_DoF6)
    except:
        success = False

    if success:
        new_outputs_na_1m['epsilon_shear'].append(outputs['epsilon_shear'][ioutput])
        new_outputs_na_1m['useful_grad_T_PSD'].append(outputs['useful_grad_T_PSD'][ioutput])
        new_outputs_na_1m['peps'].append(outputs['peps'][ioutput])
        new_outputs_na_1m['tempcor'].append(outputs['tempcor'][ioutput])
        new_outputs_na_1m['wavenumbers'].append(outputs['wavenumbers'][ioutput])
        new_outputs_na_1m['thermdiss_MSSpro'].append(outputs['thermdiss_MSSpro'][ioutput])
        new_outputs_na_1m['st_zdruzenih'].append(outputs['st_zdruzenih'][ioutput])
        new_outputs_na_1m['postaja_ruta_globina'].append(outputs['postaja_ruta_globina'][ioutput])
        new_outputs_na_1m['N2'].append(outputs['N2'][ioutput])
        new_outputs_na_1m['epsilon_grad_T_DoF6'].append(epsilon_DoF6)
        new_outputs_na_1m['epsilon_grad_T_DoF12'].append(epsilon_DoF12)
        new_outputs_na_1m['thermdiss_DoF6'].append(thermdiss_DoF6)
        new_outputs_na_1m['thermdiss_DoF12'].append(thermdiss_DoF12)
        new_outputs_na_1m['min_fitted_wn'].append(min_wn)

#pickle.dump(new_outputs_na_1m, open('grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_T' + datumska_mapa + '.p', 'wb'))