'''This file compares the results for epsilon from MSS90's data from different fitting techniques:
- MLE for temperature fluctuation gradients (from iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral from
    fitting_procedures_for_Batchelor_curve.py (already computed in dejansti_izracun_epsilon_T.py))
- MSSpro (from eallX.msb)
- MLE (from fit_Nasmyth_MLE from fitting_procedures_for_Nasmyth_curve.py)'''
from fitting_procedures_for_Nasmyth_curve import fit_Nasmyth_MLE
import matplotlib.pyplot as plt
import numpy as np
import pickle


savefigs = False # Whether it saves figures or not
mapa_za_shranjevanje_grafov = 'Saved_figures\ '[:-1]    # Where it saves the figures

datumska_mapa = 'all_2008_08_21'    # chose date ('all_2008_08_21' or 'all_2008_11_18' or'all_2009_04_20')

Slovenian_labels = False    # True for Slovenian labels, False for English labels









datum = f'{datumska_mapa[-2:]}.{datumska_mapa[-5:-3]}.{datumska_mapa[-10:-6]}'

# Load temperature data
outputs = pickle.load(open(f'grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_T{datumska_mapa}.p', 'rb'))



# COMPUTE DIFFERENCES BETWEEN EPSILON FROM BATCHELOR SPECTRUM AND SHEAR EPSILON FROM MSSpro

# Set boundaries for histograms
log10_epsilon_bounds = [-10, -1]
log10_epsilons = [[i, i+1] for i in range(log10_epsilon_bounds[0], log10_epsilon_bounds[1])]
log10_epsilons = [[-np.infty, -8], [-8, -7], [-7, -6], [-6, np.infty]]  # Comment this line to have more histograms

epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean
epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean
epsilon12_minus_epsilon6_ISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon12
epsilon12_minus_epsilon6_UNISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon12
thermdiss_my_mean_minus_thermdiss_msspro = []
too_large_peps = 0


different_epsilon_12_and_epsilon_6 = 0

for ioutput in range(len(outputs['epsilon_shear'])):
    epsilon12 = outputs['epsilon_grad_T_DoF12'][ioutput]
    epsilon6 = outputs['epsilon_grad_T_DoF6'][ioutput]
    thermdiss12 = outputs['thermdiss_DoF12'][ioutput]
    thermdiss6 = outputs['thermdiss_DoF6'][ioutput]
    thermdiss_my_mean_minus_thermdiss_msspro.append(np.mean([thermdiss6, thermdiss12]) - np.log10(outputs['thermdiss_MSSpro'][ioutput]))

    if epsilon12 != epsilon6:
        different_epsilon_12_and_epsilon_6 += 1

    epsilon_my_mean = np.mean([epsilon12, epsilon6])
    epsilon12_minus_epsilon6 = epsilon12 - epsilon6

    if epsilon_my_mean > outputs['peps'][ioutput] + 1:
        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
        nu *= 10**-6

        if 10**epsilon_my_mean / (nu * outputs['N2'][ioutput]) > 200:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon_my_mean >= log10_epsilons[jlog10_epsilon][0] and epsilon_my_mean < log10_epsilons[jlog10_epsilon][1]:
                    epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon].append(epsilon_my_mean - outputs['epsilon_shear'][ioutput])
                if epsilon12 >= log10_epsilons[jlog10_epsilon][0] and epsilon12 < log10_epsilons[jlog10_epsilon][1]:
                    epsilon12_minus_epsilon6_ISOTROPIC[jlog10_epsilon].append(epsilon12_minus_epsilon6)
        else:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon_my_mean >= log10_epsilons[jlog10_epsilon][0] and epsilon_my_mean < log10_epsilons[jlog10_epsilon][1]:
                    epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon].append(epsilon_my_mean - outputs['epsilon_shear'][ioutput])
                if epsilon12 >= log10_epsilons[jlog10_epsilon][0] and epsilon12 < log10_epsilons[jlog10_epsilon][1]:
                    epsilon12_minus_epsilon6_UNISOTROPIC[jlog10_epsilon].append(epsilon12_minus_epsilon6)
    else:
        too_large_peps += 1


print(too_large_peps, too_large_peps/len(outputs['epsilon_shear']))
print(different_epsilon_12_and_epsilon_6, different_epsilon_12_and_epsilon_6/len(outputs['epsilon_shear']))


# 1st PLOT: FINAL EPSILON FROM TEMPERATURE VS FINAL EPSILON FROM MSSpro

plt.figure(figsize=(12, 10))
labels = [r'$Re_{b} > 200$',r'$Re_{b} \leq 200$']
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon], epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]), color='k', alpha=0.85)
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--', alpha=0.85)
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--', alpha=0.85)
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':', alpha=0.85)
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':', alpha=0.85)
        if jlog10_epsilon in [1,2]:
            plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]) ,2)}', xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.77), xycoords='axes fraction')
            if len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.59), xycoords='axes fraction', color='C0')
        elif jlog10_epsilon in []:
            plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]) ,2)}', xy=(0.35, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.77), xycoords='axes fraction')
            if len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.59), xycoords='axes fraction', color='C0')

        else:
            plt.annotate(
                f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.68, 0.86), xycoords='axes fraction')
            plt.annotate(
                f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.68, 0.77), xycoords='axes fraction')
            if len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(
                    epsilon_my_mean_minus_epsilon_msspro_UNISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(
                    f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}',
                    xy=(0.68, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(
                    f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}',
                    xy=(0.68, 0.59), xycoords='axes fraction', color='C0')

    if jlog10_epsilon == 0:
        plt.legend()
    plt.title(r'$\log{(\epsilon_{\nabla T}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{\nabla T}^{MLE})} - \log{(\epsilon_{shear}^{MSSpro})}$')
    if not Slovenian_labels:
        plt.ylabel('Occurrence')
    else:
        plt.ylabel('Št. pojavitev')
    loc, labs = plt.yticks()
    maxloc = int(max(max(loc), 1))
    if maxloc < 10:
        plt.yticks([i for i in range(0, maxloc+2)])
plt.tight_layout()
if not Slovenian_labels:
    plt.suptitle(r'Difference between $\log{(\epsilon_{\nabla T}^{MLE})}$ and $\log{(\epsilon_{shear}^{MSSpro})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{\nabla T}^{MLE})}$ in $\log{(\epsilon_{shear}^{MSSpro})}$ za ' + datum)
plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_T_MLE_vs_epsilon_MSSpro_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_T_MLE_vs_epsilon_MSSpro_{datumska_mapa}_si')
plt.show()



# 2nd PLOT: MY THERMAL DISSIPATION VS MSSpro's THERMAL DISSIPATION

plt.hist(thermdiss_my_mean_minus_thermdiss_msspro, bins=30)
plt.annotate(f'mean = {round(np.mean(thermdiss_my_mean_minus_thermdiss_msspro), 2)}', xy=(0.75, 0.86), xycoords='axes fraction')
plt.annotate(f'std = {round(np.std(thermdiss_my_mean_minus_thermdiss_msspro), 2)}', xy=(0.75, 0.77), xycoords='axes fraction')
plt.axvline(np.median(thermdiss_my_mean_minus_thermdiss_msspro), color='k')
plt.axvline(np.quantile(thermdiss_my_mean_minus_thermdiss_msspro, 0.75), color='k', linestyle='--')
plt.axvline(np.quantile(thermdiss_my_mean_minus_thermdiss_msspro, 0.25), color='k', linestyle='--')
plt.axvline(np.quantile(thermdiss_my_mean_minus_thermdiss_msspro, 0.10), color='k', linestyle=':')
plt.axvline(np.quantile(thermdiss_my_mean_minus_thermdiss_msspro, 0.90), color='k', linestyle=':')
plt.xlabel(r'$\log{(\chi_T^{MLE})} - \log{(\chi_T^{MSSpro})}$')
if not Slovenian_labels:
    plt.ylabel('Occurrence')
    plt.title(r'Difference between $\log{(\chi_T^{MLE})}$ and $\log{(\chi_T^{MSSpro})}$ for ' + datum)
    if savefigs:
        plt.savefig(mapa_za_shranjevanje_grafov + f'Thermdiss_T_MLE_vs_Thermdiss_MSSpro_{datumska_mapa}')
else:
    plt.ylabel('Št. pojavitev')
    plt.title(r'Razlika med $\log{(\chi_T^{MLE})}$ in $\log{(\chi_T^{MSSpro})}$ za ' + datum)
    if savefigs:
        plt.savefig(mapa_za_shranjevanje_grafov + f'Thermdiss_T_MLE_vs_Thermdiss_MSSpro_{datumska_mapa}_si')
plt.show()





# COMPUTE DIFFERENCES BETWEEN EPSILON FROM BATCHELOR SPECTRUM AND SHEAR EPSILON FROM MLE

# Load shear data
outputs_shear = pickle.load(open(f'shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p', 'rb'))

epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_grad_T_mean
epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_grad_T_mean

samo_shear2 = False # To compare the mean epsilon from both shear sensors with the epsilon from the Batchelor spectrum,
                    # set this to False. To only compare epsilon from shear sensor 2 with the epsilon from the Batchelor
                    # spectrum, set this to True

for ioutput in range(len(outputs['epsilon_shear'])):
    print(ioutput, len(outputs['epsilon_shear']))
    prg = outputs['postaja_ruta_globina'][ioutput]
    if prg in outputs_shear['postaja_ruta_globina']:
        joutput = outputs_shear['postaja_ruta_globina'].index(prg)
        epsilon_grad_T12 = outputs['epsilon_grad_T_DoF12'][ioutput]
        epsilon_grad_T6 = outputs['epsilon_grad_T_DoF6'][ioutput]
        epsilon_shear1 = fit_Nasmyth_MLE(psd=outputs_shear['useful_PSD1'][joutput], wavenumbers=outputs_shear['wavenumbers'][joutput], min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][joutput], DoF=12)
        epsilon_shear2 = fit_Nasmyth_MLE(psd=outputs_shear['useful_PSD2'][joutput], wavenumbers=outputs_shear['wavenumbers'][joutput], min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][joutput], DoF=12)

        if samo_shear2:
            epsilon_shear_mean = epsilon_shear2
        else:
            epsilon_shear_mean = np.log10(np.mean([10**epsilon_shear1, 10**epsilon_shear2]))
        epsilon_grad_T_mean = np.mean([epsilon_grad_T12, epsilon_grad_T6])

        if epsilon_grad_T_mean > outputs['peps'][ioutput] + 1:
            if epsilon_grad_T_mean < -9:
                print(prg)
            nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
            nu *= 10 ** -6

            if 10 ** epsilon_grad_T_mean / (nu * outputs['N2'][ioutput]) > 200:
                for jlog10_epsilon in range(len(log10_epsilons)):
                    if epsilon_grad_T_mean >= log10_epsilons[jlog10_epsilon][0] and epsilon_grad_T_mean < log10_epsilons[jlog10_epsilon][1]:
                        epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon].append(epsilon_shear_mean - epsilon_grad_T_mean)
            else:
                for jlog10_epsilon in range(len(log10_epsilons)):
                    if epsilon_grad_T_mean >= log10_epsilons[jlog10_epsilon][0] and epsilon_grad_T_mean < log10_epsilons[jlog10_epsilon][1]:
                        epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon].append(epsilon_shear_mean - epsilon_grad_T_mean)

# PLOT
plt.figure(figsize=(12, 10))
labels = [r'$Re_{b} > 200$',r'$Re_{b} \leq 200$']
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon], epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--', alpha=0.85)
        plt.axvline(np.quantile(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--', alpha=0.85)
        plt.axvline(np.quantile(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':', alpha=0.85)
        plt.axvline(np.quantile(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':', alpha=0.85)
        plt.annotate('N= {}'.format(len(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] +
                               epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')
        if jlog10_epsilon in []:
            plt.annotate(f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.77), xycoords='axes fraction')
            if len(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.35, 0.59), xycoords='axes fraction', color='C0')
        elif jlog10_epsilon in []:
            plt.annotate( f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.65, 0.66), xycoords='axes fraction')
            plt.annotate( f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.65, 0.57), xycoords='axes fraction')
        elif jlog10_epsilon in []:
            plt.annotate(f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.05, 0.77), xycoords='axes fraction')
            if len(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]), 2)}',
                    xy=(0.05, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]), 2)}',
                    xy=(0.05, 0.59), xycoords='axes fraction', color='C0')
        else:
            plt.annotate(f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.68, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon] + epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.68, 0.77), xycoords='axes fraction')
            if len(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_shear_mean_minus_epsilon_grad_T_mean_UNISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]), 2)}',
                    xy=(0.68, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon_shear_mean_minus_epsilon_grad_T_mean_ISOTROPIC[jlog10_epsilon]), 2)}',
                    xy=(0.68, 0.59), xycoords='axes fraction', color='C0')

    if jlog10_epsilon == 0:
        plt.legend()
    plt.title(r'$\log{(\epsilon_{\nabla T}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear}^{MLE})} - \log{(\epsilon_{\nabla T}^{MLE})}$')
    if not Slovenian_labels:
        plt.ylabel('Occurrence')
    else:
        plt.ylabel('Št. pojavitev')
    loc, labs = plt.yticks()
    maxloc = int(max(max(loc), 1))
    if maxloc < 10:
        plt.yticks([i for i in range(0, maxloc+2)])
plt.tight_layout()
if samo_shear2:
    if not Slovenian_labels:
        plt.suptitle(r'Difference between $\log{(\epsilon_{shear2}^{MLE})}$ and $\log{(\epsilon_{\nabla T}^{MLE})}$ for ' + datum)
    else:
        plt.suptitle(r'Razlika med $\log{(\epsilon_{shear2}^{MLE})}$ in $\log{(\epsilon_{\nabla T}^{MLE})}$ za ' + datum)
    plt.subplots_adjust(top=0.92)
    if savefigs:
        if not Slovenian_labels:
            plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear2_MLE_vs_epsilon_T_MLE_{datumska_mapa}')
        else:
            plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear2_MLE_vs_epsilon_T_MLE_{datumska_mapa}_si')
else:
    if not Slovenian_labels:
        plt.suptitle(r'Difference between $\log{(\epsilon_{shear}^{MLE})}$ and $\log{(\epsilon_{\nabla T}^{MLE})}$ for ' + datum)
    else:
        plt.suptitle(r'Razlika med $\log{(\epsilon_{shear}^{MLE})}$ in $\log{(\epsilon_{\nabla T}^{MLE})}$ za ' + datum)
    plt.subplots_adjust(top=0.92)
    if savefigs:
        if not Slovenian_labels:
            plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_MLE_vs_epsilon_T_MLE_{datumska_mapa}')
        else:
            plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_MLE_vs_epsilon_T_MLE_{datumska_mapa}_si')
plt.show()
