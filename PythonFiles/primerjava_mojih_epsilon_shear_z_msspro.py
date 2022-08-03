'''This file compares the results for epsilon from MSS90's data from different fitting techniques:
- MLE (from fit_Nasmyth_MLE from fitting_procedures_for_Nasmyth_curve.py)
- IM (from fit_Nasmyth_Lueck from fitting_procedures_for_Nasmyth_curve.py)
- MSSpro (from eallX.msb)'''

from fitting_procedures_for_Nasmyth_curve import fit_Nasmyth_MLE, fit_Nasmyth_Lueck
import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy.stats

savefigs = False # Whether it saves figures or not
mapa_za_shranjevanje_grafov = 'Saved_figures\ '[:-1]    # Where it saves the figures

datumska_mapa = 'all_2008_08_21'    # chose date ('all_2008_08_21' or 'all_2008_11_18' or'all_2009_04_20')

Slovenian_labels = False    # True for Slovenian labels, False for English labels




# THE BASIC MANUAL SETTINGS END HERE





min_wn = 6  # minimum wavenumber in fit
max_wn = 35 # maximum wavenumber in fit

datum = f'{datumska_mapa[-2:]}.{datumska_mapa[-5:-3]}.{datumska_mapa[-10:-6]}'  # nicer string for date

# Load data from pickle file
outputs = pickle.load(open(f'shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p', 'rb'))


# COMPUTE DIFFERENCES

# Set boundaries for histograms
log10_epsilon_bounds = [-10, -1]
log10_epsilons = [[i, i+1] for i in range(log10_epsilon_bounds[0], log10_epsilon_bounds[1])]
log10_epsilons = [[-np.infty, -8], [-8, -7], [-7, -6], [-6, np.infty]]  # Comment this line to have more histograms

epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean
epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean
epsilon1_minus_epsilon2_ISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1
epsilon1_minus_epsilon2_ANISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1
epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean_MLE
epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean_MLE
epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1_MLE
epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1_MLE
epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean_MLE
epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC = [[] for i in log10_epsilons] # sorted by outcomes of epsilon_my_mean_MLE
too_large_peps_my_mean = 0
too_large_peps_my_mean_MLE = 0
epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1_MLE
epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1_MLE
epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1_MLE
epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC = [[] for i in log10_epsilons]  # sorted by outcomes of epsilon1_MLE

for ioutput in range(len(outputs['epsilon_shear'])):
    print(ioutput, len(outputs['epsilon_shear']))
    epsilon1 = fit_Nasmyth_Lueck(psd=outputs['useful_PSD1'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput], min_wavenumber=min_wn, max_wavenumber=max_wn)
    epsilon2 = fit_Nasmyth_Lueck(psd=outputs['useful_PSD2'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput], min_wavenumber=min_wn, max_wavenumber=max_wn)

    epsilon1_MLE = fit_Nasmyth_MLE(psd=outputs['useful_PSD1'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput], min_wavenumber=min_wn, max_wavenumber=max_wn, DoF=12)
    epsilon2_MLE = fit_Nasmyth_MLE(psd=outputs['useful_PSD2'][ioutput], wavenumbers=outputs['wavenumbers'][ioutput], temp_mean=outputs['tempcor'][ioutput], min_wavenumber=min_wn, max_wavenumber=max_wn, DoF=12)

    epsilon_my_mean = np.log10(np.mean([10**epsilon1, 10**epsilon2]))
    epsilon1_minus_epsilon2 = epsilon1 - epsilon2
    epsilon_my_mean_MLE = np.mean([epsilon1_MLE, epsilon2_MLE])
    epsilon1_MLE_minus_epsilon2_MLE = epsilon1_MLE - epsilon2_MLE
    epsilon1_MLE_minus_epsilon1_msspro = epsilon1_MLE - outputs['epsilon_shear_1_MSSpro'][ioutput]
    epsilon2_MLE_minus_epsilon2_msspro = epsilon2_MLE - outputs['epsilon_shear_2_MSSpro'][ioutput]


    if epsilon_my_mean > outputs['peps'][ioutput] + 1:
        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
        nu *= 10**-6

        if 10**epsilon_my_mean / (nu * outputs['N2'][ioutput]) > 200:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon_my_mean >= log10_epsilons[jlog10_epsilon][0] and epsilon_my_mean < log10_epsilons[jlog10_epsilon][1]:
                    epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon].append(epsilon_my_mean - outputs['epsilon_shear'][ioutput])
                if epsilon1 >= log10_epsilons[jlog10_epsilon][0] and epsilon1 < log10_epsilons[jlog10_epsilon][1]:
                    epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon].append(epsilon1_minus_epsilon2)
        else:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon_my_mean >= log10_epsilons[jlog10_epsilon][0] and epsilon_my_mean < log10_epsilons[jlog10_epsilon][1]:
                    epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon].append(epsilon_my_mean - outputs['epsilon_shear'][ioutput])
                if epsilon1 >= log10_epsilons[jlog10_epsilon][0] and epsilon1 < log10_epsilons[jlog10_epsilon][1]:
                    epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon].append(epsilon1_minus_epsilon2)
    else:
        too_large_peps_my_mean += 1

    if epsilon_my_mean_MLE > outputs['peps'][ioutput] + 1:
        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
        nu *= 10**-6

        if 10**epsilon_my_mean_MLE / (nu * outputs['N2'][ioutput]) > 200:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon_my_mean_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon_my_mean_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon].append(epsilon_my_mean_MLE - outputs['epsilon_shear'][ioutput])
                    epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon].append(epsilon_my_mean_MLE - epsilon_my_mean)
                if epsilon1_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon1_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon].append(epsilon1_MLE_minus_epsilon2_MLE)
        else:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon_my_mean_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon_my_mean_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon].append(epsilon_my_mean_MLE - outputs['epsilon_shear'][ioutput])
                    epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon].append(epsilon_my_mean_MLE - epsilon_my_mean)
                if epsilon1 >= log10_epsilons[jlog10_epsilon][0] and epsilon1 < log10_epsilons[jlog10_epsilon][1]:
                    epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon].append(epsilon1_MLE_minus_epsilon2_MLE)
    else:
        too_large_peps_my_mean_MLE += 1

    if epsilon1_MLE > outputs['peps'][ioutput] + 1:
        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
        nu *= 10**-6

        if 10**epsilon1_MLE / (nu * outputs['N2'][ioutput]) > 200:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon1_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon1_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon].append(epsilon1_MLE_minus_epsilon1_msspro)
        else:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon1_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon1_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon].append(epsilon1_MLE_minus_epsilon1_msspro)
    if epsilon2_MLE > outputs['peps'][ioutput] + 1:
        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
        nu *= 10**-6

        if 10**epsilon2_MLE / (nu * outputs['N2'][ioutput]) > 200:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon2_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon2_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon].append(epsilon2_MLE_minus_epsilon2_msspro)
        else:
            for jlog10_epsilon in range(len(log10_epsilons)):
                if epsilon2_MLE >= log10_epsilons[jlog10_epsilon][0] and epsilon2_MLE < log10_epsilons[jlog10_epsilon][1]:
                    epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon].append(epsilon2_MLE_minus_epsilon2_msspro)

print(too_large_peps_my_mean, too_large_peps_my_mean / len(outputs['epsilon_shear']))
print(too_large_peps_my_mean_MLE, too_large_peps_my_mean_MLE / len(outputs['epsilon_shear']))

# 1st PLOT: FINAL EPSILON IM vs FINAL EPSILON MSSpro
plt.figure(figsize=(12, 10))
labels = [r'$Re_{b} > 200$',r'$Re_{b} \leq 200$']
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon], epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        plt.annotate('N= {}'.format(len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] +
                               epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')

        if jlog10_epsilon not in [0]:   # Which figures to annotate on the right side
            plt.annotate(
                f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.6, 0.86), xycoords='axes fraction')
            plt.annotate(
                f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.6, 0.77), xycoords='axes fraction')
        elif jlog10_epsilon == 0:   # The first figure has a legend, so lower annotate
            plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.68), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.59), xycoords='axes fraction')
        else:   # Which figures to annotate on the left side
            plt.annotate(
                f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(
                f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}',
                xy=(0.05, 0.77), xycoords='axes fraction')
        if len(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
            plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.68), xycoords='axes fraction', color='C0')
            plt.annotate(f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.59), xycoords='axes fraction', color='C0')

    if jlog10_epsilon == 0:
        plt.legend()
    plt.title(r'$\log{(\epsilon_{shear}^{integral})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear}^{integral})} - \log{(\epsilon_{shear}^{MSSpro})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear}^{integral})}$ and $\log{(\epsilon_{shear}^{MSSpro})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear}^{integral})}$ in $\log{(\epsilon_{shear}^{MSSpro})}$ za ' + datum)

plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_integral_vs_epsilon_MSSpro_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_integral_vs_epsilon_MSSpro_{datumska_mapa}_si')
plt.show()


# 2nd PLOT: EPSILON 1 IM VS EPSILON 2 IM

plt.figure(figsize=(12, 10))
labels = [r'$Re_{b} > 200$',r'$Re_{b} \leq 200$']
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon], epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        plt.annotate('N= {}'.format(len(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] +
                               epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')
        if jlog10_epsilon in []: # Which figures to annotate on the LEFT side
            plt.annotate(f'mean = {round(np.mean(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.77), xycoords='axes fraction')
            if len(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.59), xycoords='axes fraction', color='C0')
        else: # Which figures to annotate on the RIGHT side
            plt.annotate(f'mean = {round(np.mean(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon] + epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.77), xycoords='axes fraction')
            if len(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon1_minus_epsilon2_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon1_minus_epsilon2_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.59), xycoords='axes fraction', color='C0')
    if jlog10_epsilon == 0:
        plt.legend(loc='upper right')
    plt.title(r'$\log{(\epsilon_{shear_1}^{integral})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear_1}^{integral})} - \log{(\epsilon_{shear_2}^{integral})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear_1}^{integral})}$ and $\log{(\epsilon_{shear_2}^{integral})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear_1}^{integral})}$ in $\log{(\epsilon_{shear_2}^{integral})}$ za ' + datum)

plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear1_integral_vs_epsilon_shear2_integral_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear1_integral_vs_epsilon_shear2_integral_{datumska_mapa}_si')
plt.show()



# 3rd PLOT: FINAL EPSILON MLE VS FINAL EPSILON MSSpro

plt.figure(figsize=(12, 10))
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon], epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.86), xycoords='axes fraction')
        plt.annotate(f'std = {round(np.std(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.77), xycoords='axes fraction')
        plt.annotate('N= {}'.format(len(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] +
                                        epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')

        if len(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_my_mean_MLE_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
            plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.68), xycoords='axes fraction', color='C0')
            plt.annotate(f'std = {round(np.std(epsilon_my_mean_MLE_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.59), xycoords='axes fraction', color='C0')

    if jlog10_epsilon == 0:
        plt.legend()
    plt.title(r'$\log{(\epsilon_{shear}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear}^{MLE})} - \log{(\epsilon_{shear}^{MSSpro})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear}^{MLE})}$ and $\log{(\epsilon_{shear}^{MSSpro})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear}^{MLE})}$ in $\log{(\epsilon_{shear}^{MSSpro})}$ za ' + datum)
plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_MLE_vs_epsilon_MSSpro_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_MLE_vs_epsilon_MSSpro_{datumska_mapa}_si')
plt.show()



# 4th PLOT: EPLSILON 1 MLE VS EPSILON 2 MLE

plt.figure(figsize=(12, 10))
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon], epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        plt.annotate('N= {}'.format(len(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] +
                                        epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')
        if jlog10_epsilon in []:     # Which figures to annotate on the LEFT side
            plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.77), xycoords='axes fraction')
            if len(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.59), xycoords='axes fraction', color='C0')
        else:    # Which figures to annotate on the RIGHT side
            plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.77), xycoords='axes fraction')
            if len(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon1_MLE_minus_epsilon2_MLE_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon2_MLE_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.6, 0.59), xycoords='axes fraction', color='C0')
    if jlog10_epsilon == 0:
        plt.legend(loc='upper right')
    plt.title(r'$\log{(\epsilon_{shear_1}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear_1}^{MLE})} - \log{(\epsilon_{shear_2}^{MLE})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear_1}^{MLE})}$ and $\log{(\epsilon_{shear_2}^{MLE})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear_1}^{MLE})}$ in $\log{(\epsilon_{shear_2}^{MLE})}$ za ' + datum)
plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear1_MLE_vs_epsilon_shear2_MLE_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear1_MLE_vs_epsilon_shear2_MLE_{datumska_mapa}_si')
plt.show()




# 5th PLOT: FINAL EPSILON MLE VS FINAL EPSILON IM

plt.figure(figsize=(12, 10))
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon], epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        #textstr = f'mean = {round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]) ,2)}\n' + f'std = {round(np.std(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_minus_epsilon_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}'
        plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.65, 0.86), xycoords='axes fraction')
        plt.annotate(f'std = {round(np.std(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] + epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.65, 0.77), xycoords='axes fraction')
        plt.annotate('N= {}'.format(len(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon] +
                                        epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')

        if len(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon_my_mean_MLE_minus_epsilon_my_mean_ANISOTROPIC[jlog10_epsilon]) > 0:
            #pass
            #textstr += '\n' + r'\textcolor{C0}{mean = {}}'.format(round(np.mean(epsilon_my_mean_minus_epsilon_msspro_ISOTROPIC[jlog10_epsilon]), 2))
            plt.annotate(f'mean = {round(np.mean(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.65, 0.68), xycoords='axes fraction', color='C0')
            plt.annotate(f'std = {round(np.std(epsilon_my_mean_MLE_minus_epsilon_my_mean_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.65, 0.59), xycoords='axes fraction', color='C0')

    if jlog10_epsilon == 0:
        plt.legend()
    plt.title(r'$\log{(\epsilon_{shear}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear}^{MLE})} - \log{(\epsilon_{integral}^{MSSpro})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear}^{MLE})}$ and $\log{(\epsilon_{shear}^{integral})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear}^{MLE})}$ in $\log{(\epsilon_{shear}^{integral})}$ za ' + datum)
plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_MLE_vs_epsilon_shear_integral_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear_MLE_vs_epsilon_shear_integral_{datumska_mapa}_si')
plt.show()



# 6th PLOT: EPSILON 2 MLE VS EPSILON 2 MSSpro

plt.figure(figsize=(12, 10))
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon], epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        plt.annotate('N= {}'.format(len(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] +
                                        epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')
        if jlog10_epsilon in []: # Which figures to annotate on the LEFT side
            plt.annotate(f'mean = {round(np.mean(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.77), xycoords='axes fraction')
            if len(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.59), xycoords='axes fraction', color='C0')
        else: # Which figures to annotate on the RIGHT side
            plt.annotate(f'mean = {round(np.mean(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon] + epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.77), xycoords='axes fraction')
            if len(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon2_MLE_minus_epsilon2_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon2_MLE_minus_epsilon2_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.59), xycoords='axes fraction', color='C0')
    if jlog10_epsilon == 0:
        plt.legend(loc='upper right')
    plt.title(r'$\log{(\epsilon_{shear_2}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear_2}^{MLE})} - \log{(\epsilon_{shear_2}^{MSSpro})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear_2}^{MLE})}$ and $\log{(\epsilon_{shear_2}^{MSSpro})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear_2}^{MLE})}$ in $\log{(\epsilon_{shear_2}^{MSSpro})}$ za ' + datum)
plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear2_MLE_vs_epsilon_shear2_MSSpro_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear2_MLE_vs_epsilon_shear2_MSSpro_{datumska_mapa}_si')
plt.show()



# 7th FIGURE: EPSILON 1 MSSpro VS EPSILON 1 MLE

plt.figure(figsize=(12, 10))
for jlog10_epsilon in range(len(log10_epsilons)):
    plt.subplot(3, 3, jlog10_epsilon+1)
    plt.hist([epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon], epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]], stacked=True, bins=30, label=labels)
    if len(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
        plt.axvline(np.median(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]), color='k')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon], 0.75), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon], 0.25), color='k', linestyle='--')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon], 0.10), color='k', linestyle=':')
        plt.axvline(np.quantile(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon], 0.90), color='k', linestyle=':')
        plt.annotate('N= {}'.format(len(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] +
                                        epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.86),
                     xycoords='axes fraction')
        plt.annotate('Ni= {}'.format(len(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon])), xy=(0.1, 0.66),
                     xycoords='axes fraction')
        if jlog10_epsilon in []: # Which figures to annotate on the LEFT side
            plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.77), xycoords='axes fraction')
            if len(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.05, 0.59), xycoords='axes fraction', color='C0')
        else: # Which figures to annotate on the RIGHT side
            plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.86), xycoords='axes fraction')
            plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon] + epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.77), xycoords='axes fraction')
            if len(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon]) > 0 and len(epsilon1_MLE_minus_epsilon1_msspro_ANISOTROPIC[jlog10_epsilon]) > 0:
                plt.annotate(f'mean = {round(np.mean(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.68), xycoords='axes fraction', color='C0')
                plt.annotate(f'std = {round(np.std(epsilon1_MLE_minus_epsilon1_msspro_ISOTROPIC[jlog10_epsilon]), 2)}', xy=(0.67, 0.59), xycoords='axes fraction', color='C0')
    if jlog10_epsilon == 0:
        plt.legend(loc='upper right')
    plt.title(r'$\log{(\epsilon_{shear_1}^{MLE})} \in [' + str(log10_epsilons[jlog10_epsilon][0]) + ', ' + str(log10_epsilons[jlog10_epsilon][1]) +')$')
    plt.xlabel(r'$\log{(\epsilon_{shear_1}^{MLE})} - \log{(\epsilon_{shear_1}^{MSSpro})}$')
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
    plt.suptitle(r'Difference between $\log{(\epsilon_{shear_1}^{MLE})}$ and $\log{(\epsilon_{shear_1}^{MSSpro})}$ for ' + datum)
else:
    plt.suptitle(r'Razlika med $\log{(\epsilon_{shear_1}^{MLE})}$ in $\log{(\epsilon_{shear_1}^{MSSpro})}$ za ' + datum)
plt.subplots_adjust(top=0.92)
if savefigs:
    if not Slovenian_labels:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear1_MLE_vs_epsilon_shear1_MSSpro_{datumska_mapa}')
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_shear1_MLE_vs_epsilon_shear1_MSSpro_{datumska_mapa}_si')
plt.show()