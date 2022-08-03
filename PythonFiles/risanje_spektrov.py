'''Plots the spectra at combinations station-route-depth. Shear fluctuation spectra from both shear sensors with the
best fitting Nasmyth spectrum from fit_Nasmyth_MLE from fitting_procedures_for_Nasmyth_curve.py are plotted as well as
the temperature fluctuation gradient spectrum with its best fitting Batchelor spectrum from
iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral from fitting_procedures_for_Batchelor_curve.py
(already computed in dejansti_izracun_epsilon_T.py).
If the given combination of station-route-depth is not valid in data, it will be skipped.
'''


import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
from scipy.special import erf
from fitting_procedures_for_Nasmyth_curve import fit_Nasmyth_MLE
import pickle



datumska_mapa = 'all_2008_11_18'    # chose date ('all_2008_08_21' or 'all_2008_11_18' or'all_2009_04_20')

# Set the combinations 'station_route_depth' whose spectra you want to plot (up to 6 at once)
combinations_to_plot = ['LK01_RUTA_10_14.0', 'LK01_RUTA_10_16.0', 'LK01_RUTA_03_15.0', 'LK02_RUTA_05_2.0', 'LK02_RUTA_05_9.0', 'LK03_RUTA_03_2.0']
savefigs = False    # Set to True to save the plot


# Load temperature data
outputs = pickle.load(open(f'grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_T{datumska_mapa}.p', 'rb'))
# Load shear data
outputs_shear = pickle.load(open(f'shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p', 'rb'))

def Nasmyth_navadni(k, log10epsilon, temp_mean):
    nu = 1.702747 - 0.05126103 * temp_mean + 0.0005918645 * temp_mean ** 2  # MSSpro user manual
    nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual

    epsilon = 10 ** log10epsilon
    x = k * (nu ** 3 / epsilon) ** 0.25
    return 8.05 * x ** (1 / 3) / (1 + (20.6 * x) ** 3.715) * (epsilon ** 3 / nu) ** 0.25

def Batchelor_navadni(k, log10_epsilon, temp_mean, Thermdiss=10**-6):
    nu = 1.702747 - 0.05126103 * temp_mean + 0.0005918645 * temp_mean ** 2  # MSSpro user manual
    nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual
    Dt = 1.4 * 10 ** -7  # m**2 s**-1 (MSSpro user manual)
    q = 2 * np.sqrt(3)
    epsilon = 10**log10_epsilon
    k_b = (epsilon/(nu * Dt**2))**0.25 / (2*np.pi)  # cpm

    y = k * np.sqrt(q) / k_b
    return Thermdiss * np.sqrt(q) / (Dt * k_b) * y**2 * (np.exp(-y**2)/y - np.sqrt(np.pi)*(1 - erf(y)))
def get_thermdiss(current_epsilon, final_psd, final_wavenumbers):
    '''This is originally part of iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral from
    fitting_procedures_for_Batchelor_curve.py'''
    integral = scipy.integrate.trapz(y=final_psd, x=final_wavenumbers)
    Dt = 1.4 * 10 ** -7  # m**2 s**-1 (MSSpro user manual)
    q = 2 * np.sqrt(3)
    measured_part = 6 * Dt * integral

    y_star = np.sqrt(q) * 7.9 * 10**-3
    y_min = 2 * np.pi * final_wavenumbers[0] * np.sqrt(q) * (nu * Dt ** 2 / current_epsilon)**(1/4)
    y_max = 2 * np.pi * final_wavenumbers[-1] * np.sqrt(q) * (nu * Dt ** 2 / current_epsilon)**(1/4)

    theory_0_to_inf = 9 * y_star**(4/3) / (40 * Dt**(1/6) * current_epsilon**(1/4) * nu**(1/12) * q**(1/6)) + np.exp(-y_star**2) * np.sqrt(q) * (1 - 2*y_star**2 + np.exp(y_star**2) * y_star**2 * (-3 + 2*np.sqrt(np.pi) * y_star * scipy.special.erfc(y_star))) / (6 * Dt * (current_epsilon / Dt**2 / nu)**(1/4))
    if y_star <= y_min:
        theory_min_to_max = 1 / (6 * Dt * (current_epsilon / (Dt**2 * nu))**(1/4)) * np.exp(-y_max**2 - y_min**2) * np.sqrt(q) * (np.exp(y_max**2) * (1 - 2*y_min**2) + np.exp(y_min**2)*(-1 + 2*y_max**2 - 2* np.exp(y_max**2) * np.sqrt(np.pi) * (y_max**3 * scipy.special.erfc(y_max) - y_min**3 * scipy.special.erfc(y_min))))
    elif y_star < y_max:
        theory_min_to_max = - 9 * (y_min**(4/3) - y_star**(4/3)) / (40 * Dt**(1/6) * current_epsilon**(1/4) * nu**(1/12) * q**(1/6)) \
                            + 1 / (6 * Dt * (current_epsilon / (Dt**2 * nu))**(1/4)) * np.exp(-y_max**2 - y_star**2) * np.sqrt(q) * (np.exp(y_max**2) * (1 - 2*y_star**2) + np.exp(y_star**2)*(-1 + 2*y_max**2 - 2* np.exp(y_max**2) * np.sqrt(np.pi) * (y_max**3 * scipy.special.erfc(y_max) - y_star**3 * scipy.special.erfc(y_star))))
    else:
        theory_min_to_max = 9 * (y_max**(4/3) - y_min**(4/3)) / (40 * Dt**(1/6) * current_epsilon**(1/4) * nu**(1/12) * q**(1/6))

    theoretical_fraction = theory_min_to_max / theory_0_to_inf

    return measured_part / theoretical_fraction

fig, axs = plt.subplots(4, 3, figsize=(12, 16))
def plot_stuff(axsh, axt, wnsh, psd1, psd2, es1, es2, et, T, psdt, wnt, td, prg, min_wnT=1, legende=False):
    nu = 1.702747 - 0.05126103*T + 0.0005918645*T**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual
    Dt = 1.4 * 10 ** -7  # m**2 s**-1 (MSSpro user manual)

    axsh.scatter(wnsh, psd1, color='C0', label=r'$\mathrm{PSD}_{s1}$')
    wnshs = np.linspace(wnsh[0], wnsh[-1], 500)
    axsh.plot(wnshs, Nasmyth_navadni(wnshs, es1, T), color='k', linewidth=2)
    axsh.plot(wnshs, Nasmyth_navadni(wnshs, es1, T), color='C0', linewidth=1.5)
    axsh.scatter(wnsh, psd2, color='C1', label=r'$\mathrm{PSD}_{s2}$')
    axsh.plot(wnshs, Nasmyth_navadni(wnshs, es2, T), color='k', linewidth=2)
    axsh.plot(wnshs, Nasmyth_navadni(wnshs, es2, T), color='C1', linewidth=1.5)
    axsh.set_xscale('log')
    axsh.set_xticks([1,2,3,5,10,20,30])
    axsh.set_xticklabels([1,2,3,5,10,20,30], fontsize=14)
    axsh.set_yticklabels([10**(-i) for i in range(-12, 0)], fontsize=14)
    axsh.set_xlabel(r'$\kappa$ [cpm]', fontsize=16)
    axsh.set_yscale('log')
    axsh.axvline(6, linestyle='--', linewidth=1, color='k')
    axsh.axvline(35, linestyle='--', linewidth=1, color='k')
    axsh.set_ylabel(r'PSD [$(\mathrm{s}^{-1})^{2}$/cpm]', fontsize=16)
    axsh.annotate(r'$\log{(\epsilon_{s1})}=$'+str(round(es1, 2)), xy=(0.2, 0.15), xycoords='axes fraction', color='C0', fontsize=16)
    axsh.annotate(r'$\log{(\epsilon_{s2})}=$'+str(round(es2, 2)), xy=(0.2, 0.05), xycoords='axes fraction', color='C1', fontsize=16)
    axt.scatter(wnt, psdt, color='C2', label=r'$\mathrm{PSD}_{T}$')
    k_bt = (10**et / nu / Dt ** 2) ** 0.25 / (2 * np.pi)
    k_start = 7.9*10**-3 * k_bt
    wnts = np.linspace(max([wnt[0], k_start]), wnt[-1], 500)
    axt.plot(wnts, Batchelor_navadni(wnts, et, T, 10**td), color='k', linewidth=2)
    axt.plot(wnts, Batchelor_navadni(wnts, et, T, 10**td), color='C2', linewidth=1.5)
    es = (es1 + es2)/2
    tds = get_thermdiss(10**es1, [psdt[i] for i in range(len(psdt)) if wnt[i] >= min_wnT], np.array([wn for wn in wnt if wn >= min_wnT]))
    k_bs = (10**es / nu / Dt ** 2) ** 0.25 / (2 * np.pi)
    k_stars = 7.9*10**-3 * k_bs
    print(k_stars, k_start)
    wnts = np.linspace(max([wnt[0], k_stars]), wnt[-1], 500)
    axt.plot(wnts, Batchelor_navadni(wnts, es, T, tds), color='k', linewidth=2)
    axt.plot(wnts, Batchelor_navadni(wnts, es, T, tds), color='C6', linewidth=1.5)
    axt.set_xscale('log')
    axt.set_xticks([1, 10, 100])
    axt.set_xticklabels([1, 10, 100], fontsize=14)
    axt.set_yticklabels([10**(-i) for i in range(-12, 0)], fontsize=14)
    axt.set_xlabel(r'$\kappa$ [cpm]', fontsize=16)
    axt.set_yscale('log')
    axt.axvline(min_wnT, linestyle='--', linewidth=1, color='k')
    axt.axvline(100, linestyle='--', linewidth=1, color='k')
    axt.set_ylabel(r'PSD [(°C/m$)^2$/cpm]', fontsize=16)
    axt.annotate(r'$\log{(\epsilon_{T})}=$'+str(round(et, 2)), xy=(0.2, 0.35), xycoords='axes fraction', color='C2', fontsize=16)
    axt.annotate(r'$\log{(\chi_{T,T})}=$'+str(round(td, 2)), xy=(0.2, 0.25), xycoords='axes fraction', color='C2', fontsize=16)
    axt.annotate(r'$\log{(\epsilon_{s})}=$'+str(round(es, 2)), xy=(0.2, 0.15), xycoords='axes fraction', color='C6', fontsize=16)
    axt.annotate(r'$\log{(\chi_{T,s})}=$'+str(round(np.log10(tds), 2)), xy=(0.2, 0.05), xycoords='axes fraction', color='C6', fontsize=16)
    axsh.set_title(prg[:4] + ', ruta ' + prg[10:12] + ', ' + prg[13:15] + ' dbar', fontsize=18)
    if prg[14] == '.':
        axsh.set_title(prg[:4] + ', ruta ' + prg[10:12] + ', ' + prg[13:14] + ' dbar', fontsize=18)
    plt.tight_layout()
    if legende:
        axsh.legend()
        axt.legend()

count1 = 0
count2 = 0

for ioutput in range(len(outputs['epsilon_shear'])):
    prg = outputs['postaja_ruta_globina'][ioutput]
    if prg in combinations_to_plot:
        joutput = outputs_shear['postaja_ruta_globina'].index(prg)
        epsilon_grad_T12 = outputs['epsilon_grad_T_DoF12'][ioutput]
        epsilon_shear1 = fit_Nasmyth_MLE(psd=outputs_shear['useful_PSD1'][joutput], wavenumbers=outputs_shear['wavenumbers'][joutput], min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][joutput], DoF=12)
        epsilon_shear2 = fit_Nasmyth_MLE(psd=outputs_shear['useful_PSD2'][joutput], wavenumbers=outputs_shear['wavenumbers'][joutput], min_wavenumber=6, max_wavenumber=35, temp_mean=outputs_shear['tempcor'][joutput], DoF=12)

        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2  # MSSpro user manual
        nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual, nu is viscosity
        N2 = outputs['N2'][ioutput]
        rebT = 10**epsilon_grad_T12/(nu*N2)
        rebs1 = 10**epsilon_shear1/(nu*N2)
        rebs2 = 10**epsilon_shear2/(nu*N2)
        rebs = (10**epsilon_shear1 + 10**epsilon_shear2)/2/(nu*N2)
        if count1 < 3:
            plot_stuff(axs[0, count2], axs[1, count2], psd1=outputs_shear['useful_PSD1'][joutput], psd2=outputs_shear['useful_PSD2'][joutput], wnsh=outputs_shear['wavenumbers'][joutput], T=outputs['tempcor'][ioutput], es1=epsilon_shear1, es2=epsilon_shear2,\
                   et = epsilon_grad_T12, wnt=outputs['wavenumbers'][ioutput], psdt=outputs['useful_grad_T_PSD'][ioutput], td=outputs['thermdiss_DoF12'][ioutput], prg=prg, legende=False, min_wnT=outputs['min_fitted_wn'][ioutput])
        else:
            plot_stuff(axs[2, count2-3], axs[3, count2-3], psd1=outputs_shear['useful_PSD1'][joutput], psd2=outputs_shear['useful_PSD2'][joutput], wnsh=outputs_shear['wavenumbers'][joutput], T=outputs['tempcor'][ioutput], es1=epsilon_shear1, es2=epsilon_shear2,\
                   et = epsilon_grad_T12, wnt=outputs['wavenumbers'][ioutput], psdt=outputs['useful_grad_T_PSD'][ioutput], td=outputs['thermdiss_DoF12'][ioutput], prg=prg, legende=False, min_wnT=outputs['min_fitted_wn'][ioutput])
        count1 += 1
        count2 += 1


if savefigs:
    plt.savefig(mapa_za_shranjevanje_grafov + 'primeri_spektrov.jpg', dpi=300)
plt.show()