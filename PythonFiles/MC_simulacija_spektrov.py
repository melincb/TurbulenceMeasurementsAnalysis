'''Computes Monte Carlo simulations for our fitting techniques within our range of fit (Batchelor: 1-100cpm, Nasmyth
6-35cm, wavenumber step 1cpm (both ranges and wavenumber steps can be manually changed in the beginning of the code)).
The random noise in spectrum is distributed via chi_squared distribution with 6 and 12 degrees of freedom (also can be
manually changed).

The code has three sections:
1) Plotting examples of the generated spectra
2) Computing MC simulations for the Batchelor spectrum (or just plotting their results)
3) Computing MC simulations for the Batchelor spectrum (or just plotting their results)
'''
from fitting_procedures_for_Batchelor_curve import iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral, fit_Batchelorjeve_krivulje_na_razliko_logaritmov, Batchelor_Fit_least_squares_log10
from fitting_procedures_for_Nasmyth_curve import fit_Nasmyth_least_squares_log10, fit_Nasmyth_Lueck, fit_Nasmyth_MLE
import matplotlib.pyplot as plt
import numpy as np
import scipy.special

import pickle


mapa_za_shranjevanje_grafov =  'Saved_figures\ '[:-1]    # Where it saves the figures

plot_example = False    # To see the plot example, set this to True, else False
save_figures = False    # To save the plots, set this to True, else False

compute_Batchelor_simulations = False   # If True, the entire process of MC simulations for the Batchelor spectrum is
                                        # going to run (WARNING: it may last for more than 16 hours!),
                                        # if False, only the plot of presaved simulations will appear
save_Batchelor_results_in_pickle = False   # To save the results of Batchelor simulations in .p, set to True

compute_Nasmyth_simulations = False     # If True, the entire process of MC simulations for the Nasmyth spectrum is
                                        # going to run (WARNING: it may last for a few hours!),
                                        # if False, only the plot of presaved simulations will appear
save_Nasmyth_results_in_pickle = False  # To save the results of Nasmyth simulations in .p, set to True

English_labels = True   # Set to False for Slovenian labels

min_wn_in_sh_spectrum = 6   # minimum wavenumber in the shear fluctuation spectrum
max_wn_in_sh_spectrum = 35  # maximum wavenumber in the shear fluctuation spectrum
sh_wn_step = 1              # wavenumber step in the shear fluctuation spectrum
critical_wn_in_tg_spectrum = 100    # critical wavenumber in the temperature fluctuation gradient spectrum
tg_wn_step = 1              # wavenumber step in the temperature fluctuation gradient spectrum (the minimum wavenumber
                            # in this spectrum is equal to tg_wn_step)


# THE BASIC MANUAL SETTINGS END HERE


def temperature_gradient_spectrum(wn, eps, td, temp_mean):
    '''Generates temperature fluctuation gradient spectrum'''
    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual
    Dt = 1.4 * 10 ** -7  # m**2 s**-1 (MSSpro user manual)
    q = 2 * np.sqrt(3)
    k_b = (eps / nu / Dt ** 2) ** 0.25 / (2*np.pi)  # cpm
    y = wn / k_b * np.sqrt(q)
    y_star = np.sqrt(q) * 7.9 * 10**-3
    try:    # works if y is a float or integer or np.float64, or np.float32, or ...
        if y >= y_star:
            return td * np.sqrt(q) / Dt / k_b * y ** 2 * (np.exp(-y ** 2) / y - np.sqrt(np.pi) * (1 - scipy.special.erf(y)))
        else:
            return 0.3 * td * y**(1/3) / (eps**(1/4) * q**(1/6) * nu**(1/12) * Dt**(1/6))
    except: # works if y is a np.array or a list or ...
        return_me = []
        for yy in y:
            if yy >= y_star:
                return_me.append(td * np.sqrt(q) / Dt / k_b * yy ** 2 * (np.exp(-yy ** 2) / yy - np.sqrt(np.pi) * (1 - scipy.special.erf(yy))))
            else:
                return_me.append(0.3 * td * yy**(1/3) / (eps**(1/4) * q**(1/6) * nu**(1/12) * Dt**(1/6)))
        return np.array(return_me)


def Nasmyth_navadni(k, log10epsilon, temp_mean=False):
    if temp_mean:
        nu = 1.702747 - 0.05126103 * temp_mean + 0.0005918645 * temp_mean ** 2  # MSSpro user manual
        nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual
    else:
        nu = 10**-6

    epsilon = 10 ** log10epsilon
    x = k * (nu ** 3 / epsilon) ** 0.25
    return 8.05 * x ** (1 / 3) / (1 + (20.6 * x) ** 3.715) * (epsilon ** 3 / nu) ** 0.25




if plot_example:
    plt.figure(figsize=(12, 4))
    epsiloni_slika = [-9, -8, -7, -6, -5]   # log10(epsilon) in the plot
    kji = np.linspace(tg_wn_step, critical_wn_in_tg_spectrum, int(critical_wn_in_tg_spectrum/tg_wn_step))
    for iepsilon in range(len(epsiloni_slika)):
        plt.subplot(1, 5, iepsilon+1)
        dof = 6
        epsilon = epsiloni_slika[iepsilon]
        neokrnjen_spekter = temperature_gradient_spectrum(kji, 10**epsilon, td=10 ** (-6.0), temp_mean=20)
        zasumljen_spekter_12_temp = neokrnjen_spekter * np.random.chisquare(dof, size=len(kji)) / dof
        plt.plot(kji, neokrnjen_spekter, color='k')
        plt.scatter(kji, zasumljen_spekter_12_temp, s=10, alpha=0.5, label=r'$d={}$'.format(dof))
        plt.title(r'$\log{(\epsilon)}=$'+str(epsilon), fontsize=17) # + r', $d={}$'.format(dof)
        plt.xscale('log')
        plt.yscale('log')
        plt.xticks(ticks=[1, 10, 100], labels=['1', '10', '100'], fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel(r'$k$ [cpm]', fontsize=16)
        if iepsilon == 0:
            plt.ylabel(r'PSD [(°C/m$)^2$/cpm]', fontsize=16)

        dof = 12
        zasumljen_spekter_12_temp = neokrnjen_spekter * np.random.chisquare(dof, size=len(kji)) / dof
        plt.scatter(kji, zasumljen_spekter_12_temp, s=10, alpha=0.5, label=r'$d={}$'.format(dof))
        if iepsilon == 0:
            plt.legend(fontsize=16, handletextpad=-0.5)

    plt.tight_layout()
    if save_figures:
        plt.savefig(mapa_za_shranjevanje_grafov + 'MC_generated_spectra_examples_B.jpg', dpi=400)
    plt.show()

    plt.figure(figsize=(12, 4))
    epsiloni_slika = [-9, -8, -7, -6, -5]
    kji_shear = np.linspace(min_wn_in_sh_spectrum, max_wn_in_sh_spectrum, int((max_wn_in_sh_spectrum - min_wn_in_sh_spectrum)/sh_wn_step + 1))
    for iepsilon in range(len(epsiloni_slika)):
        plt.subplot(1, 5, iepsilon+1)
        dof = 6
        epsilon = epsiloni_slika[iepsilon]
        neokrnjen_spekter = Nasmyth_navadni(kji_shear, epsilon)
        zasumljen_spekter_12_temp = neokrnjen_spekter * np.random.chisquare(dof, size=len(kji_shear)) / dof
        plt.plot(kji_shear, neokrnjen_spekter, color='k')
        plt.scatter(kji_shear, zasumljen_spekter_12_temp, s=10, alpha=0.5, label=r'$d={}$'.format(dof))
        plt.title(r'$\log{(\epsilon)}=$'+str(epsilon), fontsize=17)
        plt.xscale('log')
        plt.yscale('log')
        plt.xticks(ticks=[6, 10, 20, 30], labels=['6', '10', '20', '30'], fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel(r'$k$ [cpm]', fontsize=16)
        if iepsilon == 0:
            plt.ylabel(r'PSD [$\mathrm{s}^{-2}$/cpm]', fontsize=16)

        dof = 12
        zasumljen_spekter_12_temp = neokrnjen_spekter * np.random.chisquare(dof, size=len(kji_shear)) / dof
        plt.scatter(kji_shear, zasumljen_spekter_12_temp, s=10, alpha=0.5, label=r'$d={}$'.format(dof))

        if iepsilon == 0:
            plt.legend(fontsize=16, handletextpad=-0.5)

    plt.tight_layout()
    if save_figures:
        plt.savefig(mapa_za_shranjevanje_grafov + 'MC_generated_spectra_examples_N.jpg', dpi=400)
    plt.show()


# BATCHELOR SPECTRUM

N_samples = 1000    # Number of samples of each log10_epsilon
kji = np.linspace(tg_wn_step, critical_wn_in_tg_spectrum, int(critical_wn_in_tg_spectrum / tg_wn_step))   # wavenumbers
log10_epsilon_bounds = [-9.5, -2]   # min and max log10_epsilon for the generation of spectra
log10_epsilon_step = 0.1    # the step between adjacent log10_epsilon for the generation of spectra
# Initiate log10_epsilons for the generation of spectra
log10_epsilons = np.arange(log10_epsilon_bounds[0], log10_epsilon_bounds[1]+log10_epsilon_step, step=log10_epsilon_step)

# output errors, sorted by the true value of log10_epsilon, 12 degrees of freedom
log10_12_true_intervals_errors = [[] for i in log10_epsilons]   # for Batchelor_Fit_least_squares_log10()
log10_minus_mean_12_true_intervals_errors = [[] for i in log10_epsilons] # for fit_Batchelorjeve_krivulje_na_razliko_logaritmov()
MLE_12_true_intervals_errors = [[] for i in log10_epsilons] # for iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral() with starting log10_epsilon from fit_Batchelorjeve_krivulje_na_razliko_logaritmov()
MLE_perfect_start_12_true_intervals_errors = [[] for i in log10_epsilons] # for iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral() with starting true log10_epsilon

# output errors, sorted by output log10_epsilon, 12 degrees of freedom
log10_12_results_intervals_errors = [[] for i in log10_epsilons]
log10_minus_mean_12_results_intervals_errors = [[] for i in log10_epsilons]
MLE_12_results_intervals_errors = [[] for i in log10_epsilons]
MLE_perfect_start_12_results_intervals_errors = [[] for i in log10_epsilons]

# output errors, sorted by the true value of log10_epsilon, 6 degrees of freedom
log10_6_true_intervals_errors = [[] for i in log10_epsilons]
log10_minus_mean_6_true_intervals_errors = [[] for i in log10_epsilons]
MLE_6_true_intervals_errors = [[] for i in log10_epsilons]
MLE_perfect_start_6_true_intervals_errors = [[] for i in log10_epsilons]

# output errors, sorted by output log10_epsilon, 6 degrees of freedom
log10_6_results_intervals_errors = [[] for i in log10_epsilons]
log10_minus_mean_6_results_intervals_errors = [[] for i in log10_epsilons]
MLE_6_results_intervals_errors = [[] for i in log10_epsilons]
MLE_perfect_start_6_results_intervals_errors = [[] for i in log10_epsilons]

if compute_Batchelor_simulations:
    for ilog10_epsilon in range(len(log10_epsilons)):   # Loop through all log10_epsilons
        trenutni_log10_epsilon = log10_epsilons[ilog10_epsilon] # current log10_epsilon
        print(trenutni_log10_epsilon)
        for jspekter in range(N_samples):   # Each time generate and fit N_samples
            temperatura = np.random.uniform(5,25)   # Draw temperature (unit degreeC)
            log10_thermdiss = np.random.uniform(-10,-4) # Draw decimal logarithm of thermal dissipation (unit log10(degreeC**2/s))

            # spectrum without random noise
            neokrnjen_spekter = temperature_gradient_spectrum(kji, eps=10**trenutni_log10_epsilon, td=10**log10_thermdiss, temp_mean=temperatura)
            # spectrum with random noise with 12 degrees of freedom
            zasumljen_spekter_12_temp = neokrnjen_spekter * np.random.chisquare(12, size=len(kji)) / 12
            # spectrum with random noise with 6 degrees of freedom
            zasumljen_spekter_6_temp = neokrnjen_spekter * np.random.chisquare(6, size=len(kji)) / 6


            # Batchelor_Fit_least_squares_log10, 12 degrees of freedom
            try:
                log10_12_fit = Batchelor_Fit_least_squares_log10(zasumljen_spekter_12_temp, kji, temp_mean=temperatura, critical_wavenumber=critical_wn_in_tg_spectrum)[0]
                log10_12_true_intervals_errors[ilog10_epsilon].append(log10_12_fit - trenutni_log10_epsilon)
                for jlog10_epsilon in range(len(log10_epsilons)):
                    if log10_12_fit >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step/2 and log10_12_fit < log10_epsilons[jlog10_epsilon] + log10_epsilon_step/2:
                        log10_12_results_intervals_errors[jlog10_epsilon].append(log10_12_fit - trenutni_log10_epsilon)
            except:
                # The spectrum might be detected as falling all the time (finestructure is falsely assumed to predominate in the entire spectrum)
                # Another option is that Batchelor_Fit_least_squares_log10 has a runtime error
                pass

            # Batchelor_Fit_least_squares_log10, 6 degrees of freedom
            try:
                log10_6_fit = Batchelor_Fit_least_squares_log10(zasumljen_spekter_6_temp, kji, temp_mean=temperatura, critical_wavenumber=critical_wn_in_tg_spectrum)[0]
                log10_6_true_intervals_errors[ilog10_epsilon].append(log10_6_fit - trenutni_log10_epsilon)
                for jlog10_epsilon in range(len(log10_epsilons)):
                    if log10_6_fit >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step/2 and log10_6_fit < log10_epsilons[jlog10_epsilon] + log10_epsilon_step/2:
                        log10_6_results_intervals_errors[jlog10_epsilon].append(log10_6_fit - trenutni_log10_epsilon)
            except:
                # The spectrum might be detected as falling all the time (finestructure is falsely assumed to predominate in the entire spectrum)
                # Another option is that Batchelor_Fit_least_squares_log10 has a runtime error
                pass

            # fit_Batchelorjeve_krivulje_na_razliko_logaritmov and iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral, 12 degrees of freedom
            try:
                log10_minus_mean_12_fit = fit_Batchelorjeve_krivulje_na_razliko_logaritmov(kji, zasumljen_spekter_12_temp, temp_mean=temperatura)[0]
                log10_minus_mean_12_true_intervals_errors[ilog10_epsilon].append(log10_minus_mean_12_fit - trenutni_log10_epsilon)
                for jlog10_epsilon in range(len(log10_epsilons)):
                    if log10_minus_mean_12_fit >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step/2 and log10_minus_mean_12_fit < log10_epsilons[jlog10_epsilon] + log10_epsilon_step/2:
                        log10_minus_mean_12_results_intervals_errors[jlog10_epsilon].append(log10_minus_mean_12_fit - trenutni_log10_epsilon)

                if jspekter < 100:  # Because iterative fit is very slow, we only do 100 samples
                    MLE_fit_12 = np.log10(iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(kji, zasumljen_spekter_12_temp, temp_mean=temperatura, starting_epsilon=10**log10_minus_mean_12_fit, DoF=12, critical_wavenumber=critical_wn_in_tg_spectrum)[0])
                    MLE_12_true_intervals_errors[ilog10_epsilon].append(MLE_fit_12 - trenutni_log10_epsilon)
                    MLE_perfect_start_fit_12 = np.log10(iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(kji, zasumljen_spekter_12_temp, temp_mean=temperatura, starting_epsilon=10**trenutni_log10_epsilon, DoF=12, critical_wavenumber=critical_wn_in_tg_spectrum)[0])
                    MLE_perfect_start_12_true_intervals_errors[ilog10_epsilon].append(MLE_perfect_start_fit_12 - trenutni_log10_epsilon)

                    for jlog10_epsilon in range(len(log10_epsilons)):
                        if MLE_fit_12 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and MLE_fit_12 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                            MLE_12_results_intervals_errors[jlog10_epsilon].append(MLE_fit_12 - trenutni_log10_epsilon)
                        if MLE_perfect_start_fit_12 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and MLE_perfect_start_fit_12 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                            MLE_perfect_start_12_results_intervals_errors[jlog10_epsilon].append(MLE_perfect_start_fit_12 - trenutni_log10_epsilon)
            except:
                # The spectrum might be detected as falling all the time (finestructure is falsely assumed to predominate in the entire spectrum)
                pass

            # fit_Batchelorjeve_krivulje_na_razliko_logaritmov and iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral, 6 degrees of freedom
            try:
                log10_minus_mean_6_fit = fit_Batchelorjeve_krivulje_na_razliko_logaritmov(kji, zasumljen_spekter_6_temp, temp_mean=temperatura)[0]
                log10_minus_mean_6_true_intervals_errors[ilog10_epsilon].append(log10_minus_mean_6_fit - trenutni_log10_epsilon)
                for jlog10_epsilon in range(len(log10_epsilons)):
                    if log10_minus_mean_6_fit >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step/2 and log10_minus_mean_6_fit < log10_epsilons[jlog10_epsilon] + log10_epsilon_step/2:
                        log10_minus_mean_6_results_intervals_errors[jlog10_epsilon].append(log10_minus_mean_6_fit - trenutni_log10_epsilon)

                if jspekter < 100:  # Because iterative fit is very slow, we only do 100 samples
                    MLE_fit_6 = np.log10(iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(kji, zasumljen_spekter_6_temp, temp_mean=temperatura, starting_epsilon=10**log10_minus_mean_6_fit, DoF=6)[0])
                    MLE_6_true_intervals_errors[ilog10_epsilon].append(MLE_fit_6 - trenutni_log10_epsilon)
                    MLE_perfect_start_fit_6 = np.log10(iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(kji, zasumljen_spekter_6_temp, temp_mean=temperatura, starting_epsilon=10**trenutni_log10_epsilon, DoF=6)[0])
                    MLE_perfect_start_6_true_intervals_errors[ilog10_epsilon].append(MLE_perfect_start_fit_6 - trenutni_log10_epsilon)

                    for jlog10_epsilon in range(len(log10_epsilons)):
                        if MLE_fit_6 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and MLE_fit_6 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                            MLE_6_results_intervals_errors[jlog10_epsilon].append(MLE_fit_6 - trenutni_log10_epsilon)
                        if MLE_perfect_start_fit_6 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and MLE_perfect_start_fit_6 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                            MLE_perfect_start_6_results_intervals_errors[jlog10_epsilon].append(MLE_perfect_start_fit_6 - trenutni_log10_epsilon)
            except:
                # The spectrum might be detected as falling all the time (finestructure is falsely assumed to predominate in the entire spectrum)
                pass

else:
    log10_minus_mean_12_results_intervals_errors, log10_minus_mean_6_results_intervals_errors, log10_minus_mean_12_true_intervals_errors, log10_minus_mean_6_true_intervals_errors= pickle.load(open('log10_minus_mean_MC.p', 'rb'))
    log10_12_results_intervals_errors, log10_6_results_intervals_errors, log10_12_true_intervals_errors, log10_6_true_intervals_errors = pickle.load(open('log10_MC.p', 'rb'))
    MLE_12_results_intervals_errors, MLE_6_results_intervals_errors, MLE_12_true_intervals_errors, MLE_6_true_intervals_errors = pickle.load(open('MLE_MC.p', 'rb'))
    MLE_perfect_start_12_results_intervals_errors, MLE_perfect_start_6_results_intervals_errors, MLE_perfect_start_12_true_intervals_errors, MLE_perfect_start_6_true_intervals_errors = pickle.load(open('MLE_perfect_start_MC.p', 'rb'))


# Plot the results

plt.figure(figsize=(12, 10))
plt.subplot(2,2,1)
if not English_labels:
    plt.plot(log10_epsilons, [np.mean(log10_12_true_intervals_errors[i]) for i in range(len(log10_12_true_intervals_errors))], label=r'$d=12$, LS', linewidth=2.5)
    plt.plot(log10_epsilons, [np.mean(log10_minus_mean_12_true_intervals_errors[i]) for i in range(len(log10_minus_mean_12_true_intervals_errors))], label=r'$d=12$, LSnTD')
    plt.plot(log10_epsilons, [np.mean(MLE_12_true_intervals_errors[i]) for i in range(len(MLE_12_true_intervals_errors))], label=r'$d=12$, MLE, $\epsilon_0 = \epsilon_{\mathtt{LSnTD}}$', linewidth=2.5, color='C2')
    plt.plot(log10_epsilons, [np.mean(MLE_perfect_start_12_true_intervals_errors[i]) for i in range(len(MLE_perfect_start_12_true_intervals_errors))], label=r'$d=12$, MLE, $\epsilon_0 = \epsilon_{\mathtt{pravi}}$', color='C6')
    plt.plot(log10_epsilons, [np.mean(log10_6_true_intervals_errors[i]) for i in range(len(log10_6_true_intervals_errors))], label=r'$d=6$, LS', linewidth=2.5, linestyle='--', color='C0')
    plt.plot(log10_epsilons, [np.mean(log10_minus_mean_6_true_intervals_errors[i]) for i in range(len(log10_minus_mean_6_true_intervals_errors))], label=r'$d=6$, LSnTD', linestyle='--', color='C1')
    plt.plot(log10_epsilons, [np.mean(MLE_6_true_intervals_errors[i]) for i in range(len(MLE_6_true_intervals_errors))], label=r'$d=6$, MLE, $\epsilon_0 = \epsilon_{\mathtt{LSnTD}}$', linestyle='--', linewidth=2.5, color='C2')
    plt.plot(log10_epsilons, [np.mean(MLE_perfect_start_6_true_intervals_errors[i]) for i in range(len(MLE_perfect_start_6_true_intervals_errors))], label=r'$d=6$, MLE, $\epsilon_0 = \epsilon_{\mathtt{pravi}}$', linestyle='--', color='C6')
else:
    plt.plot(log10_epsilons, [np.mean(log10_12_true_intervals_errors[i]) for i in range(len(log10_12_true_intervals_errors))], label=r'$d=12$, LS', linewidth=2.5)
    plt.plot(log10_epsilons, [np.mean(log10_minus_mean_12_true_intervals_errors[i]) for i in range(len(log10_minus_mean_12_true_intervals_errors))], label=r'$d=12$, LSnTD')
    plt.plot(log10_epsilons, [np.mean(MLE_12_true_intervals_errors[i]) for i in range(len(MLE_12_true_intervals_errors))], label=r'$d=12$, MLE, $\epsilon_0 = \epsilon_{\mathtt{LSnTD}}$', linewidth=2.5, color='C2')
    plt.plot(log10_epsilons, [np.mean(MLE_perfect_start_12_true_intervals_errors[i]) for i in range(len(MLE_perfect_start_12_true_intervals_errors))], label=r'$d=12$, MLE, $\epsilon_0 = \epsilon_{\mathtt{true}}$', color='C6')
    plt.plot(log10_epsilons, [np.mean(log10_6_true_intervals_errors[i]) for i in range(len(log10_6_true_intervals_errors))], label=r'$d=6$, LS', linewidth=2.5, linestyle='--', color='C0')
    plt.plot(log10_epsilons, [np.mean(log10_minus_mean_6_true_intervals_errors[i]) for i in range(len(log10_minus_mean_6_true_intervals_errors))], label=r'$d=6$, LSnTD', linestyle='--', color='C1')
    plt.plot(log10_epsilons, [np.mean(MLE_6_true_intervals_errors[i]) for i in range(len(MLE_6_true_intervals_errors))], label=r'$d=6$, MLE, $\epsilon_0 = \epsilon_{\mathtt{LSnTD}}$', linestyle='--', linewidth=2.5, color='C2')
    plt.plot(log10_epsilons, [np.mean(MLE_perfect_start_6_true_intervals_errors[i]) for i in range(len(MLE_perfect_start_6_true_intervals_errors))], label=r'$d=6$, MLE, $\epsilon_0 = \epsilon_{\mathtt{true}}$', linestyle='--', color='C6')

plt.grid(linestyle=':', linewidth=0.6)
plt.axhline(0, color='k')
plt.legend(fontsize=16, loc='lower left', handletextpad=0.5, labelspacing=0.1)
if not English_labels:
    plt.xlabel(r'Pravi $\log{(\epsilon)}$', fontsize=16)
    plt.ylabel('Povp. napaka dobljenega $\log{(\epsilon)}$', fontsize=16)
else:
    plt.xlabel(r'True value of $\log{(\epsilon)}$', fontsize=16) #
    plt.ylabel(r'Mean error of $\log{(\epsilon)}$ result', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(2,2,3)
plt.plot(log10_epsilons, [np.std(log10_12_true_intervals_errors[i]) for i in range(len(log10_12_true_intervals_errors))], linewidth=2.5)
plt.plot(log10_epsilons, [np.std(log10_minus_mean_12_true_intervals_errors[i]) for i in range(len(log10_minus_mean_12_true_intervals_errors))])
plt.plot(log10_epsilons, [np.std(MLE_12_true_intervals_errors[i]) for i in range(len(MLE_12_true_intervals_errors))], color='C2')
plt.plot(log10_epsilons, [np.std(MLE_perfect_start_12_true_intervals_errors[i]) for i in range(len(MLE_perfect_start_12_true_intervals_errors))], color='C6')
plt.plot(log10_epsilons, [np.std(log10_6_true_intervals_errors[i]) for i in range(len(log10_6_true_intervals_errors))], linewidth=2.5, linestyle='--', color='C0')
plt.plot(log10_epsilons, [np.std(log10_minus_mean_6_true_intervals_errors[i]) for i in range(len(log10_minus_mean_6_true_intervals_errors))], linestyle='--', color='C1')
plt.plot(log10_epsilons, [np.std(MLE_6_true_intervals_errors[i]) for i in range(len(MLE_6_true_intervals_errors))], linestyle='--', linewidth=2.5, color='C2')
plt.plot(log10_epsilons, [np.std(MLE_perfect_start_6_true_intervals_errors[i]) for i in range(len(MLE_perfect_start_6_true_intervals_errors))], linestyle='--', color='C6')
plt.grid(linestyle=':', linewidth=0.6)
plt.ylim(bottom=0)
if not English_labels:
    plt.xlabel(r'Pravi $\log{(\epsilon)}$', fontsize=16) #r'True value of $\log{(\epsilon)}$'
    plt.ylabel(r'STD napake dobljenega $\log{(\epsilon)}$', fontsize=16) #r'Standard deviation of error of fit'
else:
    plt.xlabel(r'True value of $\log{(\epsilon)}$', fontsize=16) #
    plt.ylabel(r'STD of error of $\log{(\epsilon)}$ result', fontsize=16) #
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(2,2,2)
plt.plot(log10_epsilons, [np.mean(log10_12_results_intervals_errors[i]) for i in range(len(log10_12_results_intervals_errors))], linewidth=2.5)
plt.plot(log10_epsilons, [np.mean(log10_minus_mean_12_results_intervals_errors[i]) for i in range(len(log10_minus_mean_12_results_intervals_errors))])
plt.plot(log10_epsilons, [np.mean(MLE_12_results_intervals_errors[i]) for i in range(len(MLE_12_results_intervals_errors))])
plt.plot(log10_epsilons, [np.mean(MLE_perfect_start_12_results_intervals_errors[i]) for i in range(len(MLE_perfect_start_12_results_intervals_errors))], color='C6')
plt.plot(log10_epsilons, [np.mean(log10_6_results_intervals_errors[i]) for i in range(len(log10_6_results_intervals_errors))], linewidth=2.5, linestyle='--', color='C0')
plt.plot(log10_epsilons, [np.mean(log10_minus_mean_6_results_intervals_errors[i]) for i in range(len(log10_minus_mean_6_results_intervals_errors))], linestyle='--', color='C1')
plt.plot(log10_epsilons, [np.mean(MLE_6_results_intervals_errors[i]) for i in range(len(MLE_6_results_intervals_errors))], linewidth=2.5, color='C2')
plt.plot(log10_epsilons, [np.mean(MLE_perfect_start_6_results_intervals_errors[i]) for i in range(len(MLE_perfect_start_6_results_intervals_errors))], linestyle='--', color='C6')
plt.grid(linestyle=':', linewidth=0.6)
plt.axhline(color='k')
if not English_labels:
    plt.xlabel(r'Dobljeni $\log{(\epsilon)}$', fontsize=16) # r'$\log{(\epsilon)}$ result of fit'
    plt.ylabel(r'Povp. napaka dobljenega $\log{(\epsilon)}$', fontsize=16) # r'Mean error of fit'
else:
    plt.xlabel(r'$\log{(\epsilon)}$ result of fit', fontsize=16) #
    plt.ylabel(r'Mean error of $\log{(\epsilon)}$ result', fontsize=16) # r'Mean error of fit'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(2,2,4)
plt.plot(log10_epsilons, [np.std(log10_12_results_intervals_errors[i]) for i in range(len(log10_12_results_intervals_errors))], linewidth=2.5)
plt.plot(log10_epsilons, [np.std(log10_minus_mean_12_results_intervals_errors[i]) for i in range(len(log10_minus_mean_12_results_intervals_errors))])
plt.plot(log10_epsilons, [np.std(MLE_12_results_intervals_errors[i]) for i in range(len(MLE_12_results_intervals_errors))], linewidth=2.5, color='C2')
plt.plot(log10_epsilons, [np.std(MLE_perfect_start_12_results_intervals_errors[i]) for i in range(len(MLE_perfect_start_12_results_intervals_errors))], color='C6')
plt.plot(log10_epsilons, [np.std(log10_6_results_intervals_errors[i]) for i in range(len(log10_6_results_intervals_errors))], linewidth=2.5, linestyle='--', color='C0')
plt.plot(log10_epsilons, [np.std(log10_minus_mean_6_results_intervals_errors[i]) for i in range(len(log10_minus_mean_6_results_intervals_errors))], linestyle='--', color='C1')
plt.plot(log10_epsilons, [np.std(MLE_6_results_intervals_errors[i]) for i in range(len(MLE_6_results_intervals_errors))], linestyle='--', linewidth=2.5, color='C2')
plt.plot(log10_epsilons, [np.std(MLE_perfect_start_6_results_intervals_errors[i]) for i in range(len(MLE_perfect_start_6_results_intervals_errors))], linestyle='--', color='C6')
plt.grid(linestyle=':', linewidth=0.6)
plt.ylim(bottom=0)
if not English_labels:
    plt.xlabel(r'Dobljeni $\log{(\epsilon)}$', fontsize=16) #r'$\log{(\epsilon)}$ result of fit'
    plt.ylabel(r'STD napake dobljenega $\log{(\epsilon)}$', fontsize=16) #r'Standard deviation of error of fit'
else:
    plt.xlabel(r'$\log{(\epsilon)}$ result of fit', fontsize=16) #
    plt.ylabel(r'STD of error of $\log{(\epsilon)}$ result', fontsize=16) # r'Standard deviation of error of fit'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

if not English_labels:
    plt.suptitle(r'Prileganje Batchelorjeve krivulje generiranim podatkom', fontsize=20)
else:
    plt.suptitle(r'Batchelor fit results for generated samples', fontsize=20)

plt.tight_layout()
if save_figures:
    plt.savefig(mapa_za_shranjevanje_grafov + 'Batchelor_fit_success_comparison_DoF=12_and_6_N_samples=1000.jpg', dpi=300)# + str(dof))
if save_Batchelor_results_in_pickle:
    pickle.dump([log10_minus_mean_12_results_intervals_errors, log10_minus_mean_6_results_intervals_errors, log10_minus_mean_12_true_intervals_errors, log10_6_true_intervals_errors], open('log10_minus_mean_MC.p', 'wb'))
    pickle.dump([log10_12_results_intervals_errors, log10_6_results_intervals_errors, log10_12_true_intervals_errors, log10_6_true_intervals_errors], open('log10_MC.p', 'wb'))
    pickle.dump([MLE_12_results_intervals_errors, MLE_6_results_intervals_errors, MLE_12_true_intervals_errors, MLE_6_true_intervals_errors], open('MLE_MC.p', 'wb'))
    pickle.dump([MLE_perfect_start_12_results_intervals_errors, MLE_perfect_start_6_results_intervals_errors, MLE_perfect_start_12_true_intervals_errors, MLE_perfect_start_6_true_intervals_errors], open('MLE_perfect_start_MC.p', 'wb'))
plt.show()
plt.close()



# NASMYTH SPECTRUM

N_samples = 1000    # number of samples
kji = np.linspace(min_wn_in_sh_spectrum, max_wn_in_sh_spectrum, int((max_wn_in_sh_spectrum - min_wn_in_sh_spectrum)/sh_wn_step + 1))    # wavenumbers (in cpm)
log10_epsilon_bounds = [-9.5, -2]   # min and max log10_epsilon for the generation of spectra
log10_epsilon_step = 0.1    # the step between adjacent log10_epsilon for the generation of spectra
# Initiate log10_epsilons for the generation of spectra
log10_epsilons = np.arange(log10_epsilon_bounds[0], log10_epsilon_bounds[1]+log10_epsilon_step, step=log10_epsilon_step)

# output errors, sorted by the true value of log10_epsilon, 12 degrees of freedom
Nlog10_dof12_true_intervals_errors = [[] for i in log10_epsilons]   # for fit_Nasmyth_least_squares_log10()
NL_dof12_true_intervals_errors = [[] for i in log10_epsilons]   # for fit_Nasmyth_Lueck()
NMLE_dof12_true_intervals_errors = [[] for i in log10_epsilons] # for fit_Nasmyth_MLE()

# output errors, sorted by the true value of log10_epsilon, 6degrees of freedom
Nlog10_dof6_true_intervals_errors = [[] for i in log10_epsilons]
NL_dof6_true_intervals_errors = [[] for i in log10_epsilons]
NMLE_dof6_true_intervals_errors = [[] for i in log10_epsilons]

# output errors, sorted by output log10_epsilon, 12 degrees of freedom
Nlog10_dof12_results_intervals_errors = [[] for i in log10_epsilons]
NL_dof12_results_intervals_errors = [[] for i in log10_epsilons]
NMLE_dof12_results_intervals_errors = [[] for i in log10_epsilons]

# output errors, sorted by output log10_epsilon, 6 degrees of freedom
Nlog10_dof6_results_intervals_errors = [[] for i in log10_epsilons]
NL_dof6_results_intervals_errors = [[] for i in log10_epsilons]
NMLE_dof6_results_intervals_errors = [[] for i in log10_epsilons]

if compute_Nasmyth_simulations:
    for ilog10_epsilon in range(len(log10_epsilons)):   # Loop through all log10_epsilons
        trenutni_log10_epsilon = log10_epsilons[ilog10_epsilon] # current log10_epsilon
        print(trenutni_log10_epsilon)
        for jspekter in range(N_samples):   # Do it N_samples times
            temperatura = np.random.uniform(5,25)   # Draw temperature (unit degreeC

            # spectrum without random noise
            neokrnjen_spekter = Nasmyth_navadni(kji, trenutni_log10_epsilon, temp_mean=temperatura)
            # spectrum with random noise with 6 degrees of freedom
            zasumljen_spekter_6 = neokrnjen_spekter * np.random.chisquare(6, size=len(kji)) / 6
            # spectrum with random noise with 12 degrees of freedom
            zasumljen_spekter_12 = neokrnjen_spekter * np.random.chisquare(12, size=len(kji)) / 12

            fit_Nlog10_6 = fit_Nasmyth_least_squares_log10(zasumljen_spekter_6, kji, min_wavenumber=min_wn_in_sh_spectrum, max_wavenumber=max_wn_in_sh_spectrum, temp_mean=temperatura)[0]
            fit_Nlog10_12 = fit_Nasmyth_least_squares_log10(zasumljen_spekter_12, kji, min_wavenumber=min_wn_in_sh_spectrum, max_wavenumber=max_wn_in_sh_spectrum, temp_mean=temperatura)[0]
            fit_NL_6 = fit_Nasmyth_Lueck(psd=zasumljen_spekter_6, wavenumbers=kji, min_wavenumber=min_wn_in_sh_spectrum, max_wavenumber=max_wn_in_sh_spectrum, temp_mean=temperatura)
            fit_NL_12 = fit_Nasmyth_Lueck(psd=zasumljen_spekter_12, wavenumbers=kji, min_wavenumber=min_wn_in_sh_spectrum, max_wavenumber=max_wn_in_sh_spectrum, temp_mean=temperatura)
            fit_NMLE_6 = fit_Nasmyth_MLE(psd=zasumljen_spekter_6, wavenumbers=kji, min_wavenumber=min_wn_in_sh_spectrum, max_wavenumber=max_wn_in_sh_spectrum, temp_mean=temperatura, DoF=6)
            fit_NMLE_12 = fit_Nasmyth_MLE(psd=zasumljen_spekter_12, wavenumbers=kji, min_wavenumber=min_wn_in_sh_spectrum, max_wavenumber=max_wn_in_sh_spectrum, temp_mean=temperatura, DoF=12)


            Nlog10_dof12_true_intervals_errors[ilog10_epsilon].append(fit_Nlog10_12 - trenutni_log10_epsilon)
            Nlog10_dof6_true_intervals_errors[ilog10_epsilon].append(fit_Nlog10_6 - trenutni_log10_epsilon)
            NL_dof12_true_intervals_errors[ilog10_epsilon].append(fit_NL_12 - trenutni_log10_epsilon)
            NL_dof6_true_intervals_errors[ilog10_epsilon].append(fit_NL_6 - trenutni_log10_epsilon)
            NMLE_dof6_true_intervals_errors[ilog10_epsilon].append(fit_NMLE_6 - trenutni_log10_epsilon)
            NMLE_dof12_true_intervals_errors[ilog10_epsilon].append(fit_NMLE_12 - trenutni_log10_epsilon)

            for jlog10_epsilon in range(len(log10_epsilons)):
                if fit_Nlog10_6 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and fit_Nlog10_6 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                    Nlog10_dof6_results_intervals_errors[jlog10_epsilon].append(fit_Nlog10_6 - trenutni_log10_epsilon)
                if fit_Nlog10_12 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and fit_Nlog10_12 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                    Nlog10_dof12_results_intervals_errors[jlog10_epsilon].append(fit_Nlog10_12 - trenutni_log10_epsilon)
                if fit_NL_6 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and fit_NL_6 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                    NL_dof6_results_intervals_errors[jlog10_epsilon].append(fit_NL_6 - trenutni_log10_epsilon)
                if fit_NL_12 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and fit_NL_12 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                    NL_dof12_results_intervals_errors[jlog10_epsilon].append(fit_NL_12 - trenutni_log10_epsilon)
                if fit_NMLE_6 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and fit_NMLE_6 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                    NMLE_dof6_results_intervals_errors[jlog10_epsilon].append(fit_NMLE_6 - trenutni_log10_epsilon)
                if fit_NMLE_12 >= log10_epsilons[jlog10_epsilon] - log10_epsilon_step / 2 and fit_NMLE_12 < log10_epsilons[jlog10_epsilon] + log10_epsilon_step / 2:
                    NMLE_dof12_results_intervals_errors[jlog10_epsilon].append(fit_NMLE_12 - trenutni_log10_epsilon)

else:
    Nlog10_dof12_results_intervals_errors, Nlog10_dof6_results_intervals_errors, Nlog10_dof12_true_intervals_errors, Nlog10_dof6_true_intervals_errors = pickle.load(open('Nlog10_MC.p', 'rb'))
    NL_dof12_results_intervals_errors, NL_dof6_results_intervals_errors, NL_dof12_true_intervals_errors, NL_dof6_true_intervals_errors = pickle.load(open('NL_MC.p', 'rb'))
    NMLE_dof12_results_intervals_errors, NMLE_dof6_results_intervals_errors, NMLE_dof12_true_intervals_errors, NMLE_dof6_true_intervals_errors = pickle.load(open('NMLE_MC.p', 'rb'))

# Plot

plt.figure(figsize=(12, 10))
plt.subplot(2,2,1)
plt.plot(log10_epsilons, [np.mean(Nlog10_dof12_true_intervals_errors[i]) for i in range(len(Nlog10_dof12_true_intervals_errors))], label='$d = 12$, LS', color='C0')
plt.plot(log10_epsilons, [np.mean(NL_dof12_true_intervals_errors[i]) for i in range(len(NL_dof12_true_intervals_errors))], label='$d=12$, IM', color='C1')
plt.plot(log10_epsilons, [np.mean(NMLE_dof12_true_intervals_errors[i]) for i in range(len(NMLE_dof12_true_intervals_errors))], label='$d=12$, MLE', color='C2')
plt.plot(log10_epsilons, [np.mean(Nlog10_dof6_true_intervals_errors[i]) for i in range(len(Nlog10_dof6_true_intervals_errors))], label='$d=6$, LS', color='C0', linestyle='--')
plt.plot(log10_epsilons, [np.mean(NL_dof6_true_intervals_errors[i]) for i in range(len(NL_dof6_true_intervals_errors))], label='$d=6$, IM', color='C1', linestyle='--')
plt.plot(log10_epsilons, [np.mean(NMLE_dof6_true_intervals_errors[i]) for i in range(len(NMLE_dof6_true_intervals_errors))], label='$d=6$, MLE', color='C2', linestyle='--')
plt.grid(linestyle=':', linewidth=0.6)
plt.axhline(0, color='k')
plt.legend(fontsize=16, loc='lower left', handletextpad=0.5, labelspacing=0.1)
if not English_labels:
    plt.xlabel(r'Pravi $\log{(\epsilon)}$', fontsize=16)    #r'True value of $\log{(\epsilon)}$'
    plt.ylabel(r'Povp. napaka dobljenega $\log{(\epsilon)}$', fontsize=16)    #r'Mean error of fit'
else:
    plt.xlabel(r'True value of $\log{(\epsilon)}$', fontsize=16)    #r'True value of $\log{(\epsilon)}$'
    plt.ylabel(r'Mean error of $\log{(\epsilon)}$ result', fontsize=16)    #r'Mean error of fit'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(2,2,3)
plt.plot(log10_epsilons, [np.std(Nlog10_dof12_true_intervals_errors[i]) for i in range(len(Nlog10_dof12_true_intervals_errors))], color='C0')
plt.plot(log10_epsilons, [np.std(NL_dof12_true_intervals_errors[i]) for i in range(len(NL_dof12_true_intervals_errors))], color='C1')
plt.plot(log10_epsilons, [np.std(NMLE_dof12_true_intervals_errors[i]) for i in range(len(NMLE_dof12_true_intervals_errors))], color='C2')
plt.plot(log10_epsilons, [np.std(Nlog10_dof6_true_intervals_errors[i]) for i in range(len(Nlog10_dof6_true_intervals_errors))], color='C0', linestyle='--')
plt.plot(log10_epsilons, [np.std(NL_dof6_true_intervals_errors[i]) for i in range(len(NL_dof6_true_intervals_errors))], linestyle='--')
plt.plot(log10_epsilons, [np.std(NMLE_dof6_true_intervals_errors[i]) for i in range(len(NMLE_dof6_true_intervals_errors))], linestyle='--')
plt.grid(linestyle=':', linewidth=0.6)
plt.ylim(bottom=0)
if not English_labels:
    plt.xlabel(r'Pravi $\log{(\epsilon)}$', fontsize=16) #r'True value of $\log{(\epsilon)}$'
    plt.ylabel(r'STD napake dobljenega $\log{(\epsilon)}$', fontsize=16) #r'Standard deviation of error of fit'
else:
    plt.xlabel(r'True value of $\log{(\epsilon)}$', fontsize=16) #r'True value of $\log{(\epsilon)}$'
    plt.ylabel(r'STD of $\log{(\epsilon)}$ result', fontsize=16) #r'Standard deviation of error of fit'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(2,2,2)
plt.plot(log10_epsilons, [np.mean(Nlog10_dof12_results_intervals_errors[i]) for i in range(len(Nlog10_dof12_results_intervals_errors))], color='C0')
plt.plot(log10_epsilons, [np.mean(NL_dof12_results_intervals_errors[i]) for i in range(len(NL_dof12_results_intervals_errors))], color='C1')
plt.plot(log10_epsilons, [np.mean(NMLE_dof12_results_intervals_errors[i]) for i in range(len(NMLE_dof12_results_intervals_errors))], color='C2')
plt.plot(log10_epsilons, [np.mean(Nlog10_dof6_results_intervals_errors[i]) for i in range(len(Nlog10_dof6_results_intervals_errors))], color='C0', linestyle='--')
plt.plot(log10_epsilons, [np.mean(NL_dof6_results_intervals_errors[i]) for i in range(len(NL_dof6_results_intervals_errors))], color='C1', linestyle='--')
plt.plot(log10_epsilons, [np.mean(NMLE_dof6_results_intervals_errors[i]) for i in range(len(NMLE_dof6_results_intervals_errors))], color='C2', linestyle='--')
plt.grid(linestyle=':', linewidth=0.6)
plt.axhline(color='k')
if not English_labels:
    plt.xlabel(r'Dobljeni $\log{(\epsilon)}$', fontsize=16) #r'$\log{(\epsilon)}$ result of fit'
    plt.ylabel(r'Povp. napaka dobljenega $\log{(\epsilon)}$', fontsize=16) #r'Mean error of fit'
else:
    plt.xlabel(r'$\log{(\epsilon)}$ result of fit', fontsize=16) #r'$\log{(\epsilon)}$ result of fit'
    plt.ylabel(r'Mean error of $\log{(\epsilon)}$ result', fontsize=16) #r'Mean error of fit'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.subplot(2,2,4)
plt.plot(log10_epsilons, [np.std(Nlog10_dof12_results_intervals_errors[i]) for i in range(len(Nlog10_dof12_results_intervals_errors))], color='C0')
plt.plot(log10_epsilons, [np.std(NL_dof12_results_intervals_errors[i]) for i in range(len(NL_dof12_results_intervals_errors))], color='C1')
plt.plot(log10_epsilons, [np.std(NMLE_dof12_results_intervals_errors[i]) for i in range(len(NMLE_dof12_results_intervals_errors))], color='C2')
plt.plot(log10_epsilons, [np.std(Nlog10_dof6_results_intervals_errors[i]) for i in range(len(Nlog10_dof6_results_intervals_errors))], color='C0', linestyle='--')
plt.plot(log10_epsilons, [np.std(NL_dof6_results_intervals_errors[i]) for i in range(len(NL_dof6_results_intervals_errors))], color='C1', linestyle='--')
plt.plot(log10_epsilons, [np.std(NMLE_dof6_results_intervals_errors[i]) for i in range(len(NMLE_dof6_results_intervals_errors))], color='C2', linestyle='--')
plt.grid(linestyle=':', linewidth=0.6)
plt.ylim(bottom=0)
if not English_labels:
    plt.xlabel(r'Dobljeni $\log{(\epsilon)}$', fontsize=16) #r'$\log{(\epsilon)}$ result of fit'
    plt.ylabel(r'STD napake dobljenega $\log{(\epsilon)}$', fontsize=16) #r'Standard deviation of error of fit'
else:
    plt.xlabel(r'$\log{(\epsilon)}$ result of fit', fontsize=16) #r'$\log{(\epsilon)}$ result of fit'
    plt.ylabel(r'STD of error of $\log{(\epsilon)}$ result', fontsize=16) #r'Standard deviation of error of fit'
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

if not English_labels:
    plt.suptitle(r'Prileganje Nasmythove krivulje generiranim podatkom', fontsize=20)#r'Nasmyth fit results for generated samples'
else:
    plt.suptitle(r'Nasmyth fit results for generated samples', fontsize=20)
plt.tight_layout()
if save_figures:
    plt.savefig(mapa_za_shranjevanje_grafov + 'Nasmyth_fit_success_comparison_DoF=12_and_6_N_samples=1000.jpg', dpi=300)# + str(dof))
if save_Nasmyth_results_in_pickle:
    pickle.dump([Nlog10_dof12_results_intervals_errors, Nlog10_dof6_results_intervals_errors, Nlog10_dof12_true_intervals_errors, Nlog10_dof6_true_intervals_errors], open('Nlog10_MC.p', 'wb'))
    pickle.dump([NL_dof12_results_intervals_errors, NL_dof6_results_intervals_errors, NL_dof12_true_intervals_errors, NL_dof6_true_intervals_errors], open('NL_MC.p', 'wb'))
    pickle.dump([NMLE_dof12_results_intervals_errors, NMLE_dof6_results_intervals_errors, NMLE_dof12_true_intervals_errors, NMLE_dof6_true_intervals_errors], open('NMLE_MC.p', 'wb'))
plt.show()
