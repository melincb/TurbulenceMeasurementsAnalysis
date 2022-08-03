def iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral(wavenumbers, psd, temp_mean, starting_epsilon, DoF, n_iter=100, critical_wavenumber=100):
    ''' Fit to Batchelor curve using integration to determine thermal dissipation and Maximum Likelihood Estimate to
    calcultate epsilon. Acronym in MC_simulacija_spektrov.py: MLE.

    Inputs: (gTo = "grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m" + datumska_mapa + ".p")
        wavenumbers ... wavenumbers that correspond to values of PSD ("wavenumbers" in gTo, unit cpm),
        psd ... PSD from FFT of vertical derivatives of detrended temperature ("useful_grad_T_PSD" in gTo, unit °C**2*m**(-2)/cpm),
        temp_mean ... mean temperature ("tempcor" in gTo, unit °C),
        starting_epsilon ... epsilon with which we start our iteration (unit W/kg),
        DoF ... degrees of freedom in our PSD by which we assume the behaviour of random noise in PSD,
        n_iter ... maximum number of iterations,
        critical_wavenumber ... maximum wavenumber for fit (unit cpm)

    Outputs:
        current_epsilon ... the best estimate for epsilon (unit W/kg) (the resolution of log10(current_epsilon) is 0.002),
        current_thermdiss ... the best estimate for thermal dissipation (unit °C**2/s),
        final_wavenumbers[0] ... the lowest wavenumber at which the PSD is rising (unit cpm)
    '''

    from scipy.optimize import curve_fit
    import scipy.integrate
    import numpy as np
    from scipy.stats import chi2
    from scipy.special import erf

    # First let's calculate viscosity
    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual
    # Determine thermal diffusivity and factor q
    Dt = 1.4 * 10**-7   # m**2 s**-1 (MSSpro user manual)
    q = 2 * np.sqrt(3)

    # Take only useful PSD (i.e. wavenumbers have to be larger than 0 and lower or equal to critical wavenumber)
    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= critical_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])

    # Due to finestructure we may have descending values at the beginning of the spectrum - let's cut them out
    wn_subinterval_width = 5    # width of subinterval in which we want our PSD to be rising
    min_len_wn_for_fit = 30     # minimum width of the PSD to still make a fit
    def lin_fit(k, a, b):
        return a*k + b

    final_wavenumbers = []
    final_psd = []
    for iwn in range(len(proper_wavenumbers) - min_len_wn_for_fit):
        popt, pcov = curve_fit(lin_fit, proper_wavenumbers[iwn:iwn+wn_subinterval_width], proper_psd[iwn:iwn+wn_subinterval_width])
        if popt[0] - 2*np.sqrt(pcov[0][0]) > 0:   # a - 2*std(a) > 0 => surely we are rising!
            final_wavenumbers = proper_wavenumbers[iwn:]
            final_psd = proper_psd[iwn:]
            break

    if len(final_wavenumbers) == 0:
        print('iterative, PSD is never rising in the dedicated area!')
        return None


    def get_thermdiss(current_epsilon):
        '''Calculates thermal dissipation by integrating the PSD and assuming that the true value of epsilon is
        current_epsilon.

        Input:
            current_epsilon ... the assumed true value of epsilon (unit W/kg)

        Output:
            the estimate for thermal dissipation from integration (unit °C**2/s)
        '''
        integral = scipy.integrate.trapz(y=final_psd, x=final_wavenumbers)  # Integrate using the trapezoidal rule
        measured_part = 6 * Dt * integral   # From relation "ThermalDissipation = 6*ThermalDiffusivity*IntegralOfPSD

        y_star = np.sqrt(q) * 7.9 * 10**-3  # Where the spectrum goes from inertial-convective subrange to
                                            # viscous-convective subrange
        y_min = 2 * np.pi * wavenumbers[0] * np.sqrt(q) * (nu * Dt ** 2 / current_epsilon)**(1/4)
        y_max = 2 * np.pi * wavenumbers[-1] * np.sqrt(q) * (nu * Dt ** 2 / current_epsilon)**(1/4)

        # *************************************************************************************************************
        # The following lines contain results from analytical integrals from Wolfram Mathematica, by which we determine
        # the integral of the theoretical curves with current_epsilon between 0 and infinity and between our minimum and
        # maximum wavenumber. Both integrals linearly depend on thermal dissipation, so we divide them by thermal
        # dissipation, which won't affect their fraction.
        # *************************************************************************************************************
        theory_0_to_inf = 9 * y_star**(4/3) / (40 * Dt**(1/6) * current_epsilon**(1/4) * nu**(1/12) * q**(1/6)) + np.exp(-y_star**2) * np.sqrt(q) * (1 - 2*y_star**2 + np.exp(y_star**2) * y_star**2 * (-3 + 2*np.sqrt(np.pi) * y_star * scipy.special.erfc(y_star))) / (6 * Dt * (current_epsilon / Dt**2 / nu)**(1/4))
        if y_star <= y_min:
            theory_min_to_max = 1 / (6 * Dt * (current_epsilon / (Dt**2 * nu))**(1/4)) * np.exp(-y_max**2 - y_min**2) * np.sqrt(q) * (np.exp(y_max**2) * (1 - 2*y_min**2) + np.exp(y_min**2)*(-1 + 2*y_max**2 - 2* np.exp(y_max**2) * np.sqrt(np.pi) * (y_max**3 * scipy.special.erfc(y_max) - y_min**3 * scipy.special.erfc(y_min))))
        elif y_star < y_max:
            theory_min_to_max = - 9 * (y_min**(4/3) - y_star**(4/3)) / (40 * Dt**(1/6) * current_epsilon**(1/4) * nu**(1/12) * q**(1/6)) \
                                + 1 / (6 * Dt * (current_epsilon / (Dt**2 * nu))**(1/4)) * np.exp(-y_max**2 - y_star**2) * np.sqrt(q) * (np.exp(y_max**2) * (1 - 2*y_star**2) + np.exp(y_star**2)*(-1 + 2*y_max**2 - 2* np.exp(y_max**2) * np.sqrt(np.pi) * (y_max**3 * scipy.special.erfc(y_max) - y_star**3 * scipy.special.erfc(y_star))))
        else:
            theory_min_to_max = 9 * (y_max**(4/3) - y_min**(4/3)) / (40 * Dt**(1/6) * current_epsilon**(1/4) * nu**(1/12) * q**(1/6))
        # *************************************************************************************************************
        # End of theoretical integration
        # *************************************************************************************************************

        theoretical_fraction = theory_min_to_max / theory_0_to_inf  # Theoretical fraction of what we see vs. what we
                                                                    # would see, if we could use infinite wavenumbers

        return measured_part / theoretical_fraction # total thermal dissipation = (what we see) / (the fraction of what we see)

    def get_epsilon(current_thermdiss, DoF=DoF):
        '''Calculates epsilon using Maximum Likelihood estimate and an assumption that the true value of thermal
        dissipation is current_thermdiss.

        Inputs:
            current_thermdiss ... the assumed true value of thermal dissipation (unit °C**2/s),
            DoF ... degrees of freedom in our PSD by which we assume the behaviour of random noise in PSD

        Output:
            the best estimate of epsilon from Maximum Likelihood estimate (unit W/kg)
        '''
        def best_fit_of_trial_log10_epsilons(trial_log10_epsilons):
            '''Maximum likelihood estimate following Ruddick et al: Maximum Likelihood Spectral Fitting: The Batchelor
            Spectrum, JAOT, 2000, https://doi.org/10.1175/1520-0426(2000)017<1541:MLSFTB>2.0.CO;2

            Ruddick et all also use more complicated instrumental noise parameterization (in our case noise is
            simplified to be a Heavyside function which is 0 for wavenumbers <= 100 cpm and infinity for
            wavenumbers > 100 cpm). To change the instrumental noise to some more complicated function, refer to
            variable noise_spectrum.
            '''
            noise_spectrum = 0  # simplified case
            C11s = []   # notation of cost function from Ruddick et al
            for log10_epsilon in trial_log10_epsilons:
                epsilon = 10**log10_epsilon
                k_b = (epsilon / nu / Dt ** 2) ** 0.25 / (2*np.pi)  # Batchelor wavenumber in cpm
                k_star = 7.9 * 10**-3 * k_b # Where the spectrum goes from inertial-convective subrange to
                                            # viscous-convective subrange

                # Let's only fit at wavenumbers where the spectrum is in viscous subrange
                true_final_wavenumbers = np.array([wn for wn in final_wavenumbers if wn >= k_star])
                true_final_psd = np.array([final_psd[iwn] for iwn in range(len(final_wavenumbers)) if final_wavenumbers[iwn] in true_final_wavenumbers])
                y = true_final_wavenumbers / k_b * np.sqrt(q)
                SB = current_thermdiss * np.sqrt(q) / Dt / k_b * y ** 2 * (np.exp(-y ** 2) / y - np.sqrt(np.pi) * (1 - scipy.special.erf(y)))
                SB_plus_Sn = SB + noise_spectrum
                Sobs = true_final_psd
                C11 = sum(np.log(DoF / SB_plus_Sn * scipy.stats.chi2.pdf(DoF*Sobs / SB_plus_Sn, DoF))) / len(Sobs) # NATURAL LOGARITHM!
                C11s.append(C11)
            return trial_log10_epsilons[np.argmax(C11s)]    # The best log10_epsilon is the one with the highest cost function

        trial_log10_epsilons0 = np.linspace(-11, -1, 51)    # First a rough trial (resolution=0.2)
        best_0 = best_fit_of_trial_log10_epsilons(trial_log10_epsilons0)
        trial_log10_epsilons1 = np.linspace(best_0 - 0.3, best_0 + 0.3, 61) # Then a bit finer trial (resolution=0.01)
        best_1 = best_fit_of_trial_log10_epsilons(trial_log10_epsilons1)
        trial_log10_epsilons2 = np.linspace(best_1 - 0.015, best_1 + 0.015, 16) # Then the finest trial (resolution=0.002)
        best_2 = best_fit_of_trial_log10_epsilons(trial_log10_epsilons2)
        return 10**best_2

    # Iterate!
    current_epsilon = starting_epsilon
    for j in range(n_iter):
        old_epsilon = current_epsilon
        current_thermdiss = get_thermdiss(current_epsilon)
        current_epsilon = get_epsilon(current_thermdiss)
        if abs(np.log10(current_epsilon) - np.log10(old_epsilon)) < 10**-5:
            break   # the resolution of log10(current_epsilon) is 2*10**-3 > 10**-5

    return current_epsilon, current_thermdiss, final_wavenumbers[0]

def fit_Batchelorjeve_krivulje_na_razliko_logaritmov(wavenumbers, psd, temp_mean, critical_wavenumber=100):
    ''' Fit to Batchelor curve using the least squares fit for difference between the logarithms of PSD and the mean
    logarithm of PSD. This method is completely independent from thermal dissipation.
    Acronym in MC_simulacija_spektrov.py: LSnTD.

    Inputs: (gTo = "grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m" + datumska_mapa + ".p")
        wavenumbers ... wavenumbers that correspond to values of PSD ("wavenumbers" in gTo, unit cpm),
        psd ... PSD from FFT of vertical derivatives of detrended temperature ("useful_grad_T_PSD" in gTo, unit °C**2*m**(-2)/cpm),
        temp_mean ... mean temperature ("tempcor" in gTo, unit °C),
        critical_wavenumber ... maximum wavenumber for fit (unit cpm)

    Outputs:
        log10_epsilon ... decimal logarithm of the best estimate for epsilon (unit log10(W/kg))
        log10_epsilon_sigma ... standard deviation of the best log10_epsilon
    '''

    from scipy.optimize import curve_fit
    import numpy as np
    from functools import partial
    from scipy.special import erf


    # ***************************************************************************************************
    # Here begins section of lines which is almost identical as in
    # iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral()
    #****************************************************************************************************

    nu = 1.702747 - 0.05126103 * temp_mean + 0.0005918645 * temp_mean ** 2
    nu *= 10 ** -6
    Dt = 1.4 * 10 ** -7

    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= critical_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])
    log10_proper_psd = np.log10(proper_psd)

    wn_subinterval_width = 5
    min_len_wn_for_fit = 30

    def lin_fit(k, a, b):
        return a * k + b

    final_wavenumbers = []
    final_log10_psd = []
    for iwn in range(len(proper_wavenumbers) - min_len_wn_for_fit):
        popt, pcov = curve_fit(lin_fit, proper_wavenumbers[iwn:iwn + wn_subinterval_width],
                               log10_proper_psd[iwn:iwn + wn_subinterval_width])
        if popt[0] - 2*np.sqrt(pcov[0][0]) > 0:
            final_wavenumbers = proper_wavenumbers[iwn:]
            final_log10_psd = log10_proper_psd[iwn:]
            break

    if len(final_wavenumbers) == 0:
        print('PSD is never rising in the dedicated area!')
        return None
    # *************************************************************************************************
    # End of section of almost identical lines
    # *************************************************************************************************


    def Batchelor_log10_minus_mean(k, log10_epsilon, nu=10**-6, kmin=1, kmax=100, kstep=1):
        ''' Function that we fit to.

        kmin, kmax and kstep are in cpm.'''
        q = 2 * np.sqrt(3)
        epsilon = 10 ** log10_epsilon
        k_b = (epsilon / (nu * Dt ** 2)) ** 0.25 / (2 * np.pi)  # cpm

        y = k * np.sqrt(q) / k_b
        ys = np.arange(start=kmin, stop=kmax + kstep, step=kstep) * np.sqrt(q) / k_b

        def Snorm(y):
            return y ** 2 * (np.exp(-y ** 2) / y - np.sqrt(np.pi) * (1 - erf(y)))

        Snorms = Snorm(ys)
        log10_snorms = np.log10(Snorms)

        return np.log10(Snorm(y)) - np.mean(log10_snorms)

    # The first fit
    fitted_Batchelor_log10_minus_mean = partial(Batchelor_log10_minus_mean, nu=nu, kmin=final_wavenumbers[0], kmax=final_wavenumbers[-1], kstep=final_wavenumbers[1] - final_wavenumbers[0])
    popt, pcov = curve_fit(fitted_Batchelor_log10_minus_mean, final_wavenumbers, final_log10_psd - np.mean(final_log10_psd), bounds=([-11],[-1]))

    log10_epsilon = popt[0]
    k_b = (10**log10_epsilon / (nu * Dt ** 2)) ** 0.25 / (2 * np.pi)
    k_star = k_b * 7.9 * 10**-3 # Up to k_star the spectrum goes as k**(1/3), from k_star on it goes as Batchelor curve


    # Repeat fit from k_star on!
    wavenumbers_after_k_star = [final_wavenumbers[i] for i in range(len(final_wavenumbers)) if final_wavenumbers[i] > k_star]
    psd_after_k_star = [final_log10_psd[i] for i in range(len(final_wavenumbers)) if final_wavenumbers[i] > k_star]

    fitted_Batchelor_log10_minus_mean = partial(Batchelor_log10_minus_mean, nu=nu, kmin=wavenumbers_after_k_star[0], kmax=wavenumbers_after_k_star[-1], kstep=wavenumbers_after_k_star[1] - wavenumbers_after_k_star[0])
    popt, pcov = curve_fit(fitted_Batchelor_log10_minus_mean, wavenumbers_after_k_star, psd_after_k_star - np.mean(psd_after_k_star), bounds=([-11],[-1]))
    log10_epsilon = popt[0]
    log10_epsilon_sigma = pcov[0][0]

    return log10_epsilon, log10_epsilon_sigma



def Batchelor_Fit_least_squares_log10(psd, wavenumbers, temp_mean, critical_wavenumber):
    ''' Fit to Batchelor curve using the least squares fit for logarithm of PSD.
     Acronym in MC_simulacija_spektrov.py: LS.

    Inputs: (gTo = "grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m" + datumska_mapa + ".p")
        wavenumbers ... wavenumbers that correspond to values of PSD ("wavenumbers" in gTo, unit cpm),
        psd ... PSD from FFT of vertical derivatives of detrended temperature ("useful_grad_T_PSD" in gTo, unit °C**2*m**(-2)/cpm),
        temp_mean ... mean temperature ("tempcor" in gTo, unit °C),
        critical_wavenumber ... maximum wavenumber for fit (unit cpm)

    Outputs:
        log10_epsilon ... decimal logarithm of the best estimate for epsilon (unit log10(W/kg))
        log10_epsilon_sigma ... standard deviation of the best log10_epsilon
    '''

    import numpy as np
    from scipy.optimize import curve_fit

    # ***************************************************************************************************
    # Here begins section of lines which is almost identical as in
    # iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral()
    #****************************************************************************************************

    nu = 1.702747 - 0.05126103 * temp_mean + 0.0005918645 * temp_mean ** 2
    nu *= 10 ** -6
    Dt = 1.4 * 10 ** -7

    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= critical_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])
    log10_proper_psd = np.log10(proper_psd)

    wn_subinterval_width = 5
    min_len_wn_for_fit = 30

    def lin_fit(k, a, b):
        return a * k + b

    final_wavenumbers = []
    final_log10_psd = []
    for iwn in range(len(proper_wavenumbers) - min_len_wn_for_fit):
        popt, pcov = curve_fit(lin_fit, proper_wavenumbers[iwn:iwn + wn_subinterval_width],
                               log10_proper_psd[iwn:iwn + wn_subinterval_width])
        if popt[0] - 2*np.sqrt(pcov[0][0]) > 0:
            final_wavenumbers = proper_wavenumbers[iwn:]
            final_log10_psd = log10_proper_psd[iwn:]
            break

    if len(final_wavenumbers) == 0:
        print('PSD is never rising in the dedicated area!')
        return None
    # *************************************************************************************************
    # End of section of almost identical lines
    # *************************************************************************************************


    def log10_Batchelor(k, log10_epsilon, log10_Thermdiss_unknown):
        '''Returns decimal logarithm of Batchelor curve, so it is more useful for curve_fit. The shape of the Batchelor
        curve is from Luketina&Imberger:Determining Turbulent Kinetic Energy Dissipation from Batchelor Curve Fitting,
        JAOT, 2001, https://doi.org/10.1175/1520-0426(2001)018<0100:DTKEDF>2.0.CO;2'''
        epsilon = 10**log10_epsilon
        k_b = epsilon**0.25 / (2*np.pi * nu**0.25 * Dt**0.5)
        q = 2 * np.sqrt(3)
        alpha = np.sqrt(2*q) * k / k_b
        t = 1 / (1 + 0.2316419 * alpha) # (Luketina&Imberger)
        Q = 1 / np.sqrt(2*np.pi) * np.exp(-alpha**2 / 2) * (0.319381530*t - 0.356563782*t**2 + 1.781477837*t**3 - 1.821255978*t**4 + 1.330274429*t**5)
        S_N = alpha * (np.exp(-alpha**2 / 2) - np.sqrt(2*np.pi) * alpha * Q)
        return np.log10(10**log10_Thermdiss_unknown / (2*Dt) * np.sqrt(2*q) / (2*np.pi*k_b) * S_N)


    try:# possible "RuntimeError: Optimal parameters not found: The maximum number of function evaluations is exceeded.", occurs very rarely
        popt, pcov = curve_fit(log10_Batchelor, final_wavenumbers, final_log10_psd,
                               bounds=([-11, -11], [-1, -1]))
        log10epsilon = popt[0]

        k_b = (10 ** log10epsilon / (nu * Dt ** 2)) ** 0.25 / (2 * np.pi)
        k_star = k_b * 7.9 * 10 ** -3  # Up to k_star the spectrum goes as k**(1/3), from k_star on it goes as Batchelor curve

        # Repeat fit from k_star on!
        wavenumbers_after_k_star = [final_wavenumbers[i] for i in range(len(final_wavenumbers)) if
                                    final_wavenumbers[i] > k_star]
        psd_after_k_star = [final_log10_psd[i] for i in range(len(final_wavenumbers)) if
                            final_wavenumbers[i] > k_star]

        popt, pcov = curve_fit(log10_Batchelor, wavenumbers_after_k_star,
                               psd_after_k_star, bounds=([-11, -11], [-1, -1]))
        log10epsilon = popt[0]
        log10epsilon_sigma = np.sqrt(pcov[0][0])
        log10Thermdiss = popt[1]
        log10Thermdiss_sigma = np.sqrt(pcov[1][1])
        cov_epsilon_Thermdiss = pcov[0][1]

        return log10epsilon, log10epsilon_sigma, log10Thermdiss, log10Thermdiss_sigma, cov_epsilon_Thermdiss
    except:
        print('probably runtime error')
        return None