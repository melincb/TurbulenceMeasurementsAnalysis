def fit_Nasmyth_least_squares_log10(psd, wavenumbers, min_wavenumber, max_wavenumber, temp_mean):
    '''Least squares fit to the Nasmyth spectrum. The analytical Nasmyth spectrum is from Lueck: Calculating the Rate
    of Dissipation of Turbulent Kinetic Energy, RSI, 2013,
    https://rocklandscientific.com/wp-content/uploads/2021/12/TN-028-Calculating-the-Rate-of-Dissipation-of-Turbulent-Kinetic-Energy.pdf
    In MC_simulacija_spektrov.py it has an acronym LS.

    Inputs: (s12o is f"shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p")
        psd ... shear fluctuation PSD (unit s**(-2)/cpm, in s12o this is "useful_PSD1" or "useful_PSD2"),
        wavenumbers ... wavenumbers that correspond to PSD (unit cpm, in s12o this is "wavenumbers"),
        min_wavenumber ... the lowest wavenumber to include in the fit,
        max_wavenumber ... the highest wavenumber to include in the fit,
        temp_mean ... the mean temperature (unit degreeC, in s12o this is "tempcor")

    Outputs:
        log10_epsilon ... decimal logarithm of the best estimate for epsilon (unit log10(W/kg))
        log10_epsilon_sigma ... standard deviation of the best log10_epsilon
    '''

    import numpy as np
    from scipy.optimize import curve_fit

    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual, nu is viscosity

    proper_wavenumbers = [wn for wn in wavenumbers if wn >= min_wavenumber and wn <= max_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])
    log10_proper_psd = np.log10(proper_psd) # decimal logarithm of proper PSD

    def log10_Nasmyth(k, log10epsilon):
        '''Returns decimal logarithm of Lueck's version of the Nasmyth spectrum (k is in cpm)'''
        epsilon = 10**log10epsilon
        x = k * (nu**3 / epsilon)**0.25
        return np.log10(8.05 * x**(1/3) / (1 + (20.6 * x)**3.715) * (epsilon**3/nu)**0.25)

    # Least squares fit
    popt, pcov = curve_fit(log10_Nasmyth, proper_wavenumbers[:], log10_proper_psd[:], bounds=([-13], [-1]))
    log10epsilon = popt[0]
    log10epsilon_sigma = np.sqrt(pcov[0][0])

    return log10epsilon, log10epsilon_sigma



def fit_Nasmyth_Lueck(psd, wavenumbers, min_wavenumber, max_wavenumber, temp_mean, n_iter=20):
    '''Fit to the Nasmyth spectrum using the integral method as proposed by Lueck: Calculating the Rate
    of Dissipation of Turbulent Kinetic Energy, RSI, 2013,
    https://rocklandscientific.com/wp-content/uploads/2021/12/TN-028-Calculating-the-Rate-of-Dissipation-of-Turbulent-Kinetic-Energy.pdf,
    just the wavenumber range here is custom between min_wavenumber and max_wavenumber.
    Iteration usually converges after 3-5 iterations, n_iter=20 is just to be safe (the iterations are very fast).
    min_wavenumber and max_wavenumber are both included in the fit.
    In MC_simulacija_spektrov.py it has an acronym IM.

    Inputs: (s12o is f"shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p")
        psd ... shear fluctuation PSD (unit s**(-2)/cpm, in s12o this is "useful_PSD1" or "useful_PSD2"),
        wavenumbers ... wavenumbers that correspond to PSD (unit cpm, in s12o this is "wavenumbers"),
        min_wavenumber .. the lower boundary of the range of wavenumbers to fit (in our case 6 cpm),
        max_wavenumber .. the upper boundary of the range of wavenumbers to fit (in our case 35 cpm),
        temp_mean ... the mean temperature (unit degreeC, in s12o this is "tempcor")

    Outputs:
        np.log10(epsilon_iterated) ... decimal logarithm of the best fitting epsilon
    '''

    import numpy as np
    import scipy.integrate

    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual, nu is viscosity

    # Take wavenumbers between min_wavenumber and max_wavenumber
    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= max_wavenumber and wn >= min_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])

    # Numerical integration using the trapezoidal rule
    trapz_integration = scipy.integrate.trapz(y=proper_psd, x=proper_wavenumbers)
    # epsilon from numerical integration
    epsilon_trial = 15/2 * nu * trapz_integration

    def fraction(epsilon, nu=nu, k1=min_wavenumber, k2=max_wavenumber):
        '''Computes the analytical fraction between the integral that we compute within our wavenumber range and the
        integral that we would get, if we integrated between 0 cpm and inf cpm.'''
        def int_psi(x):
            '''Lueck 2013, Eq. (6)'''
            return np.tanh(48 * x ** (4 / 3)) - 2.9 * x ** (4 / 3) * np.exp(-22.3 * x ** (4 / 3))

        x1 = k1 * (nu ** 3 / epsilon) ** 0.25
        x2 = k2 * (nu ** 3 / epsilon) ** 0.25
        return int_psi(x2) - int_psi(x1)

    # Start iterating
    epsilon_iterated = epsilon_trial
    for iiteration in range(n_iter):
        frac = fraction(epsilon_iterated)
        epsilon_iterated = epsilon_trial / frac

    return np.log10(epsilon_iterated)

def fit_Nasmyth_MLE(psd, wavenumbers, min_wavenumber, max_wavenumber, temp_mean, DoF):
    '''Fit to the Nasmyth curve using the Maximum Likelihood Estimate.
    In MC_simulacija_spektrov.py it has an acronym MLE.

    Inputs: (s12o is f"shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m{datumska_mapa}.p")
        psd ... shear fluctuation PSD (unit s**(-2)/cpm, in s12o this is "useful_PSD1" or "useful_PSD2"),
        wavenumbers ... wavenumbers that correspond to PSD (unit cpm, in s12o this is "wavenumbers"),
        min_wavenumber .. the lower boundary of the range of wavenumbers to fit (in our case 6 cpm),
        max_wavenumber .. the upper boundary of the range of wavenumbers to fit (in our case 35 cpm),
        temp_mean ... the mean temperature (unit degreeC, in s12o this is "tempcor"),
        DoF ... degrees of freedom in our PSD by which we assume the behaviour of random noise in PSD,

    Outputs:
        best_2 ... the best estimate for the decimal logarithm of epsilon (unit lo10(W/kg)) (the resolution is 0.002)
    '''

    import numpy as np
    import scipy.stats

    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual

    proper_wavenumbers = np.array([wn for wn in wavenumbers if wn >= min_wavenumber and wn <= max_wavenumber])
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])


    def best_fit_of_trial_log10_epsilons(trial_log10_epsilons):
        '''Similar to what is proposed by Ruddick et al: Maximum Likelihood Spectral Fitting: The Batchelor
        Spectrum, JAOT, 2000, https://doi.org/10.1175/1520-0426(2000)017<1541:MLSFTB>2.0.CO;2
        for the Batchelor spectrum, just adapted for the Nasmyth spectrum using the assumption that the instrumental
        noise is negligible in our range of wavenumbers.
        To change the instrumental noise to something else, adapt this function using the hints from
        iterative_fit_of_epsilon_using_MLE_and_thermdiss_using_integral() from
        fitting_procedures_for_Batchelor_curve.py.
        '''
        C11s = []  # notation from Ruddick et al
        for log10_epsilon in trial_log10_epsilons:
            epsilon = 10 ** log10_epsilon
            x = proper_wavenumbers * (nu ** 3 / epsilon) ** 0.25
            SN = 8.05 * x ** (1 / 3) / (1 + (20.6 * x) ** 3.715) * (epsilon ** 3 / nu) ** 0.25
            Sobs = proper_psd
            C11 = sum(np.log(DoF / SN * scipy.stats.chi2.pdf(DoF * Sobs / SN, DoF))) / len(Sobs)  # NATURAL LOGARITHM!

            C11s.append(C11)
        return trial_log10_epsilons[np.argmax(C11s)]

    trial_log10_epsilons0 = np.linspace(-11, -1, 51)   # 101
    best_0 = best_fit_of_trial_log10_epsilons(trial_log10_epsilons0)
    trial_log10_epsilons1 = np.linspace(best_0 - 0.3, best_0 + 0.3, 61) # 61
    best_1 = best_fit_of_trial_log10_epsilons(trial_log10_epsilons1)
    trial_log10_epsilons2 = np.linspace(best_1 - 0.015, best_1 + 0.015, 16) # 31
    best_2 = best_fit_of_trial_log10_epsilons(trial_log10_epsilons2)
    return best_2