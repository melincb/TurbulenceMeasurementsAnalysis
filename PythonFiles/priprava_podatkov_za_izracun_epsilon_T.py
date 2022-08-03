'''This file contains three things:
 1) The function to calculate PSD from temperature profiles (and some other accompanying quantities. This is done by
    function PSD_T.
 2) An example of a plot of detrended NTCHP profiles, their vertical derivatives and PSD for all casts within the same
    combination of station-route-depth.
 3) A procedure to merge the results from all casts from PSD_T function to sets that are separated with default pressure
    step and to sets with pressure step of 1dbar. WARNING: The data that this procedure needs is not uploaded to GitHub
    due to its extreme size. You can, however, use this procedure on your own data. The results from this procedure is
    stored to pickle files. The pickle files from my data are uploaded to GitHub, so they can be used properly in other
    Python files.'''
from Branje_datotek import Branje_datotek, iskanje_minimalne_in_maksimalne_vrednosti_nekih_kolicin_za_dolocen_datum_pri_vseh_meritvah
import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob, os
from scipy.stats import gmean

mapa_za_shranjevanje_grafov =  'Saved_figures\ '[:-1]    # Where it saves the figures
mapa_s_skriptami = os.getcwd()  # directory with scripts

plot_example = True     # To see the plot example, set this to True, else False
mean_depth_of_example = 10.25    # dbar (fill free to change it to any depth between 0.75 and 15.0, 10.25 shows an example of a jump in NTC sensor)
save_figures = False    # To save the plot example, set this to True, else False

izvirna_mapa = '..' # directory of origin of data
datumska_mapa = 'all_2008_08_21'    # date of measurements ('all_2008_08_21', 'all_2008_11_18', 'all_2009_04_20)
                                    # plot example always uses datumska_mapa='all_2008_08_21', because this is the only
                                    # date with proper uploaded data


# THE BASIC MANUAL SETTINGS END HERE


def PSD_T(izvirna_mapa, datumska_mapa, postaja, ruta, tipska_mapa='cutted', tlacni_korak=0.25):
    '''
    Returns a dictionary with the following keys (the values are also dictionaries, but all in the same format):
        "all_rolled_psd" ... rolled PSD from FFT of vertical gradient of detrended temperature (unit degreeC**2*m**(-2)/cpm)
        "all_rolled_psd_without_gradient" ... rolled PSD from FFT of detrended temperautre (unit degreeC**2*m**(-2)/cpm)
        "all_max_wavenumbers" ... maximum acceptable wavenumber from FFT according to NTC sensor's response time (always more than 10-times smaller than the Nyquist wavenumner) (unit cpm)
        "all_wavenumbers" ... wavenumbers related to PSD from FFT of vertical gradient of detrended temperature (unit cpm)
        "all_wavenumbers_without_gradient" ... wavenumbers related to PSD from FFT of detrended temperautre (unit cpm)
        "detrended_temperatures" ... detrended temperature profiles (from NTCHP channel) in each bin
        "d_detrended_T_dz" ... vertical derivative of detrendet demperature
        "all_pressure_centers" ... pressure centers of bins of data that was used to calculate FFT
        "all_pressure_Thermdiss" ... Thermal dissipation from MSSpro
        "all_pressure_epsilon_shear" ... epsilon from MSSpro
        "all_pressure_epsilon_shear1_epsilon_shear2_peps" ... epsilon1, epsilon2 and pseudo epsilon from MSSpro
        "all_pressure_Tempcor_sal" ... temperature and salinity from MSSpro (averaged values from "epsilon" directory)
        "all_pressure_NTC" ... temperature (NTC, "cutted" directory) from MSSpro
        "all_pressure_N_N2" ... buoyancy frequency and its squared value

    Inputs:
        izvirna_mapa ... directory of origin of data (in my case always "..")
        datumska_mapa ... determines the date of measurements ("all_2008_08_21" or "all_2008_11_18" or "all_2009_04_20")
        tipska_mapa ... final directory of TOB files ("epsilon", "shear" or "cutted") (Beware: TOB files from different
        directories contain different quantities!)
        postaja ... station ("LK01", "LK02", ..., "LK07")
        ruta ... route ("RUTA_01", "RUTA_02", ...)
        tlacni_korak ... pressure step between bin centers [dbar]

    Bin width is always 1 dbar (can be changed manually in code, look for variable sirina_predala, however also check
    all other Python files for correct interpretation of results!).
    TOB files should contain the NTCHP channel (in my case such channels are in "cutted" directories).
    NTC_response_time is set to 15 ms (according to personal comm. with IOW this is the worst-case scenario).
    The sampling frequency of the probe is set to 1024 Hz. If any other sampling frequency is used, change all 1024 to
    your sampling frequency in Hz.
    '''

    import numpy as np
    import glob, os
    from scipy.optimize import curve_fit

    NTC_response_time = 15 * 10**-3 # seconds, worst-case scenario (personal comm. with IOW)

    mapa_s_skriptami = os.getcwd()  # directory with scripts
    bs = '\ '[0]
    mapa_s_podatki_cutted = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa    # directory with TOB files
    os.chdir(mapa_s_podatki_cutted) # Go to TOB files


    vsi_tlaki_po_predalih = []  # all binned pressures
    vse_temperature_po_predalih = []    # all binned temperatures
    vsi_counti_po_predalih = [] # all binned counts (quantity 'Count' in TOB files)
    vsi_centri_predalov = []    # all bin centers (i.e. pressure centers of bins)
    vse_detrendane_temperature = [] # all detrended temperatures
    vse_detrendane_temperature_Hann = []    # all detrended temperatures multiplied by Hann window
    vsi_dt_dz_Hann = [] # all vertical temperature gradients multiplied by Hann window
    vsi_dt_dz = []  # all vertical temperature gradients
    vsi_rolani_psdji = []   # all rolled PSD
    vsi_rolani_psdji_brez_gradienta = []    # all rolled PSD from detrended temperatures (without vertical gradient)
    vsa_valovna_stevila = []    # all wavenumbers for PSD
    vsa_valovna_stevila_brez_gradienta = [] # all wavenumbers for PSD from detrended temperatures (without vertical gradient)
    vsa_max_valovna_stevila = []    # all maximum wavenumbers
    for datoteka in glob.glob('*_' + postaja + '_' + ruta + '.tob'):    # Go through all files of interest
        print(datoteka) # Print filename
        vse_kolicine = {}   # all quantities
        enote_vseh_kolicin = {} # units of all quantities
        # ********************************************************************
        # For the following section, see Branje_datotek() in Branje_datotek.py
        # ********************************************************************
        with open(datoteka, 'r') as datoteka:   # Enter file
            datasets, units = False, False
            tokratne_kolicine = []
            count = 0
            for vrstica in datoteka:
                if not datasets:
                    count = 0
                    try:
                        if vrstica.split()[1] == 'Datasets':
                            datasets = True
                            for kolicina in vrstica.split()[2:]:
                                tokratne_kolicine.append(kolicina)
                                if kolicina not in vse_kolicine.keys():
                                    vse_kolicine[kolicina] = []
                                    enote_vseh_kolicin[kolicina] = []
                    except:
                        pass
                elif count == 1:
                    units = True
                    count_units = 0
                    for enota in vrstica.split()[1:]:
                        enote_vseh_kolicin[tokratne_kolicine[count_units]].append(enota)
                        count_units += 1
                elif count > 2:
                    count_values = 0
                    for vrednost in vrstica.split()[1:]: # this time there is no if-else, because no quantities of
                                                         # interest are logarithmic
                        vse_kolicine[tokratne_kolicine[count_values]].append(float(vrednost))
                        count_values += 1
                count += 1
        # ********************************************************************
        # End of section from Branje_datotek()
        # ********************************************************************

        # __________________________________
        # THIS TIME THERE IS NO SORTING of pressure values (we will resort to "given" pressure rather than "measured"
        # pressure

        # __________________________________
        # ATTENTION! The cutted files are supposed to be correctly cutted, however in some cases (18.11.2008) I've
        # noticed that the pressure in the top 300 rows doesn't fall at all (it actually rises a bit). Apparently the
        # probe wasn't yet falling (it is somehow slowly rising). So I decided to remove all the rows above row which
        # contains the global minimum of pressure (once the global minimum is reached, the probe allegedly falls freely).

        kolicine = list(vse_kolicine.keys())    # quantities
        nesortirane_vrednosti_kolicin = np.array(list(vse_kolicine.values()))   # unsorted values of quantities
        kateri_po_vrsti_je_tlak = kolicine.index('Press')   # index of 'Press'
        katera_po_vrsti_je_temperatura = kolicine.index('NTCHP')    # index of 'NTCHP'
        kateri_po_vrsti_je_count = kolicine.index('Count')  # index of 'Count'

        indeks_globalnega_minimuma_tlaka = np.argmin(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak])
                                            # index of measured pressure's global minimum
        nesortirane_vrednosti_kolicin = nesortirane_vrednosti_kolicin[:,indeks_globalnega_minimuma_tlaka:]
                                            # cut off data before measured pressure's global minimum

        minimalen_tlak = min(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak])    # minimum measured pressure
        maksimalen_tlak = max(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak])   # maximum measured pressure

        # ________________________________
        # BINNING to bins of width of 1 dbar, where bin center step is tlacni_korak
        # A bin is successfully made, if there exists a measured pressure which is too low for this bin and a measured
        # pressure which is too high for this bin. Thus we ensure, that the bin width is really 1 dbar.
        sirina_predala = 1 # dbar
        ta_center_predala = tlacni_korak    # pressure center of the attempted bin (the first attempt is tlacni_korak)
        centri_predalov = []    # bin centers
        temperature_po_predalih = []    # temperatures in bins
        tlaki_po_predalih = []  # pressures in bins
        counti_po_predalih = [] # Counts in bins
        povprecne_hitrosti_po_predalih = [] # probe's mean vertical velocity in bins
        koraki_v_prostoru_po_predalih = []  # spatial steps in bins
        while ta_center_predala < maksimalen_tlak:
            if minimalen_tlak < ta_center_predala - sirina_predala/2 and\
                maksimalen_tlak > ta_center_predala + sirina_predala/2: # The conditions to make the bin are satisfied.
                # The lower end of the bin is at the index at which the pressure is for the first time large enough to
                # be in the bin.
                for iindex in range(len(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak])):
                    if nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak][iindex] >= ta_center_predala - sirina_predala/2:
                        indeks_spodnjega_roba = iindex  # index of the lower end
                        break
                # The upper end of the bin is at the index at which the pressure in the inverted list of measurements is
                #  for the first time small enough to be in the bin.
                for iindex in range(1, len(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak]) - 1):
                    if nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak][-iindex] <= ta_center_predala + sirina_predala/2:
                        indeks_zgornjega_roba = -iindex # index of the upper end
                        break
                centri_predalov.append(ta_center_predala)   # Append pressure center as the bin center
                tlaki_po_predalih.append(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak]\
                            [indeks_spodnjega_roba:indeks_zgornjega_roba+1])   # Append measured pressures in this bin
                temperature_po_predalih.append(nesortirane_vrednosti_kolicin[katera_po_vrsti_je_temperatura]\
                                    [indeks_spodnjega_roba:indeks_zgornjega_roba+1]) # Append temperatures in this bin
                counti_po_predalih.append(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_count]\
                                        [indeks_spodnjega_roba:indeks_zgornjega_roba + 1]) # Append Counts in this bin

                # Calculate probe's mean vertical velocity while it was sampling for this bin
                povprecna_hitrost = sirina_predala / ((nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_count][indeks_zgornjega_roba+1] - nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_count][indeks_spodnjega_roba]) * 1 / 1024)  # sampling frequency was 1024 Hz
                # Calculate probe's mean vertical spatial step while it was sampling for this bin
                korakec_v_prostoru = sirina_predala / (nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_count][indeks_zgornjega_roba+1] - nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_count][indeks_spodnjega_roba])
                povprecne_hitrosti_po_predalih.append(povprecna_hitrost)    # Append mean velocity
                koraki_v_prostoru_po_predalih.append(korakec_v_prostoru)    # Append spatial step

            ta_center_predala += tlacni_korak   # In the following iteration try a bit larger pressure center

        # _____________________________
        # REMOVING LINEAR TREND IN TEMPERATURES (in time!), MULTIPLYING BY HANN WINDOW AND COMPUTING PSD
        def lin_funk(t, a, b):
            return a*t + b


        detrendane_temperature = [] # detrended temperatures
        detrendani_trenutni_dt_dz = []  # current derivatives of detrended temperatures
        detrendane_temperature_Hann = []    # Hann window * (detrended temperatures)
        dt_dz_Hann = [] # Hann window * (current derivatives of detrended temperatures)
        rolani_psdji = []   # rolled PSDs for cases with vertical derivation
        rolani_psdji_brez_gradienta = []    # rolled PSDs for cases in which there is no vertical derivation
        max_valovna_stevila = []    # maximum wavenumbers (in all cases)
        valovna_stevila = []    # wavenumbers for cases with vertical derivation
        valovna_stevila_brez_gradienta = [] # wavenumbers for cases in which there is no vertical derivation
        for icenter in range(len(centri_predalov)): # Go through all bins
            trenutni_counti = np.array(counti_po_predalih[icenter]) # current counts
            trenutne_temperature = np.array(temperature_po_predalih[icenter])   # current temperatures
            trenutni_korakec_v_prostoru = koraki_v_prostoru_po_predalih[icenter]    # current spatial step
            trenutna_povprecna_hitrost = povprecne_hitrosti_po_predalih[icenter]    # current mean probe's ver. velocity
            popt, pcov = curve_fit(lin_funk, trenutni_counti, trenutne_temperature) # lin. fit of current temperatures
            detrendane_temperature.append(trenutne_temperature - lin_funk(trenutni_counti, *popt))  # detrend current
                                                                                                    # temperatures
            # Calculate vertical derivatives of current detrended temperatures
            detrendan_trenutni_dt_dz = [(-detrendane_temperature[-1][i+2] + 8*detrendane_temperature[-1][i+1] - 8*detrendane_temperature[-1][i-1] + detrendane_temperature[-1][i-2]) / (12*trenutni_korakec_v_prostoru) for i in range(2, len(detrendane_temperature[-1]) - 2)]
            detrendani_trenutni_dt_dz.append(detrendan_trenutni_dt_dz.copy())   # Append vertical derivatives
            Hannovo_okno = np.hanning(len(trenutne_temperature))    # Hann window for cases in which there is no
                                                                    # vertical derivation
            Hannovo_okno_trenutni_dt_dz = np.hanning(len(detrendan_trenutni_dt_dz)) # Hann window for cases with
                                                                                    # vertical derivation
            # Multiply detrended temperatures with Hann window
            detrendane_temperature_Hann.append(Hannovo_okno*detrendane_temperature[-1])
            # Multiply vertical derivatives of detrended temperatures with Hann window
            dt_dz_Hann.append(Hannovo_okno_trenutni_dt_dz * np.array(detrendan_trenutni_dt_dz))
            N = len(dt_dz_Hann[-1])
            # Ensure that ft includes 0-frequency (i.e. we need an even list in FFT)
            if N%2 == 1:
                N -= 1
                ft = np.fft.fft(dt_dz_Hann[-1][:-1])    # FFT of vertical derivative of detrended temperature
            else:
                ft = np.fft.fft(dt_dz_Hann[-1][:])  # FFT of vertical derivative of detrended temperature
            N_brez_gradienta = len(detrendane_temperature_Hann[-1])
            # Ensure that ft_brez_gradienta includes 0-frequency (i.e. we need an even list in FFT)
            if N_brez_gradienta%2 == 1:
                N_brez_gradienta -= 1
                ft_brez_gradienta = np.fft.fft(detrendane_temperature_Hann[-1][:-1]) # FFT of detrended temperature
            else:
                ft_brez_gradienta = np.fft.fft(detrendane_temperature_Hann[-1][:])  # FFT of detrended temperature
            L = sirina_predala  # Bin width is 1 dbar
            max_valovno_stevilo = 1 / (trenutna_povprecna_hitrost * NTC_response_time)  # maximum wavenumber
            T = len(trenutni_counti) / 1024 # Measurement duration in seconds
            # Frequencies for cases with vertical derivation in Hz
            frekvence = (np.linspace(-N / (2 * T), N / (2 * T), N, endpoint=False))
            # Frequencies for cases in which there is no vertical derivation in Hz
            frekvence_brez_gradienta = (np.linspace(-N_brez_gradienta / (2 * T), N_brez_gradienta / (2 * T),\
                                                    N_brez_gradienta, endpoint=False))
            def corregeix(F, Pxx):
                '''Original for Matlab: Maria Elena Roget, University of Girona
                F... frequency in Hz
                Pxx... PSD before correction'''
                f = F
                fr = 2*np.pi*f
                T1 = np.sqrt((0.15*fr)**2+1) / (np.sqrt((0.0015*fr)**2 + 1) * np.sqrt((0.001*fr)**2 + 1))
                T2 = np.sqrt((0.012*fr)**2+1) / (np.sqrt((0.0012*fr)**2 + 1) * np.sqrt((0.001*fr)**2 + 1))
                T3 = 27 * np.pi / np.sqrt(fr**2 + (27*np.pi)**2)
                T = T1 * T2 * T3

                return Pxx / T**2

            nondeemphasized_rolled_psd = 2*np.roll(np.absolute(ft) ** 2 / sum(Hannovo_okno_trenutni_dt_dz)**2, N//2)    # ONLY POSITIVE FREQUENCIES!
            nondeemphasized_rolled_psd_brez_gradienta = 2*np.roll(np.absolute(ft_brez_gradienta) ** 2 / sum(Hannovo_okno)**2, N//2) # ONLY POSITIVE FREQUENCIES!
            # Time to deemphasize the PSD (feature of NTCHP)
            deemphasized_rolled_psd = corregeix(np.abs(frekvence), nondeemphasized_rolled_psd)  # np.abs(frekvence) allows us to use them all and doesn't affect the important (i.e. positive) frequencies
            deemphasized_rolled_psd_brez_gradienta = corregeix(np.abs(frekvence_brez_gradienta), nondeemphasized_rolled_psd_brez_gradienta)
            rolani_psdji.append(deemphasized_rolled_psd)
            rolani_psdji_brez_gradienta.append(deemphasized_rolled_psd_brez_gradienta)
            # wavenumbers in cpm
            valovna_stevila.append(np.linspace(-N / (2 * L), N / (2 * L), N, endpoint=False))
            valovna_stevila_brez_gradienta.append(np.linspace(-N_brez_gradienta / (2 * L), N_brez_gradienta / (2 * L), N_brez_gradienta, endpoint=False))
            max_valovna_stevila.append(max_valovno_stevilo)


        vsi_tlaki_po_predalih.append(tlaki_po_predalih.copy())
        vse_temperature_po_predalih.append(temperature_po_predalih.copy())
        vsi_counti_po_predalih.append(counti_po_predalih.copy())
        vsi_centri_predalov.append(centri_predalov.copy())
        vse_detrendane_temperature.append(detrendane_temperature.copy())
        vse_detrendane_temperature_Hann.append(detrendane_temperature_Hann.copy())
        vsi_dt_dz.append(detrendani_trenutni_dt_dz.copy())
        vsi_dt_dz_Hann.append(dt_dz_Hann.copy())
        vsi_rolani_psdji.append(rolani_psdji.copy())
        vsi_rolani_psdji_brez_gradienta.append(rolani_psdji_brez_gradienta.copy())
        vsa_valovna_stevila.append(valovna_stevila)
        vsa_valovna_stevila_brez_gradienta.append(valovna_stevila_brez_gradienta)
        vsa_max_valovna_stevila.append(max_valovna_stevila)

    vrni_me_sem = os.getcwd()
    os.chdir(mapa_s_skriptami)
    vsi_posamezne_kolicine = Branje_datotek(izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, tipska_mapa='epsilon', postaja=postaja, ruta=ruta, vrni_kolicine_brez_povprecenja=True)
    os.chdir(vrni_me_sem)
    vsi_tlak_Thermdiss = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['Thermdiss']]for i in range(len(vsi_posamezne_kolicine))]  # Thermdiss from MSSpro
    vsi_tlak_Tempcor_in_sal = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['Tempcor'], vsi_posamezne_kolicine[i]['sal']] for i in range(len(vsi_posamezne_kolicine))]    # Temperature and salinity from MSSpro
    vsi_tlak_epsilon_shear = None
    vsi_tlak_epsilon_shear1_epsilon_shear2_peps = None
    vsi_tlak_N_N2 = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['N'], vsi_posamezne_kolicine[i]['N^2']] for i in range(len(vsi_posamezne_kolicine))]    # Buoyancy frequency and its squared value from MSSpro

    try:
        vsi_tlak_epsilon_shear = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['epsilon']]for i in range(len(vsi_posamezne_kolicine))]    # epsilon from MSSpro
        vsi_tlak_epsilon_shear1_epsilon_shear2_peps = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['epsilon1'], vsi_posamezne_kolicine[i]['epsilon2'], vsi_posamezne_kolicine[i]['peps']]for i in range(len(vsi_posamezne_kolicine))]    # epsilon1, epsilon2 and pseudo epsilon from MSSpro
    except:
        pass


    os.chdir(mapa_s_skriptami)  # We return to our scripts
    return {'all_rolled_psd':vsi_rolani_psdji, 'all_max_wavenumbers':vsa_max_valovna_stevila, 'all_wavenumbers':vsa_valovna_stevila, 'all_pressure_centers':vsi_centri_predalov, 'all_pressure_Thermdiss':vsi_tlak_Thermdiss, 'all_pressure_epsilon_shear':vsi_tlak_epsilon_shear, 'all_pressure_epsilon_shear1_epsilon_shear2_peps':vsi_tlak_epsilon_shear1_epsilon_shear2_peps, 'all_pressure_Tempcor_sal':vsi_tlak_Tempcor_in_sal, 'detrended_temperatures':vse_detrendane_temperature, 'd_detrended_T_dz':vsi_dt_dz, 'all_pressure_NTC':vse_temperature_po_predalih, 'all_pressure_N_N2':vsi_tlak_N_N2, 'all_rolled_psd_without_gradient':vsi_rolani_psdji_brez_gradienta, 'all_wavenumbers_without_gradient':vsa_valovna_stevila_brez_gradienta}


# PLOT EXAMPLE OF GETTING SEPARATE PSD
if plot_example:
    original_izvirna_mapa = izvirna_mapa
    original_datumska_mapa = datumska_mapa
    izvirna_mapa = '..'
    datumska_mapa = 'all_2008_08_21'
    leto = datumska_mapa[6:8]
    mesec = datumska_mapa[9:11]
    dan = datumska_mapa[12:]
    tipska_mapa = 'shear'
    bs = '\ '[0]
    postaja='LK01'
    ruta='RUTA_01'
    mapa_s_podatki_cutted = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa
    globina = mean_depth_of_example


    generate_epsilon_T = PSD_T(izvirna_mapa, datumska_mapa, postaja, ruta)
    all_rolled_psd = generate_epsilon_T['all_rolled_psd']
    all_max_wavenumbers = generate_epsilon_T['all_max_wavenumbers']
    all_wavenumbers = generate_epsilon_T['all_wavenumbers']
    all_pressure_centers = generate_epsilon_T['all_pressure_centers']
    all_pressure_Thermdiss = generate_epsilon_T['all_pressure_Thermdiss']
    all_pressure_epsilon_shear = generate_epsilon_T['all_pressure_epsilon_shear']
    all_pressure_epsilon_shear1_epsilon_shear2_peps = generate_epsilon_T['all_pressure_epsilon_shear1_epsilon_shear2_peps']
    all_pressure_Tempcor_sal = generate_epsilon_T['all_pressure_Tempcor_sal']
    detrended_T = generate_epsilon_T['detrended_temperatures']
    d_detrended_T_dz = generate_epsilon_T['d_detrended_T_dz']
    all_pressure_NTC = generate_epsilon_T['all_pressure_NTC']
    all_rolled_psd_without_gradient = generate_epsilon_T['all_rolled_psd_without_gradient']
    all_wavenumbers_without_gradient = generate_epsilon_T['all_wavenumbers_without_gradient']

    min_pressure_center, max_pressure_center = min([min(i) for i in all_pressure_centers]), max([max(i) for i in all_pressure_centers])
    pressure_center_step = all_pressure_centers[0][1] - all_pressure_centers[0][0]
    all_possible_pressure_centers = np.arange(start=min_pressure_center, stop=max_pressure_center+pressure_center_step, step=pressure_center_step)
    grouped_rolled_psd = {str(pc):[] for pc in all_possible_pressure_centers}
    grouped_max_wavenumbers = {str(pc):[] for pc in all_possible_pressure_centers}
    grouped_wavenumbers = {str(pc):[] for pc in all_possible_pressure_centers}
    grouped_detrended_T = {str(pc):[] for pc in all_possible_pressure_centers}
    grouped_d_detrended_T_dz = {str(pc):[] for pc in all_possible_pressure_centers}
    grouped_pressure_NTC = {str(pc):[] for pc in all_possible_pressure_centers}
    for iprofil in range(len(all_pressure_centers)):
        if all_pressure_centers[iprofil][0] < min_pressure_center:
            for j in all_possible_pressure_centers:
                if j < all_pressure_centers[iprofil][0]:
                    grouped_rolled_psd[str(j)].append([])
                    grouped_max_wavenumbers[str(j)].append(0)
                    grouped_wavenumbers[str(j)].append([])
                    grouped_detrended_T[str(j)].append([])
                    grouped_d_detrended_T_dz[str(j)].append([])
                    grouped_pressure_NTC[str(j)].append([])

                else:
                    break
        for iglobina in range(len(all_pressure_centers[iprofil])):
            grouped_rolled_psd[str(all_pressure_centers[iprofil][iglobina])].append(all_rolled_psd[iprofil][iglobina])
            grouped_max_wavenumbers[str(all_pressure_centers[iprofil][iglobina])].append(all_max_wavenumbers[iprofil][iglobina])
            grouped_wavenumbers[str(all_pressure_centers[iprofil][iglobina])].append(all_wavenumbers[iprofil][iglobina])
            grouped_detrended_T[str(all_pressure_centers[iprofil][iglobina])].append(detrended_T[iprofil][iglobina])
            grouped_d_detrended_T_dz[str(all_pressure_centers[iprofil][iglobina])].append(d_detrended_T_dz[iprofil][iglobina])
            grouped_pressure_NTC[str(all_pressure_centers[iprofil][iglobina])].append(all_pressure_NTC[iprofil][iglobina])
        if all_pressure_centers[iprofil][-1] < max_pressure_center:
            for j in all_possible_pressure_centers:
                if j > all_pressure_centers[iprofil][-1]:
                    grouped_rolled_psd[str(j)].append([])
                    grouped_max_wavenumbers[str(j)].append(0)
                    grouped_wavenumbers[str(j)].append([])
                    grouped_detrended_T[str(j)].append([])
                    grouped_d_detrended_T_dz[str(j)].append([])
                    grouped_pressure_NTC[str(j)].append([])


    shapes = ['o', '^', 's', '*', 'x', 'D']
    linestyles = ['-', '--', ':', '-.', '-', '-']
    markers_for_lines = ['None', 'None', 'None', 'None', 'o', '^']
    plt.figure(figsize=(12,8))
    ax1 = plt.subplot(2,3,2)
    ax2 = plt.subplot(1,3,3)
    ax3 = plt.subplot(2,3,5)
    ax4 = plt.subplot(1,3,1)
    axs = [ax1, ax2, ax3, ax4]
    axs01 = axs[3].twiny()
    for i in range(len(grouped_rolled_psd[str(globina)])):
        axs[0].plot(grouped_detrended_T[str(globina)][i], np.linspace(-0.5,0.5, len(grouped_detrended_T[str(globina)][i])), color='tab:red', linestyle=linestyles[i], marker=markers_for_lines[i], markersize=2.5, zorder=1, label=str(i))
        axs[2].scatter(np.array(grouped_d_detrended_T_dz[str(globina)][i]), np.linspace(-0.5, 0.5, len(grouped_d_detrended_T_dz[str(globina)][i])), s=10, marker=shapes[i], label=str(i))
        axs[1].scatter(grouped_wavenumbers[str(globina)][i], grouped_rolled_psd[str(globina)][i], s=10, marker=shapes[i], label=str(i))
        #axs[1].plot(grouped_wavenumbers[str(globina)][i], grouped_rolled_psd[str(globina)][i], color='tab:red', linestyle=linestyles[i], marker=markers_for_lines[i], markersize=2.5, label=str(i))
        axs[1].axvline(grouped_max_wavenumbers[str(globina)][i], color=f'C{i}')
        axs[1].axvline(grouped_max_wavenumbers[str(globina)][i]*15/12, color=f'C{i}', linestyle='--')
        axs[3].plot(all_pressure_Tempcor_sal[i][1], all_pressure_Tempcor_sal[i][0], color='tab:red', linestyle=linestyles[i], marker=markers_for_lines[i], markersize=2.5)
        axs[3].plot([], [], color='k', linestyle=linestyles[i], marker=markers_for_lines[i], markersize=2.5, label=str(i))  # This is a fake plot just to make lines in the legend appear black
        axs01.plot(all_pressure_Tempcor_sal[i][2], all_pressure_Tempcor_sal[i][0], color='tab:blue', linestyle=linestyles[i], marker=markers_for_lines[i], markersize=2.5)


    axs[1].set_yscale('log')
    axs[1].set_xscale('log')
    axs[1].set_xlabel(r'$k$ [cpm]')
    axs[1].set_ylabel(r'PSD [$\degree \mathrm{C}^{2} \mathrm{m}^{-2}$/cpm]')
    axs[1].grid(linestyle=':', linewidth=0.6)
    axs[1].set_xlim(1, max([np.max(grouped_wavenumbers[str(globina)][i]) for i in range(len(grouped_rolled_psd[str(globina)]))]))
    axs[1].legend()
    axs[0].invert_yaxis()
    axs[0].legend()
    axs[0].set_ylabel(r'Distance from center [dbar]')
    axs[0].set_xlabel(r"Detrended $T$ [$\degree$C] = $T'$", color='tab:red')
    axs[0].tick_params(axis='x', labelcolor='tab:red')
    axs[0].grid(axis='y', linestyle=':', linewidth=0.6)
    axs[0].legend()
    axs[2].set_xlabel(r"$\mathrm{d}T'/\mathrm{d}z$ [$\degree$C/m]", color='tab:blue')
    axs[2].set_ylabel(r'Distance from center [dbar]')
    axs[2].tick_params(axis='x', labelcolor='tab:blue')
    axs[2].grid(axis='y', linestyle=':', linewidth=0.6)
    axs[2].invert_yaxis()
    axs[2].legend()
    axs[3].invert_yaxis()
    axs[3].legend()
    axs[3].set_ylabel(r'$p$ [dbar]')
    axs[3].set_xlabel(r'$T$ [$\degree$C]', color='tab:red')
    axs[3].tick_params(axis='x', labelcolor='tab:red')
    axs[3].grid(axis='y', linestyle=':', linewidth=0.6)
    axs01.set_xlabel(r'$S$ [PSU]', color='tab:blue')
    axs01.tick_params(axis='x', labelcolor='tab:blue')


    min_max_of_the_day = iskanje_minimalne_in_maksimalne_vrednosti_nekih_kolicin_za_dolocen_datum_pri_vseh_meritvah(['Tempcor', 'sal', 'epsilon', 'N^2', 'N', 'sig_t'], izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, tipska_mapa='epsilon')
    print('min and max of the day:', min_max_of_the_day)
    min_max_temperature_of_the_day, min_max_salinity_of_the_day = min_max_of_the_day['Tempcor'], min_max_of_the_day['sal']
    diff_in_min_max = abs((min_max_temperature_of_the_day[1] - min_max_temperature_of_the_day[0]) - (min_max_salinity_of_the_day[1] - min_max_salinity_of_the_day[0]))
    if min_max_temperature_of_the_day[1] - min_max_temperature_of_the_day[0] > min_max_salinity_of_the_day[1] - min_max_salinity_of_the_day[0]:
        temperature_xlims = min_max_temperature_of_the_day
        salinity_xlims = [min_max_salinity_of_the_day[0] - diff_in_min_max/2,\
                          min_max_salinity_of_the_day[1] + diff_in_min_max/2]
    else:
        salinity_xlims = min_max_salinity_of_the_day
        temperature_xlims = [min_max_temperature_of_the_day[0] - diff_in_min_max / 2,\
                            min_max_temperature_of_the_day[1] + diff_in_min_max / 2]

    axs[3].set_xlim(temperature_xlims[0], temperature_xlims[1])
    axs01.set_xlim(salinity_xlims[0], salinity_xlims[1])
    axs[3].fill_between(temperature_xlims, [globina + 0.5, globina + 0.5], [globina - 0.5, globina - 0.5], color='gray', alpha=0.2)

    plt.suptitle(f'{dan}.{mesec}.20{leto}, station {postaja}, route {ruta[-2:]}')
    axs[1].set_title(f'Mean depth: {globina} dbar')
    axs[0].set_title(f'Mean depth: {globina} dbar')
    axs[2].set_title(f'Mean depth: {globina} dbar')
    plt.tight_layout()
    if save_figures:
        plt.savefig(mapa_za_shranjevanje_grafov + f'determination_of_PSD_from_NTCHP_example.jpg')
    plt.show()
    izvirna_mapa = original_izvirna_mapa
    datumska_mapa = original_datumska_mapa





input("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nWARNING: the procedure of merging data is about to begin. It won't work, if you don't give it your own data.\nTo save the results of this procedure, uncomment (and perhaps slightly change) the two pickle.dump lines at the end of this file.\nPress ENTER to continue with the process, otherwise kill the running program.\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

# MERGING DATA WITH COMMON ROUTES AND STATIONS
# The rest of this file shows, how data is merged. Since the proper files were too large to be all uploaded, one can
# only run this part of the code with their own data (of course, you should specify your own stations (postaja) in the
# most outer loop and routes (ruta) in the second most outer loop). The results from this part of the code are, however,
# saved in 'grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji' + datumska_mapa + '.p' and
# 'grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m' + datumska_mapa + '.p', so other scripts can run smoothly.

# merged outputs with original pressure step
outputs = {'epsilon_shear':[], 'useful_grad_T_PSD':[], 'peps':[], 'tempcor':[], 'wavenumbers':[], 'thermdiss_MSSpro':[], 'st_zdruzenih':[], 'postaja_ruta_globina':[], 'N2':[]}
# tempcor ... temperature from final TOB files, st_zdruzenih ... number of merged measurements, postaja_ruta_globina ... station-route-depth, N2 ... squared buoyancy fraquency

# merged outputs with pressure step = 1dbar (outputs with original pressure step have their pressure centers rounded to 1dbar)
outputs_na_1m = {'epsilon_shear':[], 'useful_grad_T_PSD':[], 'peps':[], 'tempcor':[], 'wavenumbers':[], 'thermdiss_MSSpro':[], 'st_zdruzenih':[], 'postaja_ruta_globina':[], 'N2':[]}
# quantities are same as for original pressure step

for postaja in ['LK01','LK02','LK03','LK04','LK05','LK06','LK07']:  # loop over all stations
    for ruta in [f'RUTA_0{i}' for i in range(1, 10)] + [f'RUTA_1{i}' for i in range(0, 10)]:    # loop over all routes
        generirano = False
        try:
            generate_epsilon_T = PSD_T(izvirna_mapa, datumska_mapa, postaja, ruta)
            generirano = True
        except: # The data might be corrupt, which would lead to an error, but we don't want our entire loop to crash
            pass
        if generirano:
            epsilon_shears = {} # epsilons from MSSpro
            rolled_psds = {}    # rolled PSDs
            tempcori = {}   # temperatures Tempcor from TOB files of "epsilon" directory
            pepsi = {}  # pseudo epsilons from MSSpro
            thermdiss_MSSpros = {}  # thermal dissipations from MSSpro
            N2s = {}    # buoyancy frequencies from MSSpro
            epsilon_shears_1m = {}  # same quantities, just for outputs_na_1m
            rolled_psds_1m = {} # same...
            tempcori_1m = {}    # same...
            pepsi_1m = {}   # same...
            thermdiss_MSSpros_1m = {}   # same...
            N2s_1m = {} # same...
            all_rolled_psd = generate_epsilon_T['all_rolled_psd']   # outputs from PSD_T()
            all_max_wavenumbers = generate_epsilon_T['all_max_wavenumbers']
            all_wavenumbers = generate_epsilon_T['all_wavenumbers']
            all_pressure_centers = generate_epsilon_T['all_pressure_centers']
            all_pressure_Thermdiss = generate_epsilon_T['all_pressure_Thermdiss']
            all_pressure_epsilon_shear = generate_epsilon_T['all_pressure_epsilon_shear']
            all_pressure_epsilon_shear1_epsilon_shear2_peps = generate_epsilon_T['all_pressure_epsilon_shear1_epsilon_shear2_peps']
            all_pressure_Tempcor_sal = generate_epsilon_T['all_pressure_Tempcor_sal']
            detrended_T = generate_epsilon_T['detrended_temperatures']
            d_detrended_T_dz = generate_epsilon_T['d_detrended_T_dz']
            all_pressure_NTC = generate_epsilon_T['all_pressure_NTC']
            all_pressure_N_N2 = generate_epsilon_T['all_pressure_N_N2']

            min_max_pressure_centers = False
            try:
                min_pressure_center, max_pressure_center = min([min(i) for i in all_pressure_centers]), max([max(i) for i in all_pressure_centers])
                min_max_pressure_centers = True
            except:
                pass
            if min_max_pressure_centers:
                pressure_center_step = all_pressure_centers[0][1] - all_pressure_centers[0][0]
                all_possible_pressure_centers = np.arange(start=min_pressure_center, stop=max_pressure_center+pressure_center_step, step=pressure_center_step)
                grouped_rolled_psd = {str(pc):[] for pc in all_possible_pressure_centers}   # initiate dictionaries
                grouped_max_wavenumbers = {str(pc):[] for pc in all_possible_pressure_centers}
                grouped_wavenumbers = {str(pc):[] for pc in all_possible_pressure_centers}
                grouped_detrended_T = {str(pc):[] for pc in all_possible_pressure_centers}
                grouped_d_detrended_T_dz = {str(pc):[] for pc in all_possible_pressure_centers}
                grouped_pressure_NTC = {str(pc):[] for pc in all_possible_pressure_centers}

                for iprofil in range(len(all_pressure_centers)):    # Append stuff to dictionaries' lists
                    if all_pressure_centers[iprofil][0] > min_pressure_center:
                        for j in all_possible_pressure_centers:
                            if j < all_pressure_centers[iprofil][0]:
                                grouped_rolled_psd[str(j)].append([])
                                grouped_max_wavenumbers[str(j)].append(0)
                                grouped_wavenumbers[str(j)].append([])
                                grouped_detrended_T[str(j)].append([])
                                grouped_d_detrended_T_dz[str(j)].append([])
                                grouped_pressure_NTC[str(j)].append([])

                            else:
                                break
                    for iglobina in range(len(all_pressure_centers[iprofil])):
                        grouped_rolled_psd[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_rolled_psd[iprofil][iglobina])
                        grouped_max_wavenumbers[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_max_wavenumbers[iprofil][iglobina])
                        grouped_wavenumbers[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_wavenumbers[iprofil][iglobina])
                        grouped_detrended_T[str(all_pressure_centers[iprofil][iglobina])].append(
                            detrended_T[iprofil][iglobina])
                        grouped_d_detrended_T_dz[str(all_pressure_centers[iprofil][iglobina])].append(
                            d_detrended_T_dz[iprofil][iglobina])
                        grouped_pressure_NTC[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_pressure_NTC[iprofil][iglobina])
                    if all_pressure_centers[iprofil][-1] < max_pressure_center:
                        for j in all_possible_pressure_centers:
                            if j > all_pressure_centers[iprofil][-1]:
                                grouped_rolled_psd[str(j)].append([])
                                grouped_max_wavenumbers[str(j)].append(0)
                                grouped_wavenumbers[str(j)].append([])
                                grouped_detrended_T[str(j)].append([])
                                grouped_d_detrended_T_dz[str(j)].append([])
                                grouped_pressure_NTC[str(j)].append([])


                for cast_no in range(len(all_pressure_epsilon_shear)):  # Go through all casts (usually 5)
                    for i in range(len(all_pressure_epsilon_shear[cast_no][0])):
                        # depth of pressure center
                        globina = all_pressure_epsilon_shear[cast_no][0][i]
                        # depth rounded to 1dbar
                        globina_1m = round(globina + 0.001,0)  # + 0.001 due to round error (e.g. 2.5 is storred as 2.4999999999 => + 0.001 it becomes 2.50099999, which rounds to 3.0

                        # Find Tempcor, Thermdiss and N2 for this cast number
                        tale_Tempcor_in_Thermdiss_in_N2 = [0, 0, 0]
                        for j in range(len(all_pressure_Tempcor_sal[cast_no][0])):
                            if all_pressure_Tempcor_sal[cast_no][0][j] == globina:
                                tale_Tempcor_in_Thermdiss_in_N2[0] = all_pressure_Tempcor_sal[cast_no][1][j]
                                break
                        for j in range(len(all_pressure_Thermdiss[cast_no][0])):
                            if all_pressure_Thermdiss[cast_no][0][j] == globina:
                                tale_Tempcor_in_Thermdiss_in_N2[1] = all_pressure_Thermdiss[cast_no][1][j]
                                break
                        for j in range(len(all_pressure_N_N2[cast_no][0])):
                            if all_pressure_N_N2[cast_no][0][j] == globina:
                                tale_Tempcor_in_Thermdiss_in_N2[2] = all_pressure_N_N2[cast_no][2][j]


                        if globina > max_pressure_center or globina < min_pressure_center:
                            # MSSpro returns Tempcor, Thermdiss and N2 to larger depths than we can calculate PSD (we
                            # have stricter conditions for calculation of PSD)
                            pass
                        elif 0 in tale_Tempcor_in_Thermdiss_in_N2:
                            print('No Tempcor or Thermdiss or N2')
                        elif len(grouped_rolled_psd[str(globina)][cast_no]) > 0:    # This condition is not fulfilled if
                                                                                    # the combination station-route
                                                                                    # doesn't exist
                            if max(np.abs(grouped_d_detrended_T_dz[str(globina)][cast_no])) >= 500:
                                print('Jump in sensor') # Criterion to figure out whether there was a jump in the NTC
                            else:   # Finally, all conditions are fulfilled
                                epsilon_shear = all_pressure_epsilon_shear[cast_no][1][i]
                                epsilon_shear1 = all_pressure_epsilon_shear1_epsilon_shear2_peps[cast_no][1][i]
                                epsilon_shear2 = all_pressure_epsilon_shear1_epsilon_shear2_peps[cast_no][2][i]
                                peps = all_pressure_epsilon_shear1_epsilon_shear2_peps[cast_no][3][i]
                                thermdiss_MSSpro = tale_Tempcor_in_Thermdiss_in_N2[1]

                                nu = 1.702747 - 0.05126103*tale_Tempcor_in_Thermdiss_in_N2[0] + 0.0005918645*tale_Tempcor_in_Thermdiss_in_N2[0]**2    # MSSpro user manual
                                nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual (nu is viscosity)

                                if str(globina) in epsilon_shears.keys():   # This depth already has a key
                                    epsilon_shears[str(globina)].append(epsilon_shear)
                                    rolled_psds[str(globina)].append(grouped_rolled_psd[str(globina)][cast_no])
                                    tempcori[str(globina)].append(tale_Tempcor_in_Thermdiss_in_N2[0])
                                    pepsi[str(globina)].append(peps)
                                    thermdiss_MSSpros[str(globina)].append(thermdiss_MSSpro)
                                    N2s[(str(globina))].append(tale_Tempcor_in_Thermdiss_in_N2[2])
                                else:   # Iniciate a key and make its value to be a list
                                    epsilon_shears[str(globina)] = [epsilon_shear]
                                    rolled_psds[str(globina)] = [grouped_rolled_psd[str(globina)][cast_no]]
                                    tempcori[str(globina)] = [tale_Tempcor_in_Thermdiss_in_N2[0]]
                                    pepsi[str(globina)] = [peps]
                                    thermdiss_MSSpros[str(globina)] = [thermdiss_MSSpro]
                                    N2s[(str(globina))] = [tale_Tempcor_in_Thermdiss_in_N2[2]]
                                if str(globina_1m) in epsilon_shears_1m.keys():
                                    epsilon_shears_1m[str(globina_1m)].append(epsilon_shear)
                                    rolled_psds_1m[str(globina_1m)].append(grouped_rolled_psd[str(globina)][cast_no])
                                    tempcori_1m[str(globina_1m)].append(tale_Tempcor_in_Thermdiss_in_N2[0])
                                    pepsi_1m[str(globina_1m)].append(peps)
                                    thermdiss_MSSpros_1m[str(globina_1m)].append(thermdiss_MSSpro)
                                    N2s_1m[(str(globina_1m))].append(tale_Tempcor_in_Thermdiss_in_N2[2])
                                else:
                                    epsilon_shears_1m[str(globina_1m)] = [epsilon_shear]
                                    rolled_psds_1m[str(globina_1m)] = [grouped_rolled_psd[str(globina)][cast_no]]
                                    tempcori_1m[str(globina_1m)] = [tale_Tempcor_in_Thermdiss_in_N2[0]]
                                    pepsi_1m[str(globina_1m)] = [peps]
                                    thermdiss_MSSpros_1m[str(globina_1m)] = [thermdiss_MSSpro]
                                    N2s_1m[(str(globina_1m))] = [tale_Tempcor_in_Thermdiss_in_N2[2]]

            print(f'Depths with at least one measurement that fulfilled all conditions for station {postaja} and route {ruta}', epsilon_shears.keys())
            for globina in epsilon_shears.keys():
                if len(epsilon_shears[globina]) < 1:    # If you want larger number for better statistics, increase
                    print(f'For station {postaja}, route {ruta} and depth {globina} we lack statistical significance')
                else:
                    povprecni_epsilon_shear = gmean(epsilon_shears[globina])
                    #NEXT LINE ONLY WORKS WELL IF YOU USE FIRST 100 WAVENUMBERS (otherwise change 100 to proper number)
                    samo_uporabni_psds = [rolled_psds[globina][i][len(rolled_psds[globina][i])//2:len(rolled_psds[globina][i])//2+100] for i in range(len(rolled_psds[globina]))]   # Ker so razlicnih dolzin, ne moremo racunati povprecij, zato gremo ze zdaj gledati samo uporabne psd.
                    povprecni_rolled_psd = gmean(samo_uporabni_psds, axis=0)
                    povprecni_tempcor = np.mean(tempcori[globina])
                    povprecni_peps = gmean(pepsi[globina])
                    povprecni_thermdiss_MSSpro = gmean(thermdiss_MSSpros[globina])

                    # Simplified look on useful wavenumbers. They are in cpm, ONLY WORKS WELL IF BIN WIDTH IS 1 dbar
                    valovna_stevila = np.linspace(1, 100, 100)
                    outputs['epsilon_shear'].append(np.log10(povprecni_epsilon_shear))
                    outputs['peps'].append(np.log10(povprecni_peps))
                    outputs['useful_grad_T_PSD'].append(povprecni_rolled_psd)
                    outputs['tempcor'].append(povprecni_tempcor)
                    outputs['wavenumbers'].append(valovna_stevila)
                    outputs['thermdiss_MSSpro'].append(povprecni_thermdiss_MSSpro)
                    outputs['N2'].append(np.mean(N2s[globina]))
                    outputs['st_zdruzenih'].append(len(epsilon_shears[globina]))
                    outputs['postaja_ruta_globina'].append(postaja + '_' + ruta + '_' + globina)

                    print('success')

            for globina in epsilon_shears_1m.keys():
                if len(epsilon_shears_1m[globina]) < 1: # If you want larger number for better statistics, increase
                    print(f'For station {postaja}, route {ruta} and depth {globina} we lack statistical significance')
                else:
                    povprecni_epsilon_shear = gmean(epsilon_shears_1m[globina])
                    #NEXT LINE ONLY WORKS WELL, IF YOU USE FIRST 100 WAVENUMBERS (otherwise change 100 to proper number)
                    samo_uporabni_psds = [rolled_psds_1m[globina][i][len(rolled_psds_1m[globina][i])//2:len(rolled_psds_1m[globina][i])//2+100] for i in range(len(rolled_psds_1m[globina]))]   # Ker so razlicnih dolzin, ne moremo racunati povprecij, zato gremo ze zdaj gledati samo uporabne psd.
                    povprecni_rolled_psd = gmean(samo_uporabni_psds, axis=0)
                    povprecni_tempcor = np.mean(tempcori_1m[globina])
                    povprecni_peps = gmean(pepsi_1m[globina])
                    povprecni_thermdiss_MSSpro = gmean(thermdiss_MSSpros_1m[globina])

                    # Simplified look on useful wavenumbers. They are in cpm, ONLY WORKS WELL IF BIN WIDTH IS 1 dbar
                    valovna_stevila = np.linspace(1, 100, 100)
                    outputs_na_1m['epsilon_shear'].append(np.log10(povprecni_epsilon_shear))
                    outputs_na_1m['peps'].append(np.log10(povprecni_peps))
                    outputs_na_1m['useful_grad_T_PSD'].append(povprecni_rolled_psd)
                    outputs_na_1m['tempcor'].append(povprecni_tempcor)
                    outputs_na_1m['wavenumbers'].append(valovna_stevila)
                    outputs_na_1m['thermdiss_MSSpro'].append(povprecni_thermdiss_MSSpro)
                    outputs_na_1m['N2'].append(np.mean(N2s_1m[globina]))
                    outputs_na_1m['st_zdruzenih'].append(len(epsilon_shears_1m[globina]))
                    outputs_na_1m['postaja_ruta_globina'].append(postaja + '_' + ruta + '_' + globina)
                    print('success')



os.chdir(mapa_s_skriptami)  # Return to scripts (otherwise you'll dump pickle files somewhere else)
# pickle.dump(outputs, open('grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji' + datumska_mapa + '.p', 'wb'))
# pickle.dump(outputs_na_1m, open('grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m' + datumska_mapa + '.p', 'wb'))
print('DONE')