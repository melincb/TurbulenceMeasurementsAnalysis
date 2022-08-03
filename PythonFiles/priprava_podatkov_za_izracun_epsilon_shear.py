'''This file contains three things:
 1) The function to calculate PSD from shear fluctuation profiles (and some other accompanying quantities. This is done by
    function PSD_shear.
 2) A procedure to merge the results from all casts from PSD_shear function to sets that are separated with default
    pressure step and to sets with pressure step of 1dbar. WARNING: The data that this procedure needs is not uploaded
    to GitHub due to its extreme size. You can, however, use this procedure on your own data. The results from this
    procedure is stored to pickle files. The pickle files from my data are uploaded to GitHub, so they can be used
    properly in other Python files.
 3) A procedure to estimate the degrees of freedom of shear PSD. Advice: if the final histograms don't give clear output,
    use integral method for the final computations in DATOTEKA PYTHON'''
from Branje_datotek import Branje_datotek
from scipy.optimize import curve_fit
from scipy.stats.mstats import gmean
import scipy.integrate
import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob, os

mapa_za_shranjevanje_grafov =  'Saved_figures\ '[:-1]    # Where it saves the figures
mapa_s_skriptami = os.getcwd()  # directory with scripts

def PSD_shear(izvirna_mapa, datumska_mapa, postaja, ruta, tipska_mapa='shear', tlacni_korak=0.25):
    '''
    Returns a dictionary with the following keys (the values are also dictionaries, but all in the same format):
        "all_rolled_psd1" ... rolled PSD from FFT of shear fluctuations from channel "shear1" in TOB files (unit s**(-2)/cpm),
        "all_rolled_psd2" ... rolled PSD from FFT of shear fluctuations from channel "shear2" in TOB files (unit s**(-2)/cpm),
        "all_wavenumbers" ... wavenumbers related to both PSDs (unit cpm),
        "all_pressure_centers" ... pressure centers of bins of data that was used to calculate FFT (unit dbar),
        "all_pressure_Thermdiss" ... Thermal dissipation from MSSpro,
        "all_pressure_epsilon_shear" ... epsilon from MSSpro,
        "all_pressure_epsilon_shear1_epsilon_shear2_peps" ... epsilon1, epsilon2 and pseudo epsilon from MSSpro,
        "all_pressure_Tempcor_sal" ... temperature and salinity from MSSpro (averaged values from "epsilon" directory),
        "all_pressure_N_N2" ... buoyancy frequency and its squared value,
        "all_shear1" ... values of shear fluctuations from channel "shear1" in TOB files,
        "all_shear2" ... values of shear fluctuations from channel "shear2" in TOB files

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
    TOB files should contain the shear1 and shear2 channels (in my case such channels are in "shear" directories). These
    channels don't contain the true values of shear, but their fluctuations (batch job _shear_c.msb, which we use to
    make TOB files in "shear" directory, detrends the true values)!
    The sampling frequency of the probe is set to 1024 Hz. If any other sampling frequency is used, change all 1024 to
    your sampling frequency in Hz.
    '''
    import numpy as np
    import glob, os

    mapa_s_skriptami = os.getcwd()  # directory with scripts
    bs = '\ '[0]
    mapa_s_podatki_cutted = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa    # directory with TOB files
    os.chdir(mapa_s_podatki_cutted) # Go to TOB files

    vsi_tlaki_po_predalih = []  # all binned pressures
    vsi_centri_predalov = []    # all bin centers (i.e. pressure centers of bins)
    vsi_shear1_po_predalih = [] # all values of "shear1" channel in each bin
    vsi_shear2_po_predalih = [] # all values of "shear2" channel in each bin
    vsi_shear1_Hann = []    # all values of "shear1" channel in each bin multiplied by Hann window
    vsi_shear2_Hann = []    # all values of "shear2" channel in each bin multiplied by Hann window
    vsi_rolani_psdji1 = []  # all rolled PSD of FFT of elements of vsi_shear1_Hann
    vsi_rolani_psdji2 = []  # all rolled PSD of FFT of elements of vsi_shear2_Hann
    vsa_valovna_stevila = []    # all wavenumbers

    for datoteka in glob.glob('*_' + postaja + '_' + ruta + '.tob'):    # Go through all files of interest
        print(datoteka)  # Print filename
        vse_kolicine = {}  # all quantities
        enote_vseh_kolicin = {}  # units of all quantities
        # ********************************************************************
        # For the following section, see Branje_datotek() in Branje_datotek.py
        # ********************************************************************
        with open(datoteka, 'r') as datoteka:  # Enter file
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
                    for vrednost in vrstica.split()[1:]:  # this time there is no if-else, because no quantities of
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
        kateri_po_vrsti_je_shear1 = kolicine.index('shear1')    # index of 'shear1'
        kateri_po_vrsti_je_shear2 = kolicine.index('shear2')    # index of 'shear2'

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
        shear1_po_predalih = [] # shear1 in bins
        shear2_po_predalih = [] # shear2 in bins
        tlaki_po_predalih = []  # pressures in bins
        povprecne_hitrosti_po_predalih = [] # probe's mean vertical velocity in bins
        koraki_v_prostoru_po_predalih = []  # spatial steps in bins

        while ta_center_predala < maksimalen_tlak:
            if minimalen_tlak < ta_center_predala - sirina_predala/2 and\
                maksimalen_tlak > ta_center_predala + sirina_predala/2: # The conditions to make the bin are satisfied.
                # The lower end of the bin is at the index at which the pressure is for the first time large enough to
                # be in the bin.
                for iindex in range(len(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak])):
                    if nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak][iindex] >= ta_center_predala - sirina_predala/2:
                        indeks_spodnjega_roba = iindex
                        break
                # The upper end of the bin is at the index at which the pressure in the inverted list of measurements is
                #  for the first time small enough to be in the bin.
                for iindex in range(1, len(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak]) - 1):
                    if nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak][-iindex] <= ta_center_predala + sirina_predala/2:
                        indeks_zgornjega_roba = -iindex
                        break
                centri_predalov.append(ta_center_predala)   # Append pressure center as the bin center
                tlaki_po_predalih.append(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak]\
                            [indeks_spodnjega_roba:indeks_zgornjega_roba+1])   # Append measured pressures in this bin
                # Chronologically I first wrote PSD_T() in priprava_podatkov_za_epsilon_T.py and then tried to make this
                # function similar to PSD_T(). However, PSD_T() uses files from "cutted" directory, which contain Count
                # channel, whereas files in "shear" directory don't contain Count channels, so I had to make a
                # substitute for Count, which is in the following line:
                nadomestilo_count = len(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak]\
                                             [indeks_spodnjega_roba:indeks_zgornjega_roba+1])
                shear1_po_predalih.append(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_shear1]\
                                    [indeks_spodnjega_roba:indeks_zgornjega_roba+1]) # Append shear1 in this bin
                shear2_po_predalih.append(nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_shear2]\
                                    [indeks_spodnjega_roba:indeks_zgornjega_roba+1]) # Append shear1 in this bin
                # Calculate probe's mean vertical velocity while it was sampling for this bin
                povprecna_hitrost = sirina_predala / (nadomestilo_count * 1 / 1024)  # sampling frequency was 1024 Hz
                # Calculate probe's mean vertical spatial step while it was sampling for this bin
                korakec_v_prostoru = sirina_predala / (nadomestilo_count)
                povprecne_hitrosti_po_predalih.append(povprecna_hitrost)    # Append mean velocity
                koraki_v_prostoru_po_predalih.append(korakec_v_prostoru)    # Append spatial step

            ta_center_predala += tlacni_korak   # In the following iteration try a bit larger pressure center

        # _____________________________
        # MULTIPLYING BY HANN WINDOW AND COMPUTING PSD
        # In comparison to PSD_T() from priprava_podatkov_za_izracun_epsilon_T.py there is no detrending, because we
        # expect the shear1 and shear2 values to be detrended prior to entering this function (_shear_c.msb does this)

        shear1_Hann = []    # Hann window * shear1
        shear2_Hann = []    # Hann window * shear2
        rolani_psdji1 = []  # rolled PSD of (Hann window * shear 1) (unit
        rolani_psdji2 = []  # rolled PSD of (Hann window * shear 2)
        valovna_stevila = []    # wavenumbers for both rolled PSDs (unit cpm)

        for icenter in range(len(centri_predalov)): # Go through all bins
            trenutni_shear1 = np.array(shear1_po_predalih[icenter]) # current shear1
            trenutni_shear2 = np.array(shear2_po_predalih[icenter]) # current shear1
            Hannovo_okno_trenutni_shear12 = np.hanning(len(trenutni_shear1))  # Hann window (len(trenutni_shear1) =
                                                                                # = len(trenutni_shear2) always)
            # Multiply shear1 with Hann window
            shear1_Hann.append(Hannovo_okno_trenutni_shear12*trenutni_shear1)
            # Multiply shear2 with Hann window
            shear2_Hann.append(Hannovo_okno_trenutni_shear12*trenutni_shear2)
            N = len(shear1_Hann[-1])    # len(shear1_Hann[-1]) = len(shear2_Hann[-1])
            # Ensure that ft1 and ft2 include 0-frequency (i.e. we need an even list in FFT)
            if N%2 == 1:
                N -= 1
                ft1 = np.fft.fft(shear1_Hann[-1][:-1])  # FFT of Hannovo_okno_trenutni_shear12*trenutni_shear1
                ft2 = np.fft.fft(shear2_Hann[-1][:-1])  # FFT of Hannovo_okno_trenutni_shear12*trenutni_shear2
            else:
                ft1 = np.fft.fft(shear1_Hann[-1][:])    # FFT of Hannovo_okno_trenutni_shear12*trenutni_shear1
                ft2 = np.fft.fft(shear2_Hann[-1][:])    # FFT of Hannovo_okno_trenutni_shear12*trenutni_shear2
            L = sirina_predala  # Bin width is 1 dbar

            # Rolled PSD1 (USE ONLY POSITIVE FREQUENCIES!)
            rolani_psdji1.append(2*np.roll(np.absolute(ft1) ** 2 / sum(Hannovo_okno_trenutni_shear12)**2, N//2))
            # Rolled PSD2 (USE ONLY POSITIVE FREQUENCIES!)
            rolani_psdji2.append(2*np.roll(np.absolute(ft2) ** 2 / sum(Hannovo_okno_trenutni_shear12)**2, N//2))
            valovna_stevila.append(np.linspace(-N / (2 * L), N / (2 * L), N, endpoint=False))   # Append wavenumbers

        vsi_tlaki_po_predalih.append(tlaki_po_predalih.copy())
        vsi_shear1_po_predalih.append(shear1_po_predalih.copy())
        vsi_shear2_po_predalih.append(shear2_po_predalih.copy())
        vsi_centri_predalov.append(centri_predalov.copy())
        vsi_shear1_Hann.append(shear1_Hann.copy())
        vsi_shear2_Hann.append(shear1_Hann.copy())
        vsi_rolani_psdji1.append(rolani_psdji1.copy())
        vsi_rolani_psdji2.append(rolani_psdji2.copy())
        vsa_valovna_stevila.append(valovna_stevila)

    vrni_me_sem = os.getcwd()
    os.chdir(mapa_s_skriptami)
    vsi_posamezne_kolicine = Branje_datotek(izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, tipska_mapa='epsilon', postaja=postaja, ruta=ruta, vrni_kolicine_brez_povprecenja=True)
    os.chdir(vrni_me_sem)
    vsi_tlak_Thermdiss = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['Thermdiss']]for i in range(len(vsi_posamezne_kolicine))]
    vsi_tlak_Tempcor_in_sal = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['Tempcor'], vsi_posamezne_kolicine[i]['sal']] for i in range(len(vsi_posamezne_kolicine))]
    vsi_tlak_epsilon_shear = None
    vsi_tlak_epsilon_shear1_epsilon_shear2_peps = None
    vsi_tlak_N_N2 = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['N'], vsi_posamezne_kolicine[i]['N^2']] for i in range(len(vsi_posamezne_kolicine))]

    try:
        vsi_tlak_epsilon_shear = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['epsilon']]for i in range(len(vsi_posamezne_kolicine))]
        vsi_tlak_epsilon_shear1_epsilon_shear2_peps = [[vsi_posamezne_kolicine[i]['Press'], vsi_posamezne_kolicine[i]['epsilon1'], vsi_posamezne_kolicine[i]['epsilon2'], vsi_posamezne_kolicine[i]['peps']]for i in range(len(vsi_posamezne_kolicine))]
    except:
        pass


    os.chdir(mapa_s_skriptami)  # We return to our scripts
    return {'all_rolled_psd1':vsi_rolani_psdji1, 'all_rolled_psd2':vsi_rolani_psdji2, 'all_wavenumbers':vsa_valovna_stevila, 'all_pressure_centers':vsi_centri_predalov, 'all_pressure_Thermdiss':vsi_tlak_Thermdiss, 'all_pressure_epsilon_shear':vsi_tlak_epsilon_shear, 'all_pressure_epsilon_shear1_epsilon_shear2_peps':vsi_tlak_epsilon_shear1_epsilon_shear2_peps, 'all_pressure_Tempcor_sal':vsi_tlak_Tempcor_in_sal, 'all_shear1':vsi_shear1_po_predalih, 'all_shear2':vsi_shear2_po_predalih, 'all_pressure_N_N2':vsi_tlak_N_N2}



def Nasmith_fit_least_squares_log10(psd, wavenumbers, temp_mean, critical_wavenumber=np.infty, lowest_wavenumber=0):
    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual

    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= critical_wavenumber and wn >= lowest_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])
    log10_proper_psd = np.log10(proper_psd)

    def log10_Nasmyth(k, log10epsilon):
        epsilon = 10**log10epsilon
        x = k * (nu**3 / epsilon)**0.25
        return np.log10(8.05 * x**(1/3) / (1 + (20.6 * x)**3.715) * (epsilon**3/nu)**0.25)

    popt, pcov = curve_fit(log10_Nasmyth, proper_wavenumbers[:], log10_proper_psd[:], bounds=([-13], [-1]))
    log10epsilon = popt[0]
    log10epsilon_sigma = np.sqrt(pcov[0][0])

    return log10epsilon, log10epsilon_sigma

def shear_fit_k_1_3(psd, wavenumbers, temp_mean, critical_wavenumber=np.infty):
    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual

    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= critical_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])
    log10_proper_psd = np.log10(proper_psd)

    def log10_shear_k_1_3(k, log10epsilon):
        epsilon = 10**log10epsilon
        x = k * (nu**3 / epsilon)**0.25
        return np.log10(8.05 * x**(1/3) * (epsilon**3/nu)**0.25)

    popt, pcov = curve_fit(log10_shear_k_1_3, proper_wavenumbers[:], log10_proper_psd[:], bounds=([-13], [-1]))
    log10epsilon = popt[0]
    log10epsilon_sigma = np.sqrt(pcov[0][0])

    return log10epsilon, log10epsilon_sigma




input("\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nWARNING: the procedure of merging data is about to begin. It won't work, if you don't give it your own data.\nTo save the results of this procedure, uncomment (and perhaps slightly change) the two pickle.dump lines at the end of this file.\nPress ENTER to continue with the process, otherwise kill the running program.\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

# MERGING DATA WITH COMMON ROUTES AND STATIONS
# The next section of this file shows, how data is merged. Since the proper files were too large to be all uploaded, one
# can only run this part of the code with their own data (of course, you should specify your own stations (postaja) in
# the most outer loop and routes (ruta) in the second most outer loop). The results from this part of the code are,
# however, saved in ''shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji' + datumska_mapa + '.p' and
# 'shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m' + datumska_mapa + '.p', so other scripts can run smoothly.

# merged outputs with original pressure step
outputs = {'epsilon_shear':[], 'useful_PSD1':[], 'useful_PSD2':[], 'peps':[], 'tempcor':[], 'wavenumbers':[], 'st_zdruzenih':[], 'postaja_ruta_globina':[], 'epsilon_shear_1_MSSpro':[], 'epsilon_shear_2_MSSpro':[], 'N2':[]}
# tempcor ... temperature from final TOB files, st_zdruzenih ... number of merged measurements, postaja_ruta_globina ... station-route-depth, N2 ... squared buoyancy fraquency, all epsilon values are from MSSpro!
# useful PSD is PSD in range 6-35 cpm (for MSS90 probe)

# merged outputs with pressure step = 1dbar (outputs with original pressure step have their pressure centers rounded to 1dbar)
outputs_na_1m = {'epsilon_shear':[], 'useful_PSD1':[], 'useful_PSD2':[], 'peps':[], 'tempcor':[], 'wavenumbers':[], 'st_zdruzenih':[], 'postaja_ruta_globina':[], 'epsilon_shear_1_MSSpro':[], 'epsilon_shear_2_MSSpro':[], 'N2':[]}
# quantities are same as for original pressure step
for postaja in ['LK01','LK02','LK03','LK04','LK05','LK06','LK07']:
    for ruta in [f'RUTA_0{i}' for i in range(1, 10)] + [f'RUTA_1{i}' for i in range(0, 10)]:
        generirano = False
        try:
            generate_epsilon_shear = PSD_shear(izvirna_mapa, datumska_mapa, postaja, ruta)
            generirano = True
        except: # The data might be corrupt, which would lead to an error, but we don't want our entire loop to crash
            pass
        if generirano:
            rolled_psds1 = {}   # rolled PSD of shear1
            rolled_psds2 = {}   # rolled PSD of shear2
            tempcori = {}   # temperatures Tempcor from TOB files of "epsilon" directory
            pepsi = {}  # pseudo epsilons from MSSpro
            epsilon_shears = {} # epsilons from MSSpro
            epsilon_shears1 = {}    # epsilons from shear sensor 1 from MSSpro
            epsilon_shears2 = {}    # epsilons from shear sensor 1 from MSSpro
            N2s = {}    # buoyancy frequencies from MSSpro
            rolled_psds1_1m = {}    # same quantities, just for outputs_na_1m
            rolled_psds2_1m = {}    # same ...
            tempcori_1m = {}    # same ...
            pepsi_1m = {}   # same ...
            epsilon_shears_1m = {}  # same ...
            epsilon_shears1_1m = {} # same ...
            epsilon_shears2_1m = {} # same ...
            N2s_1m = {} # same ...
            all_rolled_psd1 = generate_epsilon_shear['all_rolled_psd1'] # outputs from PSD_shear()
            all_rolled_psd2 = generate_epsilon_shear['all_rolled_psd2']
            all_wavenumbers = generate_epsilon_shear['all_wavenumbers']
            all_pressure_centers = generate_epsilon_shear['all_pressure_centers']
            all_pressure_Thermdiss = generate_epsilon_shear['all_pressure_Thermdiss']
            all_pressure_epsilon_shear = generate_epsilon_shear['all_pressure_epsilon_shear']
            all_pressure_epsilon_shear1_epsilon_shear2_peps = generate_epsilon_shear[
                'all_pressure_epsilon_shear1_epsilon_shear2_peps']
            all_pressure_Tempcor_sal = generate_epsilon_shear['all_pressure_Tempcor_sal']
            all_pressure_N_N2 = generate_epsilon_shear['all_pressure_N_N2']
            all_shear1 = generate_epsilon_shear['all_shear1']
            all_shear2 = generate_epsilon_shear['all_shear2']

            min_max_pressure_centers = False
            try:
                min_pressure_center, max_pressure_center = min([min(i) for i in all_pressure_centers]), max([max(i) for i in all_pressure_centers])
                min_max_pressure_centers = True
            except:
                pass
            if min_max_pressure_centers:
                pressure_center_step = all_pressure_centers[0][1] - all_pressure_centers[0][0]
                all_possible_pressure_centers = np.arange(start=min_pressure_center,
                                                          stop=max_pressure_center + pressure_center_step,
                                                          step=pressure_center_step)

                # Initiate dictionaries
                grouped_rolled_psd1 = {str(pc): [] for pc in all_possible_pressure_centers}
                grouped_rolled_psd2 = {str(pc): [] for pc in all_possible_pressure_centers}
                grouped_wavenumbers = {str(pc): [] for pc in all_possible_pressure_centers}
                grouped_shear1 = {str(pc): [] for pc in all_possible_pressure_centers}
                grouped_shear2 = {str(pc): [] for pc in all_possible_pressure_centers}

                for iprofil in range(len(all_pressure_centers)):    # Append stuff to dictionaries' lists
                    if all_pressure_centers[iprofil][0] > min_pressure_center:
                        for j in all_possible_pressure_centers:
                            if j < all_pressure_centers[iprofil][0]:
                                grouped_rolled_psd1[str(j)].append([])
                                grouped_rolled_psd2[str(j)].append([])
                                grouped_wavenumbers[str(j)].append([])
                                grouped_shear1[str(j)].append([])
                                grouped_shear2[str(j)].append([])
                            else:
                                break
                    for iglobina in range(len(all_pressure_centers[iprofil])):
                        grouped_rolled_psd1[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_rolled_psd1[iprofil][iglobina])
                        grouped_rolled_psd2[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_rolled_psd2[iprofil][iglobina])
                        grouped_wavenumbers[str(all_pressure_centers[iprofil][iglobina])].append(
                            all_wavenumbers[iprofil][iglobina])
                        grouped_shear1[str(all_pressure_centers[iprofil][iglobina])].append(all_shear1[iprofil][iglobina])
                        grouped_shear2[str(all_pressure_centers[iprofil][iglobina])].append(all_shear2[iprofil][iglobina])
                    if all_pressure_centers[iprofil][-1] < max_pressure_center:
                        for j in all_possible_pressure_centers:
                            if j > all_pressure_centers[iprofil][-1]:
                                grouped_rolled_psd1[str(j)].append([])
                                grouped_rolled_psd2[str(j)].append([])
                                grouped_wavenumbers[str(j)].append([])
                                grouped_shear1[str(j)].append([])
                                grouped_shear2[str(j)].append([])

                for cast_no in range(len(all_pressure_epsilon_shear)):  # Go through all casts (usually 5)
                    for i in range(len(all_pressure_epsilon_shear[cast_no][0])):
                        # depth of pressure center
                        globina = all_pressure_epsilon_shear[cast_no][0][i]
                        # depth rounded to 1dbar
                        globina_1m = round(globina + 0.001, 0)  # + 0.001 due to round error (e.g. 2.5 is storred as 2.4999999999 => + 0.001 it becomes 2.50099999, which rounds to 3.0

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
                        elif len(grouped_rolled_psd1[str(globina)][cast_no]) > 0:   # This condition is not fulfilled if
                                                                                    # the combination station-route
                                                                                    # doesn't exist
                            epsilon_shear = all_pressure_epsilon_shear[cast_no][1][i]
                            epsilon_shear1 = all_pressure_epsilon_shear1_epsilon_shear2_peps[cast_no][1][i]
                            epsilon_shear2 = all_pressure_epsilon_shear1_epsilon_shear2_peps[cast_no][2][i]
                            peps = all_pressure_epsilon_shear1_epsilon_shear2_peps[cast_no][3][i]

                            nu = 1.702747 - 0.05126103*tale_Tempcor_in_Thermdiss_in_N2[0] + 0.0005918645*tale_Tempcor_in_Thermdiss_in_N2[0]**2    # MSSpro user manual
                            nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual, nu is viscosity

                            if str(globina) in epsilon_shears.keys():   # This depth already has a key
                                epsilon_shears[str(globina)].append(epsilon_shear)
                                rolled_psds1[str(globina)].append(grouped_rolled_psd1[str(globina)][cast_no])
                                rolled_psds2[str(globina)].append(grouped_rolled_psd2[str(globina)][cast_no])
                                tempcori[str(globina)].append(tale_Tempcor_in_Thermdiss_in_N2[0])
                                pepsi[str(globina)].append(peps)
                                epsilon_shears1[str(globina)].append(epsilon_shear1)
                                epsilon_shears2[str(globina)].append(epsilon_shear2)
                                N2s[(str(globina))].append(tale_Tempcor_in_Thermdiss_in_N2[2])
                            else:   # Iniciate a key and make its value to be a list
                                epsilon_shears[str(globina)] = [epsilon_shear]
                                rolled_psds1[str(globina)] = [grouped_rolled_psd1[str(globina)][cast_no]]
                                rolled_psds2[str(globina)] = [grouped_rolled_psd2[str(globina)][cast_no]]
                                tempcori[str(globina)] = [tale_Tempcor_in_Thermdiss_in_N2[0]]
                                pepsi[str(globina)] = [peps]
                                epsilon_shears1[str(globina)] = [epsilon_shear1]
                                epsilon_shears2[str(globina)] = [epsilon_shear2]
                                N2s[str(globina)] = [tale_Tempcor_in_Thermdiss_in_N2[2]]

                            if str(globina_1m) in epsilon_shears_1m.keys():
                                epsilon_shears_1m[str(globina_1m)].append(epsilon_shear)
                                rolled_psds1_1m[str(globina_1m)].append(grouped_rolled_psd1[str(globina)][cast_no])
                                rolled_psds2_1m[str(globina_1m)].append(grouped_rolled_psd2[str(globina)][cast_no])
                                tempcori_1m[str(globina_1m)].append(tale_Tempcor_in_Thermdiss_in_N2[0])
                                pepsi_1m[str(globina_1m)].append(peps)
                                epsilon_shears1_1m[str(globina_1m)].append(epsilon_shear1)
                                epsilon_shears2_1m[str(globina_1m)].append(epsilon_shear2)
                                N2s_1m[(str(globina_1m))].append(tale_Tempcor_in_Thermdiss_in_N2[2])
                            else:
                                epsilon_shears_1m[str(globina_1m)] = [epsilon_shear]
                                rolled_psds1_1m[str(globina_1m)] = [grouped_rolled_psd1[str(globina)][cast_no]]
                                rolled_psds2_1m[str(globina_1m)] = [grouped_rolled_psd2[str(globina)][cast_no]]
                                tempcori_1m[str(globina_1m)] = [tale_Tempcor_in_Thermdiss_in_N2[0]]
                                pepsi_1m[str(globina_1m)] = [peps]
                                epsilon_shears1_1m[str(globina_1m)] = [epsilon_shear1]
                                epsilon_shears2_1m[str(globina_1m)] = [epsilon_shear2]
                                N2s_1m[str(globina_1m)] = [tale_Tempcor_in_Thermdiss_in_N2[2]]

            print(f'Depths with at least one measurement that fulfilled all conditions for station {postaja} and route {ruta}', epsilon_shears.keys())
            for globina in epsilon_shears.keys():
                if len(epsilon_shears[globina]) < 1:    # If you want larger number for better statistics, increase
                    print(f'For station {postaja}, route {ruta} and depth {globina} we lack statistical significance')
                else:
                    povprecni_epsilon_shear = gmean(epsilon_shears[globina])
                    povprecni_epsilon_shear1 = gmean(epsilon_shears1[globina])
                    povprecni_epsilon_shear2 = gmean(epsilon_shears2[globina])
                    #NEXT LINE ONLY WORKS WELL IF THE INTERVAL WIDTH IS 1dbar
                    max_val_st = 40 # maximum wavenumber to be stored
                    samo_uporabni_psds1 = [rolled_psds1[globina][i][len(rolled_psds1[globina][i])//2:len(rolled_psds1[globina][i])//2+max_val_st] for i in range(len(rolled_psds1[globina]))]   # Ker so razlicnih dolzin, ne moremo racunati povprecij, zato gremo ze zdaj gledati samo uporabne psd.
                    samo_uporabni_psds2 = [rolled_psds2[globina][i][len(rolled_psds2[globina][i])//2:len(rolled_psds2[globina][i])//2+max_val_st] for i in range(len(rolled_psds2[globina]))]   # Ker so razlicnih dolzin, ne moremo racunati povprecij, zato gremo ze zdaj gledati samo uporabne psd.
                    povprecni_rolled_psd1 = gmean(samo_uporabni_psds1, axis=0)
                    povprecni_rolled_psd2 = gmean(samo_uporabni_psds2, axis=0)
                    povprecni_tempcor = np.mean(tempcori[globina])
                    povprecni_peps = gmean(pepsi[globina])

                    valovna_stevila = np.linspace(1, max_val_st, max_val_st)  # ONLY WORKS PROPERLY IF INTERVAL WIDTH IS 1dbar
                    outputs['epsilon_shear'].append(np.log10(povprecni_epsilon_shear))
                    outputs['epsilon_shear_1_MSSpro'].append(np.log10(povprecni_epsilon_shear1))
                    outputs['epsilon_shear_2_MSSpro'].append(np.log10(povprecni_epsilon_shear2))
                    outputs['peps'].append(np.log10(povprecni_peps))
                    outputs['useful_PSD1'].append(povprecni_rolled_psd1)
                    outputs['useful_PSD2'].append(povprecni_rolled_psd2)
                    outputs['tempcor'].append(povprecni_tempcor)
                    outputs['wavenumbers'].append(valovna_stevila)
                    outputs['st_zdruzenih'].append(len(epsilon_shears[globina]))
                    outputs['postaja_ruta_globina'].append(postaja + '_' + ruta + '_' + globina)
                    outputs['N2'].append(np.mean(N2s[globina]))
                    print('success')
            for globina in epsilon_shears_1m.keys():
                if len(epsilon_shears_1m[globina]) < 1:    # If you want larger number for better statistics, increase
                    print(f'For station {postaja}, route {ruta} and depth {globina} we lack statistical significance')
                else:
                    povprecni_epsilon_shear = gmean(epsilon_shears_1m[globina])
                    povprecni_epsilon_shear1 = gmean(epsilon_shears1_1m[globina])
                    povprecni_epsilon_shear2 = gmean(epsilon_shears2_1m[globina])
                    #NEXT LINE ONLY WORKS WELL IF THE INTERVAL WIDTH IS 1dbar
                    max_val_st = 40
                    samo_uporabni_psds1 = [rolled_psds1_1m[globina][i][len(rolled_psds1_1m[globina][i])//2:len(rolled_psds1_1m[globina][i])//2+max_val_st] for i in range(len(rolled_psds1_1m[globina]))]   # Ker so razlicnih dolzin, ne moremo racunati povprecij, zato gremo ze zdaj gledati samo uporabne psd.
                    samo_uporabni_psds2 = [rolled_psds2_1m[globina][i][len(rolled_psds2_1m[globina][i])//2:len(rolled_psds2_1m[globina][i])//2+max_val_st] for i in range(len(rolled_psds2_1m[globina]))]   # Ker so razlicnih dolzin, ne moremo racunati povprecij, zato gremo ze zdaj gledati samo uporabne psd.
                    povprecni_rolled_psd1 = gmean(samo_uporabni_psds1, axis=0)
                    povprecni_rolled_psd2 = gmean(samo_uporabni_psds2, axis=0)
                    povprecni_tempcor = np.mean(tempcori_1m[globina])
                    povprecni_peps = gmean(pepsi_1m[globina])

                    valovna_stevila = np.linspace(1, max_val_st, max_val_st)  # ONLY WORKS PROPERLY IF INTERVAL WIDTH IS 1dbar
                    outputs_na_1m['epsilon_shear'].append(np.log10(povprecni_epsilon_shear))
                    outputs_na_1m['epsilon_shear_1_MSSpro'].append(np.log10(povprecni_epsilon_shear1))
                    outputs_na_1m['epsilon_shear_2_MSSpro'].append(np.log10(povprecni_epsilon_shear2))
                    outputs_na_1m['peps'].append(np.log10(povprecni_peps))
                    outputs_na_1m['useful_PSD1'].append(povprecni_rolled_psd1)
                    outputs_na_1m['useful_PSD2'].append(povprecni_rolled_psd2)
                    outputs_na_1m['tempcor'].append(povprecni_tempcor)
                    outputs_na_1m['wavenumbers'].append(valovna_stevila)
                    outputs_na_1m['st_zdruzenih'].append(len(epsilon_shears_1m[globina]))
                    outputs_na_1m['postaja_ruta_globina'].append(postaja + '_' + ruta + '_' + globina)
                    outputs_na_1m['N2'].append(np.mean(N2s_1m[globina]))
                    print('success')

os.chdir(mapa_s_skriptami)  # Return to scripts (otherwise you'll dump pickle files somewhere else)
#pickle.dump(outputs, open('shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji' + datumska_mapa + '.p', 'wb'))
#pickle.dump(outputs_na_1m, open('shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m' + datumska_mapa + '.p', 'wb'))
print('DONE')



# THE FINAL PART FOCUSES ON THE DETERMINATION OF DEGREES OF FREEDOM OF PSD
# We determine the degrees of freedom (DoF) using the integral fit, which is unbiased (maximum likelihood estimate
# requires us to tell it the expected DoF to make a fit). To plot the output histograms, set plot_histograms=True
plot_histograms = False

def fit_Nasmyth_Lueck(psd, wavenumbers, temp_mean, min_wavenumber=6, max_wavenumber=35, n_iter=10):
    '''All wavenumbers should be in cpm! Iteration usually converges after 3-5 iteration, n_iter=10 is just to be safe.'''
    nu = 1.702747 - 0.05126103*temp_mean + 0.0005918645*temp_mean**2    # MSSpro user manual
    nu *= 10**-6    # m**2*s**-1, MSSpro user manual

    proper_wavenumbers = [wn for wn in wavenumbers if wn > 0 and wn <= max_wavenumber and wn >= min_wavenumber]
    proper_psd = np.array([psd[iwn] for iwn in range(len(wavenumbers)) if wavenumbers[iwn] in proper_wavenumbers])

    trapz_integration = scipy.integrate.trapz(y=proper_psd, x=proper_wavenumbers)
    epsilon_trial = 15/2 * nu * trapz_integration

    def fraction(epsilon, nu=nu, k1=min_wavenumber, k2=max_wavenumber):
        def int_psi(x):
            '''Lueck 2016, Eq. (6)'''
            return np.tanh(48 * x ** (4 / 3)) - 2.9 * x ** (4 / 3) * np.exp(-22.3 * x ** (4 / 3))

        x1 = k1 * (nu ** 3 / epsilon) ** 0.25
        x2 = k2 * (nu ** 3 / epsilon) ** 0.25
        return int_psi(x2) - int_psi(x1)

    epsilon_iterated = epsilon_trial
    #print(epsilon_iterated)
    for iiteration in range(n_iter):
        frac = fraction(epsilon_iterated)
        epsilon_iterated = epsilon_trial / frac
        #print(epsilon_iterated, frac)
    #print()
    return np.log10(epsilon_iterated)

# Load data
merged_to_1m = True # True to look at the data merged to 1m, False to look at the data merged to default pressure step
if merged_to_1m:
    outputs = pickle.load(open('shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m' + datumska_mapa + '.p', 'rb'))
else:
    outputs = pickle.load(open('shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji' + datumska_mapa + '.p', 'rb'))



ioutput = 0
sp_meja = 6
zg_meja = 35 + 1
dofs1 = {i:[] for i in range(1, max(outputs['st_zdruzenih'])+1)}
dofs2 = {i:[] for i in range(1, max(outputs['st_zdruzenih'])+1)}

for ioutput in range(len(outputs['tempcor'])):
    print(ioutput, len(outputs['tempcor']))
    shear_fit1 = fit_Nasmyth_Lueck(psd=outputs['useful_PSD1'][ioutput][sp_meja:zg_meja], wavenumbers=outputs['wavenumbers'][ioutput][sp_meja:zg_meja], temp_mean=outputs['tempcor'][ioutput])
    epsilon_shear_my1 = shear_fit1
    shear_fit2 = fit_Nasmyth_Lueck(psd=outputs['useful_PSD2'][ioutput][sp_meja:zg_meja], wavenumbers=outputs['wavenumbers'][ioutput][sp_meja:zg_meja], temp_mean=outputs['tempcor'][ioutput])
    epsilon_shear_my2 = shear_fit2
    N2 = outputs['N2'][ioutput]
    peps = outputs['peps'][ioutput]
    nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2
    nu *= 10**-6

    def log10_Nasmyth(k, log10epsilon):
        '''Theoretical curve'''
        epsilon = 10 ** log10epsilon
        nu = 1.702747 - 0.05126103 * outputs['tempcor'][ioutput] + 0.0005918645 * outputs['tempcor'][ioutput] ** 2  # MSSpro user manual
        nu *= 10 ** -6  # m**2*s**-1, MSSpro user manual
        x = k * (nu ** 3 / epsilon) ** 0.25
        return np.log10(8.05 * x ** (1 / 3) / (1 + (20.6 * x) ** 3.715) * (epsilon ** 3 / nu) ** 0.25)


    def ugotovi_dof(razmerja):
        '''Returns estimated degrees of freedom'''
        razmerja.sort()
        delezi = [i / (len(razmerja)) for i in range(1, len(razmerja) + 1)] # fractions

        def chi2cdf_x(x, dof):
            return scipy.stats.chi2.cdf(x * dof, dof)

        # Fit to chi2cdf_x
        popt, pcov = curve_fit(chi2cdf_x, razmerja, delezi, bounds=([0.1], [50]))
        pravi_dof = popt[0]


        return pravi_dof


    if 10 ** epsilon_shear_my1 / (nu * N2) > 200 and epsilon_shear_my1:
        generirani_1 = 10 ** log10_Nasmyth(outputs['wavenumbers'][ioutput], epsilon_shear_my1)
        razmerja_1 = np.array(outputs['useful_PSD1'][ioutput])[sp_meja:zg_meja] / generirani_1[sp_meja:zg_meja]
        dof_1 = ugotovi_dof(razmerja_1)
        dofs1[outputs['st_zdruzenih'][ioutput]].append(dof_1)

    if 10 ** epsilon_shear_my2 / (nu * N2) > 200 and epsilon_shear_my2:
        generirani_2 = 10 ** log10_Nasmyth(outputs['wavenumbers'][ioutput], epsilon_shear_my2)
        razmerja_2 = np.array(outputs['useful_PSD2'][ioutput])[sp_meja:zg_meja] / generirani_2[sp_meja:zg_meja]
        razmerja_2.sort()
        dof_2 = ugotovi_dof(razmerja_2)
        dofs2[outputs['st_zdruzenih'][ioutput]].append(dof_2)

if merged_to_1m:
    pickle.dump([dofs1, dofs2], open(f'dofs_na_1m_{datumska_mapa}.p', 'wb'))
else:
    pickle.dump([dofs1, dofs2], open(f'dofs_{datumska_mapa}.p', 'wb'))
print('DONE')


if plot_histograms:
    if merged_to_1m:

        dofs12 = pickle.load(open(f'dofs_na_1m{datumska_mapa}.p', 'rb'))
        dofs1 = dofs12[0]
        dofs2 = dofs12[1]
        print(len(dofs1), len(dofs2))
        labels = ['shear1', 'shear2']

        plt.hist([dofs1[len(dofs1)], dofs2[len(dofs1)]], bins=[i for i in range(40) if i%2==0], stacked=True, label=labels)
        plt.axvline(np.median(dofs1[len(dofs1)] + dofs2[len(dofs1)]), color='k')
        plt.title(f'{len(dofs1)} averaged')#(f'Št. povprečenih: {len(dofs1)}')    # (f'{len(dofs1)} averaged')
        plt.xlabel('Calculated DoF')#('Izračunane DoF')    # ('Calculated DoF')
        plt.ylabel('Occurrence')#(r'Št. pojavitev') #('Occurrence')
        plt.annotate('median = ' + str(round(np.median(dofs1[len(dofs1)] + dofs2[len(dofs1)]), 1)), xy=(0.75, 0.75), xycoords='axes fraction')
        plt.annotate('mean = ' + str(round(np.mean(dofs1[len(dofs1)] + dofs2[len(dofs1)]), 1)), xy=(0.75, 0.65), xycoords='axes fraction')
        plt.annotate('STD = ' + str(round(np.std(dofs1[len(dofs1)] + dofs2[len(dofs1)]), 1)), xy=(0.75, 0.55), xycoords='axes fraction')
        plt.annotate('q10 = ' + str(round(np.quantile(dofs1[len(dofs1)] + dofs2[len(dofs1)], 0.10), 1)), xy=(0.75, 0.45), xycoords='axes fraction')
        plt.annotate('q90 = ' + str(round(np.quantile(dofs1[len(dofs1)] + dofs2[len(dofs1)], 0.90), 1)), xy=(0.75, 0.35), xycoords='axes fraction')
        plt.legend()
        plt.savefig(mapa_za_shranjevanje_grafov + f'calculated_dof_na_1M_ONLY_ISOTROPIC_CASES_{datumska_mapa}')#_si')
        plt.show()
    else:
        dofs12 = pickle.load(open(f'dofs_{datumska_mapa}.p', 'rb'))
        dofs1 = dofs12[0]
        dofs2 = dofs12[1]
        print(len(dofs1), len(dofs2))
        labels = ['shear1', 'shear2']
        plt.figure(figsize=(12, 4))
        for j in range(len(dofs1)):
            plt.subplot(1, len(dofs1), j+1)
            plt.hist([dofs1[j+1], dofs2[j+1]], bins=[i for i in range(20)], stacked=True, label=labels)
            plt.axvline(np.median(dofs1[j+1] + dofs2[j+1]), color='k')
            plt.title(f'Št. povprečenih: {j+1}')    # (f'{j+1} averaged')
            plt.xlabel('Izračunane DoF')    # ('Calculated DoF')
            plt.ylabel(r'Št. pojavitev') #('Occurrence')
            plt.annotate('median = ' + str(round(np.median(dofs1[j+1] + dofs2[j+1]), 1)), xy=(0.4, 0.75), xycoords='axes fraction')
            plt.annotate('mean = ' + str(round(np.mean(dofs1[len(dofs1)] + dofs2[len(dofs1)]), 1)), xy=(0.4, 0.65), xycoords='axes fraction')
            plt.annotate('STD = ' + str(round(np.std(dofs1[len(dofs1)] + dofs2[len(dofs1)]), 1)), xy=(0.4, 0.55), xycoords='axes fraction')
            plt.annotate('q10 = ' + str(round(np.quantile(dofs1[len(dofs1)] + dofs2[len(dofs1)], 0.10), 1)), xy=(0.4, 0.45), xycoords='axes fraction')
            plt.annotate('q90 = ' + str(round(np.quantile(dofs1[len(dofs1)] + dofs2[len(dofs1)], 0.90), 1)), xy=(0.4, 0.35), xycoords='axes fraction')
            plt.legend()
        plt.tight_layout()
        plt.savefig(mapa_za_shranjevanje_grafov + f'calculated_dof_VSI_ONLY_ISOTROPIC_CASES_{datumska_mapa}_si')
        plt.show()
