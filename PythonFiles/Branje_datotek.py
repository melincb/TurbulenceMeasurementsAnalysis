def Branje_datotek(izvirna_mapa, datumska_mapa, tipska_mapa, postaja, ruta, tlacni_korak=0.25, vrni_kolicine_brez_povprecenja=False, ura=False):
    '''Reads TOB files (i.e. the outputs from MSSpro).

    izvirna_mapa ... directory of origin of data (in our case always "..")
    datumska_mapa ... determines date of measurements ("all_2008_08_21" or "all_2008_11_18" or "all_2009_04_20")
    tipska_mapa ... final directory of TOB files ("epsilon", "shear" or "cutted") (Beware: TOB files from different
    directories contain different quantities!)
    postaja ... station ("LK01", "LK02", ..., "LK07")
    ruta ... route ("RUTA_01", "RUTA_02", ...)
    tlacni_korak ... pressure step for merging data [dbar]
    vrni_kolicine_brez_povprecenja ... look below
    ura ... look below

    This function merges data from all files in the directory which correspond to the same station and route (if
    ura!=False (e.g. ura="0709") it looks at just one specific file with time=ura, station=postaja and route=ruta).

    If vrni_kolicine_brez_povprecenja=True it returns a list of dictionaries, where i-th dictionary corresponds to i-th
    measurement.
    If vrni_kolicine_brez_povprecenja=False it returns three separate dictionaries. The first dictionary has means of
    each quantity at each depth (all measurements included in averaging), the second one has standard deviations and the
     third one has medians. The pressure step between adjacent depths is tlacni_korak.
    The values of quantities which are used to calculate means, standard deviations and medians are never logarithmic!
    This is why the values of those quantities that are given by MSSpro in logarithmic values (i.e. epsilon, epsilon1,
    epsilon2, peps, Thermdiss) are transformed as 10**value before calculating mean, standard deviation and median. The
    final outputs are not transformed back to logarithms!
    '''
    import numpy as np
    import glob, os


    mapa_s_skriptami = os.getcwd()  # directory with scripts
    print(os.getcwd())
    bs = '\ '[0]    # an awkward way to tell Python that our symbol is backslash
    mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa   # directory with data
    print(mapa_s_podatki)
    os.chdir(mapa_s_podatki)    # Go to directory with data

    vse_kolicine = {}   # all quantities
    vse_kolicine_locene_po_izvornih_datotekah = []  # all quantities separated by file of origin
    enote_vseh_kolicin = {} # units of all quantities
    if ura:
        kje_iscem = '*_' + ura + '_' + postaja + '_' + ruta + '.tob'    # files of interest
        vrni_kolicine_brez_povprecenja = True   # If ura=True you only look at one file
    else:
        kje_iscem = '*_' + postaja + '_' + ruta + '.tob'    # files of interest

    for datoteka in glob.glob(kje_iscem):   # Go through all files of interest
        vse_kolicine_locene_po_izvornih_datotekah.append({})    # a dictionary for each cast
        with open(datoteka, 'r') as datoteka:   # Examine the file
            datasets, units = False, False
            tokratne_kolicine = []  # this cast's quantities
            count = 0   # All TOB files have the same way of ending their header. The penultimate line contains Datasets
                        # and the ultimate line contains Units. So, until we reach 'Datasets', we reset count=0, when
                        # count=1 the line states 'Units', when count=2 the line is just ';' and when count>2 we are
                        # finally dealling with the measured data.
            for vrstica in datoteka:    # Go through each line of a file
                if not datasets:
                    count = 0
                    try:
                        if vrstica.split()[1] == 'Datasets':
                            datasets = True
                            for kolicina in vrstica.split()[2:]:    # Remember each quantity
                                tokratne_kolicine.append(kolicina)
                                vse_kolicine_locene_po_izvornih_datotekah[-1][kolicina] = []
                                if kolicina not in vse_kolicine.keys(): # Check if vse_kolicine already has this quantity
                                    if len(vse_kolicine) == 0:
                                        vse_kolicine[kolicina] = []
                                        enote_vseh_kolicin[kolicina] = []
                                    else:   # This never happens in my data, however if someone else uses their data, it
                                            # could be that not all TOB files contain the same quantities. In such case,
                                            # we set all the values of the missing quantity in previous files to np.nan
                                            # This approach hasn't been properly tested, use it at your own risk! If all
                                            # your TOB files contain the same quantities, this "else" is never executed,
                                            # so no worries in this case.
                                        vse_kolicine[kolicina] = [np.nan for i in range(len(vse_kolicine[list(vse_kolicine.keys())[0]]))]
                                        enote_vseh_kolicin[kolicina] = []
                    except: # The line might be empty or it cotains just one element
                        pass
                elif count == 1: # This line contains units
                    units = True
                    count_units = 0
                    for enota in vrstica.split()[1:]:
                        enote_vseh_kolicin[tokratne_kolicine[count_units]].append(enota)
                        count_units += 1
                elif count > 2: # If count=2 line only contains ';'
                    count_values = 0    # e.g. tokratne_kolicine[0] = 'Press', tokratne_kolicine[1] = 'Tempcor', etc.
                    for vrednost in vrstica.split()[1:]:    # For each value in a line
                        if tokratne_kolicine[count_values] in ('epsilon', 'epsilon1', 'epsilon2', 'peps', 'Thermdiss'):
                            # Decimal logaritms go away
                            vse_kolicine[tokratne_kolicine[count_values]].append(10**float(vrednost))   # Append value
                            vse_kolicine_locene_po_izvornih_datotekah[-1][tokratne_kolicine[count_values]].append(10**float(vrednost))
                        else:
                            vse_kolicine[tokratne_kolicine[count_values]].append(float(vrednost))   # Append value
                            vse_kolicine_locene_po_izvornih_datotekah[-1][tokratne_kolicine[count_values]].append(float(vrednost))
                        count_values += 1
                count += 1 # While datasets=False count will always be reset to 0

    if vrni_kolicine_brez_povprecenja:  # Return values from separate TOB files!
        os.chdir(mapa_s_skriptami)  # Go back to scripts
        return vse_kolicine_locene_po_izvornih_datotekah


    # __________________________________
    # SORTING BY PRESSURE
    # WORKS ONLY IN PYTHON 3.7+ !!! (in earlier versions dictionaries' keys' order was random)
    # sorting style: https://thispointer.com/sorting-2d-numpy-array-by-column-or-row-in-python/ (24.7.2022)
    kolicine = list(vse_kolicine.keys())
    nesortirane_vrednosti_kolicin = np.array(list(vse_kolicine.values()))   # unsorted values of quantities
    kateri_po_vrsti_je_tlak = kolicine.index('Press')   # index of 'Press'

    sortirane_vrednosti_kolicin = nesortirane_vrednosti_kolicin[:, nesortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak].argsort()]    # sorted values of quantities


    # __________________________________
    # ROUNDING PRESSURE TO tlacni_korak [dbar]

    def zaokrozi_na_dp(p, dp=0.25):
        '''Round pressure, so the step is dp.'''
        return dp * np.round(p / dp)

    sortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak] = zaokrozi_na_dp(sortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak], dp=tlacni_korak)

    count_p_elements = 1
    count_p_elements_pri_tem_p = 0
    trenutni_p = sortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak, 0]
    sortirane_vrednosti_kolicin_podmatrika = [[i] for i in sortirane_vrednosti_kolicin[:, 0]]
    vse_podmatrike_sortiranih_vrednosti_kolicin = []
    while count_p_elements < np.shape(sortirane_vrednosti_kolicin)[1]:
        if trenutni_p == sortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak, count_p_elements]:
            sortirane_vrednosti_kolicin_podmatrika = sortirane_vrednosti_kolicin[:, count_p_elements_pri_tem_p:count_p_elements + 1]
        else:
            if np.shape(sortirane_vrednosti_kolicin_podmatrika)[1] > 2:
                # In other case statistics is poor (can only happen at the top or bottom of water column)
                vse_podmatrike_sortiranih_vrednosti_kolicin.append(sortirane_vrednosti_kolicin_podmatrika)
            trenutni_p = sortirane_vrednosti_kolicin[kateri_po_vrsti_je_tlak, count_p_elements]
            count_p_elements_pri_tem_p = count_p_elements
            sortirane_vrednosti_kolicin_podmatrika = [[i] for i in sortirane_vrednosti_kolicin[:, count_p_elements]]
        count_p_elements += 1
    if np.shape(sortirane_vrednosti_kolicin_podmatrika)[1] > 2:
        # In other case statistics is poor (can only happen at the top or bottom of water column)
        vse_podmatrike_sortiranih_vrednosti_kolicin.append(sortirane_vrednosti_kolicin_podmatrika)


    vse_kolicine_mean_transposed = []  # means at given pressure
    vse_kolicine_std_transposed = []   # standard deviations at given pressure
    vse_kolicine_median_transposed = [] # medians at given pressure
    for podmatrika in vse_podmatrike_sortiranih_vrednosti_kolicin:
        vse_kolicine_mean_transposed.append(np.mean(podmatrika, axis=1))
        vse_kolicine_std_transposed.append(np.std(podmatrika, axis=1, ddof=1))  # ddof=1, so denominator is n-1
        vse_kolicine_median_transposed.append(np.median(podmatrika, axis=1))
    vse_kolicine_mean = np.transpose(vse_kolicine_mean_transposed)
    vse_kolicine_std = np.transpose(vse_kolicine_std_transposed)
    vse_kolicine_median = np.transpose(vse_kolicine_median_transposed)

    slovar_kolicina_mean = {}
    slovar_kolicina_std = {}
    slovar_kolicina_median = {}

    for ikolicina in range(len(kolicine)):
        slovar_kolicina_mean[kolicine[ikolicina]] = vse_kolicine_mean[ikolicina,:]
        slovar_kolicina_std[kolicine[ikolicina]] = vse_kolicine_std[ikolicina,:]
        slovar_kolicina_median[kolicine[ikolicina]] = vse_kolicine_median[ikolicina,:]

    os.chdir(mapa_s_skriptami)  # Back to scripts
    return slovar_kolicina_mean, slovar_kolicina_std, slovar_kolicina_median



def iskanje_minimalne_in_maksimalne_vrednosti_nekih_kolicin_za_dolocen_datum_pri_vseh_meritvah(kolicine, izvirna_mapa, datumska_mapa, tipska_mapa):
    '''Finds the min and max value of given quantities on certain day.

    kolicine ... list of quantities (e.g. ['Tempcor, 'sal'])
    izvirna_mapa ... directory of origin of data (in our case always "..")
    datumska_mapa ... determines date of measurements ("all_2008_08_21" or "all_2008_11_18" or "all_2009_04_20")
    tipska_mapa ... final directory of TOB files ("epsilon", "shear" or "cutted") (Beware: TOB files from different
    directories contain different quantities!)
    '''

    import numpy as np
    import glob, os

    mapa_s_skriptami = os.getcwd()  # directory with scripts
    print(mapa_s_skriptami)
    bs = '\ '[0]
    mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa   # directory with data
    os.chdir(mapa_s_podatki)    # Go to data
    print(os.getcwd())

    minimumi_iskanih_kolicin = {kolicina:[] for kolicina in kolicine}   # minima
    maksimumi_iskanih_kolicin = {kolicina:[] for kolicina in kolicine}  # maxima
    skupno_stevilo_vrstic = 0   # total number of rows

    for datoteka in glob.glob('*.tob'):
        # ******************************************************************
        # Begining of the part, which is the same as in Branje_datotek(), except it doesn't contain complications
        # for vrni_kolicine_brez_povprecenja=True.
        # ******************************************************************

        vse_kolicine = {}
        enote_vseh_kolicin = {}
        with open(datoteka, 'r') as datoteka:
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
                                    if len(vse_kolicine) == 0:
                                        vse_kolicine[kolicina] = []
                                        enote_vseh_kolicin[kolicina] = []
                                    else:
                                        vse_kolicine[kolicina] = [np.nan for i in range(len(vse_kolicine[list(vse_kolicine.keys())[0]]))]
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
                    for vrednost in vrstica.split()[1:]:
                        if tokratne_kolicine[count_values] in ('epsilon', 'epsilon1', 'epsilon2', 'peps', 'Thermdiss'):
                            vse_kolicine[tokratne_kolicine[count_values]].append(10**float(vrednost))
                        else:
                            vse_kolicine[tokratne_kolicine[count_values]].append(float(vrednost))
                        count_values += 1
                count += 1

        # ******************************************************************
        # End the part, which is the same as in Branje_datotek().
        # ******************************************************************

        skupno_stevilo_vrstic += count - 2
        for kolicina in kolicine:
            minimumi_iskanih_kolicin[kolicina].append(np.min(vse_kolicine[kolicina]))   # Append minima from this file
            maksimumi_iskanih_kolicin[kolicina].append(np.max(vse_kolicine[kolicina]))  # Append maxima from this file

    print('Total number of rows:', skupno_stevilo_vrstic)
    # Finaly, find global minima and maxima
    dejanski_minimumi_in_maksimumi_iskanih_kolicin = {kolicina:[np.min(minimumi_iskanih_kolicin[kolicina]), np.max(maksimumi_iskanih_kolicin[kolicina])] for kolicina in kolicine}

    os.chdir(mapa_s_skriptami)  # Go back to scripts
    return dejanski_minimumi_in_maksimumi_iskanih_kolicin