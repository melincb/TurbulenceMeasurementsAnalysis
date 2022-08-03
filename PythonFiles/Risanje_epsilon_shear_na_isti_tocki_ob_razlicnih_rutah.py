'''This file plots epsilon on one station throughout the entire 25-hour field campaigns.
WARNING: The interpolation is done very badly, because the temporal distance between two routes can be from 10 minutes
to 3 hours! You can only rely on the areas, which show the true measurements (i.e. where the vertical dashed line is).
The plots from Risanje_epsilon_shear_na_vseh_tockah_ob_isti_ruti.py are much more reliable.'''

from Branje_datotek import Branje_datotek
import numpy as np
import glob, os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.tri as tri
import pickle

izvirna_mapa = '..' # directory with directories of measurements
mapa_za_shranjevanje_grafov =  'Saved_figures\ '[:-1]    # Where it saves the figures
save_figures = False    # To save the plot, set this to True, else False

datumska_mapa = 'all_2008_11_18'    # date of measurements ('all_2008_08_21', 'all_2008_11_18', 'all_2009_04_20)
postaja = 'LK02'    # Select station

# Which epsilon should we plot? Options: 'log_epsilon_MLE_mean', 'log_epsilon_MLE_1', 'log_epsilon_MLE_2',
# 'log_epsilon_IM_mean', 'log_epsilon_IM_1', 'log_epsilon_IM_2'
which_epsilon = 'log_epsilon_MLE_mean'





# THE BASIC MANUAL SETTINGS END HERE


leto = datumska_mapa[6:8]   # year
mesec = datumska_mapa[9:11] # month
dan = datumska_mapa[12:]    # day
tipska_mapa = 'epsilon' # directory, where data from each measurement is stored (only used to determine the time of each
                        # measurement, which isn't stored in f'slovar_epsilon_za_risat_NA_1m_{datumska_mapa}.p'

# Load data
izracuni_epsilon_in_Re_b = pickle.load(open(f'slovar_epsilon_za_risat_NA_1m_{datumska_mapa}.p', 'rb'))
# Pick only data for this station
izracuni_epsilon_in_Re_b_za_to_postajo = [izracuni_epsilon_in_Re_b[postaja_ruta] for postaja_ruta in izracuni_epsilon_in_Re_b.keys() if postaja_ruta[:4] == postaja]    # Predpostavka je, da so rute v pravilnem vrstnem redu!

bs = '\ '[0]
mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa   # directory with data
mapa_s_skriptami = os.getcwd()  # directory with scripts

min_ruta = 1 # Min. route index in data
max_rute = {'all_2008_08_21':14, 'all_2008_11_18':16, 'all_2009_04_20':19}  # Max. route index in data
max_ruta = max_rute[datumska_mapa]    # Max. route index for this date
moja_min_ruta = 1    # Min. station index to be plotted (WARNING - not tested for any other number than 1)
moja_max_ruta = max_ruta    # Min. station index to be plotted (WARNING - not tested for any other number than max_ruta)



# FIND TIME OF EACH STATION'S MEASUREMENTS

casi_rut = []   # times of measurements (i.e. times of routes at the chosen station)

for iruta in range(moja_min_ruta, moja_max_ruta + 1):
    ruta = 'RUTA_{:02d}'.format(iruta)
    try:    # Check if there was a measurement on this station at this route (if so, then len(tokratne_mean) > 0)
        tokratne_mean, tokratne_std, tokratne_median = Branje_datotek(izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, \
                                                                        tipska_mapa=tipska_mapa, postaja=postaja, ruta=ruta)
        if len(tokratne_mean) > 0:
            os.chdir(mapa_s_podatki)    # Go to data
            datoteke_ruta = glob.glob('*_' + postaja + '_' + ruta + '.tob') # Select all files with this combination
                                                                            # station-route
            os.chdir(mapa_s_skriptami)  # Go back to scripts
            srednja_meritev_ruta = datoteke_ruta[len(datoteke_ruta) // 2]      # Time of the medial measurement

            # Figure out time of the srednja_meritev_postaja_ruta from its name
            for iznak in range(0, len(srednja_meritev_ruta)):
                if srednja_meritev_ruta[iznak] == '_':
                    datetime_srednja = datetime(year=int('20' + srednja_meritev_ruta[iznak - 6:iznak - 4]), \
                                             month=int(srednja_meritev_ruta[iznak - 4:iznak - 2]), \
                                             day=int(srednja_meritev_ruta[iznak - 2:iznak]), \
                                             hour=int(srednja_meritev_ruta[iznak + 1:iznak + 3]), \
                                             minute=int(srednja_meritev_ruta[iznak + 3:iznak + 5]))
                    print(datetime_srednja)
                    casi_rut.append(datetime_srednja)
                    break



    except:
        pass

# SORT DATA IN 2D PLAIN AND PLOT

os.chdir(mapa_s_skriptami)  # Go to scripts
slovar_kolicina_kolicina = pickle.load(open('slovar_kolicina_kolicina.p', 'rb'))    # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_enota = pickle.load(open('slovar_kolicina_enota.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_brez_enot = pickle.load(open('slovar_kolicina_brez_enot.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py

kolicine = [which_epsilon, 'Re_b']
oznaka_kolicin_v_podatkih = [slovar_kolicina_kolicina[kolicina] for kolicina in kolicine]
izpis_kolicin_v_grafih = [slovar_kolicina_enota[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]
izpis_kolicin_brez_enot = [slovar_kolicina_brez_enot[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]


# Find time of the very first and very last cast (all stations, all routes)
os.chdir(mapa_s_podatki)  # Go to data
datoteke = glob.glob('*.tob')   # all data files
os.chdir(mapa_s_skriptami)  # Go back to scripts

prva_meritev, zadnja_meritev = datoteke[0], datoteke[-1]    # first cast, last cast
for iznak in range(0, len(prva_meritev)):
    if prva_meritev[iznak] == '_':
        datetime_prva_od_vseh_meritev = datetime(year=int('20' + prva_meritev[iznak - 6:iznak - 4]), \
                                                month=int(prva_meritev[iznak - 4:iznak - 2]), \
                                                day=int(prva_meritev[iznak - 2:iznak]), \
                                                hour=int(prva_meritev[iznak + 1:iznak + 3]), \
                                                minute=int(prva_meritev[iznak + 3:iznak + 5]))
        break

for iznak in range(0, len(zadnja_meritev)):
    if zadnja_meritev[iznak] == '_':
        datetime_zadnja_od_vseh_meritev = datetime(year=int('20' + zadnja_meritev[iznak - 6:iznak - 4]), \
                                                    month=int(zadnja_meritev[iznak - 4:iznak - 2]), \
                                                    day=int(zadnja_meritev[iznak - 2:iznak]), \
                                                    hour=int(zadnja_meritev[iznak + 1:iznak + 3]), \
                                                    minute=int(zadnja_meritev[iznak + 3:iznak + 5]))
        break

# time between the firs and the last cast
zadnja_minus_prva = int((datetime_zadnja_od_vseh_meritev - datetime_prva_od_vseh_meritev).total_seconds())//60


vse_globine = []    # all depths
for iruta_zabelezena in range(len(casi_rut)):
    vse_globine += izracuni_epsilon_in_Re_b_za_to_postajo[iruta_zabelezena]['Press'].tolist()

globalna_min_globina, globalna_max_globina = min(vse_globine), max(vse_globine) # global minimum and maximu depth
korak_globine = vse_globine[1] - vse_globine[0]   # pressure step (in my cases always 1 dbar)

# all possible depths
vse_mozne_globine = np.linspace(globalna_min_globina, globalna_max_globina, \
                                int((globalna_max_globina - globalna_min_globina)/korak_globine + 1))



x_casi = [[] for i in oznaka_kolicin_v_podatkih]    # x-axis: time
y_globine = [[] for i in oznaka_kolicin_v_podatkih] # y-axis: depth
z_vrednosti = [[] for i in oznaka_kolicin_v_podatkih]   # z-axis: value

for iruta_zabelezena in range(len(casi_rut)):
    zacetni_cas_te_meritve = casi_rut[iruta_zabelezena]
    relativni_zacetni_cas_te_meritve = int((zacetni_cas_te_meritve - datetime_prva_od_vseh_meritev).total_seconds())//60
    for iglobina in range(len(vse_mozne_globine)):
        if vse_mozne_globine[iglobina] in izracuni_epsilon_in_Re_b_za_to_postajo[iruta_zabelezena]['Press']:
            for ioznaka_kolicin_v_podatkih in range(len(oznaka_kolicin_v_podatkih)):
                x_casi[ioznaka_kolicin_v_podatkih].append(relativni_zacetni_cas_te_meritve)
                y_globine[ioznaka_kolicin_v_podatkih].append(vse_mozne_globine[iglobina])
                z_vrednosti[ioznaka_kolicin_v_podatkih].append(\
                    izracuni_epsilon_in_Re_b_za_to_postajo[iruta_zabelezena][oznaka_kolicin_v_podatkih[ioznaka_kolicin_v_podatkih]][
                        izracuni_epsilon_in_Re_b_za_to_postajo[iruta_zabelezena]['Press'] \
                            .tolist().index(vse_mozne_globine[iglobina])]
                )


# maximum depth of each cast
max_globine = [y_globine[0][i] for i in range(len(y_globine[0]) - 1) if y_globine[0][i] > y_globine[0][i+1]] + [y_globine[0][-1]]

posamezni_x_casi = []   # time of each route
for posamezen_cas in x_casi[0]:
    if posamezen_cas not in posamezni_x_casi:
        posamezni_x_casi.append(posamezen_cas)

# Initiate grid
st_izmerkov_x = moja_max_ruta - moja_min_ruta + 1
st_izmerkov_y = len(vse_mozne_globine)
multiplikator_za_mrezo_x, multiplikator_za_mrezo_y = 2.5, 2
ngridx = int(st_izmerkov_x * multiplikator_za_mrezo_x)
ngridy = int(st_izmerkov_y * multiplikator_za_mrezo_y)
konture = True

# Plot (only works if len(kolicine) > 1)
fig, axs = plt.subplots(nrows=len(kolicine), ncols=1, figsize=(6, 4 * len(kolicine)))
st_izmerkov_x = moja_max_ruta - moja_min_ruta + 1
st_izmerkov_y = len(vse_mozne_globine)

for ikolicine in range(len(kolicine)):
    xi = np.linspace(min(x_casi[0]), max(x_casi[0]), ngridx)
    yi = np.linspace(min(y_globine[0]), max(y_globine[0]), ngridy)

    triang = tri.Triangulation(x_casi[ikolicine], y_globine[ikolicine])
    if kolicine[ikolicine] == 'Re_b':
        interpolator = tri.LinearTriInterpolator(triang, np.log10(np.array(z_vrednosti[ikolicine])))
        Xi, Yi = np.meshgrid(xi, yi)
        zi = 10 ** interpolator(Xi, Yi)
    else:
        interpolator = tri.LinearTriInterpolator(triang, z_vrednosti[ikolicine])
        Xi, Yi = np.meshgrid(xi, yi)
        zi = interpolator(Xi, Yi)

    if kolicine[ikolicine] == 'Re_b':
        ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='brg', extend='both', levels=[1, 10, 100, 200],\
                                      norm=colors.LogNorm(vmin=1, vmax=200))
        axs[ikolicine].contour(xi, yi, zi, colors='k', levels=[1, 10, 100, 200],\
                                      norm=colors.LogNorm(vmin=1, vmax=200))
        fig.colorbar(ctf, label=izpis_kolicin_v_grafih[ikolicine], ax=axs[ikolicine], format='%d')
    else:   # quantity is epsilon
        ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='jet', extend='both', levels=np.arange(start=-9.0, stop=-5.01, step=0.4))
        axs[ikolicine].contour(xi, yi, zi, colors='k', levels=np.arange(start=-9.0, stop=-5.01, step=0.4))
        fig.colorbar(ctf, label=izpis_kolicin_v_grafih[ikolicine], ax=axs[ikolicine])


    axs[ikolicine].invert_yaxis()
    axs[ikolicine].set_ylabel(slovar_kolicina_enota['Press'])
    # Draw vertical lines for times of routes
    for casek in x_casi[ikolicine]:
        axs[ikolicine].axvline(casek, linestyle='--', color='grey', linewidth=0.7)

    # Fill the area without measurements with black color (interpolation interpolates gere as well, we don't want that)
    axs[ikolicine].fill_between(posamezni_x_casi, max_globine, [max(max_globine) for i in max_globine], color='k')

axs[-1].set_xlabel('Ura meritve')
axs[0].set_title('{}, {}.-{}.{}.20{}'.format(postaja, dan, str(int(dan)+1), mesec, leto))

# all possible times
vsi_mozni_casi = [str(datetime_prva_od_vseh_meritev + timedelta(minutes=1 * i))[-8:-3] \
                  for i in range(zadnja_minus_prva)]
vsi_mozni_casi_brez_str = [datetime_prva_od_vseh_meritev + timedelta(minutes=1 * i) \
                           for i in range(zadnja_minus_prva)]
vsi_mozni_casi_sekunde = [i for i in range(zadnja_minus_prva)]
oznake_x_str, oznake_x_int = [], []
for j in range(len(vsi_mozni_casi_sekunde)):
    if vsi_mozni_casi[j][1:] in ['{}:00'.format(i) for i in range(10) if i % 2 == 0]:
        if casi_rut[0] <= vsi_mozni_casi_brez_str[j] and casi_rut[-1] >= vsi_mozni_casi_brez_str[j]:
            oznake_x_int.append(j)
            oznake_x_str.append(vsi_mozni_casi[j][:2])
plt.setp(axs, xticks=oznake_x_int, xticklabels=oznake_x_str)

if save_figures:
    if which_epsilon == 'log_epsilon_MLE_mean':
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_single_station_{datumska_mapa}_{idejna_ruta}_si.jpg', dpi=200)
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_single_station_{datumska_mapa}_{idejna_ruta}_{which_epsilon}_si.jpg', dpi=200)
plt.show()