'''This file plots the epsilon profiles for the entire basin in a specific route, accompanied by the buoyancy Reynolds
number. If LK01 or LK07 wasn't measured in the specific route, it takes its values from the previous route (which was
just minutes earlier). If any other station is missing, an error occurs.'''
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

datumska_mapa = 'all_2008_08_21'    # date of measurements ('all_2008_08_21', 'all_2008_11_18', 'all_2009_04_20)
ruta = 'RUTA_01'    # select route (01-14 for 21.8.2008, 01-16 for 18.11.2008, 01-19 for 20.4.2009)

# Which epsilon should we plot? Options: 'log_epsilon_MLE_mean', 'log_epsilon_MLE_1', 'log_epsilon_MLE_2',
# 'log_epsilon_IM_mean', 'log_epsilon_IM_1', 'log_epsilon_IM_2'
which_epsilon = 'log_epsilon_MLE_mean'





# THE BASIC MANUAL SETTINGS END HERE

leto = datumska_mapa[6:8]   # year
mesec = datumska_mapa[9:11] # month
dan = datumska_mapa[12:]    # day
idejna_ruta = ruta  # the route we want
ruta_minus_1 = 'RUTA_{:02d}'.format(int(idejna_ruta[-2:]) - 1)  # the route prior to our route
tipska_mapa = 'epsilon' # directory, where data from each measurement is stored (only used to determine the time of each
                        # measurement, which isn't stored in f'slovar_epsilon_za_risat_NA_1m_{datumska_mapa}.p'
# Load data
izracuni_epsilon_in_Re_b = pickle.load(open(f'slovar_epsilon_za_risat_NA_1m_{datumska_mapa}.p', 'rb'))



bs = '\ '[0]
mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa   # directory with data (only used to determine
                                                                        # time of each measurement)
mapa_s_skriptami = os.getcwd()  # directory with scripts

min_postaja = 1         # Min. station index in data
max_postaja = 7         # Max. station index in data
moja_min_postaja = 1    # Min. station index to be plotted (WARNING - not tested for any other number than 1)
moja_max_postaja = 7    # Max. station index to be plotted (WARNING - not tested for any other number than 7)

izracuni_epsilon_in_Re_b_za_to_ruto = []    # epsilon and Re_b for this route
postaje_rute = []   # used combinations of station-route
for postaja in ['LK{:02d}'.format(ip) for ip in range(moja_min_postaja, moja_max_postaja+1)]:   # Loop through all stations
    postaja_ruta = postaja + '_' + idejna_ruta  # station-route
    postaja_ruta_minus_1 = postaja + '_' + ruta_minus_1 # station-previous route
    print(postaja_ruta, postaja_ruta_minus_1)
    if postaja_ruta in izracuni_epsilon_in_Re_b.keys():
        print('pr')
        izracuni_epsilon_in_Re_b_za_to_ruto.append(izracuni_epsilon_in_Re_b[postaja_ruta])
        postaje_rute.append(postaja_ruta)
    elif postaja_ruta_minus_1 in izracuni_epsilon_in_Re_b.keys() and postaja in ['LK01', 'LK07']:
        print('pr-1')
        izracuni_epsilon_in_Re_b_za_to_ruto.append(izracuni_epsilon_in_Re_b[postaja_ruta_minus_1])
        postaje_rute.append(postaja_ruta_minus_1)
    else:
        AssertionError # Combination station-route not found neither for the route we want nor for the previous route


# FIND TIME OF EACH STATION'S MEASUREMENTS

casi_rut = []   # Times of measurements

for pr in postaje_rute: # Loop through all used combinations of station-route
    postaja = pr[:4]    # station
    ruta = pr[5:]       # route
    print(postaja, ruta)

    os.chdir(mapa_s_podatki)    # Go to the directory with data to find all casts from this combination station-route
    datoteke_ruta = glob.glob('*_' + postaja + '_' + ruta + '.tob') # Select all files with this combination
                                                                    # station-route
    os.chdir(mapa_s_skriptami)  # Go back to scripts
    srednja_meritev_postajaruta = datoteke_ruta[len(datoteke_ruta) // 2]   # Time of the medial measurement

    # Figure out time of the srednja_meritev_postaja_ruta from its name
    for iznak in range(0, len(srednja_meritev_postajaruta)):
        if srednja_meritev_postajaruta[iznak] == '_':
            datetime_srednja = datetime(year=int('20' + srednja_meritev_postajaruta[iznak - 6:iznak - 4]), \
                                        month=int(srednja_meritev_postajaruta[iznak - 4:iznak - 2]), \
                                        day=int(srednja_meritev_postajaruta[iznak - 2:iznak]), \
                                        hour=int(srednja_meritev_postajaruta[iznak + 1:iznak + 3]), \
                                        minute=int(srednja_meritev_postajaruta[iznak + 3:iznak + 5]))
            print(datetime_srednja)
            casi_rut.append(datetime_srednja)
            break




# SORT DATA IN 2D PLAIN AND PLOT

os.chdir(mapa_s_skriptami)  # Go to scripts

slovar_kolicina_kolicina = pickle.load(open('slovar_kolicina_kolicina.p', 'rb'))    # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_enota = pickle.load(open('slovar_kolicina_enota.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_brez_enot = pickle.load(open('slovar_kolicina_brez_enot.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py


kolicine = [which_epsilon, 'Re_b']
oznaka_kolicin_v_podatkih = [slovar_kolicina_kolicina[kolicina] for kolicina in kolicine]
izpis_kolicin_v_grafih = [slovar_kolicina_enota[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]
izpis_kolicin_brez_enot = [slovar_kolicina_brez_enot[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]




vse_globine = []    # all depths with valid data
for ipostaja in range(moja_min_postaja, moja_max_postaja + 1):
    vse_globine += izracuni_epsilon_in_Re_b_za_to_ruto[ipostaja - moja_min_postaja]['Press'].tolist()

globalna_min_globina, globalna_max_globina = min(vse_globine), max(vse_globine) # global minimum and maximum of depths
korak_globine = vse_globine[1] - vse_globine[0]   # pressure step (in our case always 1 dbar)

# all possible depths
vse_mozne_globine = np.linspace(globalna_min_globina, globalna_max_globina, \
                                int((globalna_max_globina - globalna_min_globina)/korak_globine + 1))


x_postaje = [[] for i in oznaka_kolicin_v_podatkih] # x-axis: stations
y_globine = [[] for i in oznaka_kolicin_v_podatkih] # y-axis: depths
z_vrednosti = [[] for i in oznaka_kolicin_v_podatkih]   # z-axis: values

for ipostaja in range(moja_min_postaja, moja_max_postaja + 1):
    for iglobina in range(len(vse_mozne_globine)):
        if vse_mozne_globine[iglobina] in izracuni_epsilon_in_Re_b_za_to_ruto[ipostaja - moja_min_postaja]['Press']:
            for ioznaka_kolicin_v_podatkih in range(len(oznaka_kolicin_v_podatkih)):
                x_postaje[ioznaka_kolicin_v_podatkih].append(ipostaja)
                y_globine[ioznaka_kolicin_v_podatkih].append(vse_mozne_globine[iglobina])
                z_vrednosti[ioznaka_kolicin_v_podatkih].append(\
                    izracuni_epsilon_in_Re_b_za_to_ruto[ipostaja - moja_min_postaja][oznaka_kolicin_v_podatkih[ioznaka_kolicin_v_podatkih]]\
                        [izracuni_epsilon_in_Re_b_za_to_ruto[ipostaja - moja_min_postaja]['Press'].tolist().index(vse_mozne_globine[iglobina])])


# maximum depth of each cast
max_globine = [y_globine[0][i] for i in range(len(y_globine[0]) - 1) if y_globine[0][i] > y_globine[0][i+1]] + [y_globine[0][-1]]

# Initiate grid
ngridx = (moja_max_postaja - moja_min_postaja) * 3 + 1
ngridy = int((max(max_globine)*1) * 2) - 1
konture = True

# Plot
fig, axs = plt.subplots(nrows=len(kolicine), ncols=1, figsize=(6, 4 * len(kolicine)))
st_izmerkov_x = moja_max_postaja - moja_min_postaja + 1
st_izmerkov_y = len(vse_mozne_globine)

for ikolicine in range(len(kolicine)):  # len(kolicine) has to be al least 2!
    xi = np.linspace(moja_min_postaja, moja_max_postaja, ngridx)    # x-coords
    yi = np.linspace(min(y_globine[0]), max(y_globine[0]), ngridy)  # y-coords

    triang = tri.Triangulation(x_postaje[ikolicine], y_globine[ikolicine])  # triangulation

    if kolicine[ikolicine] == 'Re_b':   # The quantity is buoyancy Reynolds number
        interpolator = tri.LinearTriInterpolator(triang, np.log10(z_vrednosti[ikolicine]))
        Xi, Yi = np.meshgrid(xi, yi)
        zi = 10**interpolator(Xi, Yi)
        ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='brg', extend='both', levels=[1, 10, 100, 200], \
                                      norm=colors.LogNorm(vmin=1, vmax=200))
        ct = axs[ikolicine].contour(xi, yi, zi, colors='k', levels=[1, 10, 100, 200], \
                               norm=colors.LogNorm(vmin=1, vmax=200))
        cb = fig.colorbar(ctf, ax=axs[ikolicine], format='%d')
        cb.ax.tick_params(labelsize=12)
        cb.set_label(izpis_kolicin_v_grafih[ikolicine], fontsize=13)
        cb.add_lines(ct)
    else:   # The quantity is epsilon
        interpolator = tri.LinearTriInterpolator(triang, z_vrednosti[ikolicine])
        Xi, Yi = np.meshgrid(xi, yi)
        zi = interpolator(Xi, Yi)
        ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='jet', extend='both',
                                      levels=np.arange(start=-8.6, stop=-4.99, step=0.4))
        ct = axs[ikolicine].contour(xi, yi, zi, colors='k', levels=np.arange(start=-8.6, stop=-4.99, step=0.4))
        cb = fig.colorbar(ctf, ax=axs[ikolicine])
        cb.ax.tick_params(labelsize=12)
        cb.set_label(izpis_kolicin_v_grafih[ikolicine], fontsize=13)
        cb.add_lines(ct)

    axs[ikolicine].invert_yaxis()   # Invert y-axis
    axs[ikolicine].set_ylabel(slovar_kolicina_enota['Press'], fontsize=13)  # Set y-label
    for ip in range(moja_min_postaja, moja_max_postaja+1):
        axs[ikolicine].axvline(ip, linestyle='--', color='grey', linewidth=0.7)

    # Fill the area without measurements with black color (interpolation interpolates gere as well, we don't want that)
    axs[ikolicine].fill_between([i for i in range(moja_min_postaja, moja_max_postaja+1)], max_globine, [max(max_globine) for i in max_globine], color='k')
    # Set location of xticks
    xticks = [i for i in range(moja_min_postaja, moja_max_postaja+1)]
    axs[ikolicine].set_xticks(xticks)
    # Set labels of xticks (station name + time of the medial measurement)
    axs[ikolicine].set_xticklabels(['LK{:02d}\n{}'.format(xticks[i], str(casi_rut[i])[-8:-3]) for i in range(len(xticks))], fontsize=12)
    # Set location and labels of yticks
    yticks = [i for i in range(2, 16+1) if i%2 == 0]
    axs[ikolicine].set_yticks(yticks)
    axs[ikolicine].set_yticklabels([yt for yt in yticks], fontsize=12)

axs[-1].set_xlabel('Postaja (ura)', fontsize=13)
print(min(casi_rut))
print(str(min(casi_rut)))

# The title also tells us the direction of measurements
if casi_rut[0] < casi_rut[-1]:
    axs[0].set_title(r'Ruta {}, {}.{}.20{} {}-{}, LK{:02d}$\rightarrow$LK{:02d}'.format(ruta[-2:], str(min(casi_rut))[8:10], mesec, leto, str(min(casi_rut))[-8:-3], str(max(casi_rut))[-8:-3], moja_min_postaja, moja_max_postaja), fontsize=14.5)
else:
    axs[0].set_title(r'Ruta {}, {}.{}.20{} {}-{}, LK{:02d}$\rightarrow$LK{:02d}'.format(ruta[-2:], str(min(casi_rut))[8:10],
                                                                           mesec, leto, str(min(casi_rut))[-8:-3],
                                                                           str(max(casi_rut))[-8:-3],
                                                                           moja_max_postaja, moja_min_postaja), fontsize=14.5)


if save_figures:
    if which_epsilon == 'log_epsilon_MLE_mean':
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_crossesction_{datumska_mapa}_{idejna_ruta}_si.jpg', dpi=200)
    else:
        plt.savefig(mapa_za_shranjevanje_grafov + f'epsilon_crossesction_{datumska_mapa}_{idejna_ruta}_{which_epsilon}_si.jpg', dpi=200)
plt.show()