'''This file plots oceanographic measurements from TOB files in "epsilon" directory on one station throughout the entire
 25-hour field campaigns. WARNING: The interpolation is done very badly, because the temporal distance between two
 routes can be from 10 minutes to 3 hours! You can only rely on the areas, which show the true measurements (i.e. where
 the vertical dashed line is). The plots from Risanje_celega_kanala_na_isti_ruti.py are much more reliable and also
 much more sophisticated when it comes to contour-spacing, colormaps, etc.  To set different (perhaps English) displayed
names of quantities, go to Izdelava_pickle_slovarjev_za_risanje.py.'''

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

# Select the quantities to plot (see Izdelava_pickle_slovarjev_za_risanje.py for the names of quantities)
# Some examples: 'temperatura' ... temperature, 'slanost' ... salinity, 'gostota' ... density anomaly,
#                   'N2' ... squared buoyancy frequency, ...
kolicine = ['temperatura', 'slanost', 'gostota', 'N2']



# THE BASIC MANUAL SETTINGS END HERE


leto = datumska_mapa[6:8]   # year
mesec = datumska_mapa[9:11] # month
dan = datumska_mapa[12:]    # day
tipska_mapa = 'epsilon'  # directory with the TOB files that contain our data


# Load data

bs = '\ '[0]
mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa   # directory with data
mapa_s_skriptami = os.getcwd()  # directory with scripts

min_ruta = 1 # Min. route index in data
max_rute = {'all_2008_08_21':14, 'all_2008_11_18':16, 'all_2009_04_20':19}  # Max. route index in data
max_ruta = max_rute[datumska_mapa]    # Max. route index for this date
moja_min_ruta = 1    # Min. station index to be plotted (WARNING - not tested for any other number than 1)
moja_max_ruta = max_ruta    # Min. station index to be plotted (WARNING - not tested for any other number than max_ruta)



vse_mean = []   # all means
vse_std = []    # all standard deviations
vse_median = [] # all medians
casi_rut = []   # times of measurements (i.e. times of routes at the chosen station)

for iruta in range(moja_min_ruta, moja_max_ruta + 1): # Loop through all routes
    ruta = 'RUTA_{:02d}'.format(iruta)  # route as written in the name of a TOB file
    try:
        tokratne_mean, tokratne_std, tokratne_median = Branje_datotek(izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, \
                                                                        tipska_mapa=tipska_mapa, postaja=postaja, ruta=ruta)
        if len(tokratne_mean) > 0:  # Not all stations were measured in all routes!
            os.chdir(mapa_s_podatki)    # Go to data
            datoteke_ruta = glob.glob('*_' + postaja + '_' + ruta + '.tob') # All TOB files with this combination
                                                                            # station-route
            os.chdir(mapa_s_skriptami)  # Go back to scripts
            srednja_meritev_ruta = datoteke_ruta[len(datoteke_ruta) // 2]   # Time of the medial measurement

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

            vse_mean.append(tokratne_mean)
            vse_std.append(tokratne_std)
            vse_median.append(tokratne_median)


    except:
        pass

os.chdir(mapa_s_skriptami)  # Go to scripts
slovar_kolicina_kolicina = pickle.load(open('slovar_kolicina_kolicina.p', 'rb'))    # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_enota = pickle.load(open('slovar_kolicina_enota.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_brez_enot = pickle.load(open('slovar_kolicina_brez_enot.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py


oznaka_kolicin_v_podatkih = [slovar_kolicina_kolicina[kolicina] for kolicina in kolicine]
izpis_kolicin_v_grafih = [slovar_kolicina_enota[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]
izpis_kolicin_brez_enot = [slovar_kolicina_brez_enot[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]




# Find time of the very first and very last cast (all stations, all routes)
os.chdir(mapa_s_podatki)  # Go to data
datoteke = glob.glob('*.tob')   # all TOB files
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



# all depths
vse_globine = []
for iruta_zabelezena in range(len(casi_rut)):
    vse_globine += vse_mean[iruta_zabelezena]['Press'].tolist()

globalna_min_globina, globalna_max_globina = min(vse_globine), max(vse_globine) # global minimum and maximu depth
korak_globine = vse_globine[1] - vse_globine[0]   # pressure step (in my cases always 1 dbar)

# all possible depths
vse_mozne_globine = np.linspace(globalna_min_globina, globalna_max_globina, \
                                int((globalna_max_globina - globalna_min_globina)/korak_globine + 1))

x_casi = [[] for i in oznaka_kolicin_v_podatkih]    # x-axis: time
y_globine = [[] for i in oznaka_kolicin_v_podatkih] # y-axis: depth
z_vrednosti = [[] for i in oznaka_kolicin_v_podatkih]   # z-axis: value

for iruta_zabelezena in range(len(vse_mean)):
    zacetni_cas_te_meritve = casi_rut[iruta_zabelezena]
    relativni_zacetni_cas_te_meritve = int((zacetni_cas_te_meritve - datetime_prva_od_vseh_meritev).total_seconds())//60
    for iglobina in range(len(vse_mozne_globine)):
        if vse_mozne_globine[iglobina] in vse_mean[iruta_zabelezena]['Press']:
            for ioznaka_kolicin_v_podatkih in range(len(oznaka_kolicin_v_podatkih)):
                x_casi[ioznaka_kolicin_v_podatkih].append(relativni_zacetni_cas_te_meritve)
                y_globine[ioznaka_kolicin_v_podatkih].append(vse_mozne_globine[iglobina])
                z_vrednosti[ioznaka_kolicin_v_podatkih].append(\
                    vse_mean[iruta_zabelezena][oznaka_kolicin_v_podatkih[ioznaka_kolicin_v_podatkih]][
                        vse_mean[iruta_zabelezena]['Press'] \
                            .tolist().index(vse_mozne_globine[iglobina])]
                )

# maximum depth of each cast
max_globine = [y_globine[0][i] for i in range(len(y_globine[0]) - 1) if y_globine[0][i] > y_globine[0][i+1]] + [y_globine[0][-1]]
posamezni_x_casi = []
for posamezen_cas in x_casi[0]:
    if posamezen_cas not in posamezni_x_casi:
        posamezni_x_casi.append(posamezen_cas)

# Initiate grid
st_izmerkov_x = moja_max_ruta - moja_min_ruta + 1
st_izmerkov_y = len(vse_mozne_globine)
multiplikator_za_mrezo_x, multiplikator_za_mrezo_y = 3, 1
ngridx = int(st_izmerkov_x * multiplikator_za_mrezo_x)
ngridy = int(st_izmerkov_y * multiplikator_za_mrezo_y)
konture = True

if len(kolicine) > 1:

    fig, axs = plt.subplots(nrows=len(kolicine), ncols=1, figsize=(6, 4 * len(kolicine)))
    st_izmerkov_x = moja_max_ruta - moja_min_ruta + 1
    st_izmerkov_y = len(vse_mozne_globine)

    for ikolicine in range(len(kolicine)):
        xi = np.linspace(min(x_casi[0]), max(x_casi[0]), ngridx)
        yi = np.linspace(min(y_globine[0]), max(y_globine[0]), ngridy)
        triang = tri.Triangulation(x_casi[ikolicine], y_globine[ikolicine])

        if kolicine[ikolicine] not in ('epsilon', 'epsilon1', 'epsilon2', 'psevdo epsilon', 'N2'):
            interpolator = tri.LinearTriInterpolator(triang, z_vrednosti[ikolicine])
            Xi, Yi = np.meshgrid(xi, yi)
            zi = interpolator(Xi, Yi)
            ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='jet', extend='both')
            if konture:
                axs[ikolicine].contour(xi, yi, zi, colors='k')
        else:
            interpolator = tri.LinearTriInterpolator(triang, np.log10(z_vrednosti[ikolicine]))
            Xi, Yi = np.meshgrid(xi, yi)
            zi = 10**interpolator(Xi, Yi)
            if kolicine[ikolicine] != 'N2':
                ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='jet', extend='both',\
                                              norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi)))
                if konture:
                    axs[ikolicine].contour(xi, yi, zi, colors='k',\
                                              norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi)))
            else:
                ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='jet', extend='both',\
                                              norm=colors.LogNorm(vmin=10**-7, vmax=10**-1), levels=[10**-7, 10**-6, 10**-5, 10**-4, 10**-3, 10**-2, 10**-1])
                if konture:
                    axs[ikolicine].contour(xi, yi, zi, colors='k',\
                                              norm=colors.LogNorm(vmin=10**-7, vmax=10**-1), levels=[10**-7, 10**-6, 10**-5, 10**-4, 10**-3, 10**-2, 10**-1])

        fig.colorbar(ctf, label=izpis_kolicin_v_grafih[ikolicine], ax=axs[ikolicine])
        axs[ikolicine].invert_yaxis()
        axs[ikolicine].set_ylabel(slovar_kolicina_enota['Press'])
        axs[ikolicine].set_yticks(ticks=[i for i in range(100) if i % 2 == 0])
        axs[ikolicine].set_yticklabels([i for i in range(100) if i % 2 == 0])
        # Draw vertical lines for times of routes
        for casek in x_casi[ikolicine]:
            axs[ikolicine].axvline(casek, linestyle='--', color='grey', linewidth=0.7)

        # Cover the "sea floor"
        axs[ikolicine].fill_between(posamezni_x_casi, max_globine, [max(max_globine) for i in max_globine], color='k')

    axs[-1].set_xlabel('Ura meritve')
    axs[0].set_title('{}, {}.{}.20{}'.format(postaja, dan, mesec, leto))


    vsi_mozni_casi = [str(datetime_prva_od_vseh_meritev + timedelta(minutes=1 * i))[-8:-3] \
                      for i in range(zadnja_minus_prva)]
    vsi_mozni_casi_brez_str = [datetime_prva_od_vseh_meritev + timedelta(minutes=1 * i) \
                               for i in range(zadnja_minus_prva)]
    vsi_mozni_casi_sekunde = [i for i in range(zadnja_minus_prva)]
    oznake_x_str, oznake_x_int = [], []
    for j in range(len(vsi_mozni_casi_sekunde)):
        if vsi_mozni_casi[j][1:] in ['{}:00'.format(i) for i in range(10) if i % 2 == 0]:
            if casi_rut[0] <= vsi_mozni_casi_brez_str[j] and casi_rut[-1] >= vsi_mozni_casi_brez_str[j]:
                # zato, da nimam praznega belega polja pri straneh po nepotrebnem
                oznake_x_int.append(j)
                oznake_x_str.append(vsi_mozni_casi[j][:2])
    plt.setp(axs, xticks=oznake_x_int, xticklabels=oznake_x_str)

    if save_figures:
        zapis_kolicin_za_savefig = ''
        for kolicina in kolicine:
            zapis_kolicin_za_savefig += kolicina + '_'
        plt.savefig(mapa_za_shranjevanje_grafov + f'{zapis_kolicin_za_savefig}_{leto}_{mesec}_{dan}_{postaja}.jpg')
    plt.show()

else:
    xi = np.linspace(min(x_casi[0]), max(x_casi[0]), int(ngridx))
    yi = np.linspace(min(y_globine[0]), max(y_globine[0]), int(ngridy))

    triang = tri.Triangulation(x_casi[0], y_globine[0])
    interpolator = tri.LinearTriInterpolator(triang, z_vrednosti[0])
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    if kolicine[0] not in ('epsilon', 'epsilon1', 'epsilon2', 'psevdo epsilon'):
        ctf = plt.contourf(xi, yi, zi, cmap='jet', extend='both')
        if konture:
            plt.contour(xi, yi, zi, colors='k')
    else:
        ctf = plt.contourf(xi, yi, zi, cmap='jet', extend='both', \
                                      norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi)))
        if konture:
            plt.contour(xi, yi, zi, colors='k', \
                                      norm=colors.LogNorm(vmin=np.nanmin(zi), vmax=np.nanmax(zi)))

    plt.colorbar(ctf, label=izpis_kolicin_v_grafih[0])
    plt.gca().invert_yaxis()
    plt.ylabel(slovar_kolicina_enota['Press'])
    # Draw vertical lines for times of routes
    for casek in x_casi[0]:
        plt.axvline(casek, linestyle='--', color='grey', linewidth=0.7)

    # Cover the "sea floor"
    plt.fill_between(posamezni_x_casi, max_globine, [max(max_globine) for i in max_globine], color='k')
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
    plt.xticks(oznake_x_int, labels=oznake_x_str)
    plt.xlabel('Ura meritve')
    plt.title('{}, {}, {}.{}.20{}'.format(izpis_kolicin_brez_enot[0], postaja, dan, mesec, leto))

    if save_figures:
    zapis_kolicin_za_savefig = ''
    for kolicina in kolicine:
        zapis_kolicin_za_savefig += kolicina + '_'
    plt.savefig(mapa_za_shranjevanje_grafov + f'{zapis_kolicin_za_savefig}_{leto}_{mesec}_{dan}_{postaja}.jpg')

    plt.show()
