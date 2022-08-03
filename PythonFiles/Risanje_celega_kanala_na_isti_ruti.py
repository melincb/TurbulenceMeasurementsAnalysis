'''Plots the values of the desired quantities from the TOB files in the "epsilon" directory for each route on all
stations. If LK01 or LK07 wasn't measured in the specific route, it takes its values from the previous route (which was
just minutes earlier). If any other station is missing, an error occurs. To set different (perhaps English) displayed
names of quantities, go to Izdelava_pickle_slovarjev_za_risanje.py.'''

from Branje_datotek import Branje_datotek
import numpy as np
import glob, os
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.tri as tri
import pickle

izvirna_mapa = '..' # directory with directories of measurements
mapa_za_shranjevanje_grafov =  'Saved_figures\ '[:-1]    # Where it saves the figures
save_figures = False    # To save the plot, set this to True, else False

datumska_mapa = 'all_2008_08_21'    # date of measurements ('all_2008_08_21', 'all_2008_11_18', 'all_2009_04_20)
ruta = 'RUTA_01'    # select route (01-14 for 21.8.2008, 01-16 for 18.11.2008, 01-19 for 20.4.2009)


# Select the quantities to plot (see Izdelava_pickle_slovarjev_za_risanje.py for the names of quantities)
# Some examples: 'temperatura' ... temperature, 'slanost' ... salinity, 'gostota' ... density anomaly,
#                   'N2' ... squared buoyancy frequency, ...
kolicine = ['temperatura', 'slanost', 'gostota', 'N2']


# COMPUTE THE VALUES OF CONTOURS
dT = 0.35   # temperature step between two contours
nlevelsT = int(3.8/dT)  # number of levels of temperature
dS = 0.15   # salinity step between two contours
nlevelsS = int(1.7/dS)  # number of levels of salinity
drho = 0.2  # density step between two contours
nlevelsrho = int(2.5/drho)  # number of levels of density anomaly

# Here you can set the values of contours of your own quantities and data (the values for 'N2' are always [10**-6,
# 10**-5, 10**-4, 10**-3, 10**-2, 10**-1], you can change them further on in the code)
slovar_mej_kolicin = {'all_2008_08_21':{'temperatura':np.sort([25.1-i*dT for i in range(1, nlevelsT//2+1)]).tolist()+[25.1 + i*dT for i in range(nlevelsT//2+1)],\
                                        'slanost':np.sort([36.85-i*dS for i in range(1, nlevelsS//2+1)]).tolist()+[36.85 + i*dS for i in range(nlevelsS//2+1)],\
                                        'gostota':np.sort([24.7-i*drho for i in range(1, nlevelsrho//2+1)]).tolist()+[24.7 + i*drho for i in range(nlevelsrho//2+1)]},\
                      'all_2008_11_18':{'temperatura':np.sort([16.0-i*dT for i in range(1, nlevelsT//2+1)]).tolist()+[16.0 + i*dT for i in range(nlevelsT//2+1)],\
                                        'slanost':np.sort([36.35-i*dS for i in range(1, nlevelsS//2+1)]).tolist()+[36.35 + i*dS for i in range(nlevelsS//2+1)],\
                                        'gostota':np.sort([26.8-i*drho for i in range(1, nlevelsrho//2+1)]).tolist()+[26.8 + i*drho for i in range(nlevelsrho//2+1)]},\
                      'all_2009_04_20':{'temperatura':np.sort([13.6-i*dT for i in range(1, nlevelsT//2+1)]).tolist()+[13.6 + i*dT for i in range(nlevelsT//2+1)],\
                                        'slanost':np.sort([37-i*dS for i in range(1, nlevelsS//2+1)]).tolist()+[37 + i*dS for i in range(nlevelsS//2+1)],\
                                        'gostota':np.sort([27.8-i*drho for i in range(1, nlevelsrho//2+1)]).tolist()+[27.8 + i*drho for i in range(nlevelsrho//2+1)]}}



# THE BASIC MANUAL SETTINGS END HERE

leto = datumska_mapa[6:8]   # year
mesec = datumska_mapa[9:11] # month
dan = datumska_mapa[12:]    # day
idejna_ruta = ruta  # the route we want
ruta_minus_1 = 'RUTA_{:02d}'.format(int(idejna_ruta[-2:]) - 1)  # the route prior to our route
tipska_mapa = 'epsilon' # directory with the TOB files that contain our data


bs = '\ '[0]
mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa   # directory with data (only used to determine
                                                                        # time of each measurement)
mapa_s_skriptami = os.getcwd()  # directory with scripts

min_postaja = 1         # Min. station index in data
max_postaja = 7         # Max. station index in data
moja_min_postaja = 1    # Min. station index to be plotted (WARNING - not tested for any other number than 1)
moja_max_postaja = 7    # Max. station index to be plotted (WARNING - not tested for any other number than 7)


vse_mean = []   # all means
vse_std = []    # all standard deviations
vse_median = [] # all medians
casi_rut = []   # all times of routes

# GET DATA AND TIME OF EACH MEASUREMENT
os.chdir(mapa_s_podatki)
for ipostaja in range(moja_min_postaja, moja_max_postaja + 1):
    postaja = 'LK{:02d}'.format(ipostaja)
    try:
        tokratne_mean, tokratne_std, tokratne_median = Branje_datotek(izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, tipska_mapa=tipska_mapa, postaja=postaja, ruta=ruta)

        datoteke_ruta = glob.glob('*_' + postaja + '_' + ruta + '.tob')
        srednja_meritev_ruta = datoteke_ruta[len(datoteke_ruta) // 2]  # time of the medial route
        for iznak in range(0, len(srednja_meritev_ruta)):
            if srednja_meritev_ruta[iznak] == '_':
                datetime_srednja = datetime(year=int('20' + srednja_meritev_ruta[iznak - 6:iznak - 4]), \
                                            month=int(srednja_meritev_ruta[iznak - 4:iznak - 2]), \
                                            day=int(srednja_meritev_ruta[iznak - 2:iznak]), \
                                            hour=int(srednja_meritev_ruta[iznak + 1:iznak + 3]), \
                                            minute=int(srednja_meritev_ruta[iznak + 3:iznak + 5]))
                break
        casi_rut.append(datetime_srednja)
    except: # If 'LK01' or 'LK07' is missing, we take the data from the prior route, which ended with the missing station
        ruta_minus_1 = 'RUTA_{:02d}'.format(int(ruta[-2:]) - 1)
        tokratne_mean, tokratne_std, tokratne_median = Branje_datotek(izvirna_mapa=izvirna_mapa, datumska_mapa=datumska_mapa, tipska_mapa=tipska_mapa, postaja=postaja, ruta=ruta_minus_1)

        datoteke_ruta = glob.glob('*_' + postaja + '_' + ruta_minus_1 + '.tob')
        srednja_meritev_ruta = datoteke_ruta[len(datoteke_ruta) // 2]  # time of the medial route
        for iznak in range(0, len(srednja_meritev_ruta)):
            if srednja_meritev_ruta[iznak] == '_':
                datetime_srednja = datetime(year=int('20' + srednja_meritev_ruta[iznak - 6:iznak - 4]), \
                                            month=int(srednja_meritev_ruta[iznak - 4:iznak - 2]), \
                                            day=int(srednja_meritev_ruta[iznak - 2:iznak]), \
                                            hour=int(srednja_meritev_ruta[iznak + 1:iznak + 3]), \
                                            minute=int(srednja_meritev_ruta[iznak + 3:iznak + 5]))
                break
        casi_rut.append(datetime_srednja)


    vse_mean.append(tokratne_mean)
    vse_std.append(tokratne_std)
    vse_median.append(tokratne_median)

os.chdir(mapa_s_skriptami)

slovar_kolicina_kolicina = pickle.load(open('slovar_kolicina_kolicina.p', 'rb'))    # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_enota = pickle.load(open('slovar_kolicina_enota.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py
slovar_kolicina_brez_enot = pickle.load(open('slovar_kolicina_brez_enot.p', 'rb'))   # From Izdelava_pickle_slovarjev_za_risanje.py


oznaka_kolicin_v_podatkih = [slovar_kolicina_kolicina[kolicina] for kolicina in kolicine]
izpis_kolicin_v_grafih = [slovar_kolicina_enota[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]
izpis_kolicin_brez_enot = [slovar_kolicina_brez_enot[oznaka_kolicine_v_podatkih] for oznaka_kolicine_v_podatkih in oznaka_kolicin_v_podatkih]




# PLOT

def pripravi_lokacije_tock_in_trikotnike_za_tricontourf(vhodno_polje):
    '''This function sets the location of the grid points and the indexes of these points for each triangle.'''
    Z = np.transpose(vhodno_polje)

    lokacije_tock_na_grafu = []
    for i in range(np.shape(Z)[0]):
        for j in range(np.shape(Z)[1]):
            lokacije_tock_na_grafu.append([i,j])
            # indexes for triangulation are np.shape(Z)[1]*i + j

    def zaporedna_stevilka_lokacije_tock_na_grafu(vrstica, stolpec, Z):
        '''Computes index of a grid point in row vrstica and column stolpec.'''
        return np.shape(Z)[1] * vrstica + stolpec

    trikotniki = [] # triangles
    for j in range(np.shape(Z)[1] - 1): # -1, because there is j+1 later in the loop
        for i in range(np.shape(Z)[0] - 1): # -1, because there is i+1 later in the loop
            if Z[i+1,j] < np.infty and Z[i+1,j+1] < np.infty:
                trikotniki.append([zaporedna_stevilka_lokacije_tock_na_grafu(i, j, Z),\
                                   zaporedna_stevilka_lokacije_tock_na_grafu(i+1, j, Z),\
                                   zaporedna_stevilka_lokacije_tock_na_grafu(i, j+1, Z)])
                trikotniki.append([zaporedna_stevilka_lokacije_tock_na_grafu(i+1, j, Z),\
                                   zaporedna_stevilka_lokacije_tock_na_grafu(i+1, j+1, Z),\
                                   zaporedna_stevilka_lokacije_tock_na_grafu(i, j+1, Z)])
            elif Z[i+1, j+1] < np.infty:
                k = 0
                while Z[i+1+k, j+1] < np.infty:
                    trikotniki.append([zaporedna_stevilka_lokacije_tock_na_grafu(i, j, Z),\
                                       zaporedna_stevilka_lokacije_tock_na_grafu(i+1+k, j+1, Z),\
                                       zaporedna_stevilka_lokacije_tock_na_grafu(i+k, j+1, Z)])
                    k += 1
                    if i+1+k >= np.shape(Z)[0]:
                        break
                break
            elif Z[i+1, j] < np.infty:
                k = 0
                while Z[i+1+k, j] < np.infty:
                    trikotniki.append([zaporedna_stevilka_lokacije_tock_na_grafu(i+k, j, Z),\
                                       zaporedna_stevilka_lokacije_tock_na_grafu(i+1+k, j, Z),\
                                       zaporedna_stevilka_lokacije_tock_na_grafu(i, j+1, Z)])
                    k += 1
                    if i+1+k >= np.shape(Z)[0]:
                        break
                break
            else:
                break

    lokacije_tock_na_grafu = np.asarray(lokacije_tock_na_grafu)
    trikotniki = np.asarray(trikotniki)

    return lokacije_tock_na_grafu, trikotniki


vse_globine = []    # all depths
for ipostaja in range(moja_min_postaja, moja_max_postaja + 1):  # Loop through all stations
    vse_globine += vse_mean[ipostaja - moja_min_postaja]['Press'].tolist()

globalna_min_globina, globalna_max_globina = min(vse_globine), max(vse_globine) # global minimum and maximum depths
korak_globine = vse_globine[1] - vse_globine[0]   # pressure step (in my case always 0.25 dbar)

# all possible depths
vse_mozne_globine = np.linspace(globalna_min_globina, globalna_max_globina, \
                                int((globalna_max_globina - globalna_min_globina)/korak_globine + 1))
# list of quantities at all stations in all depths
seznam_kolicin_na_vseh_postajah_na_vseh_globinah = np.empty((len(kolicine), moja_max_postaja - moja_min_postaja + 1, \
                                                            int((globalna_max_globina - globalna_min_globina)/korak_globine + 1)))
seznam_kolicin_na_vseh_postajah_na_vseh_globinah[:, :] = np.nan # If there is no measurement, we want nan
for ipostaja in range(moja_min_postaja, moja_max_postaja + 1):  # Loop throught stations
    for iglobina in range(len(vse_mozne_globine)):  # Loop through depths
        if vse_mozne_globine[iglobina] in vse_mean[ipostaja - moja_min_postaja]['Press']:   # If this depth is in the
                                                                                            # data for this station
            for ioznaka_kolicin_v_podatkih in range(len(oznaka_kolicin_v_podatkih)):    # Append the data
                seznam_kolicin_na_vseh_postajah_na_vseh_globinah[ioznaka_kolicin_v_podatkih, ipostaja - moja_min_postaja, iglobina] =\
                    vse_mean[ipostaja - moja_min_postaja][oznaka_kolicin_v_podatkih[ioznaka_kolicin_v_podatkih]]\
                        [vse_mean[ipostaja - moja_min_postaja]['Press'].tolist().index(vse_mozne_globine[iglobina])]




x_postaje = [[] for i in oznaka_kolicin_v_podatkih] # x-axis: stations
y_globine = [[] for i in oznaka_kolicin_v_podatkih] # y-axis: depth
z_vrednosti = [[] for i in oznaka_kolicin_v_podatkih]   # z-axis: values

for ipostaja in range(moja_min_postaja, moja_max_postaja + 1):
    for iglobina in range(len(vse_mozne_globine)):
        if vse_mozne_globine[iglobina] in vse_mean[ipostaja - moja_min_postaja]['Press']:
            for ioznaka_kolicin_v_podatkih in range(len(oznaka_kolicin_v_podatkih)):
                x_postaje[ioznaka_kolicin_v_podatkih].append(ipostaja)
                y_globine[ioznaka_kolicin_v_podatkih].append(vse_mozne_globine[iglobina])
                z_vrednosti[ioznaka_kolicin_v_podatkih].append(\
                    vse_mean[ipostaja - moja_min_postaja][oznaka_kolicin_v_podatkih[ioznaka_kolicin_v_podatkih]]\
                        [vse_mean[ipostaja - moja_min_postaja]['Press'].tolist().index(vse_mozne_globine[iglobina])])

# maximum depth at each station
max_globine = [y_globine[0][i] for i in range(len(y_globine[0]) - 1) if y_globine[0][i] > y_globine[0][i+1]] + [y_globine[0][-1]]

# Initiate grid
ngridx = (moja_max_postaja - moja_min_postaja) * 3 + 1
ngridy = int((max(max_globine)*4 - 3) * 1) - 1
konture = True  # Set to False, if you don't want black contours


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    '''https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib, user unutbu'''
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
cmT = plt.get_cmap('twilight_shifted')
cmaps = {'temperatura':'plasma', 'slanost':'viridis', 'gostota':'rainbow'}

if len(kolicine) > 1:

    fig, axs = plt.subplots(nrows=len(kolicine), ncols=1, figsize=(6, 4 * len(kolicine)))

    for ikolicine in range(len(kolicine)):
        Z = np.transpose(seznam_kolicin_na_vseh_postajah_na_vseh_globinah[ikolicine, :, :])
        lokacije_tock_na_grafu, trikotniki = pripravi_lokacije_tock_in_trikotnike_za_tricontourf(np.transpose(Z))


        xi = np.linspace(np.min(lokacije_tock_na_grafu[:, 1]), np.max(lokacije_tock_na_grafu[:, 1]), ngridx)
        yi = np.linspace(np.min(lokacije_tock_na_grafu[:, 0]), np.max(lokacije_tock_na_grafu[:, 0]), ngridy)

        triang = tri.Triangulation(lokacije_tock_na_grafu[:, 1], lokacije_tock_na_grafu[:, 0], trikotniki)




        if kolicine[ikolicine] in slovar_mej_kolicin[datumska_mapa].keys():
            interpolator = tri.LinearTriInterpolator(triang, Z.flatten())
            Xi, Yi = np.meshgrid(xi, yi)
            zi = interpolator(Xi, Yi)
            ctf = axs[ikolicine].contourf(xi, yi, zi, cmap=cmaps[kolicine[ikolicine]], extend='both', levels=slovar_mej_kolicin[datumska_mapa][kolicine[ikolicine]])
            if konture:
                ct = axs[ikolicine].contour(xi, yi, zi, colors='k', levels=slovar_mej_kolicin[datumska_mapa][kolicine[ikolicine]])
        elif kolicine[ikolicine] == 'N2':
            interpolator = tri.LinearTriInterpolator(triang, np.log10(Z.flatten() + 10**-15))
            Xi, Yi = np.meshgrid(xi, yi)
            zi = 10**interpolator(Xi, Yi)
            ctf = axs[ikolicine].contourf(xi, yi, zi, norm=colors.LogNorm(vmin=10**-6, vmax=10**-1), cmap=truncate_colormap(plt.get_cmap('spring'), 0, 1, 6), extend='both', levels=[10**-6, 10**-5, 10**-4, 10**-3, 10**-2, 10**-1])
            if konture:
                ct = axs[ikolicine].contour(xi, yi, zi, colors='k', norm=colors.LogNorm(vmin=10**-6, vmax=10**-1), levels=[10**-6, 10**-5, 10**-4, 10**-3, 10**-2, 10**-1])
        else:
            interpolator = tri.LinearTriInterpolator(triang, Z.flatten())
            Xi, Yi = np.meshgrid(xi, yi)
            zi = interpolator(Xi, Yi)
            ctf = axs[ikolicine].contourf(xi, yi, zi, cmap='jet', extend='both', levels=10)
            if konture:
                ct = axs[ikolicine].contour(xi, yi, zi, colors='k', levels=10)
        cb = fig.colorbar(ctf, ax=axs[ikolicine])
        if konture:
            cb.add_lines(ct)
        cb.ax.tick_params(labelsize=12)
        cb.set_label(izpis_kolicin_v_grafih[ikolicine], fontsize=13)

        axs[ikolicine].invert_yaxis()
        axs[ikolicine].set_ylabel(slovar_kolicina_enota['Press'], fontsize=13)
        xticks = [i for i in range(moja_max_postaja - moja_min_postaja + 1)]
        axs[ikolicine].set_xticks(xticks)
        if ikolicine < len(kolicine) - 1:
            axs[ikolicine].set_xticklabels(['LK{:02d}'.format(xticks[i] + moja_min_postaja) for i in range(len(xticks))], fontsize=12)
        else:
            axs[ikolicine].set_xticklabels(['LK{:02d}\n{}'.format(xticks[i] + moja_min_postaja, str(casi_rut[i])[-8:-3]) for i in range(len(xticks))], fontsize=12)
        axs[ikolicine].set_yticks(ticks=[i for i in range(np.shape(Z)[0]) if vse_mozne_globine[i] % 2 == 0])
        axs[ikolicine].set_yticklabels([int(globina) for globina in vse_mozne_globine if globina % 2 == 0], fontsize=12)
        if ikolicine == 0:
            axs[ikolicine].set_title('Ruta {}, {}.{}.20{}'.format(ruta[-2:], dan, mesec, leto))
        # Plot the "sea floor"
        st_stolpcev = np.shape(Z)[1]
        max_globine = [np.shape(Z)[0] - 1 for stolpec in range(st_stolpcev)]
        for stolpec in range(st_stolpcev):
            ivrstica = -1
            while not Z[ivrstica, stolpec] < np.infty:
                max_globine[stolpec] = ivrstica + np.shape(Z)[0] - 1
                ivrstica -= 1
        axs[ikolicine].fill_between(range(moja_max_postaja - moja_min_postaja + 1), max_globine, [max(max_globine) for g in max_globine], color='k')
        axs[ikolicine].grid(linestyle='--', color='grey', linewidth=0.7, axis='x')

    axs[-1].set_xlabel('Postaja (ura)', fontsize=13)
    # The title also tells us the direction of the route
    if casi_rut[0] < casi_rut[-1]:
        axs[0].set_title(r'Ruta {}, {}.{}.20{} {}-{}, LK{:02d}$\rightarrow$LK{:02d}'.format(ruta[-2:], str(min(casi_rut))[8:10], mesec, leto, str(min(casi_rut))[-8:-3], str(max(casi_rut))[-8:-3], moja_min_postaja, moja_max_postaja), fontsize=14.5)
    else:
        axs[0].set_title(r'Ruta {}, {}.{}.20{} {}-{}, LK{:02d}$\rightarrow$LK{:02d}'.format(ruta[-2:], str(min(casi_rut))[8:10],
                                                                               mesec, leto, str(min(casi_rut))[-8:-3],
                                                                               str(max(casi_rut))[-8:-3],
                                                                               moja_max_postaja, moja_min_postaja), fontsize=14.5)
    if save_figures:
        if kolicine == ['temperatura', 'slanost', 'gostota', 'N2']:
            plt.savefig(mapa_za_shranjevanje_grafov + f'conditions_{datumska_mapa}_{ruta}_si.jpg', dpi=200)
        else:
            str_kolicine = ''
            for k in kolicine:
                str_kolicine += k + '_'
            plt.savefig(mapa_za_shranjevanje_grafov + f'conditions_{datumska_mapa}_{ruta}_{str_kolicine}si.jpg', dpi=200)
    plt.show()