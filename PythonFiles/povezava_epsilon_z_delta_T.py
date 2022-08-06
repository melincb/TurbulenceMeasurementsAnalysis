'''Looks at the correlation between the Temperature_water (mean of the entire water column) - Temperature_air and
epsilon values, which are measured at least 3 dbar below the sea level and at least 2 dbar above the sea floor (these
numbers can be changed in the code).
'''
import openpyxl
from pathlib import Path
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob, os
from scipy.optimize import curve_fit

savefigs = False # Whether it saves figures or not
mapa_za_shranjevanje_grafov = 'Saved_figures\ '[:-1]    # Where it saves the figures

# dates
datumske_mape = ['all_2008_08_21', 'all_2008_11_18', 'all_2009_04_20']
# sunrise and sunset per each day
dnevni_zahod_vzhod = [[datetime.datetime(2008, 8, 21, 20, 3), datetime.datetime(2008, 8, 22, 6, 13)],\
                        [datetime.datetime(2008,11,18,16,31), datetime.datetime(2008,11,19,7,10)],\
                      [datetime.datetime(2009,4,20,19,57), datetime.datetime(2009,4,21,6,9)]]

which_epsilon = 'log_epsilon_MLE_mean'  # options: 'log_epsilon_MLE_mean', 'log_epsilon_MLE_1', 'log_epsilon_MLE_2',
                                        # 'log_epsilon_IM_mean', 'log_epsilon_IM_1', 'log_epsilon_IM_2'

folder_with_xlsx_files = '..\Meritve_ARSO'

izvirna_mapa = '..'
mapa_s_skriptami = os.getcwd()     # directory with scripts


# The rest of the basic manual settings is about the excel files

filenames = ['Koper_T2m_avg2008.xlsx', 'Koper_T2m_nov2008.xlsx', 'Koper_T2m_apr2009.xlsx']
sheets = ['Koper-Luka', 'Koper-Kapitanija']

datetime_min_max = {'Koper_T2m_apr2009.xlsx':(datetime.datetime(year=2009, month=4, day=18, hour=23, minute=59), datetime.datetime(year=2009, month=4, day=21, hour=9, minute=5)),\
                    'Koper_T2m_avg2008.xlsx':(datetime.datetime(year=2008, month=8, day=19, hour=23, minute=59), datetime.datetime(year=2008, month=8, day=22, hour=9, minute=5)),\
                    'Koper_T2m_nov2008.xlsx':(datetime.datetime(year=2008, month=11, day=16, hour=23, minute=59), datetime.datetime(year=2008, month=11, day=19, hour=9, minute=5))\
                    }

LM_min_max = {'Koper_T2m_apr2009.xlsx':(datetime.datetime(year=2009, month=4, day=20, hour=7, minute=30), datetime.datetime(year=2009, month=4, day=21, hour=8, minute=14)),\
            'Koper_T2m_avg2008.xlsx':(datetime.datetime(year=2008, month=8, day=21, hour=7, minute=9), datetime.datetime(year=2008, month=8, day=22, hour=7, minute=59)),\
            'Koper_T2m_nov2008.xlsx':(datetime.datetime(year=2008, month=11, day=18, hour=6, minute=57), datetime.datetime(year=2008, month=11, day=19, hour=8, minute=22)) \
              }


# THE BASIC MANUAL SETTINGS END HERE


def najblizji_datetime_z_leve_in_desne(cas, seznam):
    '''Find the closest datetime both before and after time cas.

    cas ... time
    seznam ... list'''
    min_delta_plus = [np.infty, 0]
    min_delta_minus = [np.infty, 0]
    for dt in seznam:
        delta = (dt - cas).total_seconds()
        if delta > 0:
            if delta < min_delta_plus[0]:
                min_delta_plus[0] = delta
                min_delta_plus[1] = dt
        elif delta < 0:
            if abs(delta) < min_delta_minus[0]:
                min_delta_minus[0] = abs(delta)
                min_delta_minus[1] = dt
        else:
            min_delta_plus = [0, dt]
            min_delta_minus = [0, dt]
            break

    return min_delta_plus, min_delta_minus



used_datetimes = []
avg_Ts = []

for ifn in range(len(filenames)):
    file_name = filenames[ifn]
    xslx_file = Path(folder_with_xlsx_files, file_name)

    wb_obj = openpyxl.load_workbook(xslx_file)

    this_datetime_min_max = datetime_min_max[file_name]
    this_datetime_min_max_LK = LM_min_max[file_name]
    sheet_Luka = wb_obj[sheets[0]]
    sheet_Kapitanija = wb_obj[sheets[1]]
    used_datetime = []
    used_mean_T = []


    row = 2
    while True:
        try:
            this_datetime = sheet_Luka[f'C{row}'].value # same value at sheet_Kapitanija[f'{C{row}'].value
            if this_datetime > this_datetime_min_max[0] and this_datetime < this_datetime_min_max[1]:
                used_datetime.append(this_datetime)
                used_mean_T.append(np.mean([sheet_Luka[f'D{row}'].value, sheet_Kapitanija[f'D{row}'].value]))
        except: # EOF
            break
        row += 1
    used_datetimes.append(used_datetime)
    avg_Ts.append(used_mean_T)


def premica(x, a, b):
    '''Linear fit'''
    return a*x + b


for idatumska_mapa in range(len(datumske_mape)):    # Loop through all dates
    datumska_mapa = datumske_mape[idatumska_mapa]
    used_dt = used_datetimes[idatumska_mapa]
    used_avg_T = avg_Ts[idatumska_mapa]
    dict_used_dt_avg_T = {used_dt[i]:used_avg_T[i] for i in range(len(used_dt))}
    delta_T_epsilon_tega_dne = [[], []]
    izracuni_epsilon_in_Re_b = pickle.load(open(f'slovar_epsilon_za_risat_NA_1m_{datumska_mapa}.p', 'rb'))
    izracuni_epsilon_in_Re_b_T = pickle.load(open('grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_T' + datumska_mapa + '.p', 'rb'))
    bs = '\ '[0]
    tipska_mapa = 'epsilon'
    mapa_s_podatki = izvirna_mapa + bs + datumska_mapa + bs + tipska_mapa
    os.chdir(mapa_s_podatki)
    datoteke = glob.glob('*.tob')
    ip = 0
    plt.figure(figsize=(10, 8))
    for postaja in ['LK01','LK02','LK03','LK04','LK05','LK06','LK07']:  # Loop through all stations
        delta_T_epsilon_tega_dne = [[], []] # delta_T and epsilon for this night
        for postajaruta in izracuni_epsilon_in_Re_b.keys():
            if postajaruta[:4] == postaja:
                temperature_stolpca = []    # temperatures in this water column
                for iprg in range(len(izracuni_epsilon_in_Re_b_T['tempcor'])):
                    if izracuni_epsilon_in_Re_b_T['postaja_ruta_globina'][iprg][:len(postajaruta)] == postajaruta:
                        # if you want to limit your temperature measurements, use if float(izracuni_epsilon_in_Re_b_T['postaja_ruta_globina'][iprg][13:]) is something (float(izracuni_epsilon_in_Re_b_T['postaja_ruta_globina'][iprg][13:]) is the depth of the measurement)
                        temperature_stolpca.append(izracuni_epsilon_in_Re_b_T['tempcor'][iprg])
                T_stolpca = np.mean(temperature_stolpca)
                # Find time of the medial cast
                datoteke = glob.glob('*' + postajaruta + '.tob')
                srednja_meritev = datoteke[len(datoteke)//2]
                for iznak in range(0, len(srednja_meritev)):
                    if srednja_meritev[iznak] == '_':
                        datetime_srednja = datetime.datetime(year=int('20' + srednja_meritev[iznak - 6:iznak - 4]), \
                                                    month=int(srednja_meritev[iznak - 4:iznak - 2]), \
                                                    day=int(srednja_meritev[iznak - 2:iznak]), \
                                                    hour=int(srednja_meritev[iznak + 1:iznak + 3]), \
                                                    minute=int(srednja_meritev[iznak + 3:iznak + 5]))
                        break
                # If the medial cast was done in the night time (between sunrise and sunset), find the closest air-temperature measurements and compute the difference in temperature
                if datetime_srednja > dnevni_zahod_vzhod[idatumska_mapa][0] and datetime_srednja < dnevni_zahod_vzhod[idatumska_mapa][1]:   # NIGHT!
                    levi, desni = najblizji_datetime_z_leve_in_desne(datetime_srednja, dict_used_dt_avg_T)  # the closest datetime from the left and from the right
                    if levi[0] == 0: # in this case also desni[0]=0
                        delta_T = T_stolpca - dict_used_dt_avg_T[levi[1]]   # temperature difference (T_water - T_air)
                    else:
                        leva_T_zraka, desna_T_zraka = dict_used_dt_avg_T[levi[1]], dict_used_dt_avg_T[desni[1]]
                        leva_oddaljenost, desna_oddaljenost = levi[0], desni[0]
                        T_zraka = (leva_T_zraka * desna_oddaljenost + desna_T_zraka * leva_oddaljenost)/(leva_oddaljenost + desna_oddaljenost)
                        delta_T = T_stolpca - T_zraka   # temperature difference (T_water - T_air)
                    delta_T_epsilon_tega_dne[0].append(delta_T)

                    # HERE YOU CAN SELECT THE PREFERED DEPTH OF EPSILON MEASUREMENTS
                    delta_T_epsilon_tega_dne[1].append(np.mean(izracuni_epsilon_in_Re_b[postajaruta][which_epsilon][2:-2]))
                    # izracuni_epsilon_in_Re_b[postajaruta][which_epsilon][a:-b] contains measurements which
                    # are at least a+1 dbar below the surface and at least b dbar above the sea floor

        plt.subplot(3,3,ip+1)
        plt.scatter(delta_T_epsilon_tega_dne[0], delta_T_epsilon_tega_dne[1])
        delta_T = np.array(delta_T_epsilon_tega_dne[0])
        epsilon = np.array(delta_T_epsilon_tega_dne[1])
        popt, pcov = curve_fit(premica, delta_T, epsilon)
        y_bar = np.mean(epsilon)
        ss_tot = sum((epsilon - y_bar)**2)
        ss_res = sum((epsilon - premica(delta_T, *popt))**2)
        ss_r = sum((premica(delta_T, *popt) - y_bar)**2)
        cmat = np.cov(delta_T, epsilon)
        plt.plot(np.linspace(min(delta_T), max(delta_T)), premica(np.linspace(min(delta_T), max(delta_T)), *popt))
        plt.annotate('R = {}'.format(round(np.corrcoef(delta_T, epsilon)[0,1],2)), xy=(0.1, 0.8), xycoords='axes fraction')
        plt.title(postaja)
        plt.xlabel('delta T')
        plt.ylabel(r'$\log{(\epsilon)}$')
        ip += 1
    plt.suptitle(datumska_mapa)
    plt.tight_layout()
    if savefigs:
        plt.savefig(mapa_za_shranjevanje_grafov + 'korelacija_temperatura_epsilon_' + datumska_mapa)
    plt.show()



    os.chdir(mapa_s_skriptami)