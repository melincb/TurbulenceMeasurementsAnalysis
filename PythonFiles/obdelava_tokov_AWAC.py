'''Calculates and plots the vertical shear (flow velocity data are only available for 20.4.2009)'''

import numpy as np
import matplotlib.pyplot as plt
import datetime



mapa_s_tokovi = r'..\tokovi_AWAC\Luka_18dec2008_meritve_AWAC\ '[:-1]    # Where is the data on currents

mapa_za_shranjevanje_grafov =  'Saved_figures\ '[:-1]    # Where it saves the figures
save_figures = False    # To save the plot example, set this to True, else False
English_labels = True   # Set to False for Slovenian labels

# time of interest
datetime_min_max = (datetime.datetime(year=2009, month=4, day=18, hour=23, minute=59), datetime.datetime(year=2009, month=4, day=21, hour=9, minute=5))
# time of epsilon measurements
LK_min_max = (datetime.datetime(year=2009, month=4, day=20, hour=7, minute=30), datetime.datetime(year=2009, month=4, day=21, hour=8, minute=14))


# THE BASIC MANUAL SETTINGS END HERE



used_datetime = []
flow_v1 = []    # flow E
flow_v2 = []    # flow N
flow_v3 = []    # flow U
used_pressure_at_bottom = []
used_lines = []

with open(mapa_s_tokovi + 'LukaKP02.sen', 'r') as datoteka_sen:
    count_vrstica = 0
    for vrstica in datoteka_sen:
        a = vrstica.split()
        month, day, year, hour, minute, pressure = int(a[0]), int(a[1]), int(a[2]), int(a[3]), int(a[4]), float(a[13])
        this_datetime = datetime.datetime(year=year, month=month, day=day, hour=hour, minute=minute)
        if this_datetime > datetime_min_max[0] and this_datetime < datetime_min_max[1]:
            used_datetime.append(this_datetime)
            used_pressure_at_bottom.append(pressure)
            used_lines.append(count_vrstica)
        count_vrstica += 1

with open(mapa_s_tokovi + 'LukaKP02.v1', 'r') as datoteka_v1:
    count_vrstica = 0
    for vrstica in datoteka_v1:
        if count_vrstica in used_lines:
            flow_v1.append([float(f) for f in vrstica.split()])
        count_vrstica += 1

with open(mapa_s_tokovi + 'LukaKP02.v2', 'r') as datoteka_v2:
    count_vrstica = 0
    for vrstica in datoteka_v2:
        if count_vrstica in used_lines:
            flow_v2.append([float(f) for f in vrstica.split()])
        count_vrstica += 1


with open(mapa_s_tokovi + 'LukaKP02.v3', 'r') as datoteka_v3:
    count_vrstica = 0
    for vrstica in datoteka_v3:
        if count_vrstica in used_lines:
            flow_v3.append([float(f) for f in vrstica.split()])
        count_vrstica += 1

used_pressures = [] # depth of cell
used_shears = []    # vertical shear
used_pressures2 = []    # the same quantity as used_pressures, except it has a different shape
used_shears2 = []   # the same quantity as used_shears, except it has a different shape
used_idt2 = []  # used datetime indexes
min_pressures = []
max_pressures = []
number_of_these_used_pressures = []

for idt in range(len(used_datetime)):
    these_used_pressures = []
    these_used_shears = []
    for jf in range(len(flow_v1[idt])-1):
        press_center = used_pressure_at_bottom[idt] - 1.0 - 0.5*jf - 0.25    # Dist. to center of 1st cell was 100cm, cell width is 50cm, center between two cells is 25cm above center of the bottom cell
        min_pressure_of_upper_cell = press_center - 0.75            # 25cm to center of upper cell + 50cm to top of the triangle of upper cell
        if min_pressure_of_upper_cell > 0:  # On 20th Apr 2009 there were no surface waves
            dv1 = flow_v1[idt][jf + 1] - flow_v1[idt][jf]
            dv2 = flow_v2[idt][jf + 1] - flow_v2[idt][jf]
            dv3 = flow_v3[idt][jf + 1] - flow_v3[idt][jf]
            shear = np.sqrt(dv1**2 + dv2**2)/0.5    #  vertical shear
            these_used_pressures.append(press_center)
            these_used_shears.append(shear)
            used_shears2.append(shear)
            used_pressures2.append(press_center)
            used_idt2.append(idt)
    number_of_these_used_pressures.append(len(these_used_pressures))
    used_pressures.append(these_used_pressures)
    used_shears.append(these_used_shears)
    min_pressures.append(min(these_used_pressures))
    max_pressures.append(max(these_used_pressures))



print('min number of cells', min(number_of_these_used_pressures), 'max number of cells', max(number_of_these_used_pressures))
# Plot just those cells that are always full
for idt in range(len(used_datetime)):
    if number_of_these_used_pressures[idt] > min(number_of_these_used_pressures):
        used_pressures[idt] = used_pressures[idt][:min(number_of_these_used_pressures)]
        used_shears[idt] = used_shears[idt][:min(number_of_these_used_pressures)]
        used_pressures[idt] = used_pressures[idt][:min(number_of_these_used_pressures)]
        min_pressures[idt] = min(used_pressures[idt])


x_datetime = [[i for sh in used_shears[i]] for i in range(len(used_datetime))]
starting_LK_measurements_fraction = (LK_min_max[0] - used_datetime[0]) / (used_datetime[-1] - used_datetime[0])
ending_LK_measurements_fraction = (LK_min_max[1] - used_datetime[0]) / (used_datetime[-1] - used_datetime[0])





plt.figure(figsize=(12,4))
plt.tricontourf(np.array(used_idt2), np.array(used_pressures2), np.array(used_shears2), extend='both', cmap='gnuplot', levels=[0.06, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48])
plt.gca().invert_yaxis()
cb = plt.colorbar(pad=0.02)
if not English_labels:
    cb.set_label(label=r'Vertikalno striženje $[\mathrm{s}^{-1}]$', fontsize=16)
else:
    cb.set_label(label=r'Vertical shear $[\mathrm{s}^{-1}]$', fontsize=16)
cb.ax.tick_params(labelsize=14)
plt.fill_between([i for i in range(len(used_datetime))], min_pressures, [min(min_pressures) for mp in min_pressures], color='lightgray')
plt.fill_between([i for i in range(len(used_datetime))], max_pressures, [max(max_pressures) for mp in max_pressures], color='lightgray')
plt.xticks(ticks=[i for i in range(len(used_datetime)) if used_datetime[i].hour%4==0 and used_datetime[i].minute==0], labels=[f'{dt.hour:02d}' for dt in used_datetime if dt.hour%4 == 0 and dt.minute == 0], fontsize=14)
plt.yticks(fontsize=14)
if not English_labels:
    plt.ylabel('Tlak [dbar]', fontsize=16)
    plt.title("Vertikalno striženje iz tokomera", fontsize=18)
    plt.xlabel('Ura meritve', fontsize=16)
else:
    plt.ylabel('Pressure [dbar]', fontsize=16)
    plt.title("Vertical shear from currentmeter", fontsize=18)
    plt.xlabel('Measurement time', fontsize=16)

plt.annotate('', xy=(starting_LK_measurements_fraction, -0.002), xycoords='axes fraction', xytext=(ending_LK_measurements_fraction, -0.002), arrowprops=dict(arrowstyle="<->", color='r'))
plt.annotate('', xy=(starting_LK_measurements_fraction, 0.998), xycoords='axes fraction', xytext=(ending_LK_measurements_fraction, 0.998), arrowprops=dict(arrowstyle="<->", color='r'))
plt.axvline(starting_LK_measurements_fraction*(len(used_datetime)-1), linestyle='--', color='r')
plt.axvline(ending_LK_measurements_fraction*(len(used_datetime)-1), linestyle='--', color='r')
plt.tight_layout()
if save_figures:
    plt.savefig(mapa_za_shranjevanje_grafov + 'strizenje_AWAC.jpg', dpi=300)


plt.show()