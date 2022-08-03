'''There are three dictionaries that help us make decent labels in our plots.

The first dictionary goes as:
{'tlak':'Press', ...} (i.e. {'pressure':'Press', ...})

The second dictionary goes as:
{'Press':'Tlak [dbar]', ...} (i.e. {'Press':'Pressure [dbar]', ...})

The third dictionary goes as:
{'Press':Tlak', ...} (i.e. {'Press':'Pressure', ...})



What are these quantities (the quantities are values of the 1st dictionary and keys of the 2nd and 3rd dictionary):

a) TOB files from directory "epsilon" (see the TOB file's headers for units)
- Press: pressure
- vel: probe's vertical velocity
- epsilon1: epsilon from shear sensor 1 by batch job eallX.msb
- epsilon2: epsilon from shear sensor 2 by batch job eallX.msb
- peps: pseudo epsilon
- epsilon: epsilon from both shear sensors by batch job eallX.msb
- Tempcor: mean temperature
- Cond: conductivity
- sal: salinity
- sig_t: density anomaly
- LT-sigma: Thorpe length from sig_t
- Chl_A: chlorophyl A
- N^2: squared buoyancy frequency (haven't detected any negative values though...)
- N: suppose this is sqrt(N^2)
- Thermdiss: thermal dissipation
- LT-ntc: Thorpe length from the NTC temperature sensor (see b))
- K_rho: mass eddy diffusion coefficient

b) TOB files from directory "cutted" (see the TOB file's headers for units)
- NTC: temperature from the fast-response sensor NTC
- SHE1: shear fluctuation from sensor with the serial number 23
- SHE2: shear fluctuation from sensor with the serial number 25
- ACC: probe's (horizontal?) acceleration
- NTCHP: extra temperature channel with NTC outputs with enhanced sensitivity for large frequencies (need correction when calculating PSD)
- NTCAC: another temperature channel with NTC outputs
- Cmdc: conductivity (DC)
- Cmac: conductivity (AC)
- for others see a)

c) TOB files from directory "shear" (see the TOB file's headers for units)
- see a) and b)
'''

import pickle

prvi_slovar = {\
    'tlak':'Press',
    'hitrost padanja':'vel',\
    'epsilon1':'epsilon1',\
    'epsilon2':'epsilon2',\
    'psevdo epsilon':'peps',\
    'epsilon':'epsilon',\
    'temperatura':'Tempcor',\
    'prevodnost':'Cond',\
    'slanost':'sal',\
    'gostota':'sig_t',\
    'LT-sigma':'LT-sigma',\
    'klorofil':'Chl_A',\
    'N2':'N^2',\
    'N':'N',\
    'Thermdiss':'Thermdiss',\
    'LT-ntc':'LT-ntc',\
    'K_rho':'K_rho',\
    'ozmidov':'LO',\
    'LT-sigma/LO':'LT-sigma/LO',\
    'NTC':'NTC',\
    'SHE1':'SHE1',\
    'temperatura cutted':'Temp',\
    'SHE2':'SHE2',\
    'pospesek padanja':'ACC',\
    'NTCHP':'NTCHP',\
    'NTCAC':'NTCAC',\
    'Cmdc':'Cmdc',\
    'Cmac':'Cmac',\
    'kisik':'rawO2',\
    'Redox':'Redox',\
    'pH':'pH',\
    'cas padanja':'Time',\
    'strizenje1':'shear1',\
    'strizenje2':'shear2',\
    'psevdo strizenje':'pshear',\
    'log_epsilon_MLE_mean':'log_epsilon_MLE_mean',\
    'log_epsilon_MLE_1':'log_epsilon_MLE_1',\
    'log_epsilon_MLE_2':'log_epsilon_MLE_2',\
    'log_epsilon_IM_mean':'log_epsilon_IM_mean',\
    'log_epsilon_IM_1':'log_epsilon_IM_1',\
    'log_epsilon_IM_2':'log_epsilon_IM_2',\
    'Re_b':'Re_b',\
    'log_peps':'log_peps'
    }

pickle.dump(prvi_slovar, open('slovar_kolicina_kolicina.p', 'wb'))



drugi_slovar = {\
    'Press':'Tlak [dbar]',\
    'vel':'Hitrost padanja sonde [dbar/s]',\
    'epsilon1': r'$\epsilon_{1}$ [W/kg]',\
    'epsilon2': r'$\epsilon_{2}$ [W/kg]',\
    'peps':r'Psevdo $\epsilon$ [W/kg]',\
    'epsilon':r'$\epsilon$ [W/kg]',\
    'Tempcor':r'Temperatura [$\degree$C]',\
    'Cond':'Specifična prevodnost [mS/cm]',\
    'sal':'Slanost [PSU]',\
    'sig_t':r'Anom. gostota $[\mathrm{kg}/\mathrm{m}^{3}]$',\
    'LT-sigma':'LT-sigma [dbar]',\
    'Chl_A':r'Klorofil A $[\mu g/L]$',\
    'N^2':r'$N^2$ $[\mathrm{s}^{-2}]$',\
    'N':r'$N$ $[\mathrm{s}^{-1}]$',\
    'Thermdiss':r'Thermdiss $[K^2$/s]',\
    'LT-ntc':'LT-ntc [dbar]',\
    'K_rho':r'$K_{\rho}~[\log{(\mathrm{m}/\mathrm{s}^2)]$',\
    'LO':'LO [m]',\
    'LT-sigma/LO':'LT-sigma/LO',\
    'NTC':r'NTC $[\degree \mathrm{C}]$',\
    'SHE1':'SHE1',\
    'Temp':r'Temperatura [$\degree$C]',\
    'SHE2':'SHE2',\
    'ACC':r'Pospesek padanja sonde [m/$\mathrm{s}^2$]',\
    'NTCHP':'NTCHP [$\degree$C]',\
    'NTCAC':'NTCAC [$\degree$K]',\
    'Cmdc':'Cmdc [mS/cm]',\
    'Cmac':'Cmac [mS/cm]',\
    'rawO2':'Kisik [mV]',\
    'Redox':'Redox [mV]',\
    'pH':'pH',\
    'Time':'Čas padanja sonde [s]',\
    'shear1':r'Striženje 1 $[\mathrm{s}^{-1}]$',\
    'shear2':r'Striženje 2 $[\mathrm{s}^{-1}]$',\
    'pshear':r'Psevdo striženje $[\mathrm{s}^{-1}]$',\
    'log_epsilon_MLE_mean':r'$\log{(\epsilon\,\mathrm{[W/kg]})}$',\
    'log_epsilon_MLE_1':r'$\log{(\epsilon\,\mathrm{[W/kg]})}$',\
    'log_epsilon_MLE_2':r'$\log{(\epsilon\,\mathrm{[W/kg]})}$',\
    'log_epsilon_IM_mean':r'$\log{(\epsilon\,\mathrm{[W/kg]})}$',\
    'log_epsilon_IM_1':r'$\log{(\epsilon\,\mathrm{[W/kg]})}$',\
    'log_epsilon_IM_2':r'$\log{(\epsilon\,\mathrm{[W/kg]})}$',\
    'Re_b':r'$\mathrm{Re}_b$',\
    'log_peps':r'$\log{(\mathrm{pseudo }\epsilon\mathrm{ [W/kg]})}$'
    }

pickle.dump(drugi_slovar, open('slovar_kolicina_enota.p', 'wb'))


tretji_slovar = {\
    'Press':'Tlak',\
    'vel':'Hitrost padanja sonde',\
    'epsilon1': r'$\epsilon_{1}$',\
    'epsilon2': r'$\epsilon_{2}$',\
    'peps':r'Psevdo $\epsilon$',\
    'epsilon':r'$\epsilon$',\
    'Tempcor':r'Temperatura',\
    'Cond':'Specifična prevodnost',\
    'sal':'Slanost',\
    'sig_t':r'Anomalna gostota',\
    'LT-sigma':'LT-sigma',\
    'Chl_A':r'Klorofil A',\
    'N^2':r'$N^2$',\
    'N':r'$N$',\
    'Thermdiss':r'Thermdiss',\
    'LT-ntc':'LT-ntc',\
    'K_rho':r'$K_{\rho}',\
    'LO':'LO',\
    'LT-sigma/LO':'LT-sigma/LO',\
    'NTC':r'NTC',\
    'SHE1':'SHE1',\
    'Temp':r'Temperatura',\
    'SHE2':'SHE2',\
    'ACC':r'Pospesek padanja sonde',\
    'NTCHP':'NTCHP',\
    'NTCAC':'NTCAC',\
    'Cmdc':'Cmdc',\
    'Cmac':'Cmac',\
    'rawO2':'Kisik',\
    'Redox':'Redox',\
    'pH':'pH',\
    'Time':'Čas padanja sonde',\
    'shear1':r'Striženje 1',\
    'shear2':r'Striženje 2',\
    'pshear':r'Psevdo striženje',\
    'log_epsilon_MLE_mean':r'$\log{(\epsilon)}$',\
    'log_epsilon_MLE_1':r'$\log{(\epsilon)}$',\
    'log_epsilon_MLE_2':r'$\log{(\epsilon)}$',\
    'log_epsilon_IM_mean':r'$\log{(\epsilon)}$',\
    'log_epsilon_IM_1':r'$\log{(\epsilon)}$',\
    'log_epsilon_IM_2':r'$\log{(\epsilon)}$',\
    'Re_b':'Vzgonsko Reynoldsovo število',\
    'log_peps':r'$\log{(\mathrm{pseudo }\epsilon)}$'
    }

pickle.dump(tretji_slovar, open('slovar_kolicina_brez_enot.p', 'wb'))