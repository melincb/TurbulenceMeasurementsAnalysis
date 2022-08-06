This is a selection of the files that I used for my MSc thesis *Analysis of turbulence measurements in Basin II of Port of Koper* (University of Ljubljana, Faculty of Mathematics and Physics, 2022).

The point of these files is to determine and explain the **turbulent kinetic energy dissipation rate** (in short: *epsilon*).

The turbulence measurements were obtained by Sea&Sun's probe **MSS90**. The measurements with the velocitymeter were obtained by Nortek AS's **AWAC** (unfortunatelly, these measurements were performed around 1km away from the turbulence measurements, so they weren't very useful for our study). The measurements which weren't done by MSS90 or AWAC were provided by the Slovenian Environment Agency and are stored in the directory **Meritve_ARSO**.

All Python files are stored in **PythonFiles** directory. They work properly with Python 3.8.

The data from MSS90 probe is stored in **all_2008_08_21**, **all_2008_11_18**, and **all_2009_04_20** directories (the directory names go as *all_yyyy_mm_dd*). The MSS90 probe produces *MRD* files, which have to be manually cutted at the top and bottom of the water column (on top to cut out the waiting period, in which the sensors settle down, and on bottom to cut out the data after touching the sea floor). Due to too their too large size, only the *MRD* files for the 1st route at the station LK01 of 21.8.2008 are uploaded (see **all_2008_08_21**). The cutted *MRD* files are then storred in the *TOB* format in **all_yyyy_mm_dd/cutted** (again, only uploaded for 1st route at the station LK01 of 21.8.2008). Each *TOB* file's name goes as *prefixyymmdd_hhhh_station_route.tob*. The *TOB* files such as those from **all_yyyy_mm_dd/cutted** directories are used to calculate the temperature fluctuation gradient spectra in **priprava_podatkov_za_izracun_epsilon_T.py** (see its function *PSD_T()*). To get the *TOB* files from **all_yyyy_mm_dd/shear**, one has to run the batch job *_shear_c.msb* with suitable callibration parameters on files from **all_yyyy_mm_dd/cutted** (this running requires *msspro.exe*, which is not included). The *TOB* files from **all_yyyy_mm_dd/shear** (again, only uploaded for 1st route at the station LK01 of 21.8.2008) contain shear fluctuation (and some other) data from both shear sensors. Such files are used to calculate the shear fluctuation spectra in **priprava_podatkov_za_izracun_epsilon_shear.py** (see its function *PSD_shear()*). To get the *TOB* files from **all_yyyy_mm_dd/epsilon**, one has to run the batch job *eallX.msb* on files from **all_yyyy_mm_dd/shear** (this running also requires *msspro.exe*). These files are included for all casts and are used to determine the conditions (temperature, salinity, density, buoyancy frequency, etc., see the files' headers) by function *Branje_datotek()* (see **Branje_datotek.py**), by **Risanje_iste_tocke_ob_razlicnih_rutah.py** (plots the conditions at one route in all stations), and by **Risanje_iste_tocke_ob_razlicnih_rutah.py** (plots the conditions at the same station at all routes (WARNING: poor interpolation, **Risanje_celega_kanala_na_isti_ruti.py** is much more reliable)). These files also contain the epsilon calculations from *eallX.msb*.

## How to get *epsilon* from shear data with other methods than MSSpro:
1. Make sure your data is storred in similar files as those from **all_2008_08_21/shear**
2. Run **priprava_podatkov_za_izracun_epsilon_shear.py** (see its documentation on the top of the code, set the manual settings in the beggining of the code and perhaps some other settings in the code). This program computes PSD of shear fluctuations and merges it two ways: 
	- computes geometrical mean of all PSD from bins with the same bin centers (i.e. central pressures of the bin) and the same combination *station-route* (the output is storred in a pickle file such as *shear12_outputs_povprecenje_VSEH_po_ruti_in_postajiall_2008_11_18.p*)
	- computes geometrical mean of all PSD from bins with bin centers within the depth range of 1 dbar (e.g. all PSDs with bin centers in [1.50 dbar, 2.50 dbar) -> PSD with bin center 2 dbar,  all PSDs with bin centers in [2.50 dbar, 3.50 dbar) -> PSD with bin center 3 dbar, etc.) and the same combination *station-route* (the output is storred in a pickle file such as *shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1mall_2008_11_18.p*).
3. Run **priprava_izracunov_epsilon_shear_za_risanje.py** (see its documentation on the top of the code, set the manual settings in the beggining of the code and perhaps some other settings in the code). This program fits each PSD from files of type *shear12_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1mall_yyyy_mm_dd.p* or *shear12_outputs_povprecenje_VSEH_po_ruti_in_postajiall_yyyy_mm_dd.p* to the Nasmyth curve using *Maximum Likelihood Estimate* or *Integral Method* (see **fitting_procedures_for_Nasmyth_curve.py**). The output is a pickle file, such as *slovar_epsilon_za_risat_NA_1m_all_2008_08_21.p*, which contains:
	- *epsilon* at each combination *station-route-depth*:
		1) using Maximum Likelihood Estimate on the data from shear sensor 1,
		2) using Maximum Likelihood Estimate on the data from shear sensor 2,
		3) arithmetical average of a. and b.
		4) using Integral Method on the data from shear sensor 1,
		5) using Integral Method on the data from shear sensor 2,
		6) arithmetical average of d. and e.
	- buoyancy Reynolds number, calculated from c. in the previous item
	- some other *shear-epsilon* related stuff.
4. Plot the shear data using *slovar_epsilon_za_risat_NA_1m_all_yyyy_mm_dd.p*. Here you have two options:
	1. Run **Risanje_epsilon_shear_na_vseh_tockah_ob_isti_ruti.py** (see its documentation) (recommended). This program plots a 2-dimensional *epsilon* crossection of the basin at a fixed route.
	2. Run **Risanje_epsilon_shear_na_isti_tocki_ob_razlicnih_rutah.py** (WARNING: the interpolation makes spurious unrealistic plots). This program plots a Hovmoller diagram of *epsilon* at a fixed station at all routes.

## How to get *epsilon* from temperature data:
1) Make sure your data is storred in similar files as those from **all_2008_08_21/cutted**
2) Run **priprava_podatkov_za_izracun_epsilon_T.py** (see its documentation on the top of the code, set the manual settings in the beggining of the code and perhaps some other settings in the code). This program computes PSD of temperature fluctuation gradient and merges it two ways: 
	- computes geometrical mean of all PSD from bins with the same bin centers (i.e. central pressures of the bin) and the same combination *station-route* (the output is storred in a pickle file such as *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postajiall_2008_11_18.p*)
	- computes geometrical mean of all PSD from bins with bin centers within the depth range of 1 dbar (e.g. all PSDs with bin centers in [1.5 dbar, 2.5dbar) -> PSD with bin center 2 dbar,  all PSDs with bin centers in [2.5 dbar, 3.5dbar) -> PSD with bin center 3 dbar, etc.) and the same combination *station-route* (the output is storred in a pickle file such as *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1mall_2008_11_18.p*).
3) Run **dejanski_izracun_epsilon_T.py** (define the input pickle file). This program fits each PSD from files of type *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1mall_yyyy_mm_dd.p* or *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postajiall_yyyy_mm_dd.p* using the *iterative process with Maximum Likelihood Estimate and the integration* (see **fitting_procedures_for_Batchelor_curve.py**). The output is a pickle file, such as *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_Tall_yyyy_mm_dd.p*, which contains:
	- *shear-epsilon* from *eallX.msb*
	- *epsilon* and thermall dissipation from iterative proces with 6 degrees of freedom
	- *epsilon* and thermall dissipation from iterative proces with 12 degrees of freedom
	- the minimum wavenumber of the range of wavenumbers that was used to fit (it is independent from degrees of freedom, maximum wavenumber is always 100 cpm)
	- some other *temperature-gradient-epsilon* related stuff.
 

## Comparison of *epsilon* outputs from different techniques
### How to simultaneously plot shear fluctuation spectra from both sensors and temperature fluctuation gradient spectrum with the best fitting theoretical curves:
To do so, you need both *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_Tall_yyyy_mm_dd.p* and *slovar_epsilon_za_risat_NA_1m_all_yyyy_mm_dd.p*.
Simply run **risanje_spektrov.py** (see its documentation).
### How to make histograms to compare multiple *epsilon* outputs:
- To compare *epsilon* from shear using many different techniques (*Maximum Likelihood Estimate*, *Integral Method* and *MSSpro*'s *eallX.msb* on each sensor and their final output using both sensors) run **primerjava_mojih_epsilon_shear_z_msspro.py*.
- To compare *epsilon* from temperature (*iterative procedure using MLE*) with *epsilon* from shear (various techniques) run **primerjava_mojih_epsilon_grad_T_z_msspro.py**.


## Monte Carlo simulations
To run Monte Carlo simulations, run file **MC_simulacija_spektrov.py** (see its documentation on the top of the code and set the manual settings).

## Shear from velocitymeter
To compute shear from the velocitymeter, run **obdelava_tokov_AWAC.py**. To see the proper data, look in **tokovi_AWAC** directory.

## Correlation between *night-time epsilon* in the center of the water column and the difference between the water column temperature and the air temperature
To seek correlation between *night-time epsilon* in the center of the water column and the difference between the water column temperature and the air temperature, run **povezava_epsilon_z_delta_T.py**. To successfully run it, you need:
1) pickle file(s) with *epsilon*, such as *slovar_epsilon_za_risat_NA_1m_all_2008_08_21.p*,
2) excel file(s) with air temperature data, such as *Meritve_ARSO/Koper_T2m_avg2009.xslx*,
3) some sort of data, that includes temperature of water at each combination *station-route-depth* (e.g. I use *grad_T_outputs_povprecenje_VSEH_po_ruti_in_postaji_NA_1m_DODAN_EPSILON_GRAD_Tall_yyyy_mm_dd.p*, however you can also change the way *slovar_epsilon_za_risat_NA_1m_all_yyyy_mm_dd.p* is produced (just add this temperature to the final dictionary) so it contains this temperature),
4) select the time of interest (i.e. the night-time in your location).

## Measurements by Slovenian Environment Agency
To plot the measurements provided by Slovenian Environment Agency (the Rižana river flux in Kubed, sea level in Koper, wind and temperature), run **pretok_Rizane_in_vodostaj_morja_in_veter_in_temperatura.py** (see its documentation on the top of the code).
