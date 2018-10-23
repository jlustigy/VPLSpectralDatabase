"""
Author: Jake Lustig-Yaeger
"""

import os
import sys
import numpy as np
import readsmart as rs
import pandas as pd
from numba import jit

from spectrum_class import Spectrum, write_spectra_csv, write_spectra_meta_csv, copy_raw_spectra

"""
Set the relative path to spectra from this script
"""
path = "spectrum_files/"

"""
Below is a huge list that contains all the spectral info and data
"""
spectra = [

    #######################
    # DIRECT SPECTROSCOPY #
    #######################

    Spectrum(observation="Direct",
             star="Sun",
             planet="Mars",
             description="Present-day Mars",
             meta="Generated using SMART",
             reference="Tyler Robinson",
             path_to_file=os.path.join(path, "mars_flx_refl.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Segura Earth",
             description="Clearsky; Ocean surface",
             meta="0.1 PAL O2",
             reference="Segura et al. 2003",
             path_to_file=os.path.join(path, "smart_kasting_krelove_1e-1PAL_clr_ocean_50_80000cm_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Segura Earth",
             description="Clearsky; Ocean surface",
             meta="1 PAL O2",
             reference="Segura et al. 2003",
             path_to_file=os.path.join(path, "smart_kasting_krelove_1PAL_clr_ocean_50_80000cm_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="F2V",
             planet="Segura Earth",
             description="Clearsky; Ocean surface",
             meta="1 PAL O2",
             reference="Segura et al. 2003",
             path_to_file=os.path.join(path, "smart_kasting_krelove_F2V_1PAL_clr_ocean_50_80000cm_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="G2V",
             planet="Segura Earth",
             description="Clearsky; Ocean surface",
             meta="0.1 PAL O2; 20 ppm CH4",
             reference="Segura et al. 2003",
             path_to_file=os.path.join(path, "smart_kasting_krelove_G2V_1e-1PAL_20ppmCH4_clr_ocean_50_80000cm_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="G2V",
             planet="Segura Earth",
             description="Clearsky; Ocean surface",
             meta="0.1 PAL O2; 100 ppm CH4",
             reference="Segura et al. 2003",
             path_to_file=os.path.join(path, "smart_kasting_krelove_G2V_1e-1PAL_100ppmCH4_clr_ocean_50_80000cm_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="K2V",
             planet="Segura Earth",
             description="Clearsky; Ocean surface",
             meta="1 PAL O2",
             reference="Segura et al. 2003",
             path_to_file=os.path.join(path, "smart_kasting_krelove_K2V_1PAL_clr_ocean_50_80000cm_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="EK Draconis (G0V)",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="IR; 2 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_EKDra_2barCO2_redone_cia_clr_ocean_ir_400_3000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="EK Draconis (G0V)",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="Vis; 2 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_EKDra_2barCO2_redone_cia_clr_ocean_solar_5700_30000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="EK Draconis (G0V)",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="IR; 0.2 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_EKDra_2e-1CO2_redone_fsfm_cia_clr_ocean_ir_400_3000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="EK Draconis (G0V)",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="Vis; 0.2 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_EKDra_2e-1CO2_redone_fsfm_cia_clr_ocean_solar_5700_30000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="EK Draconis (G0V)",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="IR; 0.02 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_EKDra_2e-2CO2_redone_cia_clr_ocean_ir_400_3000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="EK Draconis (G0V)",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="Vis; 0.02 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_EKDra_2e-2CO2_redone_cia_clr_ocean_solar_5700_30000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="IR; 0.2 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_Sun_2e-1CO2_redone_clr_ocean_ir_400_3000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Abiotic O2 & O3; High CO2",
             description="Clearsky; Ocean surface",
             meta="Vis; 0.2 bar CO2",
             reference="Segura et al. 2007",
             path_to_file=os.path.join(path, "smart_spectra_highCO2_Sun_2e-1CO2_redone_clr_ocean_solar_650_30000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Serratia Marscens Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_s_marscens_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Rubrobacter Radiotolerans Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_r_radiotolerans_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Rhodopseudomonas Palustris Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_r_palustris_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Rhodobacter Sphaeroides Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_r_sphaeroides_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Phaeobacter Inhibens Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_p_inhibens_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Micrococcus Luteus Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_m_luteus_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Janthinobacterium Lividum Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_j_lividum_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Halobacterium Salinarum Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_h_salinarum_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Deinoccocus Radiodurans Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_d_radiodurans_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Chlorobium Tepidum Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_c_tepidum_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Bacterial Surface Earth",
             description="Earth with Brevibacterium Aurantiacum Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_b_aurantiacum_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Snow Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_snow_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Red Algae Water Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_redalgae_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Ocean Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_ocean_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Limestone Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_limestone_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Kaolinite Soil Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_soil_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Halophile Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_halophile_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Halite Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_halite_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Gypsum Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_gypsum_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Grass Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_grass_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Conifers Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_conifers_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Basaltic Loam Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_soil2_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Bacterial Mat Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_bacterialmat_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Material Surface Earth",
             description="Earth with Acid Mine Drainage Surface",
             meta="Modern Earth Atmosphere",
             reference="Schwieterman et al. 2015",
             path_to_file=os.path.join(path, "new_earth_clearsky_hitran08_8333_100000cm_acidminedrainage_60sza_toa.flx")
              ),

    Spectrum(observation="Direct",
             star="Sun",
             planet="Earth",
             description="Earth Spectrum at Quadrature",
             meta="True disk-integrated Earth",
             reference="Robinson et al. 2011",
             path_to_file=os.path.join(path, "earth_quadrature_radiance_refl.flx")
              ),

    Spectrum(observation="Direct",
                    star="AD Leo (M3.5V)",
                    planet="Earth",
                    description="Quiet star; No clouds; 100% Ocean",
                    meta="1 PAL O2",
                    reference="Segura et al. 2005",
                    path_to_file=os.path.join(path, "smart_Mstar_T3100_1PAL_clr_ocean_50_30000cm_r60_only_toa.flx")
              ),

    Spectrum(observation="Direct",
                    star="AD Leo (M3.5V)",
                    planet="Earth",
                    description="Active star; No clouds; 100% Ocean",
                    meta="1 PAL O2",
                    reference="Segura et al. 2005",
                    path_to_file=os.path.join(path, "smart_Mstar_ADLeo_1PAL_clr_ocean_50_30000cm_r60only_toa.flx")
              ),

    Spectrum(observation="Direct",
	         star = "Sun",
	         planet = "High CO2; Low CH4",
             description = "50% CO2 with low volcanic outgassing of reduced gases (H2)",
	         meta = "100% Ocean; No clouds",
             reference = "Domagal-Goldman et al. 2014",
             path_to_file=os.path.join(path,"Sun_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "Epsilon Eridani (K2V)",
	         planet = "High CO2; Low CH4",
             description = "50% CO2 with low volcanic outgassing of reduced gases (H2)",
	         meta = "100% Ocean; No clouds",
             reference = "Domagal-Goldman et al. 2014",
             path_to_file=os.path.join(path,"K2V_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "GJ 876 (M4V)",
	         planet = "High CO2; Low CH4",
             description = "50% CO2 with low volcanic outgassing of reduced gases (H2)",
	         meta = "100% Ocean; No clouds",
             reference = "Domagal-Goldman et al. 2014",
             path_to_file=os.path.join(path,"GJ876_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "Sigma Bootis (F2V)",
	         planet = "High CO2; Low CH4",
             description = "50% CO2 with low volcanic outgassing of reduced gases (H2)",
	         meta = "100% Ocean; No clouds",
             reference = "Domagal-Goldman et al. 2014",
             path_to_file=os.path.join(path,"F2V_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "AD Leo (M3.5V)",
	         planet = "High CO2; Low CH4",
             description = "50% CO2 with low volcanic outgassing of reduced gases (H2)",
	         meta = "100% Ocean; No clouds",
             reference = "Domagal-Goldman et al. 2014",
             path_to_file=os.path.join(path,"ADLeo_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.flx")
             ),

    # GJ 1214b uses the GJ 876 spectrum
    Spectrum(observation="Direct",
	         star = "GJ 1214",
	         planet = "GJ 1214b",
             description = "100X Solar Metallicity",
	         meta = " ",
             reference = "Charnay et al. 2015",
             path_to_file=os.path.join(path,"gj1214b_100Xsolar_hitran2012_50_100000cm_toa.flx")
             ),

    # GJ 1214b uses the GJ 876 spectrum
    Spectrum(observation="Direct",
	         star = "GJ 1214",
	         planet = "GJ 1214b",
             description = "Solar Metallicity",
	         meta = " ",
             reference = "Charnay et al. 2015",
             path_to_file=os.path.join(path,"gj1214b_solar_hitran2012_50_100000cm_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "Sun",
	         planet = "Venus",
             description = "Venus nightside",
	         meta = "Modern Solar System",
             reference = "Arney et al. 2014",
             path_to_file=os.path.join(path,"smart_venus4_30_4000_11000cm_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Sun",
             planet = "Venus",
             description = "Venus dayside",
             meta = "Modern Solar System",
             reference = "n/a",
             path_to_file=os.path.join(path,"dayside_venus_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "Sun",
             planet = "Archean Earth",
             description = "Hazy Archean Earth orbiting the modern Sun",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_modernsun.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Epsilon Eridani (K2V)",
             planet = "Archean Earth",
             description = "Hazy Archean Earth orbiting K2V star",
             meta = "1 bar; 1% CO2; 0.3% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_K_higher_CH4.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Epsilon Eridani (K2V)",
             planet = "Archean Earth",
             description = "Haze-free Archean Earth orbiting K2V star",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_K_lower_CH4.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "GJ 876 (M4V)",
             planet = "Archean Earth",
             description = "Hazy Archean Earth orbiting GJ 876",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_gj876.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Sigma Bootis (F2V)",
             planet = "Archean Earth",
             description = "Haze-free Archean Earth orbiting F2V star",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_fstar.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Hazy Archean Earth orbiting the Archean Sun",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_archeansun.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
	         star = "AD Leo (M3.5V)",
	         planet = "Archean Earth",
             description = "Haze-free Archean Earth orbiting AD Leo",
	         meta = "1 bar; 1% CO2; 0.2% CH4",
	         reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_adleo.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "AD Leo (M3.5V)",
             planet = "Archean Earth",
             description = "Hazy Archean Earth orbiting AD Leo",
             meta = "1 bar; 1% CO2; 0.9% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path,"HAZE_adleo_haze.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thin haze; Stratocumulus cloud",
             meta = "1 bar; 2% CO2; 0.32% CH4; thin haze; stratocumulus clouds",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_stratocum_2.7Ga_3.2E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
     ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thin haze; Cirrus cloud",
             meta = "1 bar; 2% CO2; 0.32% CH4; thin haze; cirrus clouds",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_cirrus_2.7Ga_3.2E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thin haze; No cloud",
             meta = "1 bar; 2% CO2; 0.32% CH4; thin haze",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_2.7Ga_3.2E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thick haze; Stratocumulus cloud",
             meta = "1 bar; 2% CO2; 0.37% CH4; thick haze; stratocumulus clouds",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_stratocum_2.7Ga_3.7E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thick haze; Cirrus cloud",
             meta = "1 bar; 2% CO2; 0.37% CH4; thick haze; cirrus clouds",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_cirrus_2.7Ga_3.7E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
             ),


    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thick haze; No cloud",
             meta = "1 bar; 2% CO2; 0.37% CH4; thick haze",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_2.7Ga_3.7E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Proterozoic Sun",
             planet = "Proterozoic Earth",
             description = "Low oxygen",
             meta = "1 bar; 1% CO2; 0.03% CH4; 1% PAL O2",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"hi_o2_proterozoic_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Proterozoic Sun",
             planet = "Proterozoic Earth",
             meta = "1 bar; 1% CO2; 0.03% CH4; 0.1% PAL O2",
             description = "High oxygen",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path, "low_o2_proterozoic_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "No haze; No cloud",
             meta = "1 bar; 2% CO2; 0.2% CH4; haze-free",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"HAZE_2.7Ga_2.00E-03ch4_rmix_2.0E-2_1.013bar_file.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "No haze; Cirrus clouds",
             meta = "1 bar; 2% CO2; 0.2% CH4; haze-free; cirrus clouds",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"cirrus_HAZE_2.7Ga_2.00E-03ch4_rmix_2.0E-2_1.013bar_1.00E-08o2_0.63gcm3.pt_toa.flx")
             ),

    Spectrum(observation="Direct",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "No haze; Stratocumulus clouds",
             meta = "1 bar; 2% CO2; 0.2% CH4; haze-free; stratocumulus clouds",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path,"stratocum_HAZE_2.7Ga_2.00E-03ch4_rmix_2.0E-2_1.013bar_1.00E-08o2_0.63gcm3.pt_toa.flx")
             ),

    #"Direct_GJ876_highO2_1bar" : "photo_clima_p01bar_3.6e-4_GJ876_dry.pt_fixed_filtered_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="High O2",
               meta="1 bar",
               reference="Schwieterman et al. 2016",
               path_to_file=os.path.join(path,"photo_clima_p01bar_3.6e-4_GJ876_dry.pt_fixed_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_GJ876_highO2_10bar" : "photo_clima_p010bar_3.6e-4CO2_GJ876_dry.pt_fixed_filtered_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="High O2",
               meta="10 bar",
               reference="Schwieterman et al. 2016",
               path_to_file=os.path.join(path,"photo_clima_p010bar_3.6e-4CO2_GJ876_dry.pt_fixed_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_GJ876_highO2_100bar" : "photo_clima_p0100bar_GJ876_dry.pt_fixed_filtered_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="High O2",
               meta="100 bar",
               reference="Schwieterman et al. 2016",
               path_to_file=os.path.join(path,"photo_clima_p0100bar_GJ876_dry.pt_fixed_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_GJ876_O2_harman" : "GJ876_worst_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="Photolysis of CO2",
               meta="1 bar",
               reference="Harman et al. 2015; Schwieterman et al. 2016",
               path_to_file=os.path.join(path,"GJ876_worst_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_ADLeo_Segura2005" : "ADLeo_Segura2005_HITRAN2012_toa.flx"
    Spectrum(observation="Direct",
               star="AD Leo (M3.5V)",
               planet="Earth",
               description="Segura",
               meta="1 bar",
               reference="Schwieterman 2016; Segura 2005",
               path_to_file=os.path.join(path,"ADLeo_Segura2005_HITRAN2012_toa.flx")
              ),

    #"Direct_Proxima Centauri_high O2_10bar_dry" : "profile_o2lb_10bar_dry.pt_filtered_hitran2012_50_100000cm_toa.flx"
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="High O2 Dry",
               meta="10 bar",
               reference="Meadows et al. 2018; Schwieterman 2016; Luger & Barnes 2015",
               path_to_file=os.path.join(path,"10bar_O2_dry.pt_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_Proxima Centauri_high O2_10bar_wet" : "profile_o2lb_10bar_h2o.pt_filtered_hitran2012_50_100000cm_toa.flx"
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="High O2 Wet",
               meta="10 bar",
               reference="Meadows et al. 2018; Schwieterman 2016; Luger & Barnes 2015",
               path_to_file=os.path.join(path,"10bar_O2_wet.pt_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_Proxima Centauri_O2 CO2_10bar" : "profile_O2_CO2_10bar_prox_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Evolved O2, CO2",
               meta="10 bar",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path,"10bar_O2_CO2_final.pt_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_Proxima Centauri_O2 CO2_90bar" : "profile_O2_CO2_90bar_prox_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Evolved O2, CO2",
               meta="90 bar",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path,"90bar_O2_CO2_profile.pt_filtered_hitran2012_50_100000cm_toa.flx")
              ),

    #"Direct_Proxima Centauri_Venus_10bar_clouds" : "fig17_smart_spectra_pandora10bar_cloudy_500_100000cm-1_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Venus-like",
               meta="10 bar; Cloudy",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path,"PCb_Venus_10bar_toa.flx")
              ),

    #"Direct_Proxima Centauri_Venus_90bar_clouds" : "fig17_smart_spectra_pandora90bar_clouds_500_100000cm-1_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Venus-like",
               meta="90 bar; Cloudy",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "PCb_Venus_90bar_toa.flx")
              ),

    #"Direct_Proxima Centauri_Gao_1bar" : "smart_gao_1bar_update_xsec_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="FP Earth",
               meta="1 bar; CO2/O2/CO",
               reference="Meadows et al. 2018; Gao et al. 2015",
               path_to_file=os.path.join(path, "smart_gao_1bar_FINAL_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; No Clouds",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_clear_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; No Clouds",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_clear_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; Cirrus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_cirrus_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; Cirrus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_cirrus_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; Stratocumulus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_stcum_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; Stratocumulus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_stcum_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; 50% No Cloud; 25% Cirrus; 25% Stratocumulus",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_combined_new_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; 50% No Cloud; 25% Cirrus; 25% Stratocumulus",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_combined_new_toa.flx")
              ),

    #"Direct_Proxima Centauri_Earth_Cirrus" : "profile_earth_prox.pt_cirrus_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Earth-like",
               meta="Cirrus Cloud",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "profile_Earth_proxb_.pt_cirrus_hitran2012_o4_noh2co_187Kstrat_toa.flx")
              ),

    #"Direct_Proxima Centauri_Earth_Clear" : "profile_earth_prox.pt_filtered_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Earth-like",
               meta="Clear",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat_toa.flx")
              ),

    #"Direct_Proxima Centauri_Earth_Stratocumulus" : "profile_earth_prox.pt_stratocum_hitran2012_50_100000cm_toa.flx",
    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Earth-like",
               meta="Stratocumulus Cloud",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "profile_Earth_proxb_.pt_stratocum_hitran2012_o4_noh2co_187Kstrat_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="Proxima Centauri",
               planet="b",
               description="Earth-like",
               meta="50% Clear; 25% Cirrus; 25% Stratocumulus",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "profile_Earth_proxb_combined_toa.flx")
              ),

    # TRAPPIST-1 c
    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="c",
               description="CO2",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_co2_93bar_c_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="c",
               description="Venus-like",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_venus1_93bar_c_smart_spectra_toa.flx")
              ),

    # TRAPPIST-1 e
    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="CO2",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_co2_10bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="CO2",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_co2_93bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="O2 outgassing",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_10bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="O2 outgassing",
               meta="100 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_100bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="O2 desiccated",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_dry_10bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="O2 desiccated",
               meta="100 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_dry_100bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="Venus-like",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_venus1_10bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="Venus-like",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_venus1_93bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="Aqua Planet",
               meta="1 bar; Clear Sky",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_h2o_ocean_1bar_e_smart_spectra_toa.flx")
              ),

    Spectrum(observation="Direct",
               star="TRAPPIST-1",
               planet="e",
               description="Aqua Planet",
               meta="1 bar; Cloudy",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_h2o_cld_ocean_1bar_e_smart_spectra_toa.flx")
              ),

    # MYSTERY

    # T-1e Clear Sky Aqua
    Spectrum(observation="Direct",
               star="Pandora",
               planet="b",
               description="1",
               meta="b1",
               reference="",
               path_to_file=os.path.join(path, "trappist1_h2o_ocean_1bar_e_smart_spectra_toa.flx")
              ),
    # T-1e Cloudy Aqua
    Spectrum(observation="Direct",
               star="Pandora",
               planet="b",
               description="2",
               meta="b2",
               reference="",
               path_to_file=os.path.join(path, "trappist1_h2o_cld_ocean_1bar_e_smart_spectra_toa.flx")
              ),

    # Epsilon Eridani Hazy Archean
    Spectrum(observation="Direct",
             star = "Pandora",
             planet = "c",
             description="1",
             meta = "c1",
             reference = "",
             path_to_file=os.path.join(path, "HAZE_K_higher_CH4.pt_toa.flx")
             ),
    # Epsilon Eridani Haze-free Archean
    Spectrum(observation="Direct",
             star = "Pandora",
             planet = "c",
             description="2",
             meta = "c2",
             reference = "",
             path_to_file=os.path.join(path, "HAZE_K_lower_CH4.pt_toa.flx")
             ),

    # Modern Earth
    Spectrum(observation="Direct",
	         star = "Pandora",
             planet = "d",
             description="1",
             meta = "d1",
             reference = "",
             path_to_file=os.path.join(path, "earth_quadrature_radiance_refl.flx")
             ),

    # Proxima b Modern Earth-like (CH4)
    Spectrum(observation="Direct",
               star="Pandora",
               planet="e",
               description="1",
               meta="e1",
               reference="",
               path_to_file=os.path.join(path, "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat_toa.flx")
              ),

    # T-1d O2 outgassing (H2O)
    Spectrum(observation="Direct",
               star="Pandora",
               planet="f",
               description="1",
               meta="f1",
               reference="",
               path_to_file=os.path.join(path, "trappist1_o2_10bar_d_smart_spectra_toa.flx")
              ),

    # Archean Earths
    Spectrum(observation="Direct",
               star="Pandora",
               planet="g",
               description="1",
               meta="g1",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_clear_new_toa.flx")
              ),
    Spectrum(observation="Direct",
               star="Pandora",
               planet="g",
               description="2",
               meta="g2",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_cirrus_new_toa.flx")
              ),
    Spectrum(observation="Direct",
               star="Pandora",
               planet="g",
               description="3",
               meta="g3",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_clear_new_toa.flx")
              ),
    Spectrum(observation="Direct",
               star="Pandora",
               planet="g",
               description="4",
               meta="g4",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_cirrus_new_toa.flx")
              ),

    # T-1c Clear Sky CO2
    Spectrum(observation="Direct",
               star="Pandora",
               planet="h",
               description="1",
               meta="h1",
               reference="",
               path_to_file=os.path.join(path, "trappist1_co2_93bar_c_smart_spectra_toa.flx")
              ),
    # T-1c Cloudy Venus
    Spectrum(observation="Direct",
               star="Pandora",
               planet="h",
               description="2",
               meta="h2",
               reference="",
               path_to_file=os.path.join(path, "trappist1_venus1_93bar_c_smart_spectra_toa.flx")
              ),


    #############################
    # TRANSMISSION SPECTROSCOPY #
    #############################

    Spectrum(observation="Transmission",
	         star = "Sun",
             planet = "Earth",
             description="Modern Earth orbiting the Sun",
             meta = "Globally averaged atmospheric structure; 1 bar; No Clouds",
             reference = "J. Lustig-Yaeger",
             path_to_file=os.path.join(path, "earth_avg_hitran2012_300_100000cm.trn")
             ),

    Spectrum(observation="Transmission",
	         star = "Sun",
             planet = "Archean Earth",
             description="Hazy Archean Earth orbiting the modern Sun",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_modernsun.ptTRAN.trn")
             ),


    Spectrum(observation="Transmission",
             star = "Epsilon Eridani (K2V)",
             planet = "Archean Earth",
             description="Hazy Archean Earth orbiting K2V star",
             meta = "1 bar; 1% CO2; 0.3% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_k2v_haze.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "Epsilon Eridani (K2V)",
             planet = "Archean Earth",
             description="Haze-free Archean Earth orbiting K2V star",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_kstar.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "GJ 876 (M4V)",
             planet = "Archean Earth",
             description="Hazy Archean Earth orbiting GJ 876",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_gj876.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "Sigma Bootis (F2V)",
             planet = "Archean Earth",
             description="Haze-free Archean Earth orbiting F2V star",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_fstar.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "Archean Sun",
             planet = "Archean Earth",
             description="Hazy Archean Earth orbiting the Archean Sun",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_archeansun.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "AD Leo (M3.5V)",
             planet = "Archean Earth",
             description="Haze-free Archean Earth orbiting AD Leo",
             meta = "1 bar; 1% CO2; 0.2% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_adleo.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "AD Leo (M3.5V)",
             planet = "Archean Earth",
             description="Hazy Archean Earth orbiting AD Leo",
             meta = "1 bar; 1% CO2; 0.9% CH4",
             reference = "Arney et al. 2017",
             path_to_file=os.path.join(path, "HAZE_adleo_haze.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "Archean Sun",
             planet = "Archean Earth",
             description="Thin haze",
             meta = "1 bar; 2% CO2; 0.32% CH4; thin haze",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path, "HAZE_2.7Ga_3.2E-03ch4_rmix_2.0E-2_1.013bar_file.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "Archean Sun",
             planet = "Archean Earth",
             description = "Thick haze",
             meta = "1 bar; 2% CO2; 0.37% CH4; thick haze",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path, "HAZE_2.7Ga_3.7E-03ch4_rmix_2.0E-2_1.013bar_file.ptTRAN.trn")
             ),

    Spectrum(observation="Transmission",
             star = "Archean Sun",
             planet = "Archean Earth",
             description="No haze",
             meta = "1 bar; 2% CO2; 0.2% CH4; haze-free",
             reference = "Arney et al. 2016",
             path_to_file=os.path.join(path, "HAZE_2.7Ga_2.00E-03ch4_rmix_2.0E-2_1.013bar_file.ptTRAN.trn")
             ),

    #"Transmission_GJ876_highO2_1bar" : "photo_clima_p01bar_3.6e-4_GJ876_dry.pt_fixed_filtered.trn",
    Spectrum(observation="Transmission",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="High O2",
               meta="1 bar",
               reference="Schwieterman et al. 2016",
               path_to_file=os.path.join(path, "photo_clima_p01bar_3.6e-4_GJ876_dry.pt_fixed_filtered.trn")
              ),

    #"Transmission_GJ876_highO2_10bar" : "photo_clima_p010bar_3.6e-4CO2_GJ876_dry.pt_fixed_filtered.trn",
    Spectrum(observation="Transmission",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="High O2",
               meta="10 bar",
               reference="Schwieterman et al. 2016",
               path_to_file=os.path.join(path, "photo_clima_p010bar_3.6e-4CO2_GJ876_dry.pt_fixed_filtered.trn")
              ),

    #"Transmission_GJ876_highO2_10bar" : "photo_clima_p0100bar_GJ876_dry.pt_fixed_filtered.trn",
    Spectrum(observation="Transmission",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="High O2",
               meta="100 bar",
               reference="Schwieterman et al. 2016",
               path_to_file=os.path.join(path, "photo_clima_p0100bar_GJ876_dry.pt_fixed_filtered.trn")
              ),

    #"Transmission_GJ876_O2_harman" : "GJ_876_highCO_hitran2012_50_100000cm.trn",
    Spectrum(observation="Transmission",
               star="GJ 876 (M4V)",
               planet="FP Earth",
               description="Photolysis of CO2",
               meta=" ",
               reference="Harman et al. 2015; Schwieterman et al. 2016",
               path_to_file=os.path.join(path, "GJ_876_highCO_hitran2012_50_100000cm.trn")
              ),

    #"Transmission_ADLeo_Segura2005" : "smartin_adleo_hyak_hitran2012_50_100000cm.trn",
    Spectrum(observation="Transmission",
               star="AD Leo (M3.5V)",
               planet="Earth",
               description="Segura",
               meta=" ",
               reference="Schwieterman 2016; Segura 2005",
               path_to_file=os.path.join(path, "smartin_adleo_hyak_hitran2012_50_100000cm.trn")
              ),

    #"Transmission_Proxima Centauri_O2_10bar_dry" : "profile_o2lb_10bar_dry.pt_filtered_transit.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="High O2 Dry",
               meta="10 bar",
               reference="Meadows et al. 2018; Schwieterman 2016; Luger & Barnes 2015",
               path_to_file=os.path.join(path, "10bar_O2_dry.pt_filtered_hitran2012_50_100000cm.trn")
              ),

    #"Transmission_Proxima Centauri_O2_10bar_wet" : "profile_o2lb_10bar_h2o.pt_filtered_transit.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="High O2 Wet",
               meta="10 bar",
               reference="Meadows et al. 2018; Schwieterman 2016; Luger & Barnes 2015",
               path_to_file=os.path.join(path, "10bar_O2_wet.pt_filtered_hitran2012_50_100000cm.trn")
              ),

    #"Transmission_Proxima Centauri_O2 CO2_10bar" : "profile_O2_CO2_10bar_prox_transit_hitran2012_50_100000cm.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Evolved O2, CO2",
               meta="10 bar",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "10bar_O2_CO2_final.pt_filtered_hitran2012_50_100000cm.trn")
              ),

    #"Transmission_Proxima Centauri_O2 CO2_90bar" : "profile_O2_CO2_90bar_prox_transit_hitran2012_50_100000cm.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Evolved O2, CO2",
               meta="90 bar",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "90bar_O2_CO2_profile.pt_filtered_hitran2012_50_100000cm.trn")
              ),

    #"Transmission_Proxima Centauri_Venus_10bar_clouds" : "fig24_tran_smart_spectra_pandora10bar_cloudy_500_100000cm-1.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Venus-like",
               meta="10 bar; Cloudy",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "PCb_Venus_10bar.trn")
              ),

    #"Transmission_Proxima Centauri_Venus_90bar_clouds" : "fig24_tran_smart_spectra_pandora90bar_clouds_500_100000cm-1.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Venus-like",
               meta="90 bar; Cloudy",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "PCb_Venus_90bar.trn")
              ),

    #"Transmission_Proxima Centauri_Gao_1bar" : "Gao2015_case3.pt_filtered_transit.trn",
    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="FP Earth",
               meta="1 bar; CO2/O2/CO",
               reference="Meadows et al. 2018; Gao et al. 2015",
               path_to_file=os.path.join(path, "smart_gao_1bar_FINAL.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; No Clouds",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_clear_new.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; No Clouds",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_clear_new.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; Cirrus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_cirrus_new.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; Cirrus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_cirrus_new.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="No Haze; Stratocumulus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_stcum_new.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Archean",
               meta="Haze; Stratocumulus Cloud",
               reference="Meadows et al. 2018; Arney et al. 2016",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_stcum_new.trn")
              ),

    Spectrum(observation="Transmission",
               star="Proxima Centauri",
               planet="b",
               description="Earth-like",
               meta="Clear",
               reference="Meadows et al. 2018",
               path_to_file=os.path.join(path, "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat.trn")
              ),

    # TRAPPIST-1 c
    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="c",
               description="CO2",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_co2_93bar_c_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="c",
               description="Venus-like",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_venus1_93bar_c_smart_spectra.trn")
              ),

    # TRAPPIST-1 e
    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="CO2",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_co2_10bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="CO2",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_co2_93bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="O2 outgassing",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_10bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="O2 outgassing",
               meta="100 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_100bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="O2 desiccated",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_dry_10bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="O2 desiccated",
               meta="100 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_o2_dry_100bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="Venus-like",
               meta="10 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_venus1_10bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="Venus-like",
               meta="92 bar",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_venus1_93bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="Aqua Planet",
               meta="1 bar; Clear Sky",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_h2o_ocean_1bar_e_smart_spectra.trn")
              ),

    Spectrum(observation="Transmission",
               star="TRAPPIST-1",
               planet="e",
               description="Aqua Planet",
               meta="1 bar; Cloudy",
               reference="Lincowski et al. 2018",
               path_to_file=os.path.join(path, "trappist1_h2o_cld_ocean_1bar_e_smart_spectra.trn")
              ),

    # MYSTERY

    # T-1e Clear Sky Aqua
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="b",
               description="1",
               meta="b1",
               reference="",
               path_to_file=os.path.join(path, "trappist1_h2o_ocean_1bar_e_smart_spectra.trn")
              ),
    # T-1e Cloudy Aqua
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="b",
               description="2",
               meta="b2",
               reference="",
               path_to_file=os.path.join(path, "trappist1_h2o_cld_ocean_1bar_e_smart_spectra.trn")
              ),

    # Epsilon Eridani Hazy Archean
    Spectrum(observation="Transmission",
             star = "Pandora",
             planet = "c",
             description="1",
             meta = "c1",
             reference = "",
             path_to_file=os.path.join(path, "HAZE_k2v_haze.ptTRAN.trn")
             ),
    # Epsilon Eridani Haze-free Archean
    Spectrum(observation="Transmission",
             star = "Pandora",
             planet = "c",
             description="2",
             meta = "c2",
             reference = "",
             path_to_file=os.path.join(path, "HAZE_kstar.ptTRAN.trn")
             ),

    # Modern Earth
    Spectrum(observation="Transmission",
	         star = "Pandora",
             planet = "d",
             description="1",
             meta = "d1",
             reference = "",
             path_to_file=os.path.join(path, "earth_avg_hitran2012_300_100000cm.trn")
             ),

    # Proxima b Modern Earth-like (CH4)
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="e",
               description="1",
               meta="e1",
               reference="",
               path_to_file=os.path.join(path, "profile_Earth_proxb_.pt_hitran2012_o4_noh2co_187Kstrat.trn")
              ),

    # T-1d O2 outgassing (H2O)
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="f",
               description="1",
               meta="f1",
               reference="",
               path_to_file=os.path.join(path, "trappist1_o2_10bar_d_smart_spectra.trn")
              ),

    # Archean Earths
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="g",
               description="1",
               meta="g1",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_clear_new.trn")
              ),
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="g",
               description="2",
               meta="g2",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.00e-02ch4_cirrus_new.trn")
              ),
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="g",
               description="3",
               meta="g3",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_clear_new.trn")
              ),
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="g",
               description="4",
               meta="g4",
               reference="",
               path_to_file=os.path.join(path, "HAZE_1.50e-02ch4_cirrus_new.trn")
              ),

    # T-1c Clear Sky CO2
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="h",
               description="1",
               meta="h1",
               reference="",
               path_to_file=os.path.join(path, "trappist1_co2_93bar_c_smart_spectra.trn")
              ),
    # T-1c Cloudy Venus
    Spectrum(observation="Transmission",
               star="Pandora",
               planet="h",
               description="2",
               meta="h2",
               reference="",
               path_to_file=os.path.join(path, "trappist1_venus1_93bar_c_smart_spectra.trn")
              ),

]

def make_new_small_database():
    """
    Make a new database at a degraded resolution so the file is relatively small
    """

    # Get current date to avoid overwritting an old database
    from datetime import datetime
    now = datetime.now()
    datetag = now.strftime("%m%d%y")
    new_name = "csv/spectra_small_%s.csv" %datetag
    new_meta = "csv/spectra_meta_%s.csv" %datetag

    """
    Write spectra to csv
    """
    #write_spectra_csv(spectra, savename="csv/spectra_small.csv", lammin=0.3, lammax=1.0)

    write_spectra_csv(spectra,
                      savename=new_name,
                      lammin=0.2, lammax=25.0,
                      dynamic_nth=True, fn=5000)
    #write_spectra_csv(spectra, savename="csv/spectra_large_061417.csv", lammin=0.2, lammax=25.0, nth=1)

    """
    Write metadata to csv
    """
    write_spectra_meta_csv(spectra, savename=new_meta)

def make_new_large_database():
    """
    Make a new database at the native spectral resolution
    """

    # Get current date to avoid overwritting an old database
    from datetime import datetime
    now = datetime.now()
    datetag = now.strftime("%m%d%y")
    new_name = "csv/spectra_large_%s.csv" %datetag
    new_meta = "csv/spectra_meta_%s.csv" %datetag

    """
    Write spectra to csv
    """
    #write_spectra_csv(spectra, savename="csv/spectra_small.csv", lammin=0.3, lammax=1.0)

    write_spectra_csv(spectra,
                      savename=new_name,
                      lammin=0.2, lammax=25.0,
                      dynamic_nth=False)
    #write_spectra_csv(spectra, savename="csv/spectra_large_061417.csv", lammin=0.2, lammax=25.0, nth=1)

    """
    Write metadata to csv
    """
    write_spectra_meta_csv(spectra, savename=new_meta)


if __name__ == "__main__":

    make_new_small_database()
    #make_new_large_database()
