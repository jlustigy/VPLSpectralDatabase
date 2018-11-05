spectra = [

    # DIRECT

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

    # TRANSMISSION

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
