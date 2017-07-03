"""
Author: Jake Lustig-Yaeger
"""

import os
import sys
import numpy as np
import readsmart as rs
import pandas as pd
from numba import jit

class Spectrum(object):
    """
    Object to hold spectrum and associated [meta]data
    """
    def __init__(self, observation=None, star=None, planet=None, description=None, meta=None, reference=None,
                 path_to_file=None, data=None):
        """
        Construct spectrum object
        """
        self.observation=observation
        self.star=star
        self.planet=planet
        self.description=description
        self.meta=meta
        self.reference=reference
        self._path_to_file=path_to_file
        self._data=None
        if (self._data is None) and (self._path_to_file is not None):
            self._data = SpectralData.from_file(self._path_to_file)
        """
        self.tag=observation.replace(" ", "")+"_"+star.replace(" ", "")+"_"+\
                 planet.replace(" ", "")+"_"+description.replace(" ", "")+"_"+\
                 meta.replace(" ", "")
        """
        self.tag=observation+"_"+star+"_"+planet+"_"+description+"_"+meta+"_"+reference

        # Echo object
        print self

    def __str__(self):
        string = self.observation+", "+self.star+", "+self.planet+", "+self.description+", "+self.meta
        return string

    @property
    def path_to_file(self):
        return self._path_to_file

    @path_to_file.setter
    def path_to_file(self, value):
        self._path_to_file = value
        if value is not None:
            self._data = SpectralData.from_file(value)
        else:
            self._data = None

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value

class SpectralData(object):
    def __init__(self, wavelength=None, wavenumber=None, toa_flux=None,
                 star_flux=None, geo_albedo=None, flux_transmission=None,
                 absorbing_radius=None, tdepth=None):
        self.wavelength=wavelength
        self.wavenumber=wavenumber
        self.toa_flux=toa_flux
        self.star_flux=star_flux
        self.geo_albedo=geo_albedo
        self.flux_transmission=flux_transmission
        self.absorbing_radius=absorbing_radius
        self.tdepth=tdepth

    @classmethod
    def from_file(cls, path):

        # Read path file
        if path.endswith(".rad"):
            # Read rad
            wavelength, wavenumber, star_flux, toa_flux, rad_streams = rs.rad(path)
            geo_albedo = toa_flux / star_flux
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            flux_transmission, absorbing_radius, tdepth = tmp, tmp, tmp
        elif path.endswith(".tran"):
            # Read tran
            wavelength, flux_transmission, absorbing_radius = rs.tran(path)
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            wavenumber, toa_flux, star_flux, geo_albedo, tdepth = tmp, tmp, tmp, tmp, tmp
        elif path.endswith(".trnst"):
            # Read trnst
            wavelength, wavenumber, absorbing_radius, tdepth = rs.trnst(path)
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            toa_flux, star_flux, geo_albedo, flux_transmission = tmp, tmp, tmp, tmp
        elif path.endswith("data.txt"):
            # Read dat file from Will/Shawn's version
            dat = np.genfromtxt(path, skip_header=8)
            # Parse
            wavelength = dat[:,0]
            wavenumber = dat[:,1]
            star_flux  = dat[:,2]
            toa_flux   = dat[:,3]
            geo_albedo = dat[:,4]
            # Fill transit arrays with empty strings
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            flux_transmission, absorbing_radius, tdepth = tmp, tmp, tmp
        elif path.endswith(".flx"):
            # Read dat file from Will/Shawn's version
            dat = np.genfromtxt(path, skip_header=8)
            # Parse
            wavelength = dat[:,0]
            wavenumber = dat[:,1]
            star_flux  = dat[:,2]
            toa_flux   = dat[:,3]
            geo_albedo = dat[:,4]
            # Fill transit arrays with empty strings
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            flux_transmission, absorbing_radius, tdepth = tmp, tmp, tmp
        elif path.endswith(".trn"):
            # Read dat file from Will/Shawn's version
            dat = np.genfromtxt(path, skip_header=8)
            # Parse
            wavelength = dat[:,0]
            absorbing_radius = dat[:,1]
            tdepth  = dat[:,2]
            # Fill transit arrays with empty strings
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            toa_flux, star_flux, geo_albedo, flux_transmission, wavenumber = tmp, tmp, tmp, tmp, tmp
        else:
            print "%s is an invalid path." %path

        # Return new class instance
        return cls(wavelength=wavelength, wavenumber=wavenumber, toa_flux=toa_flux,
                   star_flux=star_flux, geo_albedo=geo_albedo,
                   flux_transmission=flux_transmission,
                   absorbing_radius=absorbing_radius, tdepth=tdepth)

def write_spectra_meta_csv(spectra, savename="meta.csv"):
    """
    Writes the metadata for each spectrum to a csv

    Parameters
    ----------
    spectra : list
        List of `Spectrum` objects
    savename : str
        Name of the csv file to be saved
    """

    print "Writing %s..." %savename

    # Specify columns for the meta data file
    cols = [
        "Which",
        "Observation",
        "Star",
        "Planet",
        "Description",
        "Meta",
        "Reference",
        "Datafile"
    ]

    # Create an empty list for the columns
    meta = []

    # Loop over list of Spectrum objects
    i = 1
    for s in spectra:
        # Append metadata to list
        meta.append([
            i,
            s.observation,
            s.star,
            s.planet,
            s.description,
            s.meta,
            s.reference,
            s.path_to_file
        ])
        i += 1

    # Convert to np array...
    metadata = np.array(meta)

    # ...and now to pandas DataFrame
    df = pd.DataFrame(metadata, columns=cols)

    # Save
    df.to_csv(savename)

    print "Saved %s" %savename

    return



def set_nth(N, final_length=10000.):
    """
    """
    n = int(round(N / final_length))
    if n < 1:
        nth = 1
    else:
        nth = n
    return nth

def write_spectra_csv(spectra, savename="test.csv", lammin=0.1, lammax=20.0,
                      degrade=False, Res=1000, nth=1, dynamic_nth=False, fn=10000):
    """
    Write list of Spectrum objects to csv file

    Parameters
    ----------
    spectra : list
        List of `Spectrum` objects
    savename : str
        Name of the csv file to be saved
    nth : int, optional
        Use every `nth` spectral element
    """

    print "Writing %s..." %savename

    if degrade:
        import coronagraph as cg
        lam, dlam = cg.noise_routines.construct_lam(lammin, lammax, Res=Res)

    wl = []
    wn = []
    toaf = []
    starf = []
    galb = []
    ftrn = []
    absrad = []
    tdepth = []
    tag = []
    which = []

    i = 1
    for s in spectra:
        # Select wavelength range
        mask = (s.data.wavelength >= lammin) & (s.data.wavelength <= lammax)

        # Number of wavelength points
        Nlam = np.sum(mask)

        # Create meaningless
        tmp = np.chararray(Nlam, itemsize=len(s.tag))
        tmp2 = np.zeros(Nlam, dtype=int)
        tmp[:] = s.tag
        tmp2[:] = int(i)

        if dynamic_nth:
            nth = set_nth(Nlam, final_length=fn)

        tag = np.hstack([tag,tmp[0::nth]])
        wl = np.hstack([wl, s.data.wavelength[mask][0::nth]])
        wn = np.hstack([wn, s.data.wavenumber[mask][0::nth]])
        toaf = np.hstack([toaf, s.data.toa_flux[mask][0::nth]])
        starf = np.hstack([starf, s.data.star_flux[mask][0::nth]])
        galb = np.hstack([galb, s.data.geo_albedo[mask][0::nth]])
        ftrn = np.hstack([ftrn, s.data.flux_transmission[mask][0::nth]])
        absrad = np.hstack([absrad, s.data.absorbing_radius[mask][0::nth]])
        tdepth = np.hstack([tdepth, s.data.tdepth[mask][0::nth]])
        which = np.hstack([which, tmp2[0::nth]])
        i += 1

    which = np.array(which, dtype=int)

    data = np.array([tag, wl, wn, toaf, starf, galb, ftrn, absrad, tdepth, which]).T

    cols = [
        "File Name",
        "Wavelength",
        "Wavenumber",
        "Top of Atmosphere Flux",
        "Stellar Flux",
        "Geometric Albedo",
        "Flux Transmission",
        "Absorbing Radius",
        "Transit Depth",
        "Which"
    ]

    df = pd.DataFrame(data, columns=cols)

    df.to_csv(savename)

    print "Saved %s" %savename

    return

def test_init():

    # Relative path to spectra files
    path = "spectrum_files/"

    # Construct Spectrum object

    test1 = Spectrum(observation="Transmission",
                    star="Proxima Centauri",
                    planet="b",
                    description="FP Earth",
                    meta="1 bar; CO2/O2/CO",
                    reference="Meadows et al., 2017; Gao et al., 2015",
                    path_to_file=os.path.join(path, "smart_gao_1bar_update_xsec.trn")
              )

    test2 = Spectrum(observation="Direct",
                    star="AD Leo",
                    planet="Earth",
                    description="No clouds; 100% Ocean",
                    meta="1 PAL O2",
                    reference="Segura et al. 2005",
                    path_to_file=os.path.join(path, "smart_Mstar_ADLeo_1PAL_clr_ocean_50_30000cm_r60only_toa.flx")
              )

    ltest = [test1, test2]

    return ltest

def test_write_meta():
    """
    """
    spectra = test_init()

    write_spectra_meta_csv(spectra, savename="csv/test_meta.csv")

    return

def test_write_spectra():
    """
    """
    spectra = test_init()

    write_spectra_csv(spectra, savename="csv/test_spectra.csv")

    return

def copy_raw_spectra(dir="spectrum_files/",
                     offdir="jlustigy@acidrain.astro.washington.edu:/astro/store/vpl1/jlustigy/VPLSpectralDatabase/spectrum_files/",
                     exclude=[".DS_Store", ".git", ".gitignore"]):
    """
    Uses the rsync ssh copying protocol to transfer all *new* spectra to a remote directory.

    Parameters
    ----------
    """
    # define the rsync command and options
    command = "rsync -Cavz "

    # loop over exclude elements constructing the correct string
    exc = ''.join("--exclude %s " %x for x in exclude)

    # compose final command
    final_command = command+exc+dir+" "+offdir

    # use os to call final_command on the command line
    os.system(final_command)

    return

def convert_rad_to_txt(fn):
    """
    Convert .rad files to .txt files with headers

    Parameters
    ----------
    fn : str
        File name/location
    """
    # Construct new file header
    header =  """Columns:
1. Wavelength (um)
2. Wavenumber (cm-1)
3. Solar Flux incident at top of atmosphere (w/m2/um)
4. Planetary Flux (w/m2/um)
5. Albedo

Col 1        Col 2        Col 3        Col 4        Col 5"""

    # Read in rad file
    dat = np.genfromtxt(fn)
    wl = dat[:,0]
    wno = dat[:,1]
    solflx = dat[:,2]
    flx = dat[:,3]
    rads = dat[:,4:]

    # Compute reflectivity
    refl = flx/solflx

    # Construct data array
    data = np.array([wl, wno, solflx, flx, refl]).T

    # Create new file name
    nfn = fn[:-4]+"_data.txt"

    # Write new file
    np.savetxt(nfn, data, fmt='%.6e', header=header, comments="")

    return

def convert_rad_to_flx(fn):
    """
    Convert .rad files to .flx files with headers

    Parameters
    ----------
    fn : str
        File name/location
    """
    # Construct new file header
    header =  """Columns:
1. Wavelength (um)
2. Wavenumber (cm-1)
3. Solar Flux incident at top of atmosphere (w/m2/um)
4. Planetary Flux (w/m2/um)
5. Albedo

Col 1        Col 2        Col 3        Col 4        Col 5"""

    # Read in rad file
    dat = np.genfromtxt(fn)
    wl = dat[:,0]
    wno = dat[:,1]
    solflx = dat[:,2]
    flx = dat[:,3]
    rads = dat[:,4:]

    # Compute reflectivity
    refl = flx/solflx

    # Construct data array
    data = np.array([wl, wno, solflx, flx, refl]).T

    # Create new file name
    nfn = fn[:-4]+".flx"

    # Write new file
    np.savetxt(nfn, data, fmt='%.6e', header=header, comments="")

    return

def convert_txt_to_flx(fn):
    """
    Creates a copy of *_data.txt files renamed to *.flx
    """

    # Create new file name
    nfn = fn[:-9] + ".flx"

    # compose final command
    final_command = "cp %s %s" %(fn, nfn)

    # use os to call final_command on the command line
    os.system(final_command)

    return

def run_flx_conversion(spec_dir="spectrum_files"):
    """
    This (hopefully) single-use function converts all existing *_data.txt files
    to *.flx
    """

    spec_list = os.listdir(spec_dir)
    for spec in spec_list:
        fn = os.path.join(spec_dir, spec)
        if fn.endswith("_data.txt"):
            convert_txt_to_flx(fn)
        elif fn.endswith("_toa.rad"):
            convert_rad_to_flx(fn)
        else:
            pass

    return

def read_transit(path, radius=None):
    """
    Read *.tran or *.trnst output files from SMART.

    Parameters
    ----------
    path : str
        Path/Filename of *.tran file
    radius : tuple, optional
        (planet_radius, stellar_radius) for transit depth calculation [km]

    Returns
    -------
    wl : ndarray
        Wavelength [microns]
    absorbing_radius : ndarray
        Effective absorbing radius of the atmosphere (in km)
        (effective radius - solid body radius)
    tdepth : ndarray
        Transit depth

    Note
    ----
    dF = (Rp/Rs)**2 = -(Rp_solid + column3)**2/Rs**2
    """
    if path.endswith(".tran"):
        data = np.genfromtxt(path, skip_header=0)
        wl = data[:,0]
        flux_ratio = data[:,1]
        absorbing_radius = data[:,2]
        if radius is not None:
            tdepth = (radius[0] + absorbing_radius)**2/radius[1]**2
        else:
            tdepth = None
    elif path.endswith(".trnst"):
        # Read in .trnst file
        data = np.genfromtxt(path, skip_header=2)
        # Split into arrays
        wl = data[:,0]
        wn = data[:,1]
        absorbing_radius = data[:,2]
        tdepth = data[:,3]
    else:
        print "Invalid file..."
        return

    return wl, absorbing_radius, tdepth

def convert_to_trn(fn, radius=(6371., 695700.)):
    """
    Converts *.tran or *.trnst output file from SMART to *.trn file

    Parameters
    ----------
    path : str
        Path/Filename of *.tran file
    radius : tuple, optional
        (planet_radius, stellar_radius) for transit depth calculation [km]
    """

    if fn.endswith(".tran"):

        # Create new file name/loc
        nfn = fn[:-5] + ".trn"

        # Read-in data and parse
        data = np.genfromtxt(fn, skip_header=0)
        wl = data[:,0]
        flux_ratio = data[:,1]
        absorbing_radius = data[:,2]

        # Calculate transit depth using given radii
        if radius is not None:
            tdepth = (radius[0] + absorbing_radius)**2/radius[1]**2
        else:
            tdepth = None

    elif fn.endswith(".trnst"):

        # Create new file name/loc
        nfn = fn[:-6] + ".trn"

        # Read in .trnst file
        data = np.genfromtxt(fn, skip_header=2)
        # Split into arrays
        wl = data[:,0]
        wn = data[:,1]
        absorbing_radius = data[:,2]
        tdepth = data[:,3]
    else:
        print "Invalid file..."
        return

    # Construct data array
    data = np.array([wl, absorbing_radius, tdepth]).T

    # Construct new file header
    header =  """Columns:
1. Wavelength (um)
2. Effective absorbing radius of planetary atmosphere (km)
3. Transit Depth (Rp/Rs)^2

Col 1        Col 2        Col 3"""

    # Write new file
    np.savetxt(nfn, data, fmt='%.6e', header=header, comments="")

    return

def run_trn_conversion(spec_dir="spectrum_files"):
    """
    This (hopefully) single-use function converts all existing *.tran and *.trnst files
    to *.trn
    """

    spec_list = os.listdir(spec_dir)
    tran_count = 0
    for spec in spec_list:
        fn = os.path.join(spec_dir, spec)
        if fn.endswith(".trnst"):
            # .trnst files already have correct columns so radius kwarg not needed
            convert_to_trn(fn)
        elif fn.endswith(".tran"):
            tran_count += 1
        else:
            pass

    tran_list = [
        ("fig24_tran_smart_spectra_pandora10bar_cloudy_500_100000cm-1.tran", (6850., 98093.7)),
        ("fig24_tran_smart_spectra_pandora90bar_clouds_500_100000cm-1.tran", (6850., 98093.7)),
        ("fig28_HAZE_msun21_0.0Ga_1.00e-02ch4_rmix_5.0E-2__30.66fscaleTRAN.tran", (6850., 98093.7)),
        ("fig28_HAZE_msun21_0.0Ga_3.00e-02ch4_rmix_5.0E-2__30.66fscaleTRAN.tran", (6850., 98093.7)),
        ("Gao2015_case3.pt_filtered_transit.tran", (6850., 98093.7)),
        ("GJ_876_highCO_hitran2012_50_100000cm.tran", (6371., 261580.)),
        ("HAZE_2.7Ga_2.00E-03ch4_rmix_2.0E-2_1.013bar_file.ptTRAN.tran", (6371., 695508.)),
        ("HAZE_2.7Ga_3.2E-03ch4_rmix_2.0E-2_1.013bar_file.ptTRAN.tran", (6371., 695508.)),
        ("HAZE_2.7Ga_3.7E-03ch4_rmix_2.0E-2_1.013bar_file.ptTRAN.tran", (6371., 695508.)),
        ("HAZE_adleo_haze.ptTRAN.tran", (6371., 271248.)),
        ("HAZE_adleo.ptTRAN.tran", (6371., 271248.)),
        ("HAZE_archeansun.ptTRAN.tran", (6371., 695508.)),
        ("HAZE_fstar.ptTRAN.tran", (6371., 995271.)),
        ("HAZE_gj876.ptTRAN.tran", (6371., 261580.)),
        ("HAZE_k2v_haze.ptTRAN.tran", (6371., 511198.)),
        ("HAZE_kstar.ptTRAN.tran", (6371., 511198.)),
        ("HAZE_modernsun.ptTRAN.tran", (6371., 695508.)),
        ("photo_clima_p01bar_3.6e-4_GJ876_dry.pt_fixed_filtered.tran", (6371., 261580.)),
        ("photo_clima_p010bar_3.6e-4CO2_GJ876_dry.pt_fixed_filtered.tran", (6371., 261580.)),
        ("photo_clima_p0100bar_GJ876_dry.pt_fixed_filtered.tran", (6371., 261580.)),
        ("profile_earth_prox.pt_filtered_transit_hitran2012_50_100000cm.tran", (6850., 98093.7)),
        ("profile_O2_CO2_10bar_prox_transit_hitran2012_50_100000cm.tran", (6850., 98093.7)),
        ("profile_O2_CO2_90bar_prox_transit_hitran2012_50_100000cm.tran", (6850., 98093.7)),
        ("profile_o2lb_10bar_dry.pt_filtered_transit.tran", (6850., 98093.7)),
        ("profile_o2lb_10bar_h2o.pt_filtered_transit.tran", (6850., 98093.7)),
        ("smartin_adleo_hyak_hitran2012_50_100000cm.tran", (6371., 271248.))
        ]

    # Convert tran files to trn
    for tup in tran_list:
        fn = os.path.join(spec_dir, tup[0])
        convert_to_trn(fn, radius=tup[1])

    return

################
# MOLECULES
################

class Molecule(object):
    """
    Object to hold molecular absorption information
    """
    def __init__(self, name=None, texname=None, bandcenters=None):
        """
        Parameters
        ----------
        name : str
        texname : str
        bandcenters : array-like
        """
        self.name=name
        self.texname=texname
        self.bandcenters=bandcenters

def write_molecular_csv(molecules, savename="test2.csv"):
    """
    Write list of Molecule objects to csv file

    Parameters
    ----------
    molecules : list of molecules
        List of Molecule objects
    savename : str
        Name of the csv file to be saved
    """
    name = []
    tex = []
    bandcenters = []

    for m in molecules:

        # get name
        tmp = np.chararray(len(m.bandcenters), itemsize=len(m.name))
        tmp[:] = m.name
        name = np.hstack([name,tmp])

        # get tex name
        tmp = np.chararray(len(m.bandcenters), itemsize=len(m.texname))
        tmp[:] = m.texname
        tex = np.hstack([tex,tmp])

        # get bandcenters
        bandcenters = np.hstack([bandcenters, m.bandcenters])

    data = np.array([name, tex, bandcenters]).T

    cols = [
        "Molecule Name",
        "LaTeX Name",
        "Bandcenter [um]",
    ]

    df = pd.DataFrame(data, columns=cols)

    df.to_csv(savename)

    print "Saved %s" %savename

    return

def write_molecular_csv_hr(molecules, lammin=0.1, lammax=20.0, dwno=1.0, savename="test2.csv"):
    """
    Write list of Molecule objects to csv file

    Parameters
    ----------
    molecules : list of molecules
        List of Molecule objects
    savename : str
        Name of the csv file to be saved
    """

    # Generate wavelength grid for molecules
    wno = np.arange(1e4/lammax, 1e4/lammin+dwno, dwno)
    lam = 1e4/wno[::-1]

    #
    vals = np.chararray(len(lam), itemsize=10)
    vals[:] = ""

    for m in molecules:

        for i in range(len(m.bandcenters)):
            # get wavelength grid index closest to bandcenter value
            index = np.argmin(np.fabs(lam - m.bandcenters[i]))
            # Set molecule name
            vals[index] = m.name

    data = np.array([lam, vals]).T

    cols = [
        "Wavelength",
        "Molecule",
    ]

    df = pd.DataFrame(data, columns=cols)

    df.to_csv(savename)

    print "Saved %s" %savename

    return

def test_molecules():

    # Create list of molecules and their features
    molecules = [
        Molecule("O2", "O_2", [0.2, 0.69, 0.76, 1.27]),
        Molecule("CO", "CO", [1.6, 2.35])
    ]

    # Create test csv file
    write_molecular_csv(molecules, savename="csv/test_mols.csv")

    # Create test csv file
    write_molecular_csv_hr(molecules, lammin=0.2, lammax=20.0, dwno=1.0, savename="csv/test_mols_hr.csv")

if __name__ == "__main__":
    test_init()
    test_molecules()
    test_write_spectra()
    test_write_meta()
