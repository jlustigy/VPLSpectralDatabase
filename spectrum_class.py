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
                 absorbing_radius=None):
        self.wavelength=wavelength
        self.wavenumber=wavenumber
        self.toa_flux=toa_flux
        self.star_flux=star_flux
        self.geo_albedo=geo_albedo
        self.flux_transmission=flux_transmission
        self.absorbing_radius=absorbing_radius

    @classmethod
    def from_file(cls, path):

        # Read path file
        if path.endswith(".rad"):
            # Read rad
            wavelength, wavenumber, star_flux, toa_flux, rad_streams = rs.rad(path)
            geo_albedo = toa_flux / star_flux
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            flux_transmission, absorbing_radius = tmp, tmp
        elif path.endswith(".tran"):
            # Read tran
            wavelength, flux_transmission, absorbing_radius = rs.tran(path)
            tmp = np.chararray(len(wavelength), itemsize=1)
            tmp[:] = ""
            wavenumber, toa_flux, star_flux, geo_albedo = tmp, tmp, tmp, tmp
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
            flux_transmission, absorbing_radius = tmp, tmp
        else:
            print "%s is an invalid path." %path

        # Return new class instance
        return cls(wavelength=wavelength, wavenumber=wavenumber, toa_flux=toa_flux,
                   star_flux=star_flux, geo_albedo=geo_albedo,
                   flux_transmission=flux_transmission,
                   absorbing_radius=absorbing_radius)

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
    for s in spectra:
        # Append metadata to list
        meta.append([
            s.observation,
            s.star,
            s.planet,
            s.description,
            s.meta,
            s.reference,
            s.path_to_file
        ])

    # Convert to np array...
    metadata = np.array(meta)

    # ...and now to pandas DataFrame
    df = pd.DataFrame(metadata, columns=cols)

    # Save
    df.to_csv(savename)

    print "Saved %s" %savename

    return



def write_spectra_csv(spectra, savename="test.csv", lammin=0.1, lammax=20.0,
                      degrade=False, Res=1000):
    """
    Write list of Spectrum objects to csv file

    Parameters
    ----------
    spectra : list
        List of `Spectrum` objects
    savename : str
        Name of the csv file to be saved
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
    tag = []

    for s in spectra:
        mask = (s.data.wavelength >= lammin) & (s.data.wavelength <= lammax)
        tmp = np.chararray(np.sum(mask), itemsize=len(s.tag))
        tmp[:] = s.tag
        tag = np.hstack([tag,tmp])
        wl = np.hstack([wl, s.data.wavelength[mask]])
        wn = np.hstack([wn, s.data.wavenumber[mask]])
        toaf = np.hstack([toaf, s.data.toa_flux[mask]])
        starf = np.hstack([starf, s.data.star_flux[mask]])
        galb = np.hstack([galb, s.data.geo_albedo[mask]])
        ftrn = np.hstack([ftrn, s.data.flux_transmission[mask]])
        absrad = np.hstack([absrad, s.data.absorbing_radius[mask]])

    data = np.array([tag, wl, wn, toaf, starf, galb, ftrn, absrad]).T

    cols = [
        "File Name",
        "Wavelength",
        "Wavenumber",
        "Top of Atmosphere Flux",
        "Stellar Flux",
        "Geometric Albedo",
        "Flux Transmission",
        "Absorbing Radius"
    ]

    df = pd.DataFrame(data, columns=cols)

    df.to_csv(savename)

    print "Saved %s" %savename

    return

def test_init():

    # Relative path to spectra files
    path = "spectra_files/"

    # Construct Spectrum object

    test1 = Spectrum(observation="Transmission",
                    star="Proxima Centauri",
                    planet="b",
                    description="FP Earth",
                    meta="1 bar; CO2/O2/CO",
                    reference="Meadows et al., 2017; Gao et al., 2015",
                    path_to_file=os.path.join(path, "smart_gao_1bar_update_xsec.trnst")
              )

    test2 = Spectrum(observation="Direct",
                    star="AD Leo",
                    planet="Earth",
                    description="No clouds; 100% Ocean",
                    meta="1 PAL O2",
                    reference="Segura et al. 2005",
                    path_to_file=os.path.join(path, "smart_Mstar_ADLeo_1PAL_clr_ocean_50_30000cm_r60only_toa_data.txt")
              )

    ltest = [test1, test2]

    return ltest

def test_write_meta():
    """
    """
    spectra = test_init()

    write_spectra_meta_csv(spectra, savename="test_meta.csv")

    return

def test_write_spectra():
    """
    """
    spectra = test_init()

    write_spectra_csv(spectra, savename="test_spectra.csv")

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

def test_molecules():

    # Create list of molecules and their features
    molecules = [
        Molecule("O2", "O_2", [0.2, 0.69, 0.76, 1.27]),
        Molecule("CO", "CO", [1.6, 2.35])
    ]

    # Create test csv file
    write_molecular_csv(molecules, savename="test_mols.csv")

if __name__ == "__main__":
    #test_init()
    #test_molecules()
    test_write_spectra()
    test_write_meta()
