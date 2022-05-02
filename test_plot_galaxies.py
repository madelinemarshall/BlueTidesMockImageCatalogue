import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from astropy.io import fits
import pickle
import mock_functions



filterList=[('Euclid','NISP','Y'),('HST','WFC3','f105w'),
                   ('JWST','NIRCAM','F090W'),
                   ('Roman','Roman','F106'),('Subaru','HSC','y'),('VISTA','VIRCAM','Y')]
        

def test_invalid_filter():
    z=11
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2,9.5))#z11

    out = cat.plot_galaxies('JWST','NIRCam','F444M')
    assert out == 'ERROR: Specified filter "F444M" is not in filter list, aborted'
    

def test_invalid_JWST_instrument():
    z=11
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2,9.5))#z11
    out = cat.plot_galaxies('JWST','MN','F444W')
    assert out =='ERROR: Specified JWST instrument "MN" is not NIRCAM or MIRI, aborted'
    
    
def test_invalid_instrument():
    z=11
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2,9.5))#z11
    cat.plot_galaxies('Euclid','INVALID','Y') 
    #Should correct itself and run as normal
    return 

def test_invalid_telescope():
    z=11
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2,9.5))#z11
    out = cat.plot_galaxies('Spitzer','NIRCam','F444M')
    assert out == "ERROR: Telescope 'Spitzer' doesn't match those available ('HST','JWST','Euclid','Roman','VISTA','Subaru')"

    
def test_z12_two_galaxies_one_filter():
    z=12
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1,2])
    cat.plot_galaxies('VISTA','VIRCAM','Y')
    plt.close('all')
    return     

    
def test_z12_all_telescopes_all_filters():
    z=12
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1,2])
    cat.plot_galaxies('ALL','ALL','ALL')
    plt.close('all')
    return   
 
    
def test_z12_all_telescopes_Y():
    z=12
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1,2])
    cat.plot_galaxies('ALL','ALL','Y')
    plt.close('all')
    return   
    
    
def test_z12_one_galaxy():
    z=12
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('VISTA','VIRCAM','ALL')
    cat.plot_galaxies('Euclid','NISP','ALL')
    cat.plot_galaxies('HST','WFC3','ALL')
    cat.plot_galaxies('JWST','NIRCAM','ALL')
    cat.plot_galaxies('JWST','MIRI','ALL')
    cat.plot_galaxies('Roman','WFI','ALL')
    cat.plot_galaxies('Subaru','HSC','ALL')
    cat.plot_galaxies('ALL','ALL','ALL')
    cat.plot_galaxies('ALL','ALL','Y')
    plt.close('all')
    return 



def test_z11():
    z=11
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('VISTA','VIRCAM','ALL')
    cat.plot_galaxies('Euclid','NISP','ALL')
    cat.plot_galaxies('HST','WFC3','ALL')
    cat.plot_galaxies('JWST','NIRCAM','ALL')
    cat.plot_galaxies('JWST','WFI','ALL')
    cat.plot_galaxies('Roman','Roman','ALL')
    cat.plot_galaxies('Subaru','HSC','ALL')
    cat.plot_galaxies('ALL','ALL','ALL')
    cat.plot_galaxies('ALL','ALL','Y')
    return 


def test_z10():
    z=10
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('VISTA','VIRCAM','ALL')
    cat.plot_galaxies('Euclid','NISP','ALL')
    cat.plot_galaxies('HST','WFC3','ALL')
    cat.plot_galaxies('JWST','NIRCAM','ALL')
    cat.plot_galaxies('JWST','MIRI','ALL')
    cat.plot_galaxies('Roman','WFI','ALL')
    cat.plot_galaxies('Subaru','HSC','ALL')
    plt.close('all')
    return 
    
    
def test_z9():
    z=9
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('VISTA','VIRCAM','ALL')
    cat.plot_galaxies('Euclid','NISP','ALL')
    cat.plot_galaxies('HST','WFC3','ALL')
    cat.plot_galaxies('JWST','NIRCAM','ALL')
    cat.plot_galaxies('JWST','MIRI','ALL')
    cat.plot_galaxies('Roman','WFI','ALL')
    cat.plot_galaxies('Subaru','HSC','ALL')
    cat.plot_galaxies('ALL','ALL','ALL')
    cat.plot_galaxies('ALL','ALL','Y')
    plt.close('all')
    return 
    
    
def test_z8():
    z=8
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('VISTA','VIRCAM','ALL')
    cat.plot_galaxies('Euclid','NISP','ALL')
    cat.plot_galaxies('HST','WFC3','ALL')
    cat.plot_galaxies('JWST','NIRCAM','ALL')
    cat.plot_galaxies('JWST','MIRI','ALL')
    cat.plot_galaxies('Roman','WFI','ALL')
    cat.plot_galaxies('Subaru','HSC','ALL')
    cat.plot_galaxies('ALL','ALL','ALL')
    cat.plot_galaxies('ALL','ALL','Y')
    plt.close('all')
    return 



def test_z7():
    z=7
    cat=mock_functions.Catalogue(z=z,catalogue_path='../')
    cat.select_galaxies_with_indices([1])
    # Probably over-kill to test each of the options like this, but just to be sure there are no issues
    cat.plot_galaxies('ALL','ALL','ALL')
    cat.plot_galaxies('ALL','ALL','Y')
    cat.plot_galaxies('VISTA','VIRCAM','ALL')
    cat.plot_galaxies('Euclid','NISP','ALL')
    cat.plot_galaxies('HST','WFC3','ALL')
    cat.plot_galaxies('JWST','NIRCAM','ALL')
    cat.plot_galaxies('JWST','MIRI','ALL')
    cat.plot_galaxies('Roman','WFI','ALL')
    cat.plot_galaxies('Subaru','HSC','ALL')
    plt.close('all')
    return 