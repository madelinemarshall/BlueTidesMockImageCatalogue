import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import bluetides_mock_catalogue_functions as mock_functions
import matplotlib
import astropy
import sys
import add_noise

version = ".".join(map(str, sys.version_info[:3]))
print('python=={}'.format(version))
print('numpy=={}'.format(np.__version__))
print('pandas=={}'.format(pd.__version__))
print('matplotlib=={}'.format(matplotlib.__version__))
print('astropy=={}'.format(astropy.__version__))

filterList = [('euclid', 'nisp', 'y'), ('hst', 'wfc3', 'f105w'),
              ('jwst', 'nircam', 'f090w'),
              ('roman', 'roman', 'f106'), ('subaru', 'hsc', 'y'), ('vista', 'vircam', 'y')]

image_path='/home/mmarshal/FinalImages/FullCatalogue/' #Directory containing the catalogue fits images
jwst_filepath = '/home/mmarshal/FinalImages/' #Directory containing the jwst_nircam_depths.csv file


def test_get_aperture_limits_invalid_tele():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])
 
    tele = 'invalid'
    errorstr = "ERROR: Telescope {} is not in available list: ['jwst', 'hst', 'roman', 'subaru', 'vista']".\
                   format(tele)
    image = add_noise.Image(cat,tele,'nircam','f150w') 
    out = image.get_aperture_limits()
    assert out == errorstr

def test_get_aperture_limits_tele():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])
    for tele in ['hst', 'roman', 'subaru', 'vista']:
        image = add_noise.Image(cat,tele,'nircam','f150w') 
        image.get_aperture_limits()
    return

def test_get_aperture_limits_jwst():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])
 
    image = add_noise.Image(cat,'jwst','nircam','f277w')
    image.get_aperture_limits(jwst_filepath)
    return

def test_get_aperture_limits_jwst_nofilepath():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])
 
    image = add_noise.Image(cat,'jwst','nircam','f277w')
    image.get_aperture_limits() #no filepath specified
    return

def test_add_noise():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])
 
    image = add_noise.Image(cat,'jwst','nircam','f277w')
    image.get_image_params()
    image.open_fits()

    image.get_aperture_limits(jwst_filepath)
    image.add_noise(10000)
    return

def test_add_noise_invalid_jwst_exposure():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])

    image = add_noise.Image(cat,'jwst','nircam','f277w')
    image.get_image_params()
    image.open_fits()

    image.get_aperture_limits(jwst_filepath)
    image.add_noise(20000)
    return

def test_add_noise_invalid_vista_exposureTime():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    cat.select_galaxies_with_indices([1,2])
    cat.catalogueSelected
 
    image = add_noise.Image(cat,'vista','vircam','y')
    image.get_image_params()
    image.open_fits()

    image.get_aperture_limits()
    image.add_noise(20000)
    return

def test_add_noise_invalid_vista_filter():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    
    cat.select_galaxies_with_indices([1,2])

    cat.catalogueSelected
 
    image = add_noise.Image(cat,'vista','vircam','z')
    image.get_image_params()
    image.open_fits()

    image.get_aperture_limits()
    image.add_noise(10000)
    return

def test_add_noise_invalid_vista_filter():
    z=12
    cat=mock_functions.Catalogue(z=z,image_path=image_path)
    
    cat.select_galaxies_with_indices([1,2])

    cat.catalogueSelected
 
    image = add_noise.Image(cat,'vista','vircam','y')
    image.get_image_params()
    image.open_fits()

    image.get_aperture_limits()
    image.add_noise(10000)
    return



