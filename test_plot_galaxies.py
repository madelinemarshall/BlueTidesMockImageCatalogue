import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import bluetides_mock_catalogue_functions as mock_functions
import matplotlib
import astropy
import sys

version = ".".join(map(str, sys.version_info[:3]))
print('python=={}'.format(version))
print('numpy=={}'.format(np.__version__))
print('pandas=={}'.format(pd.__version__))
print('matplotlib=={}'.format(matplotlib.__version__))
print('astropy=={}'.format(astropy.__version__))

filterList = [('euclid', 'nisp', 'y'), ('hst', 'wfc3', 'f105w'),
              ('jwst', 'nircam', 'f090w'),
              ('roman', 'roman', 'f106'), ('subaru', 'hsc', 'y'), ('vista', 'vircam', 'y')]


def test_invalid_filter():
    z = 11
    # cat=mock_functions.Catalogue(z=z,image_path='../')
    # NOTE: if running in a directory other than where the catalogue is saved, will need to specify an image_path
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2, 9.5))

    out = cat.plot_galaxies('jwst', 'NIRCam', 'f444m')
    assert out == 'ERROR: Specified filter "f444m" is not in filter list, aborted'


def test_invalid_jwst_instrument():
    z = 11
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2, 9.5))
    out = cat.plot_galaxies('jwst', 'mn', 'F444W')
    assert out == 'ERROR: Specified jwst instrument "mn" is not nircam, miri, or all, aborted'


def test_invalid_instrument():
    z = 11
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2, 9.5))
    cat.plot_galaxies('euclid', 'INVALID', 'y')
    # Should correct itself and run as normal
    return


def test_invalid_telescope():
    z = 11
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_constraints(stellar_mass_constraint=(9.2, 9.5))
    out = cat.plot_galaxies('Spitzer', 'NIRCam', 'F444M')
    assert out == "ERROR: Telescope 'spitzer' doesn't match those available " +\
                  "('hst', 'jwst', 'euclid', 'roman', 'vista', 'subaru')"


def test_z12_two_galaxies_one_filter():
    z = 12
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1, 2])
    cat.plot_galaxies('vista', 'vircam', 'y')
    plt.close('all')
    return


def test_z12_all_telescopes_all_filters():
    z = 12
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1, 2])
    cat.plot_galaxies('all', 'all', 'all')
    plt.close('all')
    return


def test_z12_all_telescopes_Y():
    z = 12
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1, 2])
    cat.plot_galaxies('all', 'all', 'y')
    plt.close('all')
    return


def test_z12_one_galaxy():
    z = 12
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    # Probably over-kill to test each of the options like this, but just to be sure there are no issues
    cat.plot_galaxies('vista', 'vircam', 'all')
    cat.plot_galaxies('euclid', 'nisp', 'all')
    cat.plot_galaxies('hst', 'wfc3', 'all')
    cat.plot_galaxies('jwst', 'nircam', 'all')
    cat.plot_galaxies('jwst', 'miri', 'all')
    cat.plot_galaxies('jwst', 'all', 'all')
    cat.plot_galaxies('roman', 'wfi', 'all')
    cat.plot_galaxies('subaru', 'hsc', 'all')
    cat.plot_galaxies('all', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'y')
    plt.close('all')
    return


def test_z12_one_galaxy_colorbar():
    z = 12
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    # Probably over-kill to test each of the options like this, but just to be sure there are no issues
    cat.plot_galaxies('vista', 'vircam', 'all', colorbar=True)
    cat.plot_galaxies('euclid', 'nisp', 'all', colorbar=True)
    cat.plot_galaxies('hst', 'wfc3', 'all', colorbar=True)
    cat.plot_galaxies('jwst', 'nircam', 'all', colorbar=True)
    cat.plot_galaxies('jwst', 'miri', 'all', colorbar=True)
    cat.plot_galaxies('jwst', 'all', 'all', colorbar=True)
    cat.plot_galaxies('roman', 'wfi', 'all', colorbar=True)
    cat.plot_galaxies('subaru', 'hsc', 'all', colorbar=True)
    cat.plot_galaxies('all', 'all', 'all', colorbar=True)
    cat.plot_galaxies('all', 'all', 'y', colorbar=True)
    plt.close('all')
    return


def test_z11():
    z = 11
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('vista', 'vircam', 'all')
    cat.plot_galaxies('jwst', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'y')
    plt.close('all')
    return


def test_z11_multiple_galaxies_colorbar():
    z = 11
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1, 2, 3])
    cat.plot_galaxies('vista', 'vircam', 'all', colorbar=True)
    cat.plot_galaxies('jwst', 'all', 'all', colorbar=True)
    cat.plot_galaxies('all', 'all', 'all', colorbar=True)
    cat.plot_galaxies('all', 'all', 'y', colorbar=True)
    plt.close('all')
    return


def test_z10():
    z = 10
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('euclid', 'nisp', 'all')
    cat.plot_galaxies('jwst', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'y')
    plt.close('all')
    return


def test_z9():
    z = 9
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('jwst', 'all', 'all')
    cat.plot_galaxies('subaru', 'hsc', 'all')
    cat.plot_galaxies('all', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'y')
    plt.close('all')
    return


def test_z8():
    z = 8
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('hst', 'wfc3', 'all')
    cat.plot_galaxies('jwst', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'y')
    plt.close('all')
    return


def test_z7():
    z = 7
    cat = mock_functions.Catalogue(z=z)
    cat.select_galaxies_with_indices([1])
    cat.plot_galaxies('all', 'all', 'all')
    cat.plot_galaxies('all', 'all', 'y')
    cat.plot_galaxies('vista', 'vircam', 'all')
    cat.plot_galaxies('jwst', 'all', 'all')
    plt.close('all')
    return
