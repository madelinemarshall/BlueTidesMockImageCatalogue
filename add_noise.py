import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt


class empty:
    pass


class Background():
    """
    A class to represent a background, for adding to a noiseless image.

    Author: Jussi Kuusisto
    Original repository: https://github.com/jussikuusisto/background_noise_gen

    """

    def __init__(self, pixel_scale, aperture_f_limit, aperture_significance, aperture_radius, verbose=False):
        '''
        :param pixel_scale:             Pixel scale (arcsec/pixel)
        :param aperture_f_limit:        Aperture flux limit (nJy) (The flux detected at aperture_significance)
        :param aperture_significance:   Aperture significance (S/N)
        :param aperture_radius:         Aperture radius (in arcsec)
        :param verbose:                 Verbose mode toggle
        '''

        self.pixel_scale = pixel_scale

        self.aperture = empty()
        self.aperture.flux_limit = aperture_f_limit
        self.aperture.radius = aperture_radius / self.pixel_scale  # aperture radius in pixels
        self.aperture.significance = aperture_significance
        self.aperture.noise = self.aperture.flux_limit/self.aperture.significance  # nJy
        self.aperture.background = self.aperture.noise**2
        self.aperture.area = np.pi * self.aperture.radius**2

        self.pixel = empty()
        self.pixel.background = self.aperture.background/self.aperture.area
        self.pixel.noise = np.sqrt(self.pixel.background)  # nJy

        if verbose:
            print('assumed aperture radius: {0:.2f} pix'.format(self.aperture.radius))
            print('noise in aperture: {0:.2f} nJy'.format(self.aperture.noise))
            print('noise in pixel: {0:.2f} nJy'.format(self.pixel.noise))

    def create_background_image(self, CutoutWidth):
        '''
        :param CutoutWidth:             Cutout image size in pixels

        returns:                        image object: img.bkg gives the background noise map
        '''

        img = empty()
        img.pixel_scale = self.pixel_scale

        img.noise = self.pixel.noise * np.ones((CutoutWidth, CutoutWidth))
        img.wht = 1./img.noise**2

        random_map = np.random.randn(*img.noise.shape)
        img.bkg = self.pixel.noise * random_map  # nJy

        return img


def mag_to_flux(mag):
    """Convert AB magnitude to flux in nJy"""
    flux = 10**(-0.4*(mag-8.90))*1e9
    return flux


class Image():
    """
    A class to represent an image.

    Attributes:
    catalogue -- Catalogue: The catalogue of selected galaxies
    telescope -- Str: Telescope of the mock image. Current options are 'jwst', 'hst', 'roman', 'subaru' and 'vista'
    instrument -- Str: Instrument of the mock image
    filt -- Str: Filter of the mock image

    Methods:
    get_image_params -- Load in additional image parameters (FOV in kpc, width in arcsec, super_samp the rate at which
                        the native pixels are sub-sampled, and the native_pixel_scale in arcsec.
    open_fits -- Open the fits files for the galaxies in the catalogue, for the specified telescope/filter.
    add_noise -- Add noise to the noiseless images (from fits) based on the exposure time and telescope sensitivity.
    """

    def __init__(self, cat, telescope, instrument, filt):
        """
        Parameters:
        catalogue -- Catalogue: The catalogue of selected galaxies
        telescope -- Str: Telescope of the mock image
        instrument -- Str: Instrument of the mock image
        filt -- Str: Filter of the mock image
        """
        self.catalogue = cat
        self.telescope = telescope
        self.instrument = instrument
        self.filt = filt

    def get_image_params(self):
        """Read in the relevant image parameters from the fits header.

        Additional attributes:
        FOV -- Float: Field of view of the image in pkpc
        width -- Float: Width of the image in arcsec
        super_samp -- Float: The rate at which the native pixels are sub-sampled
        native_pixel_scale -- Float: The native (non-sub-sampled)  pixel scale of the instrument, in arcsec/pixel.
        """
        if self.catalogue.z == 7:
            ff = fits.open(self.catalogue.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                           format(self.telescope, self.instrument, self.catalogue.z, 1, self.filt))
        else:
            ff = fits.open(self.catalogue.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.\
                           format(self.telescope, self.instrument, self.catalogue.z, self.filt))

        self.FOV = ff[0].header['FIELDOFVIEW_PKPC']
        self.width = ff[0].header['FIELDOFVIEW_ARCSEC']
        self.super_samp = ff[0].header['SUPERSAMPLINGRATE']
        self.native_pixel_scale = ff[0].header['RESOLUTION_ARCSEC']*ff[0].header['SUPERSAMPLINGRATE']

        ff.close()
        return

    def open_fits(self):
        """Open the fits files for the galaxies in the catalogue, for the specified telescope/filter.

        Additional class attributes:
        data -- Dict: Dictionary containing the 2D image arrays for each of the galaxies in the selected catalogue.
        """
        data = {}

        if self.catalogue.z == 7:
            ff = {}
            for fn in np.unique(self.catalogue.catalogueSelected['fileNumber'].values):
                ff[fn] = fits.open(self.catalogue.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                   format(self.telescope, self.instrument, self.catalogue.z, fn, self.filt))
            for ii in range(0, self.catalogue.nGalsSelected):
                data[ii] = ff[self.catalogue.catalogueSelected['fileNumber'].values[ii]][
                              self.catalogue.catalogueSelected['extensionNumber'].values[ii]].data

            for fn in np.unique(self.catalogue.catalogueSelected['fileNumber'].values):
                ff[fn].close()

        else:
            ff = fits.open(self.catalogue.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.\
                           format(self.telescope, self.instrument, self.catalogue.z, self.filt))
            for ii in range(0, self.catalogue.nGalsSelected):
                data[ii] = ff[self.catalogue.catalogueSelected['extensionNumber'].values[ii]].data

            ff.close()

        self.data = data

        return

    def get_aperture_limits(self, jwst_filepath=False):
        """Read in the aperture parameters relating to the noise in each image.

        Parameters:
        jwst_filepath -- Str: Directory containing the jwst_nircam_depths.csv file.
                         Recommended if self.telescope='jwst' (default False)

        Additional attributes:
        ap_f_limits -- Dict: The aperture flux limits for the given telescope, as a dictionary with keys = filters
                         and values = dictionary, where the second dictionary has keys = available exposure times
                         and values = corresponding flux limit in nJy. i.e. ap_f_limits[FILT][EXP_TIME] = LIM nJy
        ap_sig -- Int: The significance of the aperture flux limit, in multiples of sigma. Generally 5 or 10 sigma.
        ap_radius -- Float: The radius of the aperture used to calculate the aperture flux limit. In arcsec. This is a
                        dictionary for JWST as there are 3 pixel scales, for the different filters.

        Notes:
        Only a select number of exposure times are currently available. These are generally those used in the Mock
        Catalogue release paper, see for details. Additional exposure times can be hard-coded into to this function
        as necessary.
        """
        tele = self.telescope.lower()

        if tele == 'jwst':
            if jwst_filepath is not False:
                ap_f_limits_NIRCAM = pd.read_csv(jwst_filepath+'jwst_nircam_depths.csv', index_col=0)
            else:
                ap_f_limits_NIRCAM = pd.DataFrame(columns=['jwst.nircam.f090w', 'jwst.nircam.f115w',
                                                           'jwst.nircam.f150w', 'jwst.nircam.f200w',
                                                           'jwst.nircam.f277w', 'jwst.nircam.f356w',
                                                           'jwst.nircam.f410w', 'jwst.nircam.f444w'])
                print('WARNING: tele is jwst, but jwst_filepath is not specified. ' +
                      'Most NIRCAM exposure times will not have a specified depth.')
            # Limits for NIRCam wide-band plus F410M filters, from 0.55 to 10ks exposures, spaced by 0.05ks
            # Extracted from Figure 1 in:
            # https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-predicted-performance/nircam-imaging-sensitivity

            # COSMOS-Web depths
            cosmosdata = {'jwst.nircam.f090w': np.nan, 'jwst.nircam.f115w': mag_to_flux(27.36)*2,
                          'jwst.nircam.f150w': mag_to_flux(27.59)*2, 'jwst.nircam.f200w': np.nan,
                          'jwst.nircam.f277w': mag_to_flux(28.05)*2, 'jwst.nircam.f356w': np.nan,
                          'jwst.nircam.f410m': np.nan, 'jwst.nircam.f444w': mag_to_flux(27.72)*2}
            # *2 because they give 5sigma values. From proposal pdf, https://www.stsci.edu/jwst/phase2-public/1727.pdf

            # JADES-Medium depths
            jadesdata = {'jwst.nircam.f090w': mag_to_flux(28.6), 'jwst.nircam.f115w': mag_to_flux(28.8),
                         'jwst.nircam.f150w': mag_to_flux(28.9), 'jwst.nircam.f200w': mag_to_flux(29),
                         'jwst.nircam.f277w': mag_to_flux(28.6), 'jwst.nircam.f356w': mag_to_flux(28.6),
                         'jwst.nircam.f410m': mag_to_flux(28.1), 'jwst.nircam.f444w': mag_to_flux(28.3)}
            # From Rieke+19, doi:10.1017/s1743921319008950

            ap_f_limits_NIRCAM.loc[774] = cosmosdata
            ap_f_limits_NIRCAM.loc[12000] = jadesdata

            ap_f_limits_MIRI = pd.DataFrame(columns=['jwst.miri.f560w', 'jwst.miri.f770w'], index=[10000])
            ap_f_limits_MIRI['jwst.miri.f560w'][10000] = 130
            ap_f_limits_MIRI['jwst.miri.f770w'][10000] = 240
            # From https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-predicted-performance/miri-sensitivity

            ap_f_limits = pd.concat([ap_f_limits_NIRCAM, ap_f_limits_MIRI], axis=1)
            ap_sig = 10

            pixel_scale = {}
            ap_radius = {}

            for filt in ['jwst.nircam.f090w', 'jwst.nircam.f115w', 'jwst.nircam.f150w', 'jwst.nircam.f200w']:
                pixel_scale[filt] = 0.031
                ap_radius[filt] = 2.58*pixel_scale[filt]         # aperture radius in arcsec
            for filt in ['jwst.nircam.f277w', 'jwst.nircam.f356w', 'jwst.nircam.f410m', 'jwst.nircam.f444w']:
                pixel_scale[filt] = 0.063
                ap_radius[filt] = 2.58*pixel_scale[filt]         # aperture radius in arcsec
            for filt in ['jwst.miri.f560w', 'jwst.miri.f770w']:
                pixel_scale[filt] = 0.11
                ap_radius[filt] = 2.58*pixel_scale[filt]         # aperture radius in arcsec

        elif tele == 'hst':
            filt_list = ['hst.wfc3.f125w', 'hst.wfc3.f160w', 'hst.wfc3.f105w']
            ap_f_limits = pd.DataFrame(columns=filt_list, index=[1900, 3600, 10000])
            ap_f_limits['hst.wfc3.f160w'][3600] = mag_to_flux(27.0)  # CANDELS EGS, from Grogin+11
            ap_f_limits['hst.wfc3.f125w'][1900] = mag_to_flux(27.2)  # CANDELS EGS, from Grogin+11
            ap_f_limits['hst.wfc3.f105w'][10000] = mag_to_flux(27.93)  # From the exposure time calculator
            ap_sig = 5
            ap_radius = 0.2         # aperture radius in arcsec

        elif tele == 'roman':
            ap_f_limits = pd.DataFrame(columns=['roman.wfi.f106'], index=[10000])
            ap_f_limits['roman.wfi.f106'][10000] = mag_to_flux(26.83)  # From the exposure time calculator
            ap_sig = 10
            ap_radius = 0.286         # aperture radius in arcsec

        elif tele == 'subaru':
            ap_f_limits = pd.DataFrame(columns=['subaru.hsc.y'], index=[10000])
            ap_f_limits['subaru.hsc.y'][10000] = mag_to_flux(23.94)  # From the exposure time calculator
            ap_sig = 10
            ap_radius = 2         # aperture radius in arcsec

        elif tele == 'vista':
            ap_f_limits = pd.DataFrame(columns=['vista.vircam.y'], index=[10000])
            ap_f_limits['vista.vircam.y'][10000] = mag_to_flux(23.50)  # From the exposure time calculator
            ap_sig = 10
            ap_radius = 2         # aperture radius in arcsec

        else:
            errorstr = "ERROR: Telescope {} is not in available list: ['jwst', 'hst', 'roman', 'subaru', 'vista']".\
                       format(tele)
            print(errorstr)
            return errorstr

        self.ap_f_limits = ap_f_limits
        self.ap_sig = ap_sig
        self.ap_radius = ap_radius
        return

    def add_noise(self, exp_time):
        """Add noise to the noiseless images (from fits) based on the exposure time and telescope sensitivity.

        Parameters:
        exp_time -- Float: the exposure time of the image, in seconds.

        Additional attributes:
        noisy_data -- Dict: Dictionary containing the 2D image arrays for each of the galaxies in the selected catalogue
                            Noise has been added, on top of the noiseless images saved in self.data
        background_level -- Float: The noise level (per pixel, not per aperture) added to the images, in nJy
        exposure_time -- Float: Exposure time used in the noisy_data images, in seconds
        """
        if exp_time in self.ap_f_limits[self.telescope+'.'+self.instrument+'.'+self.filt].keys():
            self.ap_f_limit = self.ap_f_limits[self.telescope+'.'+self.instrument+'.'+self.filt][exp_time]
            if self.telescope == 'jwst':
                self.ap_radius = self.ap_radius[self.telescope+'.'+self.instrument+'.'+self.filt]
        else:
            print('ERROR: exposure time {} is unavailable. Choose another, '.format(exp_time) +
                  'or code in the depth at this exposure time in add_noise.py, get_aperture_limits')
            print('Available times: ', list(self.ap_f_limits[self.telescope+'.'+self.instrument+'.'+self.filt].keys()))
            return

        np.random.seed(0)

        self.noisy_data = {}

        for ii in range(0, self.catalogue.nGalsSelected):
            # Add shot noise
            full_img = self.data[ii] * exp_time
            full_img[full_img < 0] = 0
            noisy_full_img = np.random.poisson(full_img)
            img_data = noisy_full_img / exp_time

            # create background image object (cutoutwidth in pixels)
            background_object = Background(pixel_scale=self.native_pixel_scale/self.super_samp,
                                           aperture_f_limit=self.ap_f_limit, aperture_significance=self.ap_sig,
                                           aperture_radius=self.ap_radius, verbose=False)

            img_bkg = background_object.create_background_image(np.shape(img_data)[0])
            img_bkg_data = img_bkg.bkg  # nJy

            self.noisy_data[ii] = img_data+img_bkg_data

        self.background_level = background_object.pixel.noise
        self.exposure_time = exp_time

        return

    def plot_noisy_images(self):
        """Plot the images with added noise"""

        fig, ax = plt.subplots(1, self.catalogue.nGalsSelected, figsize=(8, 3))
        for ii in range(0, self.catalogue.nGalsSelected):
            ax[ii].imshow(self.noisy_data[ii])

            res = (self.native_pixel_scale/self.super_samp)
            ax[ii].set_yticks([0, 1/res])
            ax[ii].set_yticklabels(['0', '1'])
            ax[ii].set_ylim([-0.5, (self.width-res)/res-0.5])
            ax[ii].set_xticks([0, 1/res])
            ax[ii].set_xticklabels(['0', '1'])
            ax[ii].set_xlim([-0.5, (self.width-res)/res-0.5])
            ax[ii].set_xlabel('arcsec')
            ax[ii].set_ylabel('arcsec')
        return fig

    def plot_noiseless_images(self):
        """Plot the noise-less images"""

        fig, ax = plt.subplots(1, self.catalogue.nGalsSelected, figsize=(8, 3))
        for ii in range(0, self.catalogue.nGalsSelected):
            ax[ii].imshow(self.data[ii])

            res = (self.native_pixel_scale/self.super_samp)  # pixel scale in arcsec
            ax[ii].set_yticks([0, 1/res])  # plot a tick at 0 and 1 arcsec
            ax[ii].set_yticklabels(['0', '1'])
            ax[ii].set_ylim([-0.5, (self.width-res)/res-0.5])
            ax[ii].set_xticks([0, 1/res])
            ax[ii].set_xticklabels(['0', '1'])
            ax[ii].set_xlim([-0.5, (self.width-res)/res-0.5])
            ax[ii].set_xlabel('arcsec')
            ax[ii].set_ylabel('arcsec')
        return fig
