import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator, FormatStrFormatter


rcParams['font.size'] = (9)
rcParams['figure.figsize'] = (4, 3.4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


class Catalogue:
    """
    A class to represent a catalogue.

    Attributes:
    image_path -- Str: Path to BlueTides Mock Catalogue .csv and .fits files
    z -- Int: redshift of interest
    catalogue -- Table: Properties of galaxies at redshift z as read from .csv file
    stellarMass -- Array: Stellar mass of all galaxies at redshift z, in solar masses
    BHmass -- Array: Black hole mass of all galaxies at redshift z, in solar masses
    haloMass -- Array: Halo mass of all galaxies at redshift z, in solar masses
    nGals -- Int: Number of galaxies in the catalgue at redshift z

    Methods:
    print_constraint_options -- Print a list of catalogue columns with their minimum and maximum values
    select_galaxies_with_constraints -- Select galaxies that satisfy the requested property constraints
    select_galaxies_with_indices -- Select galaxies from the catalogue with specified indices
    plot_galaxies -- Plot the galaxies in the selected sub-catalogue
    """

    def __init__(self, z=7, image_path='/home/mmarshal/FinalImages/FullCatalogue/'):
        """Load in catalogue from the csv file and save columns as attributes for ease of access.

        Keyword arguments:
        z -- Redshift of snapshot. Available options are 7, 8, 9, 10, 11 or 12 (default 7)
        image_path -- Directory containing the hlsp fits and csv files
                      (default MM's directory - can change to your own here)
        """
        cat = pd.read_csv(image_path+'hlsp_bluetides_multi_multi_all_multi_v1_sim.csv', index_col=0)
        self.image_path = image_path
        self.z = z
        cat = cat.loc[(cat["redshift"] == z)]
        self.catalogue = cat
        self.stellarMass = np.array(self.catalogue['stellarMass'])
        self.BHmass = np.array(self.catalogue['BHmass'])
        self.haloMass = np.array(self.catalogue['haloMass'])
        self.nGals = len(self.stellarMass)
        return

    def print_constraint_options(self, lum_constraint_type='lum', flux_constraint_type='flux'):
        """Print a list of columns with their minimum and maximum values.

        Keyword arguments:
        lum_constraint_type -- Str: 'lum' for UV constraint specified in log luminosity, in erg/s/Hz,
                                  or 'mag' for UV constraint specified in magnitudes.
        flux_constraint_type -- Str: 'flux' for telescope constraint specified in log flux, in nJy,
                                  or 'mag' for telescope constraint specified in magnitudes.
        """
        cols = self.catalogue.columns.values
        for ii, col in enumerate(cols):
            width = 32
            if ii < 3:
                print('{: <{}}:  (Min, Max) = {:g}, {:g}'.format(col, width, np.min(self.catalogue[col]),
                                                                 np.max(self.catalogue[col])))

            elif col[0:4] == 'flux':
                if flux_constraint_type == 'flux':
                    print('{: <{}}:  (Min, Max) = {:0.1f}, {:0.1f}'.format(col+' [nJy]', width,
                                                                           np.min(self.catalogue[col]),
                                                                           np.max(self.catalogue[col])))
                elif flux_constraint_type == 'mag':
                    max_mag = -2.5*np.log10(np.min(self.catalogue[col])*1e-9)+8.90
                    min_mag = -2.5*np.log10(np.max(self.catalogue[col])*1e-9)+8.90
                    print('{: <{}}:  (Min, Max) = {:0.1f}, {:0.1f}'.format(col+' [AB mag]', width, min_mag, max_mag))

                else:
                    errorstr = "ERROR: Specified constraint type '{}' is not 'flux' or 'mag', aborted".\
                               format(flux_constraint_type)
                    print(errorstr)
                    return

            elif col == 'lum_FUV':
                if lum_constraint_type == 'lum':
                    print('{: <{}}:  (Min, Max) = {:0.2e}, {:0.2e}'.format(col + ' [erg/s/Hz]', width,
                                                                           np.min(self.catalogue[col]),
                                                                           np.max(self.catalogue[col])))
                elif lum_constraint_type == 'mag':
                    max_mag = lum_to_mag(np.max(self.catalogue[col]), self.z)
                    min_mag = lum_to_mag(np.min(self.catalogue[col]), self.z)
                    print('{: <{}}:  (Min, Max) = {:0.1f}, {:0.1f}'.format(col + ' [AB mag]', width, min_mag, max_mag))

                else:
                    errorstr = "ERROR: Specified constraint type '{}' is not 'lum' or 'mag', aborted".\
                               format(flux_constraint_type)
                    print(errorstr)
                    return

            elif col[0:3] == 'rad':
                    print('{: <{}}:  (Min, Max) = {:0.2f}, {:0.2f}'.format(col+' [pkpc]', width,
                                                                           np.min(self.catalogue[col]),
                                                                           np.max(self.catalogue[col])))
            elif col == 'BHluminosity':
                print('{: <{}}:  (Min, Max) = {:0.2e} = {:0.2e}'.format(col+' [erg/s]', width,
                                                                        np.min(self.catalogue[col]),
                                                                        np.max(self.catalogue[col])))

            else:
                print('{: <{}}:  (Min, Max) = {:0.2e} = {:0.2e}'.format(col+' [Msolar]', width,
                                                                        np.min(self.catalogue[col]),
                                                                        np.max(self.catalogue[col])))
        return

    def select_galaxies_with_constraints(self, stellar_mass_constraint=None, bh_mass_constraint=None,
                                         halo_mass_constraint=None, uv_lum_constraint=None, uv_constraint_type='lum',
                                         flux_constraint=None, flux_constraint_type='flux', radius_constraint=None):
        """Select galaxies from the catalogue that satisfy the requested property constraints.

        Keyword arguments:
        stellar_mass_constraint -- Tuple: (min stellar mass, max stellar mass). Specified in log space (default None)
        bh_mass_constraint -- Tuple: (min black hole mass, max black hole mass). Specified in log space (default None)
        halo_mass_constraint -- Tuple: (min halo mass, max halo mass). Specified in log space (default None)
        uv_lum_constraint -- Tuple: (min FUV (1500A) luminosity/magnitude, max FUV luminosity/magnitude) (default None)
        uv_constraint_type -- Str: 'lum' for constraint specified in log luminosity, in log(erg/s/Hz) (default),
                                or 'mag' for constraint specified in magnitudes.
        flux_constraint -- List of tuples: [(min flux filter 1, max flux in filter 1), ...,
                                            (min flux in filter N, max flux in filter N)] (default None)
        flux_constraint_type -- Str: 'flux' for constraint specified in flux, in nJy (default),
                                  or 'mag' for constraint specified in magnitudes.
        radius_constraint -- List of tuples: [(min radius filter 1, max radius in filter 1), ...,
                                              (min flux in filter N, max flux in filter N)], in pkpc (default None)

        Additional class attributes:
        catalogueSelected -- Table: Properties of galaxies in the selected sub-catalogue
        stellarMassSelected -- Array: Stellar mass of selected galaxies, in solar masses
        BHmassSelected -- Array: Black hole mass of selected galaxies, in solar masses
        haloMassSelected -- Array: Halo mass of selected galaxies, in solar masses
        nGalsSelected -- Int: Number of galaxies in the selected sub-catalogue
        """

        # Check which galaxies match the specified constraints
        if stellar_mass_constraint:
            stellar_mass_constraint = (self.stellarMass > 10**stellar_mass_constraint[0]) &\
                                      (self.stellarMass < 10**stellar_mass_constraint[1])
        else:
            stellar_mass_constraint = np.ones(self.nGals, dtype='Bool')

        if bh_mass_constraint:
            bh_mass_constraint = (self.BHmass > 10**bh_mass_constraint[0]) & (self.BHmass < 10**bh_mass_constraint[1])
        else:
            bh_mass_constraint = np.ones(self.nGals, dtype='Bool')

        if halo_mass_constraint:
            halo_mass_constraint = (self.haloMass > 10**halo_mass_constraint[0]) &\
                                   (self.haloMass < 10**halo_mass_constraint[1])
        else:
            halo_mass_constraint = np.ones(self.nGals, dtype='Bool')

        if uv_lum_constraint:
            if uv_constraint_type == 'lum':
                uv_constraint = (self.catalogue['lum_FUV'] > 10**uv_lum_constraint[0]) &\
                                         (self.catalogue['lum_FUV'] < 10**uv_lum_constraint[1])
            elif uv_constraint_type == 'mag':
                uv_lum_constraint_0 = mag_to_lum(uv_lum_constraint[0], self.z)
                uv_lum_constraint_1 = mag_to_lum(uv_lum_constraint[1], self.z)
                uv_constraint = (self.catalogue['lum_FUV'] > uv_lum_constraint_0) &\
                                (self.catalogue['lum_FUV'] < uv_lum_constraint_1)

            else:
                errorstr = "ERROR: Specified constraint type '{}' is not 'lum' or 'mag', aborted".\
                           format(flux_constraint_type)
                print(errorstr)
                self.stellarMassSelected = []
                self.BHmassSelected = []
                self.haloMassSelected = []
                self.nGalsSelected = 0
                self.catalogueSelected = []
                return errorstr

        else:
            uv_constraint = np.ones(self.nGals, dtype='Bool')

        if flux_constraint:
            nFilters = len(flux_constraint)
            flux_constraint_allfilt = np.ones(self.nGals, dtype='Bool')
            for ff in range(0, nFilters):
                filt = flux_constraint[ff][0]
                if flux_constraint_type == 'flux':
                    flux_constraint_ff = (self.catalogue['flux_'+filt.lower()] > flux_constraint[ff][1]) &\
                                         (self.catalogue['flux_'+filt.lower()] < flux_constraint[ff][2])
                elif flux_constraint_type == 'mag':
                    mag_ff = -2.5*np.log10(self.catalogue['flux_'+filt.lower()]*1e-9)+8.90
                    if flux_constraint[ff][2] > flux_constraint[ff][1]:
                        flux_constraint_ff = (mag_ff > flux_constraint[ff][1]) &\
                                             (mag_ff < flux_constraint[ff][2])
                    elif flux_constraint[ff][2] < flux_constraint[ff][1]:
                        # In case magnitudes are provided in reverse order (i.e. min brightness instead of min value)
                        flux_constraint_ff = (mag_ff < flux_constraint[ff][1]) &\
                                             (mag_ff > flux_constraint[ff][2])
                    else:
                        errorstr = "ERROR: need one magnitude constraint to be larger than the other, aborted"
                        print(errorstr)
                        self.stellarMassSelected = []
                        self.BHmassSelected = []
                        self.haloMassSelected = []
                        self.nGalsSelected = 0
                        self.catalogueSelected = []
                        return errorstr
                else:
                    errorstr = "ERROR: Specified constraint type '{}' is not 'flux' or 'mag', aborted".\
                               format(flux_constraint_type)
                    print(errorstr)
                    self.stellarMassSelected = []
                    self.BHmassSelected = []
                    self.haloMassSelected = []
                    self.nGalsSelected = 0
                    self.catalogueSelected = []
                    return errorstr

                flux_constraint_allfilt = np.logical_and(flux_constraint_ff, flux_constraint_allfilt)
            flux_constraint = flux_constraint_allfilt
        else:
            flux_constraint = np.ones(self.nGals, dtype='Bool')

        if radius_constraint:
            nFilters = len(radius_constraint)
            radius_constraint_allfilt = np.ones(self.nGals, dtype='Bool')
            for ff in range(0, nFilters):
                filt = radius_constraint[ff][0]
                radius_constraint_ff = (self.catalogue['radius_'+filt.lower()] > radius_constraint[ff][1]) &\
                                       (self.catalogue['radius_'+filt.lower()] < radius_constraint[ff][2])
                radius_constraint_allfilt = np.logical_and(radius_constraint_ff, radius_constraint_allfilt)
            radius_constraint = radius_constraint_allfilt
        else:
            radius_constraint = np.ones(self.nGals, dtype='Bool')

        selection = stellar_mass_constraint & bh_mass_constraint & halo_mass_constraint & uv_constraint &\
            flux_constraint & radius_constraint

        if len(selection[selection]) == 0:
            print('Error: No galaxies satisfy these criteria')
            self.stellarMassSelected = []
            self.BHmassSelected = []
            self.haloMassSelected = []
            self.nGalsSelected = 0
            self.catalogueSelected = []

        else:
            self.stellarMassSelected = np.array(self.catalogue['stellarMass'][selection])
            self.BHmassSelected = np.array(self.catalogue['BHmass'][selection])
            self.haloMassSelected = np.array(self.catalogue['haloMass'][selection])
            self.nGalsSelected = len(self.stellarMass[selection])
            self.catalogueSelected = self.catalogue[selection]
        return

    def select_galaxies_with_indices(self, indices):
        """Select galaxies from the catalogue with specified indices.

        Arguments:
        indices -- list of required indices

        Additional class attributes:
        catalogueSelected -- Table: Properties of galaxies in the selected sub-catalogue
        stellarMassSelected -- Array: Stellar mass of selected galaxies, in solar masses
        BHmassSelected -- Array: Black hole mass of selected galaxies, in solar masses
        haloMassSelected -- Array: Halo mass of selected galaxies, in solar masses
        nGalsSelected -- Int: Number of galaxies in the selected sub-catalogue
        """
        full_inds = self.catalogue.index
        selection = np.in1d(full_inds, indices)

        self.stellarMassSelected = np.array(self.catalogue['stellarMass'][selection])
        self.BHmassSelected = np.array(self.catalogue['BHmass'][selection])
        self.haloMassSelected = np.array(self.catalogue['haloMass'][selection])
        self.nGalsSelected = len(self.stellarMass[selection])
        self.catalogueSelected = self.catalogue[selection]
        return

    def plot_galaxies(self, telescope, instrument, filt, colorbar=False):
        """Plot the galaxies in the selected sub-catalogue

        Arguments:
        telescope -- Str: Requested telescope. Options are the telescope names or 'all'
        instrument -- Str: Requested instrument. Options are the specific instrument, 'all' if telescope is 'all',
                      or 'all' for both MIRI and NIRCam for JWST
        filt -- Str: Requested filter. Options are the specific filter, 'all', or 'y' if telescope is 'all' to plot
                all y-band equivalents

        Returns:
        figure(s)
        """
        # Telescope, instrument and filters all lowercase
        telescope = telescope.lower()
        instrument = instrument.lower()
        filt = filt.lower()

        ##### All telescopes with Y band #####
        if telescope == 'all' and filt == 'y':
            ## Plot first set of telescopes, all have 6x6kpc FOV ##
            filterList = [('jwst', 'nircam', 'f115w'), ('hst', 'wfc3', 'f105w'), ('roman', 'wfi', 'f106')]

            # Plot for one galaxy
            if self.nGalsSelected == 1:
                if colorbar is True:
                    fig, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                           figsize=(len(filterList)*2, self.nGalsSelected*2),
                                           gridspec_kw={'hspace': 0.15, 'wspace': 0.35})
                else:
                    fig, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                           figsize=(len(filterList)*2, self.nGalsSelected*2),
                                           gridspec_kw={'hspace': 0.1, 'wspace': 0.1})

                for jj in range(0, len(filterList)):
                    ax[jj].set_xlabel('{} {} {}'.format(filterList[jj][0].upper(), filterList[jj][1].upper(),
                                                        filterList[jj][2].upper()))
                    if self.z == 7:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                       format(filterList[jj][0], filterList[jj][1], self.z,
                                              self.catalogueSelected['fileNumber'].values[0], filterList[jj][2]))
                    else:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                       format(filterList[jj][0], filterList[jj][1], self.z, filterList[jj][2]))

                    im = ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data, cmap='Greys')
                    if colorbar is True:
                        divider = make_axes_locatable(ax[jj])
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                        cbar = plt.colorbar(im, cax=cax)
                        cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
                        cbar.ax.set_title('nJy', loc='left', pad=-5)

                    res = ff[0].header['RESOLUTION_PKPC']
                    ax[jj].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xlabel('pkpc')
                    ff.close()

            # Plot for multiple galaxies
            else:
                if colorbar is True:
                    fig, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                           figsize=(len(filterList)*1.8, self.nGalsSelected*1.5),
                                           gridspec_kw={'hspace': 0.15, 'wspace': 0.35})
                else:
                    fig, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                           figsize=(len(filterList)*1.5, self.nGalsSelected*1.5),
                                           gridspec_kw={'hspace': 0.1, 'wspace': 0.1})

                for jj in range(0, len(filterList)):
                    ax[0, jj].set_title('{}\n{} {}'.format(filterList[jj][0].upper(), filterList[jj][1].upper(),
                                                           filterList[jj][2].upper()), pad=10)
                    for ii in range(0, self.nGalsSelected):
                        if self.z == 7:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                           format(filterList[jj][0], filterList[jj][1], self.z,
                                                  self.catalogueSelected['fileNumber'].values[ii], filterList[jj][2]))
                        else:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                           format(filterList[jj][0], filterList[jj][1], self.z, filterList[jj][2]))

                        im = ax[ii, jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,
                                               cmap='Greys')
                        if colorbar is True:
                            #plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
                            divider = make_axes_locatable(ax[ii, jj])
                            cax = divider.append_axes('right', size='5%', pad=0.05)
                            cbar = plt.colorbar(im, cax=cax)
                            cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
                            cbar.ax.set_title('nJy', loc='left', pad=-5)

                        res = ff[0].header['RESOLUTION_PKPC']
                        if jj == 0:
                            ax[ii, jj].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii, jj].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            ax[ii, jj].set_ylabel('pkpc')
                        else:
                            ax[ii, jj].set_yticks([])

                        if ii < self.nGalsSelected-1:
                            ax[ii, jj].set_xticks([])
                        else:
                            ax[ii, jj].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii, jj].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii, jj].set_xlim(0, ff[0].header['FIELDOFVIEW_PKPC']/res)
                        ax[ii, jj].set_ylim(0, ff[0].header['FIELDOFVIEW_PKPC']/res)

                        ax[ii, jj].xaxis.set_label_coords(.5, -.08)
                        ax[ii, jj].yaxis.set_label_coords(-0.08, 0.5)

                        if jj == 0:
                            ax[ii, jj].set_ylabel('pkpc')

                        ff.close()

                    ax[-1, jj].set_xlabel('pkpc')

            ## Plot second set of telescopes, all have 10x10kpc FOV ##
            filterList = [('euclid', 'nisp', 'y'), ('vista', 'vircam', 'y'), ('subaru', 'hsc', 'y')]

            # Plot for one galaxy
            if self.nGalsSelected == 1:
                fig2, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                        figsize=(len(filterList)*2, self.nGalsSelected*2),
                                        gridspec_kw={'hspace': 0.1, 'wspace': 0.1})

                for jj in range(0, len(filterList)):
                    ax[jj].set_xlabel('{} {} {}'.format(filterList[jj][0].upper(), filterList[jj][1].upper(),
                                      filterList[jj][2].upper()))
                    if self.z == 7:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                       format(filterList[jj][0], filterList[jj][1], self.z,
                                              self.catalogueSelected['fileNumber'].values[0], filterList[jj][2]))
                    else:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                       format(filterList[jj][0], filterList[jj][1], self.z, filterList[jj][2]))

                    im = ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data, cmap='Greys')
                    if colorbar is True:
                        #plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
                        divider = make_axes_locatable(ax[jj])
                        cax = divider.append_axes('right', size='5%', pad=0.05)
                        cbar = plt.colorbar(im, cax=cax)
                        cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
                        cbar.ax.set_title('nJy', loc='left', pad=-5)

                    res = ff[0].header['RESOLUTION_PKPC']
                    ax[jj].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xlabel('pkpc')
                    ff.close()

            # Plot for multiple galaxies
            else:
                if colorbar is True:
                    fig2, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                            figsize=(len(filterList)*1.8, self.nGalsSelected*1.5),
                                            gridspec_kw={'hspace': 0.15, 'wspace': 0.35})
                else:
                    fig2, ax = plt.subplots(self.nGalsSelected, len(filterList),
                                            figsize=(len(filterList)*1.5, self.nGalsSelected*1.5),
                                            gridspec_kw={'hspace': 0.1, 'wspace': 0.1})

                for jj in range(0, len(filterList)):
                    ax[0, jj].set_title('{}\n{} {}'.format(filterList[jj][0].upper(), filterList[jj][1].upper(),
                                                           filterList[jj][2].upper()), pad=10)
                    for ii in range(0, self.nGalsSelected):
                        if self.z == 7:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                           format(filterList[jj][0], filterList[jj][1], self.z,
                                                  self.catalogueSelected['fileNumber'].values[ii], filterList[jj][2]))
                        else:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                           format(filterList[jj][0], filterList[jj][1], self.z, filterList[jj][2]))
                        im = ax[ii, jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,
                                               cmap='Greys')
                        if colorbar is True:
                            #plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
                            divider = make_axes_locatable(ax[ii, jj])
                            cax = divider.append_axes('right', size='5%', pad=0.05)
                            cbar = plt.colorbar(im, cax=cax)
                            cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
                            cbar.ax.set_title('nJy', loc='left', pad=-5)
                        res = ff[0].header['RESOLUTION_PKPC']

                        if jj == 0:
                            ax[ii, jj].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii, jj].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            ax[ii, jj].set_ylabel('pkpc')
                        else:
                            ax[ii, jj].set_yticks([])

                        if ii < self.nGalsSelected-1:
                            ax[ii, jj].set_xticks([])
                        else:
                            ax[ii, jj].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii, jj].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii, jj].set_xlim(0, ff[0].header['FIELDOFVIEW_PKPC']/res)
                        ax[ii, jj].set_ylim(0, ff[0].header['FIELDOFVIEW_PKPC']/res)

                        ax[ii, jj].xaxis.set_label_coords(.5, -.08)
                        ax[ii, jj].yaxis.set_label_coords(-0.08, 0.5)

                        if jj == 0:
                            ax[ii, jj].set_ylabel('pkpc')

                        ff.close()

                    ax[-1, jj].set_xlabel('pkpc')

            return fig, fig2

        ##### All telescopes with all filters #####
        elif telescope == 'all' and filt == 'all':
            self.plot_galaxies('jwst', 'nircam', filt='all')
            self.plot_galaxies('jwst', 'miri', filt='all')
            for telescope in ['hst', 'euclid', 'roman', 'vista', 'subaru']:
                self.plot_galaxies(telescope, instrument, filt='all')
            return

        ##### A single telescope #####
        else:
            # Check that telescope is found.
            if telescope not in ['hst', 'jwst', 'euclid', 'roman', 'vista', 'subaru']:
                errorstr = "ERROR: Telescope '{}' doesn't match those available ".format(telescope) +\
                           "('hst', 'jwst', 'euclid', 'roman', 'vista', 'subaru')"
                print(errorstr)
                return errorstr

            # Check selected an available instrument
            if telescope == 'euclid':
                if instrument != 'nisp':
                    print('WARNING: Specified euclid instrument "{}" is not nisp, converted to nisp'.format(instrument))
                    instrument = 'nisp'

            elif telescope == 'hst':
                if instrument != 'wfc3':
                    print('WARNING: Specified hst instrument "{}" is not wfc3, converted to wfc3'.format(instrument))
                    instrument = 'wfc3'

            elif telescope == 'roman':
                if instrument != 'wfi':
                    print('WARNING: Specified roman instrument "{}" is not wfi, converted to wfi'.format(instrument))
                    instrument = 'wfi'

            elif telescope == 'vista':
                if instrument != 'vircam':
                    print('WARNING: Specified vista instrument "{}" is not vircam, converted to vircam'.
                          format(instrument))
                    instrument = 'vircam'

            elif telescope == 'subaru':
                if instrument != 'hsc':
                    print('WARNING: Specified subaru instrument "{}" is not hsc, converted to hsc'.format(instrument))
                    instrument = 'hsc'

            elif telescope == 'jwst':
                if instrument not in ['nircam', 'miri', 'all']:
                    errorstr = 'ERROR: Specified jwst instrument "{}" is not nircam, miri, or all, aborted'.\
                               format(instrument)
                    print(errorstr)
                    return errorstr

            filterList = {'nisp': ['y', 'j', 'h'],
                          'wfc3': ['f105w', 'f125w', 'f140w', 'f160w'],
                          'nircam': ['f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w'],
                          'miri': ['f560w', 'f770w'],
                          'all': ['f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w', 'f560w', 'f770w'],
                          'wfi': ['f087', 'f106', 'f129', 'f146', 'f158', 'f184'],
                          'hsc': ['z', 'y'],
                          'vircam': ['z', 'y', 'j', 'h', 'ks']}

            ## All available filters ##
            if filt == 'all':

                # Plot for one galaxy
                if self.nGalsSelected == 1:
                    fig, ax = plt.subplots(self.nGalsSelected, len(filterList[instrument]),
                                           figsize=(len(filterList[instrument])*1.8, self.nGalsSelected*1.5),
                                           gridspec_kw={'bottom': 0.05, 'top': 0.95})
                    plt.suptitle('{} {}'.format(telescope.upper(), instrument.upper()))

                    for jj, filt in enumerate(filterList[instrument]):
                        ax[jj].set_xlabel(filt.upper())

                        if self.z == 7:
                            if telescope == 'jwst' and instrument == 'all':
                                if jj < 7:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                                   format(telescope, 'nircam', self.z,
                                                          self.catalogueSelected['fileNumber'].values[0], filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                                   format(telescope, 'miri', self.z,
                                                          self.catalogueSelected['fileNumber'].values[0], filt))
                            else:
                                ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                               format(telescope, instrument, self.z,
                                                      self.catalogueSelected['fileNumber'].values[0], filt))

                        else:
                            if telescope == 'jwst' and instrument == 'all':
                                if jj < 7:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                                   format(telescope, 'nircam', self.z, filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                                   format(telescope, 'miri', self.z, filt))
                            else:
                                ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                               format(telescope, instrument, self.z, filt))

                        im = ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data, cmap='Greys')
                        if colorbar is True:
                            #plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
                            divider = make_axes_locatable(ax[jj])
                            cax = divider.append_axes('right', size='5%', pad=0.05)
                            cbar = plt.colorbar(im, cax=cax)
                            cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
                            cbar.ax.set_title('nJy', loc='left', pad=-5)
                        res = ff[0].header['RESOLUTION_PKPC']
                        ax[jj].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj].set_xlabel('pkpc')
                        ff.close()

                    plt.suptitle('{} {}'.format(telescope.upper(), instrument.upper()))

                # Plot for multiple galaxies
                else:
                    if colorbar is True:
                        fig, ax = plt.subplots(self.nGalsSelected, len(filterList[instrument]),
                                               figsize=(len(filterList[instrument])*1.5, self.nGalsSelected*1.3),
                                               gridspec_kw={'wspace': 0.32, 'bottom': 0.05, 'top': 0.95, 'left': 0.04,
                                                            'right': 0.96})
                    else:
                        fig, ax = plt.subplots(self.nGalsSelected, len(filterList[instrument]),
                                               figsize=(len(filterList[instrument])*1.5, self.nGalsSelected*1.5),
                                               gridspec_kw={'hspace': 0.1, 'wspace': 0.1})

                    for jj, filt in enumerate(filterList[instrument]):
                        ax[0, jj].set_title(filt.upper())
                        for ii in range(0, self.nGalsSelected):
                            if self.z == 7:
                                if telescope == 'jwst' and instrument == 'all':
                                    if jj < 7:
                                        ff = fits.open(self.image_path +
                                                       'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                                       format(telescope, 'nircam', self.z,
                                                              self.catalogueSelected['fileNumber'].values[ii], filt))
                                    else:
                                        ff = fits.open(self.image_path +
                                                       'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                                       format(telescope, 'miri', self.z,
                                                              self.catalogueSelected['fileNumber'].values[ii], filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                                   format(telescope, instrument, self.z,
                                                          self.catalogueSelected['fileNumber'].values[ii], filt))

                            else:
                                if telescope == 'jwst' and instrument == 'all':
                                    if jj < 7:
                                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                                       format(telescope, 'nircam', self.z, filt))
                                    else:
                                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                                       format(telescope, 'miri', self.z, filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                                   format(telescope, instrument, self.z, filt))

                            im = ax[ii, jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,
                                                   cmap='Greys')
                            if colorbar is True:
                                #plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
                                divider = make_axes_locatable(ax[ii, jj])
                                cax = divider.append_axes('right', size='5%', pad=0.05)
                                cbar = plt.colorbar(im, cax=cax)
                                cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4, steps=[1, 2, 3, 4, 5]))
                                cbar.ax.set_title('nJy', loc='left', pad=-7)
                            res = ff[0].header['RESOLUTION_PKPC']

                            if jj == 0:
                                ax[ii, jj].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                                ax[ii, jj].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                                ax[ii, jj].set_ylabel('pkpc')
                            else:
                                ax[ii, jj].set_yticks([])

                            if ii < self.nGalsSelected-1:
                                ax[ii, jj].set_xticks([])
                            else:
                                ax[ii, jj].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                                ax[ii, jj].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])

                            ax[ii, jj].set_xlim(0, ff[0].header['FIELDOFVIEW_PKPC']/res)
                            ax[ii, jj].set_ylim(0, ff[0].header['FIELDOFVIEW_PKPC']/res)

                            ax[ii, jj].xaxis.set_label_coords(.5, -.08)
                            ax[ii, jj].yaxis.set_label_coords(-0.06, 0.5)

                            ff.close()

                        ax[-1, jj].set_xlabel('pkpc')

            ## Single filter ##
            else:

                if filt not in filterList[instrument]:
                    errorstr = 'ERROR: Specified filter "{}" is not in filter list, aborted'.format(filt)
                    print(errorstr)
                    return errorstr

                else:
                    if colorbar is True:
                        fig, ax = plt.subplots(1, self.nGalsSelected, figsize=(self.nGalsSelected*2, 2),
                                               gridspec_kw={'wspace': 0.35})
                    else:
                        fig, ax = plt.subplots(1, self.nGalsSelected, figsize=(self.nGalsSelected*2, 2))

                    plt.suptitle('{} {} {}'.format(telescope.upper(), instrument.upper(), filt.upper()))
                    for ii in range(0, self.nGalsSelected):
                        if self.z == 7:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.
                                           format(telescope, instrument, self.z,
                                                  self.catalogueSelected['fileNumber'].values[ii], filt))
                        else:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.
                                           format(telescope, instrument, self.z, filt))
                        if self.nGalsSelected == 1:
                            im = ax.imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data, cmap='Greys')
                        else:
                            im = ax[ii].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,
                                               cmap='Greys')
                        if colorbar is True:
                            #plt.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)
                            divider = make_axes_locatable(ax[ii])
                            cax = divider.append_axes('right', size='5%', pad=0.05)
                            cbar = plt.colorbar(im, cax=cax)
                            cax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
                            cbar.ax.set_title('nJy', loc='left', pad=-5)

                        res = ff[0].header['RESOLUTION_PKPC']
                        ax[ii].set_yticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[ii].set_yticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii].set_xticks([0, ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[ii].set_xticklabels(['0', str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii].set_xlabel('pkpc')
                        ax[ii].xaxis.set_label_coords(.5, -.08)
                        ff.close()
            return fig


def lum_to_mag(L, z):
    """ Converts luminosity to magnitude, used particularly for the rest-frame FUV 1500A band (Fontanot+12, see Ni+19)

    Inputs:
    L -- Float: Luminosity (erg/s/Hz)
    z -- Int: Redshift, either 7, 8, 9, 10, 11, or 12.

    Returns:
    mag -- Flaot: Apparent magnitude (AB mag)
    """
    # Values taken from Gnedin's cosmological calculator for BlueTides cosmology, assuming Omega0=0.3278,
    # H0=69.7km/s/Mpc
    if z == 7:
        dlMpc = 70878.29  # Luminosity distance
        dm = 49.25  # Distance modulus
        dmK = 46.99  # Distance modulus + k-correction
    elif z == 8:
        dlMpc = 82688.79
        dm = 49.59
        dmK = 47.2
    elif z == 9:
        dlMpc = 94637.09
        dm = 49.88
        dmK = 47.38
    elif z == 10:
        dlMpc = 106715.22
        dm = 50.14
        dmK = 47.54
    elif z == 11:
        dlMpc = 118912.91
        dm = 50.38
        dmK = 47.68
    elif z == 12:
        dlMpc = 131198
        dm = 50.59
        dmK = dm-2.78
    else:
        errorstr = 'ERROR: Specified redshift "{}" is not in [7, 8, 9, 10, 11, 12], aborted'.format(z)
        print(errorstr)
        return errorstr

    dl = (dlMpc*1e6*3.086e+18)  # cm
    mag = -2.5 * np.log10(L/(4*np.pi*dl**2)) - 48.60-dm
    # No k-correction required for LUV -> MUV
    return mag


def mag_to_lum(mag, z):
    """ Converts magnitude to luminosity, used particularly for the rest-frame FUV 1500A band (Fontanot+12, see Ni+19)

    Inputs:
    mag -- Float: Apparent magnitude (AB mag)
    z -- Int: Redshift, either 7, 8, 9, 10, 11, or 12.

    Returns:
    L -- Float: Luminosity (erg/s/Hz)
    """
    # Values taken from Gnedin's cosmological calculator for BlueTides cosmology, assuming Omega0=0.3278,
    # H0=69.7km/s/Mpc
    if z == 7:
        dlMpc = 70878.29  # Luminosity distance
        dm = 49.25  # Distance modulus
        dmK = 46.99  # Distance modulus + k-correction
    elif z == 8:
        dlMpc = 82688.79
        dm = 49.59
        dmK = 47.2
    elif z == 9:
        dlMpc = 94637.09
        dm = 49.88
        dmK = 47.38
    elif z == 10:
        dlMpc = 106715.22
        dm = 50.14
        dmK = 47.54
    elif z == 11:
        dlMpc = 118912.91
        dm = 50.38
        dmK = 47.68
    elif z == 12:
        dlMpc = 131198
        dm = 50.59
        dmK = dm-2.78
    else:
        errorstr = 'ERROR: Specified redshift "{}" is not in [7, 8, 9, 10, 11, 12], aborted'.format(z)
        print(errorstr)
        return errorstr

    dl = (dlMpc*1e6*3.086e+18)  # cm
    L = 10**(-0.4*(mag+dm+48.60))*(4*np.pi*dl**2)
    # No k-correction required for LUV -> MUV
    return L
