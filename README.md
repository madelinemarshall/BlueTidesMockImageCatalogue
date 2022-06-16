# Code in this Directory
Here we provide code to help access and visualise the galaxy images contained in the BlueTides Mock Image Catalogue.
- bluetides_mock_catalogue_functions.py: helper code to read, access and select galaxies within the catalogue, and plot their images
- plot_selected_galaxies.ipynb: example notebook walking through how to use the functions in bluetides_mock_catalogue_functions.py. Recreates Figures 2 and 3 from the Mock Catalogue Release paper.
- test_plot_galaxies.py: pytest file for testing that bluetides_mock_catalogue_functions.py works as intended.



# BlueTides Mock Image Catalogue

### Details
HLSP Contributor: Madeline Marshall, NRC Herzberg Astronomy & Astrophysics
(Madeline.Marshall at nrc-cnrc.gc.ca) <br>
HLSP Authors: Katelyn Watts, Stephen Wilkins, Tiziana Di Matteo, Jussi K. Kuusisto, William J. Roper, Aswin P. Vijayan, 
Yueying Ni, Yu Feng, Rupert A.C. Croft <br>
Date: TBD <br>
Version: 1 <br>

### Overview
This catalogue contains mock images of ~100,000 galaxies from the BlueTides
hydrodynamical simulation at z=7,8,9,10,11 and 12. This includes mock images 
of these galaxies with the James Webb, Hubble, Roman, and Euclid Space 
Telescopes, as well as Subaru and VISTA.

This v1.0 release of the BlueTides Mock Image Catalogue contains:
- 540 fits image files, containing images of all of the galaxies for each of the redshift, telescope, and PSF-option combinations
- 1 csv table containing assorted physical and observed properties of each BlueTides galaxy in the catalogue
- This "README" file

The primary catalogue files have names of the form:
```
hlsp_bluetides_TELESCOPE_INSTRUMENT_zREDSHIFT_FILTER_vVERSION_sim-[PSF,NOPSF].fits
```
as discussed in detail below. The file format is FITS with multiple image extensions.

A table showing the available telescope configurations, alongside key image properties, 
is shown below.

| Telescope | Instrument | Filters | Pixel Scale ('') | Image FOV (pkpc)  |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| JWST | NIRCam (SW) | F090W, F115W, F150W, F200W | 0.0155 | 6  |
| | NIRCam (LW) | F277W, F356W, F410M, F444W | 0.0315 | 6  |
| | MIRI | F560W, F770W | 0.055 | 6  |
| HST | WFC3 | F105W, F125W, F140W, F160W | 0.065 | 6  |
| Roman | WFI | F087, F106, F129, F146, F158, F184 | 0.055 | 6  |
| Euclid | NISP | Y, J, H | 0.15 | 10  |
| VISTA | VIRCam | Z, Y, J, H, Ks | 0.17 | 10  |
| Subaru | HSC | z, y | 0.085 | 10  |


### Catalogue Creation

These observations are created with the [SynthObs](https://github.com/stephenmwilkins/SynthObs) package for producing synthetic
observations. The image creation
is described in detail in Marshall et al. (submitted).
Briefly, each star particle in a galaxy in BlueTides is assigned a spectral energy
distribution (SED) using BPASS (version 2.2.1; Stanway & Eldridge
2018) and CLOUDY (Ferland et al. 2017). To create the images the SED from each
star particle in the FOV is convolved with the instrument filter transmission curve.
The images of the galaxies are made in the 'face-on' direction, defined by the
angular momentum of the galaxy particles.
The light from each star particle is smoothed adaptively based on the distance
to its nearest neighbours.

We consider the FOV of 6x6 kpc around each galaxy, except for Euclid, Subaru, and VISTA,
which have much wider PSFs and so require a larger FOV (10x10 kpc) to contain
the galaxy emission. The images are binned onto a pixel scale of 0.5 times the
native pixel scale. This assumes that given sufficient dithering in the
observations, the final image can sub-sample the original pixels by a factor of 2.
The resulting image is convolved with the PSF of the instrument. The PSFs are 
obtained via TinyTim for HST, and WebbPSF for JWST and Roman. Euclid NISP is 
assumed to have a Gaussian PSF with FWHM of 0.175'' in the Y-band, 0.24'' in J, 
and 0.28'' in H. Ground-based Subaru and VISTA are assumed to have Gaussian PSFs 
with FWHM of 0.6'' and 0.66'' respectively, corresponding to the typical site 
seeing. These images are marked as sim-psf in the collection. We also provide 
images that have _not_ been convolved with a model PSF, marked as sim-nopsf.

The catalogue images are made without noise, which can be added in post-processing.
Code to add noise is available via [GitHub](https://github.com/madelinemarshall/BlueTidesMockImageCatalogue).

### Catalogue & Image Format

#### Available fits files:

Each fits file contains all of the galaxies at the specified redshift for the defined
telescope/instrument/filter combination.

| Redshift | Number of Galaxies |
| ------------- | ------------- |
| z=7 |  71,052* |
| z=8 | 22,144 |
| z=9 | 5,606 |
| z=10 | 1,279 |
| z=11 | 244 |
| z=12 | 31 |

*Note that to reduce the large file sizes for the z=7 snapshot, we split the z=7
catalogue into 4 separate files, containing 17764, 17763, 17763 and 17762 galaxies. 
These files are marked as z7-file1, z7-file2, z7-file3 and z7-file4 in the collection.

The filenames show the telescope, instrument, filter, and redshift:
```
hlsp_bluetides_TELESCOPE_INSTRUMENT_TARGET_FILTER_vVERSION_sim-[PSF,NOPSF].fits
```
where TARGET = zREDSHIFT if z>7, and TARGET = z7-fileN for N in 1,2,3,4 for z=7.

Example:
PSF convolution:
```
hlsp_bluetides_jwst_nircam_z7-file1_f150w_v1_sim-psf.fits
hlsp_bluetides_jwst_nircam_z8_f150w_v1_sim-psf.fits
```
No PSF convolution:
```
hlsp_bluetides_jwst_nircam_z7-file1_f150w_v1_sim-nopsf.fits
hlsp_bluetides_jwst_nircam_z8_f150w_v1_sim-nopsf.fits
```

These images contain a header file, and then a separate image extension for each galaxy.


##### Header file:

- Provides key information about the simulated images in this file
- Contains a dummy array of 1s, with the same shape as the galaxy images.

Example header file:

```
SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                  -64 / number of bits per data pixel
NAXIS   =                    2 / number of data axes
NAXIS1  =                   46 / length of data axis 1
NAXIS2  =                   46 / length of data axis 2
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
HIERARCH TELESCOPE = 'JWST    '              * The simulated telescope
HIERARCH INSTRUMENT = 'NIRCAM  '             * The simulated INSTRUMENT
HIERARCH FILTER = 'F410M   '                 * The simulated FILTER
HIERARCH FIELDOFVIEW_PKPC =  6               * The width and height of the image in pkpc
HIERARCH FIELDOFVIEW_ARCSEC = 1.50107        * The width and height of the image in arcsec (assuming BlueTides cosmology)
HIERARCH RESOLUTION_PKPC = 0.12591028        * The pixel scale of the image in pkpc
HIERARCH RESOLUTION_ARCSEC = 0.0315          * The pixel scale of the image in arcsec
HIERARCH SUPERSAMPLINGRATE = 2               * The sub-sampling rate of the native pixel scale
HIERARCH REDSHIFT =          7               * The redshift of the snapshot
DOI     = '10.17909/er09-4527'               * Digital Object Identifier for the HLSP data collection                       
HLSPID  = 'BlueTides'                        * The identifier (acronym) for this HLSP collection                         
HLSPNAME= 'BlueTides Mock Image Catalogue'   * Title for HLSP project, long form
SIMULATD=                    T               * This is simulated data                                   
HLSPLEAD= 'Madeline Marshall'                * Full name of HLSP project lead                                   
HLSPVER =                    1               * Version identifier for this HLSP product                                   
LICENSE = 'CC BY 4.0'                        * License for use of these data                                   
LICENURL= 'https://creativecommons.org/licenses/by/4.0/'       * Data license URL                 
HLSPTARG= 'z7-file1'                         * zREDSHIFT if z>7, and z7-fileN for N=1, 2, 3, or 4 if z=7.                                 
OBSERVAT= 'hst     '                         * Simulated observatory for the mock images                                 
TELESCOP= 'hst     '                         * Simulated telescope for the mock images                                             
INSTRUME= 'wfc3    '                         * Simulated instrument for the mock images                                             
FILTER  = 'f160w   '                         * Name of filter used for the mock images                                   
PRODTYPE= 'sim-psf '                         
* Type of data. sim = simulated data. psf = convolved with model PSF. nopsf =  not convolved with a PSF.                    
BUNIT   = 'nJy     '                         * Brightness unit for array values
END
```

##### Image extensions

- Each extension contains an image of a galaxy in BlueTides at that redshift <br>
- The flux is in units of nJy <br>
- The header gives identifiable and useful information about the galaxy <br>

Example image header for the 9th fits extension (the 8th image extension):

```
XTENSION= 'IMAGE   '           / IMAGE extension
BITPIX  =                  -64 / number of bits per data pixel
NAXIS   =                    2 / number of data axes
NAXIS1  =                   46 / length of data axis 1
NAXIS2  =                   46 / length of data axis 2
PCOUNT  =                    0 / required keyword; must = 0
GCOUNT  =                    1 / required keyword; must = 1
HIERARCH GALAXYINDEX =       7                  * The 8th galaxy (indexes from 0)
HIERARCH MASS_SOLARMASSES = '1.177723e+09'      * The solar mass measured from BlueTides, in solar masses
HIERARCH FLUX_NJY = '41.73964'                  * The total flux in this image, in nJy
END
```

#### Accompanying catalogue csv file:

Alongside the fits files containing the galaxy images, we provide a table of galaxy properties for use when accessing galaxies within the catalogue, 'hlsp_bluetides_multi_multi_all_multi_v1_sim.csv'. This one table contains all of the relevant properties for the galaxies at each redshift, and in each telescope/instrument/filter combination.

The table columns are:
```
redshift                    * The redshift snapshot that the galaxy is in
fileNumber                  * If z=7, gives the file number (1,2,3,4) in which the galaxy is located. If z>7, this is set to 0 as this is redundant.
extensionNumber             * The fits file extension which corresponds to this galaxy in the relevant fits file. For z>7, this is 1 greater than the index, as the first fits extension is the header. This is particularly important for z=7 where the galaxies are split into 4 files.
stellarMass                 * Galaxy stellar mass in solar masses 
BHmass                      * Central black hole mass in solar masses 
BHluminosity                * Intrinsic AGN luminosity, calculated from the black hole accretion rate, in  erg/s
haloMass                    * Friends-of-friends halo mass in solar masses
lum_FUV                     * FUV (1500 Angstrom) luminosity in erg/s/Hz 
flux_jwst.nircam.f090w      * Total flux in the image, for the specified telescope.instrument.filter (all lower-case), in nJy
radius_jwst.nircam.f090w    * Effective half-light radius of the galaxy as measured from the image, for the specified telescope.instrument.filter (all lower-case), in pkpc
...                         * [all available telescope.instrument.filter combinations]
flux_spitzer.irac.ch1       * Total flux in a 10x10 kpc FOV image with Spitzer IRCA Channel 1, in nJy (these images are not provided in the catalogue as all galaxies are unresolved)
flux_spitzer.irac.ch2       * Total flux in a 10x10 kpc FOV image with Spitzer IRCA Channel 2, in nJy (these images are not provided in the catalogue as all galaxies are unresolved)
```

For details on the calculation of these properties, see the Mock Catalogue Release paper as well as [Marshall et al. 2022](https://doi.org/10.1093/mnras/stac380) particularly for the galaxy size calculation.

Example code using this table to select galaxies from the catalogue fits files is available via [GitHub](https://github.com/madelinemarshall/BlueTidesMockImageCatalogue).
