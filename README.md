# BlueTides Mock Image Catalogue

### Details
Contributor: Madeline Marshall, NRC Herzberg Astronomy & Astrophysics
(Madeline.Marshall@nrc-cnrc.gc.ca) <br>
Date: TBD <br>
Version: 1 <br>

### Overview
This catalogue contains mock images of ~100,000 galaxies from the BlueTides
hydrodynamical simulation at z=7,8,9,10,11 and 12. This includes mock images 
of these galaxies with the James Webb, Hubble, Roman, and Euclid Space 
Telescopes, as well as Subaru and VISTA.

This v1.0 release of the BlueTides Mock Image Catalogue contains:
- 540 image files, one for each of the telescope, redshift, and PSF-option combinations
- This "README" file

The primary catalog files have names of the form:
```
hlsp_bluetides_TELESCOPE_INSTRUMENT_zREDSHIFT_FILTER_vVERSION_sim-[PSF,NOPSF].fits
```
as discussed in detail below. The file format is FITS with multiple image extensions.

A table showing the available telescope configurations, alongside key image properties, is shown below.

| Telescope | Instrument | Filters | Pixel Scale ('') | Image FOV (pkpc)  |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| JWST | NIRCam (SW) | F090W, F115W, F150W, F200W | 0.0155 | 6  |
| | NIRCam (LW) | F277W, F356W, F410M, F444W | 0.0315 | 6  |
| | MIRI | F560W, F770W | 0.055 | 6  |
| HST | WFC3 | F105W, F125W, F140W, F160W | 0.065 | 6  |
| Euclid | NISP | Y, J, H | 0.15 | 10  |
| Roman | WFI | F087, F106, F129, F146, F158, F184 | 0.055 | 6  |
| VISTA | VIRCam | Z, Y, J, H, Ks | 0.17 | 10  |
| Subaru | HSC | z, y | 0.085 | 6  |


### Catalogue Creation

These observations are created with the SynthObs package for producing synthetic
observations (https://github.com/stephenmwilkins/SynthObs). The image creation
is described in detail in Marshall et al. (in prep).
Briefly, each star particle in a galaxy in BlueTides is assigned a spectral energy
distribution (SED) using BPASS (version 2.2.1; Stanway & Eldridge
2018) and CLOUDY (Ferland et al. 2017). To create the images the SED from each
star particle in the FOV is convolved with the instrument filter transmission curve.
The images of the galaxies are made in the 'face-on' direction, defined by the
angular momentum of the galaxy particles.
The light from each star particle is smoothed adaptively based on the distance
to its nearest neighbours.

We consider the FOV of 6x6 kpc around each galaxy, except for Euclid and VISTA,
which have much wider PSFs and so require a larger FOV (10x10 kpc) to contain
the galaxy emission. The images are binned onto a pixel scale of 0.5 times the
native pixel scale. This assumes that given sufficient dithering in the
observations, the final image can sub-sample the original pixels by a factor of 2.
The resulting image is convolved with the PSF of the instrument, obtained
via TinyTim for HST, WebbPSF for JWST and Roman, **Euclid from wherever Steve got that from??**,
and ground-based Subaru and VISTA from a simple Gaussian with FWHM of 0.6'' and
0.66'' respectively corresponding to the typical site seeing. These images are
marked as sim-psf in the collection. We also provide images that have _not_ been
convolved with a model PSF, marked as sim-nopsf.

The catalogue images are made without noise, which can be added in post-processing.
Code to add noise is available at https://github.com/madelinemarshall/BlueTidesMockImageCatalogue

### Catalogue & Image Format

#### Available files:

Each file contains all of the galaxies at the specified redshift for the defined
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
catalogue into 4 separate files, each containing 17,763 galaxies. These files are
marked as z7-file1, z7-file2, z7-file3 and z7-file4 in the collection.

The filenames show the telescope, instrument, filter, and redshift:
```
hlsp_bluetides_TELESCOPE_INSTRUMENT_zREDSHIFT_FILTER_vVERSION_sim-[PSF,NOPSF].fits
```

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


#### Header file:

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
END
```

#### Image extensions

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

Example code for using the catalogue and and plotting these images is given at
https://github.com/madelinemarshall/BlueTidesMockImageCatalogue
