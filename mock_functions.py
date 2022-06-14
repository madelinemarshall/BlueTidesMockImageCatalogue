import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from astropy.io import fits
import pickle

import matplotlib
matplotlib.rcParams['font.size'] = (9)
matplotlib.rcParams['figure.figsize'] = (4,3.4)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

class Catalogue:

    def __init__(self,z=7,image_path='/home/mmarshal/FinalImages/FullCatalogue/'):
        # Load in catalogue and save columns as variables for ease of access
        cat = pd.read_csv(image_path+'BlueTidesMockImageCatalogue.csv',index_col=0)
        self.image_path = image_path
        self.z = z
        cat=cat.loc[(cat["redshift"]==z)]
        self.catalogue = cat
        self.stellarMass = np.array(self.catalogue['stellarMass'])
        self.BHmass = np.array(self.catalogue['BHmass'])
        self.haloMass = np.array(self.catalogue['haloMass'])
        self.nGals = len(self.stellarMass)
        return


    def print_constraint_options(self):
        cols = self.catalogue.columns.values
        for col in cols:
            width = 25
            print('{: <{}}:  Min = {:0.2e},  Max = {:0.2e}'.format(col,width,np.min(self.catalogue[col]),np.max(self.catalogue[col])))
        return


    def select_galaxies_with_constraints(self,stellar_mass_constraint = None,bh_mass_constraint = None,halo_mass_constraint = None,
                  lum_constraint = None,rad_constraint = None):
        # All in log space except rad_FUV
    
    
        # Check which galaxies match the specified constraints
        if stellar_mass_constraint:
            stellar_mass_constraint = (self.stellarMass > 10**stellar_mass_constraint[0]) & (self.stellarMass < 10**stellar_mass_constraint[1])
        else:
            stellar_mass_constraint = np.ones(self.nGals,dtype='Bool')
        
      
        if bh_mass_constraint:
            bh_mass_constraint = (self.BHmass > 10**bh_mass_constraint[0]) & (self.BHmass < 10**bh_mass_constraint[1])
        else:
            bh_mass_constraint = np.ones(self.nGals,dtype='Bool')
        
        
        if halo_mass_constraint:
            halo_mass_constraint = (self.haloMass > 10**halo_mass_constraint[0]) & (self.haloMass < 10**halo_mass_constraint[1])
        else:
            halo_mass_constraint = np.ones(self.nGals,dtype='Bool')
       
        
        if lum_constraint:
            nFilters = len(lum_constraint)
            lum_constraint_allfilt = np.ones(self.nGals,dtype='Bool')
            for ff in range(0,nFilters):
                filt = lum_constraint[ff][0]
                lum_constraint_ff = (self.catalogue['lum_'+filt.lower()] > 10**lum_constraint[ff][1]) & (self.catalogue['lum_'+filt.lower()] < 10**lum_constraint[ff][2])
                lum_constraint_allfilt = np.logical_and(lum_constraint_ff,lum_constraint_allfilt)
            lum_constraint = lum_constraint_allfilt
        else:
            lum_constraint = np.ones(self.nGals,dtype='Bool')
      
    
        if rad_constraint:
            rad_constraint_allfilt = np.ones(self.nGals,dtype='Bool')
            for ff in range(0,nFilters):
                filt = rad_constraint[ff][0]
                rad_constraint_ff = (self.catalogue['rad_'+filt.lower()] > rad_constraint[ff][1]) & (self.catalogue['rad_'+filt.lower()] < rad_constraint[ff][2])
                rad_constraint_allfilt = np.logical_and(rad_constraint_ff,rad_constraint_allfilt)
            rad_constraint = rad_constraint_allfilt
        else:
            rad_constraint = np.ones(self.nGals,dtype='Bool')
          
        selection = stellar_mass_constraint & bh_mass_constraint & halo_mass_constraint & \
            lum_constraint & rad_constraint
    
        if len(selection[selection])==0:
            print('Error: No galaxies satisfy these criteria')
   
        else: 
            self.stellarMassSelected = np.array(self.catalogue['stellarMass'][selection])
            self.BHmassSelected = np.array(self.catalogue['BHmass'][selection])
            self.haloMassSelected = np.array(self.catalogue['haloMass'][selection])
            self.nGalsSelected = len(self.stellarMass[selection])
            self.catalogueSelected = self.catalogue[selection]
        return 
    
    
    def select_galaxies_with_indices(self,indices):
        full_inds = self.catalogue.index
        selection = np.in1d(full_inds,indices)

        self.stellarMassSelected = np.array(self.catalogue['stellarMass'][selection])
        self.BHmassSelected = np.array(self.catalogue['BHmass'][selection])
        self.haloMassSelected = np.array(self.catalogue['haloMass'][selection])
        self.nGalsSelected = len(self.stellarMass[selection])
        self.catalogueSelected = self.catalogue[selection]
        return 



    def plot_galaxies(self,telescope,instrument,filt):
        # Telescope, instrument and filters all lowercase
        telescope=telescope.lower()
        instrument=instrument.lower()
        filt=filt.lower()
        
        ##### All telescopes with Y band #####
        if telescope == 'all' and filt == 'y':
            ## Plot first set of telescopes, all have 6x6kpc FOV ##
            filterList = [('jwst','nircam','f115w'),('hst','wfc3','f105w'),('roman','wfi','f106')]
        
            # Plot for one galaxy
            if self.nGalsSelected == 1:
                fig,ax = plt.subplots(self.nGalsSelected,len(filterList),
                    figsize = (len(filterList)*2,self.nGalsSelected*2),
                    gridspec_kw = {'hspace':0.1,'wspace':0.1})

                for jj in range(0,len(filterList)):
                    ax[jj].set_xlabel('{} {} {}'.format(filterList[jj][0].upper(),filterList[jj][1].upper(),filterList[jj][2].upper()))
                    if self.z==7:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,self.catalogueSelected['fileNumber'].values[0],filterList[jj][2]))
                    else:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,filterList[jj][2]))

                    ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data,cmap='Greys')
                    res=ff[0].header['RESOLUTION_PKPC']
                    ax[jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xlabel('pkpc')
                    ff.close()
                #plt.suptitle('Y-band images')
                
            # Plot for multiple galaxies
            else:
                fig,ax = plt.subplots(self.nGalsSelected,len(filterList),
                            figsize = (len(filterList)*1.5,self.nGalsSelected*1.5),
                            gridspec_kw = {'hspace':0.1,'wspace':0.1})

                for jj in range(0,len(filterList)):
                    ax[0,jj].set_title('{}\n{} {}'.format(filterList[jj][0].upper(),filterList[jj][1].upper(),filterList[jj][2].upper()))
                    for ii in range(0,self.nGalsSelected):
                        if self.z==7:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,self.catalogueSelected['fileNumber'].values[ii],filterList[jj][2]))
                        else:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,filterList[jj][2]))
                        ax[ii,jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                        res=ff[0].header['RESOLUTION_PKPC']

                        
                        if jj==0:
                            ax[ii,jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii,jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            ax[ii,jj].set_ylabel('pkpc')
                        else:
                            ax[ii,jj].set_yticks([])

                        if ii<self.nGalsSelected-1:
                            ax[ii,jj].set_xticks([])
                        else:
                            ax[ii,jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii,jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii,jj].set_xlim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                        ax[ii,jj].set_ylim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                        
                        ax[ii,jj].xaxis.set_label_coords(.5, -.08)
                        ax[ii,jj].yaxis.set_label_coords(-0.08, 0.5) 
                        
                        if jj==0:
                            ax[ii,jj].set_ylabel('pkpc')
    
                        ff.close()
        
                    ax[-1,jj].set_xlabel('pkpc')
            
            
            ## Plot second set of telescopes, all have 10x10kpc FOV ##
            
            filterList =[('euclid','nisp','y'),
                          ('vista','vircam','y'),('subaru','hsc','y')]
            # Plot for one galaxy
                                 
            if self.nGalsSelected == 1:
                fig2,ax = plt.subplots(self.nGalsSelected,len(filterList),
                    figsize = (len(filterList)*2,self.nGalsSelected*2),
                    gridspec_kw = {'hspace':0.1,'wspace':0.1})

                for jj in range(0,len(filterList)):
                    ax[jj].set_xlabel('{} {} {}'.format(filterList[jj][0].upper(),filterList[jj][1].upper(),filterList[jj][2].upper()))
                    if self.z==7:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,self.catalogueSelected['fileNumber'].values[0],filterList[jj][2]))
                    else:
                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,filterList[jj][2]))

                    ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data,cmap='Greys')
                    res=ff[0].header['RESOLUTION_PKPC']
                    ax[jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xlabel('pkpc')
                    ff.close()
                #plt.suptitle('Y-band images')

                
            # Plot for multiple galaxies
            else:
                fig2,ax = plt.subplots(self.nGalsSelected,len(filterList),
                            figsize = (len(filterList)*1.5,self.nGalsSelected*1.5),
                            gridspec_kw = {'hspace':0.1,'wspace':0.1})

                for jj in range(0,len(filterList)):
                    ax[0,jj].set_title('{}\n{} {}'.format(filterList[jj][0].upper(),filterList[jj][1].upper(),filterList[jj][2].upper()))
                    for ii in range(0,self.nGalsSelected):
                        if self.z==7:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,self.catalogueSelected['fileNumber'].values[ii],filterList[jj][2]))
                        else:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(filterList[jj][0],filterList[jj][1],self.z,filterList[jj][2]))
                        ax[ii,jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                        res=ff[0].header['RESOLUTION_PKPC']

                        
                        if jj==0:
                            ax[ii,jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii,jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            ax[ii,jj].set_ylabel('pkpc')
                        else:
                            ax[ii,jj].set_yticks([])

                        if ii<self.nGalsSelected-1:
                            ax[ii,jj].set_xticks([])
                        else:
                            ax[ii,jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[ii,jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii,jj].set_xlim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                        ax[ii,jj].set_ylim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                        
                        ax[ii,jj].xaxis.set_label_coords(.5, -.08)
                        ax[ii,jj].yaxis.set_label_coords(-0.08, 0.5) 
                        
                        if jj==0:
                            ax[ii,jj].set_ylabel('pkpc')
    
                        ff.close()
        
                    ax[-1,jj].set_xlabel('pkpc')
                
            return fig,fig2
     
    
    
        ##### All telescopes with all filters #####
        elif telescope == 'all' and filt == 'all': 
            self.plot_galaxies('jwst','nircam',filt = 'all')
            self.plot_galaxies('jwst','miri',filt='all')  
            for telescope in ['hst','euclid','roman','vista','subaru']:
                self.plot_galaxies(telescope,instrument,filt = 'all')
            return

        
        
        ##### A single telescope ##### 
        else:
            # Check that telescope is found. 
            if telescope not in ['hst','jwst','euclid','roman','vista','subaru']:
                    errorstr="ERROR: Telescope '{}' doesn't match those available ('hst','jwst','euclid','roman','vista','subaru')".format(telescope)
                    print(errorstr)
                    return errorstr

            # Check selected an available instrument
            if telescope=='euclid':
                if instrument!='nisp':
                    print('WARNING: Specified euclid instrument "{}" is not nisp, converted to nisp'.format(instrument)) 
                    instrument='nisp'

            elif telescope=='hst':
                if instrument!='wfc3':
                    print('WARNING: Specified hst instrument "{}" is not wfc3, converted to wfc3'.format(instrument)) 
                    instrument='wfc3'

            elif telescope=='roman':
                if instrument!='wfi':
                    print('WARNING: Specified roman instrument "{}" is not wfi, converted to wfi'.format(instrument)) 
                    instrument='wfi'

            elif telescope=='vista':
                if instrument!='vircam':
                    print('WARNING: Specified vista instrument "{}" is not vircam, converted to vircam'.format(instrument)) 
                    instrument='vircam'        

            elif telescope=='subaru':
                if instrument!='hsc':
                    print('WARNING: Specified subaru instrument "{}" is not hsc, converted to hsc'.format(instrument)) 
                    instrument='hsc'

            elif telescope=='jwst':
                if instrument not in ['nircam','miri','all']:
                    errorstr = 'ERROR: Specified jwst instrument "{}" is not nircam, miri, or all, aborted'.format(instrument)
                    print(errorstr)
                    return errorstr

            filterList={'nisp':['y','j','h'],'wfc3':['f105w','f125w','f140w','f160w'],
                   'nircam':['f090w','f115w','f150w','f200w','f277w','f356w','f444w'],'miri':['f560w','f770w'],
                   'all':['f090w','f115w','f150w','f200w','f277w','f356w','f444w','f560w','f770w'],
                   'wfi':['f087','f106','f129','f146','f158','f184'],'hsc':['z','y'],'vircam':['z','y','j','h','ks']}
                        
                    
            ## All available filters ##
            if filt=='all':

                # Plot for one galaxy
                if self.nGalsSelected==1:
                    fig,ax=plt.subplots(self.nGalsSelected,len(filterList[instrument]),
                        figsize=(len(filterList[instrument])*2,self.nGalsSelected*2),
                        gridspec_kw={'bottom':0.05,'top':0.95})
                    plt.suptitle('{} {}'.format(telescope.upper(),instrument.upper()))
                  
                    for jj,filt in enumerate(filterList[instrument]):
                        ax[jj].set_xlabel(filt.upper())
                        
                        if self.z==7:
                            if telescope=='jwst' and instrument=='all':
                                if jj<7:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,'nircam',self.z,self.catalogueSelected['fileNumber'].values[0],filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,'miri',self.z,self.catalogueSelected['fileNumber'].values[0],filt))
                            else:
                                ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,instrument,self.z,self.catalogueSelected['fileNumber'].values[0],filt))
                                
                        else:
                            if telescope=='jwst' and instrument=='all':
                                if jj<7:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,'nircam',self.z,filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,'miri',self.z,filt))
                            else:
                                ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,instrument,self.z,filt))
                            
                        ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data,cmap='Greys')
                        res=ff[0].header['RESOLUTION_PKPC']
                        ax[jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj].set_xlabel('pkpc')
                        ff.close()
                        
                    plt.suptitle('{} {}'.format(telescope.upper(),instrument.upper()))
                    
                    
                # Plot for multiple galaxies
                else:
                    fig,ax=plt.subplots(self.nGalsSelected,len(filterList[instrument]),
                                figsize=(len(filterList[instrument])*1.5,self.nGalsSelected*1.5),
                                gridspec_kw = {'hspace':0.1,'wspace':0.1})

                    for jj,filt in enumerate(filterList[instrument]):
                        ax[0,jj].set_title(filt.upper())
                        for ii in range(0,self.nGalsSelected):
                            if self.z==7:
                                if telescope=='jwst' and instrument=='all':
                                    if jj<7:
                                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,'nircam',self.z,self.catalogueSelected['fileNumber'].values[ii],filt))
                                    else:
                                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,'miri',self.z,self.catalogueSelected['fileNumber'].values[ii],filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,instrument,self.z,self.catalogueSelected['fileNumber'].values[ii],filt))

                            else:
                                if telescope=='jwst' and instrument=='all':
                                    if jj<7:
                                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,'nircam',self.z,filt))
                                    else:
                                        ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,'miri',self.z,filt))
                                else:
                                    ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,instrument,self.z,filt))

                            ax[ii,jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                            res=ff[0].header['RESOLUTION_PKPC']

                            if jj==0:
                                ax[ii,jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                                ax[ii,jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                                ax[ii,jj].set_ylabel('pkpc')
                            else:
                                ax[ii,jj].set_yticks([])
                                
                            if ii<self.nGalsSelected-1:
                                ax[ii,jj].set_xticks([])
                            else:
                                ax[ii,jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                                ax[ii,jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            

                            ax[ii,jj].set_xlim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                            ax[ii,jj].set_ylim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)

                            ax[ii,jj].xaxis.set_label_coords(.5, -.08)
                            ax[ii,jj].yaxis.set_label_coords(-0.06, 0.5)              
                      

                            ff.close()

                        ax[-1,jj].set_xlabel('pkpc')
            
            
            ## Single filter ##
            else:
                    
                if filt not in filterList[instrument]:
                    errorstr='ERROR: Specified filter "{}" is not in filter list, aborted'.format(filt)
                    print(errorstr)
                    return errorstr
                    
                else:
                    fig,ax=plt.subplots(1,self.nGalsSelected,figsize=(self.nGalsSelected*2,2))
                    plt.suptitle('{} {} {}'.format(telescope.upper(),instrument.upper(),filt.upper()))
                    for ii in range(0,self.nGalsSelected):
                        if self.z==7:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}-file{}_{}_v1_sim-psf.fits'.format(telescope,instrument,self.z,self.catalogueSelected['fileNumber'].values[ii],filt))
                        else:
                            ff = fits.open(self.image_path+'hlsp_bluetides_{}_{}_z{}_{}_v1_sim-psf.fits'.format(telescope,instrument,self.z,filt))
                        if self.nGalsSelected==1:
                            ax.imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                        else:
                            ax[ii].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                    
                        res=ff[0].header['RESOLUTION_PKPC']
                        ax[ii].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[ii].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[ii].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii].set_xlabel('pkpc')
                        ax[ii].xaxis.set_label_coords(.5, -.08)
                        ff.close()
            return fig
