import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from astropy.io import fits
import pickle


class Catalogue:

    def __init__(self,z=7,catalogue_path='./',image_path='/home/mmarshal/BLUETIDES/kwatts/Images/Noiseless_FITS/MaddieUpdate/'):
        # Load in catalogue and save columns as variables for ease of access
        cat = pd.read_csv(catalogue_path+'TestCatalogue.csv',index_col=0)
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
                print(ff,lum_constraint[ff][0],lum_constraint[ff][1],lum_constraint[ff][2])
                filt = lum_constraint[ff][0]
                lum_constraint_ff = (self.catalogue['lum_'+filt] > 10**lum_constraint[ff][1]) & (self.catalogue['lum_'+filt] < 10**lum_constraint[ff][2])
                lum_constraint_allfilt = np.logical_and(lum_constraint_ff,lum_constraint_allfilt)
            lum_constraint = lum_constraint_allfilt
        else:
            lum_constraint = np.ones(self.nGals,dtype='Bool')
      
    
        if rad_constraint:
            rad_constraint_allfilt = np.ones(self.nGals,dtype='Bool')
            for ff in range(0,nFilters):
                print(ff,rad_constraint[ff][0],rad_constraint[ff][1],rad_constraint[ff][2])
                filt = rad_constraint[ff][0]
                rad_constraint_ff = (self.catalogue['rad_'+filt] > rad_constraint[ff][1]) & (self.catalogue['rad_'+filt] < rad_constraint[ff][2])
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
        
        # If wanting all telescopes with Y band
        if telescope.upper() == 'ALL' and filt.upper() == 'Y':
            filterList = [('JWST','NIRCAM','F115W'),('HST','WFC3','f105w'),('Euclid','NISP','Y'),('Roman','WFI','F106'),
                          ('VISTA','VIRCAM','Y'),('Subaru','HSC','y')]
        
            # Plot for one galaxy
            if self.nGalsSelected == 1:
                fig,ax = plt.subplots(self.nGalsSelected,len(filterList),
                    figsize = (len(filterList)*2,self.nGalsSelected*2),
                    gridspec_kw = {'bottom':0.15,'top':0.95})

                for jj in range(0,len(filterList)):
                    ax[jj].set_xlabel('{} {} {}'.format(filterList[jj][0],filterList[jj][1],filterList[jj][2]))
                    print(self.catalogueSelected['extensionNumber'].values)
                    if self.z==7:
                        print(self.catalogueSelected['fileNumber'].values)
                        ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}_{}.fits'.format(filterList[jj][0],
                            self.z,filterList[jj][0],filterList[jj][1],filterList[jj][2],self.catalogueSelected['fileNumber'].values[0]))
                    else:
                        ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}.fits'.format(filterList[jj][0],
                            self.z,filterList[jj][0],filterList[jj][1],filterList[jj][2]))

                    ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data,cmap='Greys')
                    res=ff[0].header['RESOLUTION_PKPC']
                    ax[jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                    ax[jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                    ax[jj].set_xlabel('pkpc')
                    # Need to plot ind+1 as there is a dummy Extension 0 holding only the complete header
                    ff.close()
                plt.suptitle('Y-band images')
                
            # Plot for multiple galaxies
            else:
                fig,ax = plt.subplots(len(filterList),self.nGalsSelected,
                            figsize = (self.nGalsSelected*2,len(filterList)*2),
                            gridspec_kw = {'bottom':0.05,'top':0.95,'hspace':0.07})

                for jj in range(0,len(filterList)):
                    ax[jj,0].set_ylabel('{} {} {}'.format(filterList[jj][0],filterList[jj][1],filterList[jj][2]))
                    print(self.catalogueSelected['extensionNumber'].values)
                    for ii in range(0,self.nGalsSelected):
                        if self.z==7:
                            print(self.catalogueSelected['fileNumber'].values[ii])
                            ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}_{}.fits'.format(filterList[jj][0],
                            self.z,filterList[jj][0],filterList[jj][1],filterList[jj][2],self.catalogueSelected['fileNumber'].values[0]))
                        else:
                            ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}.fits'.format(filterList[jj][0],
                            self.z,filterList[jj][0],filterList[jj][1],filterList[jj][2],self.catalogueSelected['fileNumber'].values[0]))
                        ax[jj,ii].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                        res=ff[0].header['RESOLUTION_PKPC']
                        if ii==0:
                            ax[jj,ii].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[jj,ii].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        else:
                            ax[jj,ii].set_yticks([])
                        ax[jj,ii].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj,ii].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj,ii].set_xlim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                        ax[jj,ii].set_ylim(0,ff[0].header['FIELDOFVIEW_PKPC']/res)
                        ax[jj,ii].set_xlabel('pkpc')
                        ax[jj,ii].xaxis.set_label_coords(.5, -.08)
                        # Need to plot ind+1 as there is a dummy Extension 0 holding only the complete header
                        plt.suptitle('Y-band images')
                        ff.close()
            return fig
     
    
        # If wanting all telescopes with all filters
        elif telescope.upper() == 'ALL' and filt.upper() == 'ALL':
            self.plot_galaxies('JWST','NIRCAM',filt = 'ALL')
            self.plot_galaxies('JWST','MIRI',filt='ALL')  
            for telescope in ['HST','Euclid','Roman','VISTA','Subaru']:
                self.plot_galaxies(telescope,instrument,filt = 'ALL')
            return

        # If wanting a single telescope
        else:
            # Check that telescope is found. 
            if telescope not in ['HST','JWST','Euclid','Roman','VISTA','Subaru']:
                # Convert mismatched case for HST and JWST
                if telescope.lower() in ['HST'.lower(),'JWST'.lower(),'VISTA'.lower()]:
                    telescope = telescope.upper()
                else:
                    errorstr="ERROR: Telescope '{}' doesn't match those available ('HST','JWST','Euclid','Roman','VISTA','Subaru')".format(telescope)
                    print(errorstr)
                    return errorstr


            # All instruments are uppercase
            instrument = instrument.upper()

            # Check selected an available instrument
            if telescope=='Euclid':
                if instrument!='NISP':
                    print('WARNING: Specified Euclid instrument "{}" is not NISP, converted to NISP'.format(instrument)) 
                    instrument='NISP'

            elif telescope=='HST':
                if instrument!='WFC3':
                    print('WARNING: Specified HST instrument "{}" is not WFC3, converted to WFC3'.format(instrument)) 
                    instrument='WFC3'

            elif telescope=='Roman':
                if instrument!='WFI':
                    print('WARNING: Specified Roman instrument "{}" is not WFI, converted to WFI'.format(instrument)) 
                    instrument='WFI'

            elif telescope=='VISTA':
                if instrument!='VIRCAM':
                    print('WARNING: Specified VISTA instrument "{}" is not VIRCAM, converted to VIRCAM'.format(instrument)) 
                    instrument='VIRCAM'        

            elif telescope=='Subaru':
                if instrument!='HSC':
                    print('WARNING: Specified Subaru instrument "{}" is not HSC, converted to HSC'.format(instrument)) 
                    instrument='HSC'

            elif telescope=='JWST':
                if instrument not in ['NIRCAM','MIRI']:
                    errorstr = 'ERROR: Specified JWST instrument "{}" is not NIRCAM or MIRI, aborted'.format(instrument)
                    print(errorstr)
                    return errorstr

            filterList={'NISP':['Y','J','H'],'WFC3':['f105w','f125w','f140w','f160w'],
                   'NIRCAM':['F090W','F115W','F150W','F200W','F277W','F356W','F444W'],'MIRI':['F560W','F770W'],
                   'WFI':['F087','F106','F129','F146','F158','F184'],'HSC':['z','y'],'VIRCAM':['Z','Y','J','H','Ks']}
                        
                    
            # Plot all available filters
            if filt.upper()=='ALL':

                # Plot for one galaxy
                if self.nGalsSelected==1:
                    fig,ax=plt.subplots(self.nGalsSelected,len(filterList[instrument]),
                        figsize=(len(filterList[instrument])*2,self.nGalsSelected*2),
                        gridspec_kw={'bottom':0.05,'top':0.95})
                    plt.suptitle('{} {}'.format(telescope,instrument))
                    for jj,filt in enumerate(filterList[instrument]):
                        ax[jj].set_xlabel(filt)
                        if self.z==7:
                            print(self.catalogueSelected['fileNumber'].values)
                            ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}_{}.fits'.format(telescope,self.z,telescope,instrument,filt,self.catalogueSelected['fileNumber'].values[0]))
                        else:
                            ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}.fits'.format(telescope,self.z,telescope,instrument,filt))
                        ax[jj].imshow(ff[self.catalogueSelected['extensionNumber'].values[0]].data,cmap='Greys')
                        res=ff[0].header['RESOLUTION_PKPC']
                        ax[jj].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[jj].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[jj].set_xlabel('pkpc')
                        ff.close()
                        # Need to plot ind+1 as there is a dummy Extension 0 holding only the complete header

                    plt.suptitle('{} {}'.format(telescope,instrument))
                    
                # Plot for multiple galaxies
                else:
                    fig,ax=plt.subplots(len(filterList[instrument]),self.nGalsSelected,
                                figsize=(self.nGalsSelected*2,len(filterList[instrument])*2),
                                gridspec_kw={'bottom':0.05,'top':0.95})

                    for jj,filt in enumerate(filterList[instrument]):
                        ax[jj,0].set_ylabel(filt)
                        for ii in range(0,self.nGalsSelected):
                            #print(ii,jj)
                            if self.z==7:
                                print(self.catalogueSelected['fileNumber'].values[ii])
                                ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}_{}.fits'.format(telescope,self.z,telescope,instrument,filt,self.catalogueSelected['fileNumber'].values[ii]))
                            else:
                                ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}.fits'.format(telescope,self.z,telescope,instrument,filt))
                            print(self.catalogueSelected['extensionNumber'].values[ii])
                            ax[jj,ii].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                            res=ff[0].header['RESOLUTION_PKPC']
                            ax[jj,ii].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[jj,ii].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            ax[jj,ii].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                            ax[jj,ii].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                            ax[jj,ii].set_xlabel('pkpc')
                            # Need to plot ind+1 as there is a dummy Extension 0 holding only the complete header
                            ff.close()
            # Plot single filter
            else:
                # Filters uppercase for JWST and Euclid, lowercase for HST and Spitzer
                if telescope in ['JWST','Euclid','Roman']:
                    filt=filt.upper()
                elif telescope in ['HST','Subaru']:
                    filt=filt.lower() 
                    
                if filt not in filterList[instrument]:
                    errorstr='ERROR: Specified filter "{}" is not in filter list, aborted'.format(filt)
                    print(errorstr)
                    return errorstr
                    
                else:
                    fig,ax=plt.subplots(1,self.nGalsSelected,figsize=(self.nGalsSelected*2,2))
                    plt.suptitle('{} {} {}'.format(telescope,instrument,filt))
                    for ii in range(0,self.nGalsSelected):
                        if self.z==7:
                            print(self.catalogueSelected['fileNumber'].values[ii])
                            ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}_{}.fits'.format(telescope,self.z,telescope,instrument,filt,self.catalogueSelected['fileNumber'].values[ii]))
                        else:
                            ff = fits.open(self.image_path+'{}/z{}/{}.{}.{}.fits'.format(telescope,self.z,telescope,instrument,filt))
                        if self.nGalsSelected==1:
                            ax.imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                        else:
                            ax[ii].imshow(ff[self.catalogueSelected['extensionNumber'].values[ii]].data,cmap='Greys')
                    # Need to plot ind+1 as there is a dummy Extension 0 holding only the complete header
                        res=ff[0].header['RESOLUTION_PKPC']
                        ax[ii].set_yticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[ii].set_yticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii].set_xticks([0,ff[0].header['FIELDOFVIEW_PKPC']/res])
                        ax[ii].set_xticklabels(['0',str(ff[0].header['FIELDOFVIEW_PKPC'])])
                        ax[ii].set_xlabel('pkpc')
                        ff.close()
            return fig
