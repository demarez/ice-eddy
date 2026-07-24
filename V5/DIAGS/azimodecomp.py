###author Charly de Marez
### (very) inspired from Thomas Meunier's anamodfun.m

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import sys
sys.path.append('/home2/datahome/cdemarez/MODULES_PY') 
from tools import *


def AZIMUTHAL_MODES(x,y,VAR,nr,nth,modes):

    ####################################################################
    ####### routine to decompose a field into azimuthal modes ##########
    ####################################################################

    #we decompose the quantity into azimuthal modes
    #for each r, the phase and amplitude is different


    #interpolate on a (r,theta) grid
    r=np.linspace(0,np.max(np.sqrt(x**2+y**2)),nr)
    theta=np.linspace(0,2*np.pi,nth)
    dth=np.diff(theta)[0]
    r,theta=np.meshgrid(r,theta)
    r=r.T
    theta=theta.T
    x_polar=r*np.cos(theta)
    y_polar=r*np.sin(theta)

    #the field in r,theta coordinates
    art=griddata((np.ravel(x),np.ravel(y)), np.ravel(VAR), (x_polar,y_polar), method='cubic')


    cj=np.zeros((len(modes),nr))
    sj=np.zeros((len(modes),nr))
    ampl=np.zeros((len(modes),nr))
    phaz=np.zeros((len(modes),nr))
    psi=np.zeros((len(modes),nr,nth))
    psi_xy=np.zeros((len(modes),x.shape[0],x.shape[1]))

    #loop over modes
    for imod in range(len(modes)):
        #the mode on which we decompose
        mode=modes[imod]
        #integration of VAR*cos(theta) dtheta
        cj[imod,:]=np.nansum(art[:,:]*np.cos(theta[:,:]*mode),axis=1)*dth
        #integration of VAR*sin(theta) dtheta
        sj[imod,:]=np.nansum(art[:,:]*np.sin(theta[:,:]*mode),axis=1)*dth

        if mode !=0:
            cj[imod,:]=cj[imod,:]/np.pi
            sj[imod,:]=-sj[imod,:]/np.pi
        else :
            cj[imod,:]=cj[imod,:]/(2*np.pi)
            sj[imod,:]=-sj[imod,:]/(2*np.pi)
        #the amplitude and the phase of the mode depending on the radial coordinate
        ampl[imod,:] = cj[imod,:]*cj[imod,:]+sj[imod,:]*sj[imod,:]
        phaz[imod,:] = np.arctan2(sj[imod,:],cj[imod,:])
        phaz[imod,:]=np.where(phaz[imod,:]<0,phaz[imod,:]+2.0*np.pi,phaz[imod,:])

        #the quantity projected on the mode
        psi[imod:,:]=np.sqrt(ampl[imod,:])[:,None]*np.cos(mode*theta[:,:]+phaz[imod,:][:,None])
        psi_xy[imod,:,:]=griddata((np.ravel(x_polar),np.ravel(y_polar)),
                                  np.ravel(psi[imod,:,:]), (x,y), method='cubic')


    return r,theta,x_polar,y_polar,art,ampl,phaz,psi,psi_xy


def AZIMUTHAL_MODES_fast(x,y,x_polar,y_polar,r,theta,nr,nth,dth,vtx_x_to_r,wts_x_to_r,vtx_r_to_x, wts_r_to_x,VAR,modes):

    ####################################################################
    ####### routine to decompose a field into azimuthal modes ##########
    ####################################################################

    #we decompose the quantity into azimuthal modes
    #for each r, the phase and amplitude is different


    
    #the field in r,theta coordinates
    art=interpolate(VAR.flatten(), vtx_x_to_r,wts_x_to_r)
    art=art.reshape(x_polar.shape[0],x_polar.shape[1])
    


    cj=np.zeros((len(modes),nr))
    sj=np.zeros((len(modes),nr))
    ampl=np.zeros((len(modes),nr))
    phaz=np.zeros((len(modes),nr))
    psi=np.zeros((len(modes),nr,nth))
    psi_xy=np.zeros((len(modes),x.shape[0],x.shape[1]))

    #loop over modes
    for imod in range(len(modes)):
        #the mode on which we decompose
        mode=modes[imod]
        #integration of VAR*cos(theta) dtheta
        cj[imod,:]=np.nansum(art[:,:]*np.cos(theta[:,:]*mode),axis=1)*dth
        #integration of VAR*sin(theta) dtheta
        sj[imod,:]=np.nansum(art[:,:]*np.sin(theta[:,:]*mode),axis=1)*dth

        if mode !=0:
            cj[imod,:]=cj[imod,:]/np.pi
            sj[imod,:]=-sj[imod,:]/np.pi
        else :
            cj[imod,:]=cj[imod,:]/(2*np.pi)
            sj[imod,:]=-sj[imod,:]/(2*np.pi)
        #the amplitude and the phase of the mode depending on the radial coordinate
        ampl[imod,:] = cj[imod,:]*cj[imod,:]+sj[imod,:]*sj[imod,:]
        phaz[imod,:] = np.arctan2(sj[imod,:],cj[imod,:])
        phaz[imod,:]=np.where(phaz[imod,:]<0,phaz[imod,:]+2.0*np.pi,phaz[imod,:])

        #the quantity projected on the mode
        psi[imod:,:]=np.sqrt(ampl[imod,:])[:,None]*np.cos(mode*theta[:,:]+phaz[imod,:][:,None])
        
        psi_xy_tmp=interpolate(psi[imod,:,:].flatten(),vtx_r_to_x,wts_r_to_x)
        psi_xy_tmp[imod,:,:]=psi_xy_tmp.reshape(x.shape[0],x.shape[1])
        
        
        


    return art,ampl,phaz,psi,psi_xy


def AZIMUTHAL_MODES_fast_v2(x,y,r_x,rmax,x_polar,y_polar,r,theta,nr,nth,dth,vtx_x_to_r,wts_x_to_r,vtx_r_to_x, wts_r_to_x,VAR,modes):

    ####################################################################
    ####### routine to decompose a field into azimuthal modes ##########
    ####################################################################

    #we decompose the quantity into azimuthal modes
    #for each r, the phase and amplitude is different


    
    #the field in r,theta coordinates
    art=interpolate(VAR.flatten(), vtx_x_to_r,wts_x_to_r)
    art=art.reshape(x_polar.shape[0],x_polar.shape[1])
    


    cj=np.zeros((len(modes),nr))
    sj=np.zeros((len(modes),nr))
    ampl=np.zeros((len(modes),nr))
    phaz=np.zeros((len(modes),nr))
    psi=np.zeros((len(modes),nr,nth))
    psi_xy=np.zeros((len(modes),x.shape[0],x.shape[1]))

    #loop over modes
    for imod in range(len(modes)):
        #the mode on which we decompose
        mode=modes[imod]
        #integration of VAR*cos(theta) dtheta
        cj[imod,:]=np.nansum(art[:,:]*np.cos(theta[:,:]*mode),axis=1)*dth
        #integration of VAR*sin(theta) dtheta
        sj[imod,:]=np.nansum(art[:,:]*np.sin(theta[:,:]*mode),axis=1)*dth

        if mode !=0:
            cj[imod,:]=cj[imod,:]/np.pi
            sj[imod,:]=-sj[imod,:]/np.pi
        else :
            cj[imod,:]=cj[imod,:]/(2*np.pi)
            sj[imod,:]=-sj[imod,:]/(2*np.pi)
        #the amplitude and the phase of the mode depending on the radial coordinate
        ampl[imod,:] = cj[imod,:]*cj[imod,:]+sj[imod,:]*sj[imod,:]
        phaz[imod,:] = np.arctan2(sj[imod,:],cj[imod,:])
        phaz[imod,:]=np.where(phaz[imod,:]<0,phaz[imod,:]+2.0*np.pi,phaz[imod,:])

        #the quantity projected on the mode
        psi[imod:,:]=np.sqrt(ampl[imod,:])[:,None]*np.cos(mode*theta[:,:]+phaz[imod,:][:,None])
        
        psi_xy_tmp=interpolate(psi[imod,:,:].flatten(),vtx_r_to_x,wts_r_to_x)
        psi_xy_tmp=psi_xy_tmp.reshape(x.shape[0],x.shape[1])
        
        
        
        
        
        psi_xy[imod,:,:]=np.where(r_x>(rmax-2),np.nan,psi_xy_tmp) #the -2 is needed to remove problem of interpolationat the edge !!
        

    return art,ampl,phaz,psi,psi_xy




