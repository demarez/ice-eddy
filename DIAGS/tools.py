#!/usr/bin/python
# Filename: romstools.py


#Netcdf IO module
from netCDF4 import Dataset

#module for numerics
import numpy as np

#module for copy
from copy import copy

""
from itertools import product


import math








##################
#####CHARLY
##################

########################### done for amedargo


def zlevs_old(h,zeta, hc, Cs_r, Cs_w):
    N=Cs_r.shape[0]
    z_r=np.zeros((N,h.shape[0],h.shape[1]))
    
    for k in range(N):
        cff= hc*((k+1-N)-0.5)/N
        z_r[k,:,:]=zeta+(zeta+h)*(cff+Cs_r[k]*h)/(hc+h)
    
    z_w=np.zeros((N+1,h.shape[0],h.shape[1]))
    z_w[0,:,:] = - h
    z_w[-1,:,:] = zeta
    
    for k in range(1,N+1):
        cff= hc*((k-N))/N
        z_w[k,:,:]=zeta+(zeta+h)*(cff+Cs_w[k]*h)/(hc+h)
    
    return z_r,z_w

def get_name_large_simu(itime,Nsteps):
    name='croco_his.'
    itmp=itime%Nsteps
    ifile=(itime//Nsteps)*Nsteps
    name=name+str(ifile).zfill(5)+'.nc'
    return name,itmp
    



##############@functions to do the interpolation as in griddata
### cf https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids


import scipy.interpolate as spint
import scipy.spatial.qhull as qhull
import numpy as np

def interp_weights(xy, uv,d=2):
    tri = qhull.Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

def interpolate(values, vtx, wts):
    return np.einsum('nj,nj->n', np.take(values, vtx), wts)



def do_interp_weights_calc(X,Y,Xi,Yi):           #####beware, it takes y and x !!!!
    xy=np.zeros([X.shape[0]*X.shape[1],2])
    xy[:,0]=Y.flatten()
    xy[:,1]=X.flatten()
    uv=np.zeros([Xi.shape[0]*Xi.shape[1],2])
    uv[:,0]=Yi.flatten()
    uv[:,1]=Xi.flatten()
    vtx, wts = interp_weights(xy, uv)
    return vtx,wts


def do_interp_weights_calc_flat(X,Y,Xi,Yi):           #####version which takes flat X and Y
    xy=np.zeros([X.shape[0],2])
    xy[:,0]=Y
    xy[:,1]=X
    uv=np.zeros([Xi.shape[0]*Xi.shape[1],2])
    uv[:,0]=Yi.flatten()
    uv[:,1]=Xi.flatten()
    vtx, wts = interp_weights(xy, uv)
    return vtx,wts



#compute ur and utheta

## Routine to decompose the velocity into an azimuthhal and a radial component
def extract_u_proj_azim(x,y,U_take,V_take):
    #this version takes lon and lat wrt the center
    Unorm_take=np.sqrt(U_take**2+V_take**2)
    # now do the projection on the tilted axis between the float and the center oof the eddy
    alph=np.arctan2(y,x)
    beta=np.arctan2(V_take,U_take)
    gam=beta-alph
    #compute componant on tilted axis
    Ur=Unorm_take*np.cos(gam)
    Utheta=Unorm_take*np.sin(gam)
    #print("U theta=",Utheta)
    #print("U r=",Ur)
    return Ur,Utheta






############################done for leewa

def choose_zone(field,lon,lat,corners):
    tmp1=lon>corners[0]
    tmp2=lon<corners[2]
    tmp3=lat>corners[1]
    tmp4=lat<corners[3]

    zone=tmp1&tmp2&tmp3&tmp4
    return field[zone]

def print_date(sec):
	#date a chercher dans le fichier ini :
	#CHABU:	
	#year_ini=2013
	#month_ini=7
	#day_ini=18

	#CHABAM:	
	year_ini=2016
	month_ini=2
	day_ini=25
    
	#LEEWA v1 ficher ini=chabam_his.00100.nc:	
	year_ini=2016
	month_ini=3
	day_ini=20

	sec_ini=(month_ini-1)*30*24*3600+(day_ini-1)*24*3600

	sec_tot=sec_ini+sec
	year=sec_tot//(365*24*3600)
	tmp1=sec_tot%(365*24*3600)	
	month=tmp1//(30*24*3600)
	#print(month)
	tmp2=tmp1%(30*24*3600)
	day=tmp2//(24*3600)
	#print (year,month,day)

	year=year_ini+year
	if month==12 : month=1
	else: month=month+1
	day=day+1

	return np.array([year,month,day])







def sph_distance(lat1, lon1,lat2, lon2):
	radius = 6371. # km

	dlat = math.radians(lat2-lat1)
	dlon = math.radians(lon2-lon1)
	a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1))* math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
	c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
	d = radius * c
	return d



def peakdet_org(v, delta,aff=None, x=None):
    
    maxtab = []
    mintab = []
    if x is None:
        x = np.arange(len(v))
    v = np.asarray(v)
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN
    lookformax = True
    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        if lookformax:
            if this < mx - delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
                
    maxta=np.array(maxtab)
    minta=np.array(mintab)
    mins_x=minta#[:, 0]
    maxs_x=maxta#[:, 0]
    
    return mins_x,maxs_x




def H_bottom_layer (step_file,z,rho,rho_choice):

    deltaH=np.zeros(z[0,:].shape)
    for idx in range (len(z[0,:])):
        tmp,=np.where(rho[step_file,:,idx]>rho_choice)
        tmp2=z[tmp,idx]
        deltaH[idx]=np.abs(tmp2[-1]-tmp2[0])
    return deltaH


def H_bottom_mixed_layer (step_file,z,rho,h,dif_rho):
    deltaH=np.zeros(z[0,:].shape)

    for idx in range (len(z[0,:])):
        rho_max_bottom=np.max(rho[step_file,-1,idx])
        tmp,=np.where(rho[step_file,:,idx]>(rho_max_bottom-dif_rho))
        tmp2=z[tmp,idx]
        deltaH[idx]=np.abs(tmp2[-1]-tmp2[0])

    z_bml=-h+deltaH
    return deltaH,z_bml


def H_bottom_mixed_layer_3D (z_r,rho,h,dif_rho):
    print('Computing hbbl')
    rho_max_bottom=rho[0,:,:]

    deltaH=np.zeros(rho_max_bottom.shape)

    for idx in range(len(rho[0,0,:])):
        #print(int(idx/len(rho[0,0,:])*100))
        for idy in range(len(rho[0,:,0])):

            tmp,=np.where(rho[:,idy,idx]>(rho_max_bottom[idy,idx]-dif_rho))
            tmp2=z_r[tmp,idy,idx]
            deltaH[idy,idx]=np.abs(tmp2[-1]-tmp2[0])

    z_bml=-h+deltaH
    return deltaH,z_bml   


def rotate_uv(u,v,angle):
    urot = u*np.cos(np.deg2rad(angle)) - v*np.sin(np.deg2rad(angle))
    vrot = u*np.sin(np.deg2rad(angle)) + v*np.cos(np.deg2rad(angle))
    return urot,vrot




###################
#####PAULINE--modif
###################

from numpy import arange,max,abs,where,ceil,zeros,amin,diag,tile,nan
import sys
from scipy.interpolate import interp1d,interp2d,griddata	


def get_section(variable,lon,lat,f,mask,h,zeta, hc, Cs_r, Cs_w,lonsec,latsec,N,Nb_sigma):
        ""
        def Forder(var):
            return np.asfortranarray(var.T,dtype=np.float64)
        ""
        [M,L]	= lon.shape
        # dl:largest interval between 2 grid cells == resolution (in deg)
        dl      = 1*max((max(max(abs(lon[1:M,:]-lon[0:M-1,:]))), 
                                max(max(abs(lon[:,1:L]-lon[:,0:L-1]))),
                                max(max(abs(lat[1:M,:]-lat[0:M-1,:]))),
                                max(max(abs(lat[:,1:L]-lat[:,0:L-1])))))	



        [M,L]  	= lon.shape				# may have changed if u or v variable	
        minlon	= np.min(lonsec)			# find minimum subgrid (new square domain dilimited by the seciton) limits
        minlat	= np.min(latsec)
        maxlon	= np.max(lonsec)
        maxlat	= np.max(latsec)
        D = abs( (latsec[-1]-latsec[0])*lon - (lonsec[-1]-lonsec[0])*lat + lonsec[-1]*latsec[0] - latsec[-1]*lonsec[0]) / np.sqrt( (latsec[-1]-latsec[0])**2 + (lonsec[-1]-lonsec[0])**2)

        isub = np.where( (D<2*dl) & (lon > min(lonsec)-2*dl) & (lon < max(lonsec)+2*dl) & (lat > min(latsec) -2*dl) & (lat < max(latsec)+2*dl) )


        londata = lon[isub[0],isub[1]] 
        latdata = lat[isub[0],isub[1]]
        #mask    = mask[isub[0],isub[1]]  
        '''
        imin    = np.min(np.where(lon[1,:]>=minlon))#[1][0]               # 'where' is a 2D array, ([j-coord] [i-coord]),
        imax    = np.max(np.where(lon[1,:]<=maxlon))#[1][0]               # imin is the 1st occurence of the 2nd tabular
        jmin    = np.min(np.where(lat[:,1]>=minlat))#[0][1]
        jmax    = np.max(np.where(lat[:,1]<=maxlat))#[0][1]
        print 'indices of the section : imin,imax,jmin,jmax = ',imin,imax,jmin,jmax
        lon     = lon[jmin:jmax,imin:imax]              # extract the subgrid, overwrite lon, lat, mask
        lat     = lat[jmin:jmax,imin:imax]
        mask    = mask[jmin:jmax,imin:imax]
        '''    
        Npts	= 2*ceil(max((abs(lonsec[1]-lonsec[0])/dl,abs(latsec[1]-latsec[0])/dl)))
        #print 'Npts=',Npts
        lonsec	= arange(lonsec[0],lonsec[1]+0.5*(lonsec[1]-lonsec[0])/Npts,(lonsec[1]-lonsec[0])/Npts)
        latsec	= arange(latsec[0],latsec[1]+0.5*(latsec[1]-latsec[0])/Npts,(latsec[1]-latsec[0])/Npts)
        #londata = lon[1,:]
        #latdata = lat[:,1]		
        Npts	= lonsec.shape[0]
        VAR	= np.zeros((N,Npts))				

        #to get the bottom line (zri) !!!!!! fancy plot
        '''
        for k in arange(0,Nb_sigma,1):
            f_interp =  interp2d(londata, latdata,zr[k,jmin:jmax,imin:imax], kind='cubic')
            zri[:,Nb_sigma-1-k] = [float(f_interp(XX,YY)) for XX,YY in zip(lonsec,latsec)]
        '''
        #londata_mat = lon
        #latdata_mat = lat

        lonsec_mat = tile(lonsec,(Npts,1))
        latsec_mat = tile(latsec,(Npts,1)) 
        (elem,coef) =  get_tri_coef(londata,latdata,lonsec_mat,latsec_mat) 
        elem = elem[0,:,:]
        coef = coef[0,:,:]




        (zr,zw)	= zlevs(h,zeta, hc, Cs_r, Cs_w)
        #print (zr.shape)
        #zr = Forder(zr)
        #zw = Forder(zw)
        #print (zr.shape)

        zri = np.zeros((Nb_sigma,Npts))
        for k in arange(0,Nb_sigma,1):   
            zr_input = zr[k,isub[0],isub[1]]
            zri[Nb_sigma-1-k,:] = np.sum(coef*zr_input[elem],1)
        if N>1:
            for k in arange(0,N,1): 
                var_input= variable[k,isub[0],isub[1]]
                VAR[N-1-k,:]  = np.sum(coef*var_input[elem],1)
                #print k
        else :
            for k in arange(0,N,1): 
                var_input      = variable[isub[0],isub[1]]
                VAR[N-1-k,:]  = np.sum(coef*var_input[elem],1)

        dist= zeros((Npts,))



        for i in arange(Npts):	
            dist[i]= sph_distance(latsec[0],lonsec[0],latsec[i],lonsec[i])#/1000 #distance in km

        dist= tile(dist,(N,1))# duplicate the distance to facilitate the plot management

        return dist,VAR,zri,isub




############
#############



def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof




####################
# ####ORIGINAL-TOOLS
# ###################

def zlevs(h,zeta, hc, Cs_r, Cs_w):
    N=Cs_r.shape[0]
    z_r=np.zeros((N,h.shape[0],h.shape[1]))
    
    for k in range(N):
        cff= hc*((k+1-N)-0.5)/N
        z_r[k,:,:]=zeta+(zeta+h)*(cff+Cs_r[k]*h)/(hc+h)
    
    z_w=np.zeros((N+1,h.shape[0],h.shape[1]))
    z_w[0,:,:] = - h
    z_w[-1,:,:] = zeta
    
    for k in range(1,N+1):
        cff= hc*((k-N))/N
        z_w[k,:,:]=zeta+(zeta+h)*(cff+Cs_w[k]*h)/(hc+h)
    
    return z_r,z_w

#######################################################
# interpolate a 3D variable on a horizontal level of constant depth
# ######################################################

def vinterp_V2(var,z,depth,topo=None,cubic=0):
    [N,Mp,Lp]=z.shape
    #depth=0, just take surface field
    if depth>0: 
        varz = np.nan
    #Simple linear interpolation
    elif cubic==0:
        levs2=sum(z<depth)-1
        levs2[levs2==N-1]=N-2
        levs2[levs2==-1]=0
        levs1=levs2+1
        X,Y=np.meshgrid(np.arange(0,Lp),np.arange(0,Mp))
        pos1=levs1,Y,X
        pos2=levs2,Y,X
        z1=z[pos1]
        z2=z[pos2]
        v1=var[pos1]
        v2=var[pos2]
        varz = (((v1-v2)*depth+v2*z1-v1*z2)/(z1-z2))
        if np.any(topo!=None): 			
            varz[depth<-1*topo]=np.nan
    #Cubic interpolation (see ShchepetkinMcWilliams08.pdf)
    elif cubic==1:

        #print ('cubic interpolation')

        #find the closest level BELOW depth
        levs2=sum(z<depth)-1
        levs1=copy(levs2)
        #levs2[levs1==N-1]=N-2

        #cubic interpolation will use 4 values of var and z in the vertical (2 below and 2 above depth)
        Nlev = 4 

        #prepare arrays for intermediate variables:
        X,Y=np.meshgrid(np.arange(0,Lp),np.arange(0,Mp))
        levs=np.zeros((Nlev,Mp,Lp),int); Xlev=np.zeros((Nlev,Mp,Lp),int); Ylev=np.zeros((Nlev,Mp,Lp),int)
        for ilev in range(Nlev):
            levs[ilev,:,:]=levs2+ilev-1
            Xlev[ilev,:,:]=X
            Ylev[ilev,:,:]=Y


        levs[levs>N-1]=N-1
        levs[levs<0]=0

        pos=levs,Y,X; zz=z[pos]; vark=var[pos]


        #######################################################

        test0=np.zeros((Mp,Lp)); test0[levs2==-1]=1; 
        test1=np.zeros((Mp,Lp)); test1[levs2==0]=1;
        testN1=np.zeros((Mp,Lp)); testN1[levs2==N-2]=1; 
        testN=np.zeros((Mp,Lp)); testN[levs2==N-1]=1;
 
        #######################################################

        zz[1:-1,:,:] = testN * zz[:-2,:,:] + test0 * zz[2:,:,:] + (1 - test0 - testN) * zz[1:-1,:,:]

        dzz = zz[1:,:,:]- zz[:-1,:,:]; 
        dzz[-1,:,:] = testN1 * dzz[1,:,:] + (1-testN1)* dzz[-1,:,:]
        dzz[0,:,:] = test1 * dzz[1,:,:] + (1-test1)* dzz[0,:,:]

        vark[1:-1,:,:] = testN * vark[:-2,:,:] + test0 * vark[2:,:,:] + (1 - test0 - testN) * vark[1:-1,:,:]

        dvark = vark[1:,:,:]-vark[:-1,:,:]; 
        dvark[-1,:,:] = testN1 * dvark[1,:,:] + (1-testN1)* dvark[-1,:,:]
        dvark[0,:,:] = test1 * dvark[1,:,:] + (1-test1)* dvark[0,:,:]


        FC0 = (dvark[1:,:,:]+dvark[:-1,:,:])*dzz[1:,:,:]*dzz[:-1,:,:]
        FC1 = (dzz[1:,:,:]+dzz[:-1,:,:])*dvark[1:,:,:]*dvark[:-1,:,:]
        val=dvark[1:,:,:]*dvark[:-1,:,:]

        FC0[val<=0]=1; FC1[val<=0]=0; FC = FC1/FC0


        #######################################################
        cff = 1/dzz[1,:,:]; p=depth-zz[1,:,:]; q=zz[2,:,:]-depth

        varz = cff*(q*vark[1,:,:]+p*vark[2,:,:]- (1-test0-testN) * cff*p*q*(cff*(q-p)*dvark[1,:,:]+p*FC[1,:,:]-q*FC[0,:,:]))     


        #######################################################

        varz[depth<-1*topo]=np.nan

   
    return varz





#######################################################
# interpolate a 3D variable on horizontal levels of constant depths 
# ######################################################

def vinterps(var,z,depths,topo, cubic=0):

    [N,Mp,Lp]=var.shape
    Nz=len(depths)

    #if var not on rho-grid: interpolate z and topo to the same grid than var (u,v,or psi)
    if var.shape!=z.shape:
        if (var.shape[1]==z.shape[1]-1) and (var.shape[2]==z.shape[2]-1):
            z = rho2psi(z); topo = rho2psi(topo)
        elif (var.shape[1]==z.shape[1]-1):
            z = rho2v(z); topo = rho2v(topo)
        elif (var.shape[2]==z.shape[2]-1):
            z = rho2u(z); topo = rho2u(topo)

    if len(depths)==1:
        vnew=vinterp(var,z,depths[0],topo,cubic)

    else:
        [N,Mp,Lp]=var.shape; Nz=len(depths); vnew=np.zeros((Nz, Mp,Lp))
        for iz in range(0, Nz, 1):
            vnew[iz,:,:]=vinterp(var,z,depths[iz],topo,cubic)

    return vnew




#######################################################
# Transfert a field at psi points to rho points
# ######################################################

def psi2rho(var_psi):

    if np.rank(var_psi)<3:
        var_rho = psi2rho_2d(var_psi)
    else:
        var_rho = psi2rho_3d(var_psi)

    return var_rho


""
def psi2rho_2d(var_psi):

    [M,L]=var_psi.shape
    Mp=M+1
    Lp=L+1
    Mm=M-1
    Lm=L-1

    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,1:L]=0.25*(var_psi[0:Mm,0:Lm]+var_psi[0:Mm,1:L]+var_psi[1:M,0:Lm]+var_psi[1:M,1:L])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]

    return var_rho

""
def psi2rho_3d(var_psi):


    [Nz,Mz,Lz]=var_psi.shape
    var_rho=np.zeros((Nz,Mz+1,Lz+1))

    for iz in range(0, Nz, 1):    
        var_rho[iz,:,:]=psi2rho_2d(var_psi[iz,:,:])


    return var_rho

#######################################################
# Transfert a field at rho points to psi points
# ######################################################

def rho2psi(var_rho):

    if np.rank(var_rho)<3:
        var_psi = rho2psi_2d(var_rho)
    else:
        var_psi = rho2psi_3d(var_rho)

    return var_psi


""
def rho2psi_2d(var_rho):

    var_psi = 0.25*(var_rho[1:,1:]+var_rho[1:,:-1]+var_rho[:-1,:-1]+var_rho[:-1,1:])

    return var_psi

""
def rho2psi_3d(var_rho):

    var_psi = 0.25*(var_rho[:,1:,1:]+var_rho[:,1:,:-1]+var_rho[:,:-1,:-1]+var_rho[:,:-1,1:])

    return var_psi


#######################################################
# Transfert a 2 or 3-D field at rho points to u points
# ######################################################

def rho2u(var_rho):

    if np.rank(var_rho)<3:
        var_u = rho2u_2d(var_rho)
    else:
        var_u = rho2u_3d(var_rho)

    return var_u

def rho2u_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,0:L]+var_rho[:,1:Lp])

    return var_u


def rho2u_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    L=Lp-1
    var_u=0.5*(var_rho[:,:,0:L]+var_rho[:,:,1:Lp])

    return var_u

#######################################################
# Transfert a 3-D field at rho points to v points
# ######################################################

def rho2v(var_rho):

    if np.rank(var_rho)<3:
        var_v = rho2v_2d(var_rho)
    else:
        var_v = rho2v_3d(var_rho)

    return var_v

""
def rho2v_2d(var_rho):

    [Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[0:M,:]+var_rho[1:Mp,:]);

    return var_v

""
def rho2v_3d(var_rho):

    [N,Mp,Lp]=var_rho.shape
    M=Mp-1
    var_v=0.5*(var_rho[:,0:M,:]+var_rho[:,1:Mp,:]);

    return var_v

#######################################################
# Transfert a 2-D field at u points to the rho points
# ######################################################

def u2rho(var_u):


    if np.rank(var_u)<3:
        var_rho = u2rho_2d(var_u)
    else:
        var_rho = u2rho_3d(var_u)

    return var_rho

""
def u2rho_2d(var_u):

    [Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[:,1:L]=0.5*(var_u[:,0:Lm]+var_u[:,1:L])
    var_rho[:,0]=var_rho[:,1]
    var_rho[:,Lp-1]=var_rho[:,L-1]
    return var_rho

""
def u2rho_3d(var_u):

    [N,Mp,L]=var_u.shape
    Lp=L+1
    Lm=L-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,:,1:L]=0.5*(var_u[:,:,0:Lm]+var_u[:,:,1:L])
    var_rho[:,:,0]=var_rho[:,:,1]
    var_rho[:,:,Lp-1]=var_rho[:,:,L-1]
    return var_rho


#######################################################
# Transfert a 2 or 2-D field at v points to the rho points
# ######################################################

def v2rho(var_v):

    if np.rank(var_v)<3:
        var_rho = v2rho_2d(var_v)
    else:
        var_rho = v2rho_3d(var_v)

    return var_rho

""
def v2rho_2d(var_v):

    [M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((Mp,Lp))
    var_rho[1:M,:]=0.5*(var_v[0:Mm,:]+var_v[1:M,:])
    var_rho[0,:]=var_rho[1,:]
    var_rho[Mp-1,:]=var_rho[M-1,:]

    return var_rho

""
def v2rho_3d(var_v):

    [N,M,Lp]=var_v.shape
    Mp=M+1
    Mm=M-1
    var_rho=np.zeros((N,Mp,Lp))
    var_rho[:,1:M,:]=0.5*(var_v[:,0:Mm,:]+var_v[:,1:M,:])
    var_rho[:,0,:]=var_rho[:,1,:]
    var_rho[:,Mp-1,:]=var_rho[:,M-1,:]

    return var_rho

#######################################################
# Compute vorticity of a 2-D field
# ######################################################
"""
    vrt = dv/dx - du/dy

    u,v on originals u and v grids, respectively
    vrt is outputed on psi grid 

"""


def vort(u,v,pm,pn):

    if np.rank(u)<3:
        vrt = vort_2d(u,v,pm,pn)
    else:
        vrt = vort_3d(u,v,pm,pn)

    return vrt

""
def vort_2d(u,v,pm,pn):

    [Mp,Lp]=pm.shape
    L=Lp-1
    M=Mp-1

    dm_u=2*u/(pm[:,0:L]+pm[:,1:Lp])
    dn_v=2*v/(pn[0:M,:]+pn[1:Mp,:])

    iA_q=0.0625*(pm[0:M,0:L]+pm[0:M,1:Lp]+pm[1:Mp,1:Lp]+pm[1:Mp,0:L])\
        *(pn[0:M,0:L]+pn[0:M,1:Lp]+pn[1:Mp,1:Lp]+pn[1:Mp,0:L])
    
    vrt=iA_q*(dn_v[:,1:Lp]-dn_v[:,0:L]-dm_u[1:Mp,:]+dm_u[0:M,:])

    return vrt


""
def vort_3d(u,v,pm,pn):

    [Nz,Mz,Lz]=u.shape
    [Mp,Lp]=pm.shape

    vrt = np.zeros((Nz,Mp-1,Lp-1))

    for iz in range(0, Nz, 1):    
        vrt[iz,:,:]=vort_2d(u[iz,:,:],v[iz,:,:],pm,pn)

    return vrt


#######################################################
# Compute divergence of a 2-D field
# ######################################################
'''
div = dudx + dvdy

div is computed at rho points (first and last points as nan)
'''



def div(u,v,pm,pn):

    [Mp,Lp]=pm.shape
    L=Lp-1
    M=Mp-1
    Lm=L-1
    Mm=M-1


    dudx = np.zeros((Mp,Lp))*np.nan
    dudx[:,1:L] = (u[:,1:L]-u[:,0:Lm])*pm[:,1:L]
    #dudx[:,0] = dudx[:,1]
    #dudx[:,L] = dudx[:,Lm]

    dvdy = np.zeros((Mp,Lp))*np.nan
    dvdy[1:M,:] = (v[1:M,:]-v[0:Mm,:])*pn[1:M,:]
    #dvdy[0,:] = dvdy[1,:]
    #dvdy[M,:] = dvdy[Mm,:]

    var = dudx + dvdy

    return var


def divs(u,v,pm,pn):

    if np.rank(u)>2:

        [Nz,Mz,Lz]=u.shape
        [Mp,Lp]=pm.shape

        divs = np.zeros((Nz,Mp,Lp))*np.nan

        for iz in range(0, Nz, 1):
            divs[iz,:,:]=div(u[iz,:,:],v[iz,:,:],pm,pn)

    else:
    
        divs=div(u,v,pm,pn)

    return divs




#######################################################
# Compute density and Brunt-Vaissala frequency
# ######################################################
'''

compute rho from equation of state

rho on rho (vert and hor) grid

bvf is computed only if z_w is not None
bvf computed on rho-w grid (first and last levels set to 0) 


'''


def rho_eos(Tt,Ts,z_r,g,rho0,z_w=None):

    print('rho calculation...')
    
    if np.ndim(Tt)==2:
        [M,L]=Tt.shape
    else:
        [N,M,L]=Tt.shape

    A00=+19092.56;A01=+209.8925;
    A02=-3.041638;A03=-1.852732e-3;A04=-1.361629e-5;A10=104.4077;
    A11=-6.500517;A12=+0.1553190;A13=2.326469e-4;AS0=-5.587545;
    AS1=+0.7390729;AS2=-1.909078e-2;B00=+4.721788e-1;B01=+1.028859e-2;
    B02=-2.512549e-4;B03=-5.939910e-7;B10=-1.571896e-2;B11=-2.598241e-4;
    B12=+7.267926e-6;BS1=+2.042967e-3;E00=+1.045941e-5;E01=-5.782165e-10;
    E02=+1.296821e-7;E10=-2.595994e-7;E11=-1.248266e-9;E12=-3.508914e-9;

    QR=+999.842594;Q01=+6.793952e-2;Q02=-9.095290e-3;
    Q03=+1.001685e-4;Q04=-1.120083e-6;Q05=+6.536332e-9;Q10=+0.824493;
    Q11=-4.08990e-3;Q12=+7.64380e-5;Q13=-8.24670e-7;Q14=+5.38750e-9;
    QS0=-5.72466e-3;QS1=+1.02270e-4;QS2=-1.65460e-6;Q20=+4.8314e-4;

    sqrtTs=Ts ** 0.5;
    
    K0=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))\
    +Ts*(A10+Tt*(A11+Tt*(A12+Tt*A13))\
    +sqrtTs*(AS0+Tt*(AS1+Tt*AS2)));
    
    K1=B00+Tt*(B01+Tt*(B02+Tt*B03))\
    +Ts*(B10+Tt*(B11+Tt*B12)+sqrtTs*BS1);
    
    K2=E00+Tt*(E01+Tt*E02)\
    +Ts*(E10+Tt*(E11+Tt*E12));
    
    rho1=QR+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))\
    +Ts*(Q10+Tt*(Q11+Tt*(Q12+Tt*(Q13+Tt*Q14)))\
    +sqrtTs*(QS0+Tt*(QS1+Tt*QS2))+Ts*Q20);

    rho=rho1/(1+0.1*z_r/(K0-z_r*(K1-z_r*K2)));

    #######################################################

    if np.any(z_w)!=None:

        print('bvf calculation...')
        bvf=0.*z_w;
        cff=g/rho0;

        bvf[1:N,:]=-cff*(rho1[1:N,:,:]/\
        (1.+0.1*z_w[1:N,:,:]/\
        ( K0[1:N,:,:]-z_w[1:N,:,:]*(K1[1:N,:,:]-z_w[1:N,:,:]*K2[1:N,:,:])))\
        -rho1[0:N-1,:,:]/( 1.+0.1*z_w[1:N,:,:]/\
        ( K0[0:N-1,:,:]-z_w[1:N,:,:]*(K1[0:N-1,:,:]-z_w[1:N,:,:]*K2[0:N-1,:,:]))))\
        /(z_r[1:N,:,:]-z_r[0:N-1,:,:]);


        return [rho,bvf]

    else:

        return rho

#######################################################
# Compute thermal expansion and saline contraction coefficients
# ######################################################

def alphabeta(Tt,Ts,rho0):

    Q01=6.793952E-2; Q02=-9.095290E-3;
    Q03=+1.001685E-4; Q04=-1.120083E-6; Q05=+6.536332E-9;
    U00=+0.824493; U01=-4.08990E-3; U02=+7.64380E-5 ;
    U03=-8.24670E-7; U04=+5.38750E-9; V00=-5.72466E-3 ;
    V01=+1.02270E-4; V02=-1.65460E-6; W00=+4.8314E-4;


    sqrtTs=Ts ** 0.5;
    cff=1/rho0

    alpha=-cff*( Q01+Tt*( 2.*Q02+Tt*( 3.*Q03+Tt*(4.*Q04 +Tt*5.*Q05 )))\
    +Ts*( U01+Tt*( 2.*U02+Tt*(3.*U03 +Tt*4.*U04 ))+sqrtTs*( V01+Tt*2.*V02)))

    beta= cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))\
    +1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2.*W00*Ts )


    return [alpha, beta]








#######################################################
# z-derivative (messy, needs to be rewritten)
# ######################################################

def diffz(var,z):

    if np.rank(var)<3:
        dvardz = diffz_2d(var,z)
    else:
        dvardz = diffz_3d(var,z)

    return dvardz


    #######################

def diffz_2d(var,z):

    [N,M]=var.shape
    dvardz = var*np.nan

    if np.rank(z)==2:
        dvardz[1:,:] = (var[1:,:]-var[:-1,:])/(z[1:,:]-z[:-1,:])
        dvardz[0,:] =  dvardz[1,:]
    else:
        for iz in range(1,N):
            dvardz[iz,:] = (var[iz,:]-var[iz-1,:])/(z[iz]-z[iz-1])
        dvardz[0,:] =  dvardz[1,:]

    return dvardz

    #######################

def diffz_3d(var,z):

    [N,M,L]=var.shape
    dvardz = var*np.nan

    if np.rank(z)==3:
        dvardz[1:,:,:] = (var[1:,:,:]-var[:-1,:,:])/(z[1:,:,:]-z[:-1,:,:])
        dvardz[0,:,:] =  dvardz[1,:,:]
    else:
        for iz in range(1,N):
            dvardz[iz,:,:] = (var[iz,:,:]-var[iz-1,:,:])/(z[iz]-z[iz-1])
        dvardz[0,:,:] =  dvardz[1,:,:]

    return dvardz

#######################################################
# x-derivative from rho-grid to u-grid
# ######################################################

def diffx(var,pm):

    if np.rank(var)<3:
        dvardx = diffx_2d(var,pm)
    else:
        dvardx = diffx_3d(var,pm)

    return dvardx

""
def diffx_3d(var,pm):

    [N,M,L]=var.shape

    dvardx = np.zeros((N,M,L-1))

    for iz in range(0, N):    
        dvardx[iz,:,:]=diffx_2d(var[iz,:,:],pm)

    return dvardx

""
def diffx_2d(var,pm):


    if np.rank(pm)==2: dvardx = (var[:,1:]-var[:,:-1])*0.5*(pm[:,1:]+pm[:,:-1])

    else: dvardx = (var[:,1:]-var[:,:-1])*pm

    return dvardx


#######################################################
# y-derivative from rho-grid to v-grid
# ######################################################

def diffy(var,pn):

    if np.rank(var)<3: dvardy = diffy_2d(var,pn)
    else: dvardy = diffy_3d(var,pn)

    return dvardy

    #######################

def diffy_3d(var,pn):

    [N,M,L]=var.shape
    dvardy = np.zeros((N,M-1,L))
    for iz in range(0, N): dvardy[iz,:,:]=diffy_2d(var[iz,:,:],pn)

    return dvardy

    #######################


def diffy_2d(var,pn):

    if np.rank(pn)==2: dvardy = (var[1:,:]-var[:-1,:])*0.5*(pn[1:,:]+pn[:-1,:])
    else: dvardy = (var[1:,:]-var[:-1,:])*pn

    return dvardy

""
def diffysmooth(var,pn):

    if np.rank(var)<3: dvardy = diffy_2dsmooth(var,pn)
    else: dvardy = diffy_3dsmooth(var,pn)

    return dvardy

    #######################

def diffy_3dsmooth(var,pn):

    [N,M,L]=var.shape
    dvardy = np.zeros((N,M,L))

    for iz in range(0, N): dvardy[iz,:,:]=diffy_2dsmooth(var[iz,:,:],pn)

    return dvardy

    #######################

def diffy_2dsmooth(var,pn):

    [M,L]=var.shape
    dvardy = np.zeros((M,L))

    dvardy[1:-1,:] = (var[2:,:]-var[:-2,:])/2*(pn[1:-1,:])
    dvardy[0,:] = dvardy[1,:]
    dvardy[-1,:] = dvardy[-2,:]

    return dvardy





#######################################################
# Compute horizontal derivatives on sigma-levels (1st order)
# ######################################################
'''
var on rho-rho grid
dvardxi on psi-rho grid
'''

def diffxi(var,pm,z_r,z_w=None,newz=None,mask=None):


    if z_r.shape[0]<=2:
        dvardxi = diffxi_2d(var,pm,z_r,z_w,newz,mask)
    else:
        dvardxi = diffxi_3d(var,pm,z_r,z_w,newz,mask)

    ##############################################

    return dvardxi

#######################################################
# ######################################################


def diffxi_3d(var,pm,z_r,z_w=None,newz=None,mask=None):


    if newz==None: newz = (z_r[:,:,:-1] + z_r[:,:,1:])/2
    else: newz = rho2u(newz)

    dvardxi = np.zeros((var.shape[0],var.shape[1],var.shape[2]-1))

    ##############################################

    varzp = vinterpsF(var[:,:,1:],newz,z_r[:,:,1:],z_w[:,:,1:])
    varzm = vinterpsF(var[:,:,:-1],newz,z_r[:,:,:-1],z_w[:,:,:-1])

    dvardxi = (varzp - varzm )*0.5*(pm[:,:-1]+pm[:,1:])  

    ##############################################

    return dvardxi


#######################################################
# ######################################################


def diffxi_2d(var,pm,z_r,z_w=None,newz=None,mask=None):

    dvardxi = np.zeros((z_r.shape[1],z_r.shape[2]-1))

    ##############################################

    if newz==None: newz = (z_r[0,:,:-1] + z_r[0,:,1:])/2
    else: newz = rho2u(newz)

    dz0 = (z_r[0,:,1:]-newz)
    dz1 = (newz-z_r[1,:,1:])
    varzp = (dz1*var[0,:,1:] + dz0*var[1,:,1:])/(z_r[0,:,1:]-z_r[1,:,1:])

    dz0 = (z_r[0,:,:-1]-newz)
    dz1 = (newz-z_r[1,:,:-1])
    varzm = (dz1*var[0,:,:-1] + dz0*var[1,:,:-1])/(z_r[0,:,:-1]-z_r[1,:,:-1])

    dvardxi = (varzp - varzm )*0.5*(pm[:,:-1]+pm[:,1:])
    ##############################################

    return dvardxi









#######################################################
# Compute horizontal derivatives on sigma-levels (1st order)
# ######################################################

'''
var on rho-rho grid
dvardxi on psi-rho grid
'''

def diffeta(var,pn,z_r,z_w=None,newz=None,mask=None):


    if z_r.shape[0]<=2:
        dvardeta = diffeta_2d(var,pn,z_r,z_w,newz,mask)
    else:
        dvardeta = diffeta_3d(var,pn,z_r,z_w,newz,mask)

    ##############################################

    return dvardeta


#######################################################
# ######################################################


def diffeta_3d(var,pn,z_r,z_w=None,newz=None,mask=None):


    if newz==None: newz = (z_r[:,:-1,:] + z_r[:,1:,:])/2
    else: newz = rho2v(newz)

    dvardeta = np.zeros((var.shape[0],var.shape[1]-1,var.shape[2]))

    ##############################################

    varzp = vinterpsF(var[:,1:,:],newz,z_r[:,1:,:],z_w[:,1:,:])
    varzm = vinterpsF(var[:,:-1,:],newz,z_r[:,:-1,:],z_w[:,:-1,:])


    dvardeta = (varzp - varzm )*0.5*(pn[:-1,:]+pn[1:,:])

    ##############################################


    return dvardeta



#######################################################
# Compute horizontal derivatives on sigma-levels (1st order)
# ######################################################



def diffeta_2d(var,pn,z_r,z_w=None,newz=None,mask=None):

    dvardeta = np.zeros((z_r.shape[1]-1,z_r.shape[2]))

    ##############################################

    if newz==None: newz = (z_r[0,:-1,:] + z_r[0,1:,:])/2
    else: newz = rho2v(newz)

    dz0 = (z_r[0,1:,:]-newz)
    dz1 = (newz-z_r[1,1:,:])
    varzp = (dz1*var[0,1:,:] + dz0*var[1,1:,:])/(z_r[0,1:,:]-z_r[1,1:,:])

    dz0 = (z_r[0,:-1,:]-newz)
    dz1 = (newz-z_r[1,:-1,:])
    varzm = (dz1*var[0,:-1,:] + dz0*var[1,:-1,:])/(z_r[0,:-1,:]-z_r[1,:-1,:])

    dvardeta = (varzp - varzm )*0.5*(pn[:-1,:]+pn[1:,:])

    ##############################################


    return dvardeta





#######################################################
# operator grad.
# ######################################################

def grad(var,pm=1.,pn=1.):


    amplitude = u2rho(diffx(var,pm)**2) + v2rho(diffy(var,pn)**2)

    ##############################################


    return amplitude






#######################################################
# Compute vertical velocity (see Wvlcty.F)
# ######################################################

def getw(u,v,pm,pn,z_r,z_w):

    
    #rho grid = Nrho,Mrho,Lrho
    #u grid = Nrho,Mrho,Lu
    #v grid = Nrho,Mv,Lrho
    #w grid = Nw,Mrho,Lrho

    [Mrho,Lrho]=pm.shape
    Mv=Mrho-1
    Lu=Lrho-1

    if len(u.shape)==3: Nrho=u.shape[0]
    else: Nrho = 1

    Nw=Nrho+1

    ###########################

    flxu = (rho2u(z_w[1:,:,:]-z_w[:-1,:,:]))/(0.5*(pn[:,1:]+pn[:,:-1]))*u
    flxv = (rho2v(z_w[1:,:,:]-z_w[:-1,:,:]))/(0.5*(pm[1:,:]+pm[:-1,:]))*v

    ###########################

    wrk = np.zeros((Nw,Mrho,Lrho))*np.nan
    wvlc = np.zeros((Nw,Mrho,Lrho))*np.nan

    wrk[0,1:-1,1:-1]=0
    wvlc[0,1:-1,1:-1]=0
    for iz in range(1,Nw):
        wrk[iz,1:-1,1:-1] = wrk[iz-1,1:-1,1:-1] - pm[1:-1,1:-1] * pn[1:-1,1:-1] * (flxu[iz-1,1:-1,1:] - flxu[iz-1,1:-1,:-1] + flxv[iz-1,1:,1:-1] - flxv[iz-1,:-1,1:-1])

    if Nw==2:
        wvlc[1,:,:] = wrk[1,:,:]

    else:
        #move to vertical rho points
        wvlc[1,:,:] = -0.125*wrk[2,:,:] + 0.75*wrk[1,:,:] + 0.375 * wrk[0,:,:]
        wvlc[-1,:,:] = 0.375*wrk[-1,:,:] + 0.75*wrk[-2,:,:] - 0.125 * wrk[-3,:,:]
        for iz in range(2,Nw-1):
            wvlc[iz,:,:] = 0.5625*(wrk[iz-1,:,:]+wrk[iz,:,:])-0.0625*(wrk[iz-2,:,:]+wrk[iz+1,:,:])


    #add contributions due to S-coord slopes (u*dz/dx and v*dz/dy)
    Wx = u*(z_r[:,:,1:]-z_r[:,:,:-1])*(pm[:,1:]+pm[:,:-1])
    Wy = v*(z_r[:,1:,:]-z_r[:,:-1,:])*(pn[1:,:]+pn[:-1,:])


    wvlc[1:,1:-1,1:-1] = wvlc[1:,1:-1,1:-1] + 0.25 * (Wy[:,1:,1:-1]+Wy[:,:-1,1:-1]+Wx[:,1:-1,1:]+Wx[:,1:-1,:-1])

    return [wvlc]




###################################################################################
# Bunch of fucntions for interpolation
# ##################################################################################


def transform_grid(jmean,angle,old_Ly,Lya, Lyb=0.):
    '''
    Define transformation matrix for local rotated axes
    (equivalent to interp.griddata but much faster)
    
    In this version  x-dimensions are supposed to be the same
    '''
    
    if Lyb==0.: Lyb=Lya

    xrot=np.arange(len(jmean)); yrot= np.arange(-Lya,Lyb+1,1.)

    Xrot,Yrot = np.meshgrid(xrot,yrot)
    Xrot,Yrot = Xrot.T,Yrot.T

    for ix in range(Xrot.shape[0]): Yrot[ix,:] = Yrot[ix,:] + jmean[ix]

    newXrot,newYrot = 0.*Xrot,0.*Yrot

    for ix in range(Xrot.shape[0]):
        xycenter = [xrot[ix],jmean[ix]]
        for iy in range(len(yrot)):
            [newXrot[ix,iy],newYrot[ix,iy]] = rotate_2d([Xrot[ix,iy],Yrot[ix,iy]],np.array(xycenter),-1*angle[ix])

    x=np.arange(len(jmean)); y=np.arange(old_Ly); X,Y = np.meshgrid(x,y);  X,Y =X.T,Y.T

    [elem,coef] = get_tri_coef(X,Y,newXrot,newYrot)

    return [elem,coef,newXrot,newYrot]




""
from scipy.spatial import Delaunay
import time as tm


def get_tri_coef(X,Y,newX,newY,verbose=0):

    """
    Inputs:
        origin lon and lat 2d arrays (X,Y)
        child lon and lat 2d arrays (newX,newY)

    Ouputs:
        elem - pointers to 2d gridded data (at lonp,latp locations) from
            which the interpolation is computed (3 for each child point)
        coef - linear interpolation coefficients
    Use:
        To subsequently interpolate data from Fp to Fc, the following
        will work:      Fc  = sum(coef.*Fp(elem),3);  This line  should come in place of all
        griddata calls. Since it avoids repeated triangulations and tsearches (that are done
        with every call to griddata) it should be much faster.
    """
    

    Xp = np.array([X.ravel(), Y.ravel()]).T
    Xc = np.array([newX.ravel(), newY.ravel()]).T


    #Compute Delaunay triangulation
    if verbose==1: tstart = tm.time()
    tri = Delaunay(Xp)
    if verbose==1: print ('Delaunay Triangulation', tm.time()-tstart)

    #Compute enclosing simplex and barycentric coordinate (similar to tsearchn in MATLAB)
    npts = Xc.shape[0]
    p = np.zeros((npts,3))

    points = tri.points[tri.vertices[tri.find_simplex(Xc)]]

    if verbose==1: tstart = tm.time()
    for i in range(npts):

        if verbose==1: print (np.float(i)/npts)

        if tri.find_simplex(Xc[i])==-1:  #Point outside triangulation
             p[i,:] = p[i,:] * np.nan

        else:

            if verbose==1: tstart = tm.time()
            A = np.append(np.ones((3,1)),points[i] ,axis=1)
            if verbose==1: print ('append A', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            B = np.append(1., Xc[i])
            if verbose==1: print ('append B', tm.time()-tstart)

            if verbose==1: tstart = tm.time()
            p[i,:] = np.linalg.lstsq(A.T,B.T)[0]
            if verbose==1: print ('solve', tm.time()-tstart)




    if verbose==1: print ('Coef. computation 1', tm.time()-tstart)

    if verbose==1: tstart = tm.time()
    elem = np.reshape(tri.vertices[tri.find_simplex(Xc)],(newX.shape[0],newY.shape[1],3))
    coef = np.reshape(p,(newX.shape[0],newY.shape[1],3))
    if verbose==1: print ('Coef. computation 2', tm.time()-tstart)

    return [elem,coef]


""
def rotate_2d(pts,cnt,ang=np.pi/4):
    '''pts = {} Rotates points(nx2) about center cnt(2) by angle ang(1) in radian'''
    return np.dot(pts-cnt,np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]]))+cnt


""
def rotuv(angle,u,v,**kwargs):
    '''
    Rotate winds or u,v to lat,lon coord -> result on psi grid
    '''

    'rotate vectors by geometric angle'
    urot = u*np.cos(angle) - v*np.sin(angle)
    vrot = u*np.sin(angle) + v*np.cos(angle)

    return [urot,vrot]



""

def local_angle(jmean,imean=None,smoothing=1):
    '''   
    Get angle for local stream following axe
    '''
    
    angle = np.zeros(len(jmean))
    if imean==None: imean = np.arange(len(jmean))
    for i in range(1,len(jmean)-1):
        angle[i] = -1*np.arctan((jmean[i+1]-jmean[i-1])/(imean[i+1]-imean[i-1]))

    angle[0] =angle[1]; angle[-1] = angle[-2]
    if smoothing==1: angle = smooth(angle,min(20,jmean.shape[0]-1),'flat')
    
    return angle
        













# End of romstools.py
