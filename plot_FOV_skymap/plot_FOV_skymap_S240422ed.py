#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 13:22:10 2019

@author: ducoin
"""

"""
This file allow to plot a given skymap and square representing telescopes FOV to compare the size
"""

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import time
import mocmanip
from matplotlib import cm
from matplotlib import colors
import matplotlib
import argparse


def load_cylon():
    __all__ = ()
    
    for name in ['cylon']:
        # Read in color map RGB data.
    
        print("Building the dedicated colormap :",name)
    #    with pkg_resources.resource_stream(__name__,'/home/ducoin/Desktop/MXT_counterpart_research/cylon.csv') as f:
        data = np.loadtxt('cylon.csv', delimiter=',')
    
        # Create color map.
        cmap = colors.LinearSegmentedColormap.from_list(name, data)
        # Assign in module.
        locals().update({name: cmap})
        # Register with Matplotlib.
        cm.register_cmap(cmap=cmap)
        cmap.set_under('w')
        # Generate reversed color map.
        name += '_r'
        data = data[::-1]
        cmap = colors.LinearSegmentedColormap.from_list(name, data)
        # Assign in module.
        locals().update({name: cmap})
        # Register with Matplotlib.
        cm.register_cmap(cmap=cmap)
        cmap.set_under('w')
    return

def getSquarePixels(ra_pointing, dec_pointing, tileSide, nside, alpha = 0.6, color='#FFFFFF',fill=True,facecolor=None, edgecolor='#FFFFFF', return_corners=False, label=None, zorder=3,hatch=None):

    area = tileSide*tileSide

    decCorners = (dec_pointing - tileSide / 2.0, dec_pointing + tileSide / 2.0)
#    print("decCorners =", decCorners)
    radecs = []
    
    #security for the periodic limit conditions
    for d in decCorners:
        if d > 90.:
            d = 180. - d
        elif d < -90.:
            d = -180 - d
        
#        print("decCorners =", decCorners)
        raCorners = (ra_pointing - (tileSide / 2.0) / np.cos(np.deg2rad(d)) , ra_pointing + (tileSide / 2.0) / np.cos(np.deg2rad(d)))
#        print("raCorners =", raCorners)

        #personal proposition
        for r in raCorners:
            if r > 360.:
                r = r - 360.
            elif r < 0.:
                r = 360. + r
            radecs.append([r,d])
        
        
    radecs = np.array(radecs)
    idx1 = np.where(radecs[:,0]>=180.0)[0] 
    idx2 = np.where(radecs[:,0]<180.0)[0]
    idx3 = np.where(radecs[:,0]>300.0)[0]
    idx4 = np.where(radecs[:,0]<60.0)[0]
    if (len(idx1)>0 and len(idx2)>0) and not (len(idx3)>0 and len(idx4)>0):
        alpha = 0.0
    

#    print("radecs =", radecs)    
    
    #personal security for the degenerescence
    a = radecs[0].tolist()
    b = radecs[1].tolist()
    c = radecs[2].tolist()
    d = radecs[3].tolist()
    if a==b or a ==c or a==d or b==c or b==d or c==d:
        radecs[2,0] = radecs[2,0]+180
        if radecs[2,0] > 360:
            radecs[2,0] = 360 - radecs[2,0]
        radecs[3,0] = radecs[3,0]+180
        if radecs[3,0] > 360:
            radecs[3,0] = 360 - radecs[2,0]
            
        if return_corners:
            
            return [],[],[],[],[]
        else:
                
            return [], [], [], []
#    idx1 = np.where(np.abs(radecs[:,1])>=87.0)[0] 
#    if len(idx1) == 4:
#        return [], [], [], []

    idx1 = np.where((radecs[:,1]>=87.0) | (radecs[:,1]<=-87.0))[0]
    if len(idx1)>0:
        radecs = np.delete(radecs, idx1[0], 0)

    if return_corners:
            
        corners = radecs
    
    xyz = []
    for r, d in radecs:
        xyz.append(hp.ang2vec(r, d, lonlat=True))



    npts, junk = radecs.shape
    if npts == 4:
        xyz = [xyz[0], xyz[1],xyz[3], xyz[2]]

#        ipix = hp.query_polygon(nside, np.array(xyz))
        ipix = hp.query_polygon(nside, np.array(xyz), inclusive=True)


    else:

        ipix = hp.query_polygon(nside, np.array(xyz))

    #idx1 = np.where((radecs[:,1]>=70.0) | (radecs[:,1]<=-70.0))[0]
    #idx2 = np.where((radecs[:,0]>300.0) | (radecs[:,0]<60.0))[0]
    #if (len(idx1) == 0) or (len(idx2) > 0):
    #    return [], [], [], []

    xyz = np.array(xyz)
    proj = hp.projector.MollweideProj(rot=None, coord=None) 
    x,y = proj.vec2xy(xyz[:,0],xyz[:,1],xyz[:,2])
    xy = np.zeros(radecs.shape)
    xy[:,0] = x
    xy[:,1] = y
    path = matplotlib.path.Path(xy)
    patch = matplotlib.patches.PathPatch(path, alpha=alpha, color=color,facecolor=facecolor, fill=fill, zorder=zorder, edgecolor=edgecolor, linewidth=1, label=label,hatch=hatch,)

    if return_corners:
        
        return ipix, radecs, patch, area, corners
        
    return ipix, radecs, patch, area

def getSquarePixels1(ra_pointing, dec_pointing, tileSide, tileSide1, nside, alpha = 0.6, color='#FFFFFF',fill=True,facecolor=None, edgecolor='#FFFFFF', return_corners=False, label=None, zorder=3,hatch=None):

    area = tileSide*tileSide1

    decCorners = (dec_pointing - tileSide1 / 2.0, dec_pointing + tileSide1 / 2.0)
#    print("decCorners =", decCorners)
    radecs = []
    
    #security for the periodic limit conditions
    for d in decCorners:
        if d > 90.:
            d = 180. - d
        elif d < -90.:
            d = -180 - d
        
#        print("decCorners =", decCorners)
        raCorners = (ra_pointing - (tileSide / 2.0) / np.cos(np.deg2rad(d)) , ra_pointing + (tileSide / 2.0) / np.cos(np.deg2rad(d)))
#        print("raCorners =", raCorners)

        #personal proposition
        for r in raCorners:
            if r > 360.:
                r = r - 360.
            elif r < 0.:
                r = 360. + r
            radecs.append([r,d])
        
        
    radecs = np.array(radecs)
    idx1 = np.where(radecs[:,0]>=180.0)[0] 
    idx2 = np.where(radecs[:,0]<180.0)[0]
    idx3 = np.where(radecs[:,0]>300.0)[0]
    idx4 = np.where(radecs[:,0]<60.0)[0]
    if (len(idx1)>0 and len(idx2)>0) and not (len(idx3)>0 and len(idx4)>0):
        alpha = 0.0
    

#    print("radecs =", radecs)    
    
    #personal security for the degenerescence
    a = radecs[0].tolist()
    b = radecs[1].tolist()
    c = radecs[2].tolist()
    d = radecs[3].tolist()
    if a==b or a ==c or a==d or b==c or b==d or c==d:
        radecs[2,0] = radecs[2,0]+180
        if radecs[2,0] > 360:
            radecs[2,0] = 360 - radecs[2,0]
        radecs[3,0] = radecs[3,0]+180
        if radecs[3,0] > 360:
            radecs[3,0] = 360 - radecs[2,0]
            
        if return_corners:
            
            return [],[],[],[],[]
        else:
                
            return [], [], [], []
#    idx1 = np.where(np.abs(radecs[:,1])>=87.0)[0] 
#    if len(idx1) == 4:
#        return [], [], [], []

    idx1 = np.where((radecs[:,1]>=87.0) | (radecs[:,1]<=-87.0))[0]
    if len(idx1)>0:
        radecs = np.delete(radecs, idx1[0], 0)

    if return_corners:
            
        corners = radecs
    
    xyz = []
    for r, d in radecs:
        xyz.append(hp.ang2vec(r, d, lonlat=True))



    npts, junk = radecs.shape
    if npts == 4:
        xyz = [xyz[0], xyz[1],xyz[3], xyz[2]]

#        ipix = hp.query_polygon(nside, np.array(xyz))
        ipix = hp.query_polygon(nside, np.array(xyz), inclusive=True)


    else:

        ipix = hp.query_polygon(nside, np.array(xyz))

    #idx1 = np.where((radecs[:,1]>=70.0) | (radecs[:,1]<=-70.0))[0]
    #idx2 = np.where((radecs[:,0]>300.0) | (radecs[:,0]<60.0))[0]
    #if (len(idx1) == 0) or (len(idx2) > 0):
    #    return [], [], [], []

    xyz = np.array(xyz)
    proj = hp.projector.MollweideProj(rot=None, coord=None) 
    x,y = proj.vec2xy(xyz[:,0],xyz[:,1],xyz[:,2])
    xy = np.zeros(radecs.shape)
    xy[:,0] = x
    xy[:,1] = y
    path = matplotlib.path.Path(xy)
    patch = matplotlib.patches.PathPatch(path, alpha=alpha, color=color,facecolor=facecolor, fill=fill, zorder=zorder, edgecolor=edgecolor, linewidth=1, label=label,hatch=hatch,)

    if return_corners:
        
        return ipix, radecs, patch, area, corners
        
    return ipix, radecs, patch, area


    
def ipixs_in_percentage(skymap, percentage):
    """
    Finding ipix indices confined in a given percentage.
        
    Input parameters
    ----------------
    skymap : list
        given skymap from a GW, hp.read_map() 
    percentage : float
        float inside [0,1], percentage of skymap probability
        
    """
    
    start= time.time()
    
    sorted_skymap = skymap[skymap.argsort()]
    index_list = np.array([i for i in range(0,len(skymap))])
    index_list = index_list[skymap.argsort()]        
    
    sorted_skymap = sorted_skymap.tolist()
    index_list = index_list.tolist()
    
    sorted_skymap.reverse()
    index_list.reverse()
    
    if percentage == 1:
        
        return index_list
    
    else:
    
        sum_proba = 0   
        ipixs_in_percentage = []    
        
        i=0
        
        while sum_proba < percentage :
        
            ipixs_in_percentage += [index_list[i]]
            
            sum_proba += sorted_skymap[i]
            
            i += 1
        end = time.time()
        print("running time ipix_in_percentage =", end-start)
        return ipixs_in_percentage


def healpix_skymap(file_skymap, proba):
    """
    This function open a skymap using healpix and convert it to an list of pixels.    
    Return the maps as a list, the proba as a list, the nside, distmu as a list and distsigma as a list.
  
    Input parameters
    ----------------
    file_skymap : str
        path of the skymap .fits file    
    proba : float
        the probability if the skymap (ex: 0.9 = 90% skymap opened)
    """
    #initialize internal object
    moc_map = mocmanip.MOC_confidence_region()
    
    #load the skymap probability
    print("running skymap : {}".format(file_skymap))
    skymap, header = hp.read_map(file_skymap, h=True, verbose=False)
    print(header)
    #load also the skymap distance info distmu, distsigma and distnorm (read the read_me;txt file for more informations)
    skymap_distmu,_ = hp.read_map(file_skymap, field=1, h=True, verbose=False)
    skymap_distsigma,_ = hp.read_map(file_skymap, field=2, h=True, verbose=False)
    skymap_distnorm,_ = hp.read_map(file_skymap, field=3, h=True, verbose=False)

    """
    if this bug because no 3D map. Use something like 
    
    
    if is3D:
                    healpix_data, header = hp.read_map(filename, field=(0,1,2,3), verbose=False,h=True)
                    distmu_data = healpix_data[1]
                    distsigma_data = healpix_data[2]
                    prob_data = healpix_data[0]
                    norm_data = healpix_data[3]
                    map_struct["distmu"] = distmu_data / params["DScale"]
                    map_struct["distsigma"] = distsigma_data / params["DScale"]
                    map_struct["prob"] = prob_data
                    map_struct["distnorm"] = norm_data
                else:
                    prob_data, header = hp.read_map(filename, field=0, verbose=False,h=True)
                    prob_data = prob_data / np.sum(prob_data)
                    map_struct["prob"] = prob_data
    """

    
    #initialize internal objects
    moc_map.hpx = skymap
    moc_map.nside = hp.npix2nside(len(moc_map.hpx))

    #find pixels inside the prob_percentage contour
    #return pixels ordered by probability (from max to min)
    moc_map.ipixs = ipixs_in_percentage(skymap, proba)

    #we have to select and order skymap_distmu, skymap_distsigma and skymap_distnorm
    #in the same way of the result of ipixs_in_percentage()
    ordered_skymap_distmu = []
    ordered_skymap_distsigma = []
    ordered_skymap_distnorm = []
    for i in range(0,len(moc_map.ipixs)):
    
        ordered_skymap_distmu += [skymap_distmu[moc_map.ipixs[i]]]
        ordered_skymap_distsigma += [skymap_distsigma[moc_map.ipixs[i]]]
        ordered_skymap_distnorm += [skymap_distnorm[moc_map.ipixs[i]]]

    return moc_map.ipixs, moc_map.hpx, moc_map.nside, ordered_skymap_distmu, ordered_skymap_distsigma, ordered_skymap_distnorm


def add_edges(rotation=[]):

    hp.graticule(verbose=False)
    plt.grid(True)
    lons = np.arange(0.0,360,30.0)
    lats = np.zeros(lons.shape)
    
    if rotation != []:
        
        written_lons = lons - rotation[1]

        for i in range(len(written_lons)):
            
            if written_lons[i]>360:
                written_lons[i] = written_lons[i]-360
            elif written_lons[i]<0:
                written_lons[i] = 360 + written_lons[i]
    
    l=0
    for lon, lat in zip(lons,lats):
        if rotation != []:

            hp.projtext(lon-1,lat+1,"%.0f"%written_lons[l],lonlat=True,fontsize=14)
            l+=1
        else:
            hp.projtext(lon,lat,"%.0f"%lon,lonlat=True,fontsize=14)
    lats = np.arange(-60.0,90,30.0)
    lons = np.zeros(lons.shape)
    for lon, lat in zip(lons,lats):
        if lat !=0:
            hp.projtext(lon-1,lat+1,"%.0f"%lat,lonlat=True,fontsize=14)
    plt.text(-2.24,0,'RA (deg)',fontsize=14)
    plt.text(0,1.01,'Dec (deg)',fontsize=14)


def rotate_map(hmap, rot_theta, rot_phi):
    """
    Take hmap (a healpix map array) and return another healpix map array 
    which is ordered such that it has been rotated in (theta, phi) by the 
    amounts given.
    """
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t,p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside))) #theta, phi

    # Define a rotator
    r = hp.Rotator(deg=True, rot=[rot_phi,rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t,p)

    # Interpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot)

    return rot_map


def plot_FOV(file_skymap,title,center_square,fov,output_path,rotation_ra=170,color='lightgreen'):
    """
    This function produce the plot of observed tiles for a given skymap.
    You can apply a rotation in RA of the skymap, by default it fit with the 
    GraceDB convention.
  
    Input parameters
    ----------------
    file_skymap : str
        path of the skymap .fits file    
    title : str
        Title of the plot, can be None if you don't want title
    center_square : list
        list of the tiles center, like : [[ra1,dec1],[ra2,dec2]...]
    fov : float
        field of view of the tiles, it is the size of a side of the square,
        if your field of view is 1.4x1.4 put here 1.4
    output_path : str
        complet path of the output file, ex: /home/dupont/Desktop/super_plot.png
    rotation_ra : float
        rotation you want to apply to the skymap, default fit with the GraceDB
        convention.
    """
    print("ploting the observation skymap...")
    #load the csv containing the color params
    load_cylon()
    
    #load skymap
    skymap, skymap_proba, nside, _, _, _ = healpix_skymap(file_skymap, proba=0.9)
    
    #load the rotation on ra axis
    rotation_angles = [0,rotation_ra]# BE CAREFUL ONLY ROTATION IN RA WILL WORK WITH THE CURRENT VERSION (RA is on the right)
    
    #plot the proba skymap with the rotation
    hp.mollview(rotate_map(skymap_proba, rotation_angles[0],rotation_angles[1]),title=None,  cbar=False, cmap='cylon', unit='skymap probability', hold=True)
    ax = plt.gca()
    #add edges with rotation taked into account
    add_edges(rotation_angles)
    
    if title is not None:
        plt.title(title, pad=30, fontsize=18)
    
    #define the tile side
    tileSide = [fov] * len(center_square)
    
    #define the color for the tiles as lightgreen
    colors = [color] * len(center_square)
    
    for i in range(len(tileSide)):
        
        _,_,patch,_ = getSquarePixels(center_square[i][0]+rotation_angles[1], center_square[i][1], tileSide[i],\
                                      nside, alpha = 0.8, color=colors[i], edgecolor='red',\
                                      return_corners=False,)
   
        hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)    

    #need to plot forcing the resolution to be device independent    
    fig = plt.gcf() # get current figure

    
    DPI = fig.get_dpi()
    fig.set_size_inches(1920.0/float(DPI),1080.0/float(DPI))

    plt.savefig(output_path)
 
    return


################################### MAIN ######################################
file_skymap = '/home/ducoin/gwemopt-mxt.git/Bilby_S240422ed.gz,0'
title = None
output_path = '/home/ducoin/Desktop/S240422ed_ITM_ASTEP.svg'
#TCA
center_square_1 = [[21.137  , 57.548], [17.772  , 57.548], [23.450  , 60.784],\
                   [13.258  , 53.837], [12.289  , 50.125], [17.483  , 50.125],\
                   [20.656  , 51.981], [23.665  , 53.837], [25.656  , 59.404], [11.859  , 48.270],\
                   [25.966  , 57.548], [14.408  , 57.548], [16.335  , 53.837], [20.236  , 55.692],\
                   [19.412  , 53.837], [18.597  , 59.404], [17.021  , 55.692], [14.586  , 48.270],\
                   [14.754  , 51.981], [24.764  , 55.692]]



#TCH
center_square_2 = [[28.310  , 14.586]]



#TRE

center_square_3 = [[41.019  , 34.827], [238.671 , -34.718], [181.050 , -84.859], [238.950 , -38.809],\
                   [182.100 , -38.809], [187.555 , -38.809], [235.866 , -26.536], [34.827  , 30.736],\
                   [36.154  , 34.827], [38.100  , 26.645], [37.222  , 22.555]]



#FRAM A  FRAM-Auger
center_square_6 = [[31.705  , 18.000], [31.886  , 18.973], [32.914  , 18.973], [35.581  , 21.892],\
                   [36.416  , 20.919], [33.943  , 18.973], [34.971  , 18.973], [32.727  , 18.000],\
                   [35.376  , 20.919], [31.525  , 17.027], [32.542  , 17.027], [35.172  , 19.946],\
                   [36.628  , 21.892], [28.156  , 15.081], [30.168  , 15.081]]


#FRAM CTA
center_square_7 = [[16.074  , 54.837],[17.005  , 55.273],[19.382  , 57.021],\
                   [16.829  , 54.837],[18.926  , 56.147],[15.009  , 53.963],\
                   [18.584  , 57.021],[17.584  , 54.837],[15.749  , 53.963],\
                   [18.147  , 56.147],[16.488  , 53.963],[17.368  , 56.147],\
                   [16.242  , 55.273],[19.706  , 56.147],[16.589  , 56.147],\
                   [17.767  , 55.273],[17.227  , 53.963],[15.319  , 54.612],\
                   [18.338  , 54.837],[18.530  , 55.273]]

#OAJ
center_square_4 = [[16.471  , 54.545],[17.027  , 55.909],[20.140  , 57.273],\
                   [19.459  , 55.909],[15.949  , 53.182],[18.823  , 54.545],\
                   [17.622  , 57.273],[14.118  , 54.545],[13.671  , 53.182],\
                   [18.228  , 53.182],[15.460  , 51.818],[21.022  , 58.636],\
                   [22.657  , 57.273]]    
    
    
fov_1 = 1.9
fov_2 = 1.9
fov_3 = 4.2
fov_4 = 1.4
#fov_5 = 0.27
fov_6 = 1
fov_7 = 0.45



#load the csv containing the color params
load_cylon()

#load skymap
skymap, skymap_proba, nside, _, _, _ = healpix_skymap(file_skymap, proba=0.9)

#load the rotation on ra axis
rotation_angles = [0,170]# BE CAREFUL ONLY ROTATION IN RA WILL WORK WITH THE CURRENT VERSION (RA is on the right)

#plot the proba skymap with the rotation
hp.mollview(rotate_map(skymap_proba, rotation_angles[0],rotation_angles[1]),title=None,  cbar=False, cmap='cylon', unit='skymap probability', hold=True)
ax = plt.gca()
#add edges with rotation taked into account
add_edges(rotation_angles)

if title is not None:
    plt.title(title, pad=30, fontsize=18)

#define the tile side
tileSide_1 = [fov_1] * len(center_square_1)
tileSide_2 = [fov_2] * len(center_square_2)
tileSide_3 = [fov_3] * len(center_square_3)
tileSide_4 = [fov_4] * len(center_square_4)
#tileSide_5 = [fov_5] * len(center_square_4)
tileSide_6 = [fov_6] * len(center_square_6)
tileSide_7 = [fov_7] * len(center_square_7)

#define the color for the tiles as lightgreen
colors_tiles_1 = ["green"] * len(center_square_1)
colors_tiles_2 = ["deepskyblue"] * len(center_square_2)
colors_tiles_3 = ["yellow"] * len(center_square_3)
colors_tiles_4 = ["lime"] * len(center_square_4)
colors_tiles_6 = ["purple"] * len(center_square_6)
colors_tiles_7 = ["lightgrey"] * len(center_square_7)


"""
#TRE
for i in range(len(tileSide_3)):

    if i==0:
        label='TAROT Reunion'
    else:
        label=None
        
    _,_,patch,_ = getSquarePixels(center_square_3[i][0]+rotation_angles[1], center_square_3[i][1], tileSide_3[i],\
                                  nside, alpha = 0.8, color=colors_tiles_3[i], edgecolor=colors_tiles_3[i],\
                                  return_corners=False,label=label,zorder=1,facecolor=None)
   
    hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)   


#TCH
for i in range(len(tileSide_2)):
    
    
    if i==0:
        label='TAROT Chile'
    else:
        label=None

    _,_,patch,_ = getSquarePixels(center_square_2[i][0]+rotation_angles[1], center_square_2[i][1], tileSide_2[i],\
                                  nside, alpha = 0.8, color=colors_tiles_2[i], edgecolor=colors_tiles_2[i],\
                                  return_corners=False,label=label)
   
    hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)   
   


#TCA
for i in range(len(tileSide_1)):
    
    if i==1:
        label='TAROT Calern'
    else:
        label=None
    _,_,patch,_ = getSquarePixels(center_square_1[i][0]+rotation_angles[1], center_square_1[i][1], tileSide_1[i],\
                                  nside, alpha = 0.8, color=colors_tiles_1[i],label=label, edgecolor=colors_tiles_1[i],\
                                  return_corners=False,fill=True,)
   
    hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)    

#OAJ
for i in range(len(tileSide_4)):
    if i==0:
        label='OAJ'
    else:
        label=None
        
    _,_,patch,_ = getSquarePixels1(center_square_4[i][0]+rotation_angles[1], center_square_4[i][1], tileSide_4[i], tileSide_4[i], \
                                  nside, alpha = 0.8, color=colors_tiles_4[i], edgecolor=colors_tiles_4[i],\
                                  return_corners=False,label=label,facecolor=None)
   
    hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)   

#FRAM A  FZU-Auger
for i in range(len(tileSide_6)):
    
    if i==1:
        label='FRAM-Auger'
    else:
        label=None
    _,_,patch,_ = getSquarePixels(center_square_6[i][0]+rotation_angles[1], center_square_6[i][1], tileSide_6[i],\
                                  nside, alpha = 0.8, color=colors_tiles_6[i],label=label, edgecolor=colors_tiles_6[i],\
                                  return_corners=False,fill=True,)
   
    hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)    


#FRAM-CTA-N
for i in range(len(tileSide_7)):
    
    if i==1:
        label='FRAM-CTA-N'
    else:
        label=None
    _,_,patch,_ = getSquarePixels(center_square_7[i][0]+rotation_angles[1], center_square_7[i][1], tileSide_7[i],\
                                  nside, alpha = 0.8, color=colors_tiles_7[i],label=label, edgecolor=colors_tiles_7[i],\
                                  return_corners=False,fill=True,)
   
    hp.projaxes.HpxMollweideAxes.add_patch(ax,patch)    


"""





#ITM
ra_decs1 = [[121.806 , -31.958],\
           [122.729 , -30.858],[124.189 , -20.603],[121.461 , -25.980],[123.465, -29.125],\
           [122.773 , -30.875],[124.117 , -27.806],[124.879 , -19.457],[124.731 , -30.241],\
           [120.772 , -26.200]]

    
#ASTEP
ra_decs2 = [[122.830 , -25.683],[125.606 ,-27.077],[124.765 ,-37.066],\
           [126.580 , -20.848],[123.144 , -20.582],[124.243 , -25.870],[123.580, -32.388],\
           [123.504 , -24.646],[123.960 , -22.806],[122.261 , -24.451]]
    
#dummie
ra_decs3 = [[17.63 ,55.57],[15.20 ,  53.98],[14.73 ,  52.61]]
    
    

#dummie
ra_decs4 = [[14.865   , 52.3671],[52.591 , 143.31],[52.598 , 176.99],\
           [14.291 , 52.368],[14.755 , 52.709],[14.121 , 52.408],[14.746 , 52.649],\
           [14.726 , 52.610],[52.6097, 141.38],[14.750 , 52.587],[14.073 , 52.515],\
           [14.746 , 52.649],[19.863 , 58.521],[19.858 , 58.407],[19.717 , 58.478],\
           [19.534 , 58.567],[19.954 , 58.628],[15.442 , 55.949],[15.506 , 55.974],\
           [22.541 , 57.846],[14.865 , 52.367],[20.290 , 57.658],[17.654 , 55.603],\
           [21.124 , 56.182],[21.005 , 56.165],[22.157 , 57.924],[26.311 , 59.817],\
           [14.865 , 52.367],[17.907 , 55.439]]
    


 


for i in range(len(ra_decs1)):
    
    if i ==0:
        label='ITM'
    else:
        label = None
    hp.projscatter(ra_decs1[i][0]+rotation_angles[1], ra_decs1[i][1], lonlat=True, c='fuchsia', cmap='cylon', marker='*', alpha=1, s=100,zorder=4,label=label)

 



for i in range(len(ra_decs2)):
    
    if i ==0:
        label='ASTEP'
    else:
        label = None
    hp.projscatter(ra_decs2[i][0]+rotation_angles[1], ra_decs2[i][1], lonlat=True, c='w', cmap='cylon', marker='*', alpha=1, s=100,label=label,zorder=4)

"""
for i in range(len(ra_decs3)):
    
    if i ==0:
        label='VIRT'
    else:
        label = None
    hp.projscatter(ra_decs3[i][0]+rotation_angles[1], ra_decs3[i][1], lonlat=True, c="dodgerblue", cmap='cylon', marker='*', alpha=1, s=100,label=label,zorder=4)
"""
"""
for i in range(len(ra_decs4)):
    
    if i ==0:
        label='kilonova-catcher'
    else:
        label = None
    hp.projscatter(ra_decs4[i][0]+rotation_angles[1], ra_decs4[i][1], lonlat=True, c="olive", cmap='cylon', marker='*', alpha=1, s=100,label=label,zorder=4)
"""


#VIRT OT
ra_decs_OT1 = [[18.34020833, 49.64797222], [27.18966667, 51.43008333]]

"""
for i in range(len(ra_decs_OT1)):
    
    if i ==0:
        label='VIRT OT'
    else:
        label = None
    hp.projscatter(ra_decs_OT1[i][0]+rotation_angles[1], ra_decs_OT1[i][1], lonlat=True,label=label, c="dodgerblue", cmap='cylon', marker='P', alpha=1, s=80,zorder=4)
"""

"""
#Lijiang 2.4-m
ra_decs_OT2 = [[262.791487,-8.450722],[258.341454,-9.964466]]

for i in range(len(ra_decs_OT2)):
    
    if i ==0:
        label='Lijiang 2.4m'
    else:
        label = None
    hp.projscatter(ra_decs_OT2[i][0]+rotation_angles[1], ra_decs_OT2[i][1], lonlat=True, c='fuchsia', cmap='cylon', marker='P', alpha=1, s=80,zorder=4)


"""
plt.legend(facecolor='coral', prop={'size': 20})
#plt.xlim(-0.3,0.1)
#plt.ylim(-0.55,-0.18)
#need to plot forcing the resolution to be device independent    
fig = plt.gcf() # get current figure


DPI = fig.get_dpi()
fig.set_size_inches(1920.0/float(DPI),1080.0/float(DPI))
#fig.set_size_inches(1320.0/float(DPI),1080.0/float(DPI))
plt.savefig(output_path)









"""
center_square = [[0.0,84.5454], [18.9474,84.5454], [14.4,82.7273], [0.0,86.3636],[311.7526,60.9091], [310.5882,59.0909],\
                 [308.1356,53.6364],[307.3171,51.8182],[314.1177,59.0909],[306.5625,50.0],[305.8647,48.1818],\
                 [308.5714,48.1818],[306.7606,44.5454],[309.375,50.0],[309.375,50.0]]
"""
"""
plot_FOV(file_skymap,title,center_square,fov,output_path,rotation_ra=170,color='blue')
"""
"""
0.0    84.5454    0.0108
18.9474    84.5454    0.0094
14.4    82.7273    0.0079
0.0    86.3636    0.0066
311.7526    60.9091    0.0057
310.5882    59.0909    0.0056
308.1356    53.6364    0.0054
307.3171    51.8182    0.0054
314.1177    59.0909    0.0054
306.5625    50.0    0.0053
305.8647    48.1818    0.0053
308.5714    48.1818    0.0052
306.7606    44.5454    0.0051
309.375    50.0    0.005
0.0    82.7273    0.005
"""


"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='GRANDMA/SVOM Circular Reporting System.')

    parser.add_argument('--skymap',
                        dest='skymap',
                        required=True,
                        metavar='skymap_path',
                        type=str,
                        help='path to the skymap you want to plot')

    parser.add_argument('--title',
                        dest='title',
                        required=False,
                        metavar='title',
                        type=str,
                        help="Title of the plot, can be None if you don't want\
                         title")

    parser.add_argument('--output_path',
                        dest='output_path',
                        required=True,
                        metavar='output_path',
                        type=str,
                        help='complet path of the output file, \
                        ex: /home/dupont/Desktop/super_plot.png')

    parser.add_argument('--center_square',
                        dest='center_square',
                        required=True,
                        metavar='center_square',
                        type=float,
                        nargs='+',
                        action='append',
#                        default=[],
                        help='RA-Dec of the squares center, the format after\
                        parse should be [[ra1,dec1],[ra2,dec2]...] so use \
                        something like --center_square 0 5 --center_square 0 170')

    parser.add_argument('--fov',
                        dest='fov',
                        required=True,
                        metavar='fov',
                        type=float,
                        help='field of view of the tiles, it is the size \
                        of a side of the square, if your field of view is \
                        1.4x1.4 put here 1.4')
    
    parser.add_argument('--rotation_ra',
                        dest='rotation_ra',
                        required=False,
                        metavar='rotation_ra',
                        type=float,
                        help='rotation you want to apply to the skymap, \
                        default fit with the GraceDB convention.')


    args = parser.parse_args()
    
    if args.rotation_ra is None:
    
        plot_FOV(args.skymap,args.title,args.center_square,args.fov,args.output_path)
    else :
        
        plot_FOV(args.skymap,args.title,args.center_square,args.fov,args.output_path,args.rotation_ra)
"""  
