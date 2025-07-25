#!/usr/bin/python

import math
import h5py
import pickle
from shapely.geometry import Polygon
import matplotlib
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
matplotlib.rcParams.update({'font.size': 14})
from scipy.interpolate import make_smoothing_spline
from scipy.interpolate import splrep, BSpline, splev
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interpolate as interp
import os, sys
import glob
from matplotlib import ticker as mticker
import matplotlib
from matplotlib.patches import Rectangle
#matplotlib.rcParams.update({'font.size': 10})
import matplotlib as mpl
import requests
import json
import pandas as pd
from astropy.time import Time
from scipy import interpolate
import astropy_healpix as ah
import pandas as pd
import healpy as hp
import ligo.skymap.distance as ligodist
import ligo.skymap.plot
import lxml.etree
import numpy as np
import requests
import regions
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from ligo.gracedb.rest import GraceDb
from ligo.skymap import moc
from ligo.skymap.bayestar import rasterize
from ligo.skymap.io import read_sky_map
from matplotlib import pyplot as plt
from mocpy import MOC
from scipy.interpolate import PchipInterpolator
from scipy.stats import norm
from astropy.table import Table
from astropy.coordinates import Angle
from gwemopt.chipgaps import get_decam_quadrant_moc, get_ztf_quadrant_moc
from regions import PolygonSkyRegion

class S240422ed_observations():
    def __init__(self, light_curves, tmin, tmax, skymap):
        self.tmin = tmin
        self.tmax = tmax
        self.light_curves = light_curves
        self.skymap = skymap
        self.T0 = '2024-04-22T21:35:13.411201'
        self.event = 'S240422ed'
        self.angle = [0.00, 25.84, 36.87, 45.57, 53.13, 60.00, 66.42, 72.54, 78.46, 84.26, 90.00]
        self.moc_order = 10

        BASE = "https://treasuremap.space/api/v1"
        TARGET = "pointings"
        API_TOKEN = "3AkLCduZwOAaW_1-5uFbMa_j4ULIXNwnFyQZsA"

        json_params = {
            "api_token": API_TOKEN,
            "status": "completed",
            "graceid": self.event,
            "depth_unit": "ab_mag"
        }

        url = "{}/{}".format(BASE, TARGET)

        r = requests.get(url=url, json=json_params)
        t = []
        depth = []
        ra = []
        dec = []
        instrumentID = []
        filters = []

        for i in range(len(json.loads(r.text))):
            pointing = json.loads(r.text)[i]
            #print(pointing)
            filters.append(pointing['band'])
            t.append(pointing['time'])
            depth.append(pointing['depth'])
            instrumentID.append(int(pointing['instrumentid']))
            ra.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[1]))
            dec.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[2]))

        t = Time(t, format='isot', scale='utc')
        t_mjd = t.to_value('mjd')
        T02 = Time(self.T0, format='isot', scale='utc')
        T0_mjd = T02.to_value('mjd')
        depth = np.array(depth)
        t_from_T0 = np.array(t_mjd - T0_mjd)
        ra = np.array(ra)
        dec = np.array(dec)
        instrumentID = np.array(instrumentID)
        filters = np.array(filters)

        i1 = np.flatnonzero(instrumentID != 38)
        i2 = np.flatnonzero(filters != 'r')
        idx = np.intersect1d(i1,i2)

        #print(len(filters)-len(idx))
        t_from_T0 = t_from_T0[idx]
        depth = depth[idx]
        instrumentID = instrumentID[idx]
        ra = ra[idx]
        dec = dec[idx]
        filters = filters[idx]

        i3 = np.flatnonzero(np.array(instrumentID) != 38)
        i4 = np.flatnonzero(filters != 'z')
        idx2 = np.intersect1d(i3, i4)
        #print(len(idx2))
        t_from_T0 = t_from_T0[idx2]
        depth = depth[idx2]
        instrumentID = np.array(instrumentID)[idx2]
        ra = np.array(ra)[idx2]
        dec = np.array(dec)[idx2]
        filters = filters[idx2]

        df = pd.read_csv('atlas_exposures.csv', index_col=None)
        times_ATLAS = df['mjd_t0'].to_numpy()
        filter_ATLAS = df['filter'].to_numpy()
        ra_ATLAS = df['raDeg'].to_numpy()
        dec_ATLAS = df['decDeg'].to_numpy()
        mag_ATLAS = df['mag5sig'].to_numpy()
        instrumentID_ATLAS = 6000 * np.ones(len(mag_ATLAS))

        ra = np.append(ra, ra_ATLAS)
        dec = np.append(dec, dec_ATLAS)
        t_from_T0 = np.append(t_from_T0, times_ATLAS)
        depth = np.append(depth, mag_ATLAS)
        instrumentID = np.append(instrumentID, instrumentID_ATLAS)
        filters = np.append(filters, filter_ATLAS)

        observations_Decam = '/Users/Virgo/Desktop/results_decam.csv'
        f = pd.read_csv(observations_Decam, index_col=None)

        t_decam = []
        times = list(f['date_obs'])
        for i in range(len(times)):
            time_decam = times[i]
            obs_t = time_decam.split('+')[0]
            t_decam.append(obs_t)
        depth_decam = f['depth'].to_numpy()
        ra_decam = f['ra'].to_numpy()
        dec_decam = f['dec'].to_numpy()
        band_decam = list(f['ifilter'])

        t_decam = Time(t_decam, format='isot', scale='utc')
        t_mjd_decam = t_decam.to_value('mjd')
        #depth = np.append()
        t_from_T0_decam = np.array(t_mjd_decam - T0_mjd)

        indices_notnan = []
        for i in range(len(depth_decam)):
            if np.isnan(depth_decam[i]) == False:
                indices_notnan.append(i)

        t_from_T0_decam = t_from_T0_decam[indices_notnan]
        depth_decam = depth_decam[indices_notnan]
        ra_decam = ra_decam[indices_notnan]
        dec_decam = dec_decam[indices_notnan]
        instrumentID_decam = 38 * np.ones(len(ra_decam))
        print(np.max(depth_decam))

        filter_decam = np.array([])
        for i in range(len(band_decam)):
            list_band = band_decam[i].split(' ')
            if list_band[0] == 'r':
                filter_decam = np.append(filter_decam, 'r')
            elif list_band[0] == 'z':
                #if 'z DECam' in band_decam[i]:
                filter_decam = np.append(filter_decam, 'z')
            else:
                print('problem filter not covered: ', list_band)

        ra = np.append(ra, ra_decam)
        dec = np.append(dec, dec_decam)
        t_from_T0 = np.append(t_from_T0, t_from_T0_decam)
        depth = np.append(depth, depth_decam)
        instrumentID = np.append(instrumentID, instrumentID_decam)
        filters = np.append(filters, filter_decam)


        df = pd.read_csv('other-Obs-S240422ed.csv', index_col=None)
        time_other = df['Time'].to_numpy()
        filters_other = df['Filter'].to_numpy()
        instrument_other = df['Instrument'].to_numpy()
        ra_other = df['RA'].to_numpy()
        dec_other = df['DEC'].to_numpy()
        mag_other = df['Mag'].to_numpy()

        idx1 = np.flatnonzero(filters_other == 'sdssr')
        filters_other[idx1] = 'r'
        idx2 =  np.flatnonzero(filters == 'sdssg')
        filters_other[idx2] = 'g'

        t_from_T0_other = np.array([])
        for t in time_other:
            if ':' in t:
                t_other = t.replace(' ', 'T')
                t_other = Time(t_other, format='isot', scale='utc')
                t_mjd_other = t_other.to_value('mjd')
                t_from_T0_other = np.append(t_from_T0_other, t_mjd_other - T0_mjd)
            else:
                t_from_T0_other = np.append(t_from_T0_other, float(t))

        c_0 = SkyCoord(ra_other[0].replace(' ', '') + ' ' + dec_other[0].replace(' ', ''), unit=(u.hourangle, u.deg))
        ra_other[0] = c_0.ra.degree
        dec_other[0] = c_0.dec.degree

        c_1 = SkyCoord(ra_other[1].replace(' ', '') + ' ' + dec_other[1].replace(' ', ''), unit=(u.hourangle, u.deg))
        ra_other[1] = c_1.ra.degree
        dec_other[1] = c_1.dec.degree

        name_inst = ['Master-NET', 'KNC-BBO', 'TRT', 'ASTEP', 'FRAM-Auger', 'KAO', 'Les Makes-T60', 'KNC-T30']
        idx_inst = [5000, 18000, 17000, 16000, 15000, 14000, 13000, 19000]

        for i in range(len(name_inst)):
            idx3 = np.flatnonzero(instrument_other == name_inst[i])
            instrument_other[idx3] = idx_inst[i]

        ra = np.append(ra, ra_other)
        dec = np.append(dec, dec_other)
        t_from_T0 = np.append(t_from_T0, t_from_T0_other)
        depth = np.append(depth, mag_other)
        instrumentID = np.append(instrumentID, instrument_other)
        filters = np.append(filters, filters_other)

        idxs = np.intersect1d(np.flatnonzero(t_from_T0 >= self.tmin),
                              np.flatnonzero(t_from_T0 < self.tmax))
        
        self.ra = ra[idxs]
        self.dec = dec[idxs]
        self.filters = filters[idxs]
        self.instrumentID = instrumentID[idxs]
        self.depth = depth[idxs]
        self.t_from_T0 = t_from_T0[idxs]

        print(len(self.ra))

    def mocs(self):
        skmp = read_sky_map(self.skymap, moc=True, distances=True)
        instrument = []
        mocs = []
        median_distances = []
        fraction = []
        dic_compatible = {}
        for i in range(len(self.instrumentID)):
            if self.instrumentID[i] == 0:
                moc = self.moc_TESS
            else:
                telescope, FOV1, FOV2, reg = find_instrument(self.instrumentID[i])
                instrument.append(telescope)
                center = SkyCoord(self.ra[i], self.dec[i], unit="deg", frame="icrs")
                if telescope == 'ZTF':
                    moc = get_ztf_quadrant_moc(
                        self.ra[i], self.dec[i], max_depth=self.moc_order
                    )
                elif telescope == "DeCAM":
                    moc = get_decam_quadrant_moc(
                        self.ra[i], self.dec[i], max_depth=self.moc_order
                    )
                else:
                    if reg == 'circle':
                        region = regions.CircleSkyRegion(center, FOV1 * u.deg)
                    else:
                        region = regions.RectangleSkyRegion(center, FOV1 * u.deg, FOV2 * u.deg)

                    moc = MOC.from_astropy_regions(region,
                                                   max_depth=self.moc_order)  # changed from 18
            mocs.append(moc)

        return mocs

class S230518h_observations():
    def __init__(self, light_curves, tmin, tmax, skymap):
        self.tmin = tmin
        self.tmax = tmax
        self.light_curves = light_curves
        self.skymap = skymap
        self.T0 = '2023-05-18T12:59:08.151155'
        self.event = 'S230518h'
        self.angle = [0.00, 25.84, 36.87, 45.57, 53.13, 60.00, 66.42, 72.54, 78.46, 84.26, 90.00]
        self.moc_order = 9

        t = np.array([])
        depth = np.array([])
        ra = np.array([])
        dec = np.array([])
        instrumentID = np.array([])
        filters = np.array([])
        
        t = Time(t, format='isot', scale='utc')
        t_mjd = t.to_value('mjd')
        T02 = Time(self.T0, format='isot', scale='utc')                                                                                                            
        T0_mjd = T02.to_value('mjd')
        t_from_T0 = np.array(t_mjd - T0_mjd)

        df = pd.read_csv('KMTNet_obs_S230518h.csv', index_col=None)
        times2 = df['Delta_t'].to_numpy()
        r_a2 = df['RA'].to_numpy()
        decl2 = df['Dec'].to_numpy()
        m2 = df['Depth'].to_numpy()
        filters2 = df['Filter'].to_numpy()

        t_from_T0 = np.append(t_from_T0, times2)
        depth = np.append(depth, m2)
        filters = np.append(filters, filters2)
        instrumentID = np.append(instrumentID, np.ones(len(decl2)) * 95)

        ra_C2 = np.array([])
        dec_C2 = np.array([])

        for i in range(len(r_a2)):
            ra_list = r_a2[i].split(':')
            c2ra = Angle(ra_list[0]+'h'+ra_list[1]+'m'+ra_list[2]+'s').deg
            c2dec = Angle(decl2[i] + ' degrees').to_value()
            ra_C2 = np.append(ra_C2,c2ra)
            dec_C2 = np.append(dec_C2, c2dec)
            
        ra = np.append(ra, ra_C2)
        dec = np.append(dec, dec_C2)

        df = pd.read_csv('atlas_exposures.csv', index_col=None)
        times_ATLAS = df['mjd_t0'].to_numpy()
        filter_ATLAS = df['filter'].to_numpy()
        ra_ATLAS = df['raDeg'].to_numpy()
        dec_ATLAS = df['decDeg'].to_numpy()
        mag_ATLAS = df['mag5sig'].to_numpy()
        instrumentID_ATLAS = 6000 * np.ones(len(mag_ATLAS))

        ra = np.append(ra, ra_ATLAS)
        dec = np.append(dec, dec_ATLAS)
        t_from_T0 = np.append(t_from_T0, times_ATLAS)
        depth = np.append(depth, mag_ATLAS)
        instrumentID = np.append(instrumentID, instrumentID_ATLAS)
        filters = np.append(filters, filter_ATLAS)


        df = pd.read_csv('TESS_obs_S230518h_copie.csv', index_col=None)

        t_T0_TESS = df['t-from-T0'].to_numpy()
        indices = np.flatnonzero(t_T0_TESS <= 6)
        
        t_T0_TESS = t_T0_TESS[indices]
        mag_TESS = df['Mag'].to_numpy()[indices]

        filters_TESS = np.array([])
        RA_TESS = np.array([])
        DEC_TESS = np.array([])
        ID_TESS = np.array([])

        for i in range(len(t_T0_TESS)):
            RA_TESS = np.append(RA_TESS, 'polygon')
            DEC_TESS = np.append(DEC_TESS, 'polygon')
            ID_TESS = np.append(ID_TESS, 0)
            filters_TESS = np.append(filters_TESS, 'tess')

        t_from_T0 = np.append(t_from_T0, t_T0_TESS)                                                                                
        ra = np.append(ra, RA_TESS)                                                                                                 
        dec = np.append(dec, DEC_TESS)                                                                                              
        filters = np.append(filters, filters_TESS)                                                                            
        depth = np.append(depth, mag_TESS)                                                                                        
        instrumentID = np.append(instrumentID, ID_TESS)

        sector = 65
        with open('corners.pkl', 'rb') as f:
            loaded_corners = pickle.load(f)
        polygons = {}
        for cam in range(1, 5):
            polygons[cam] = {}
            for ccd in range(1, 5):
                polygons[cam][ccd] = [SkyCoord(polygon_edges, unit='deg', frame='icrs')
                                      for polygon_edges in loaded_corners[sector][cam][ccd]]
                for vertices in polygons[cam][ccd]:
                    region_sky = PolygonSkyRegion(vertices=vertices)
                    if cam == 1 and ccd == 1:
                        moc_TESS = MOC.from_astropy_regions(region_sky, self.moc_order)
                    else:
                        moc_TESS = moc_TESS.union(MOC.from_astropy_regions(region_sky, self.moc_order))

        idxs = np.intersect1d(np.flatnonzero(t_from_T0 >= self.tmin),
                              np.flatnonzero(t_from_T0 < self.tmax))

        self.ra = ra[idxs]
        self.dec = dec[idxs]
        self.filters = filters[idxs]
        self.instrumentID = instrumentID[idxs]
        self.depth = depth[idxs]
        self.t_from_T0 = t_from_T0[idxs]
        self.moc_TESS = moc_TESS

        print(len(self.ra))

    def mocs(self):
        skmp = read_sky_map(self.skymap, moc=True, distances=True)
        instrument = []
        mocs = []
        median_distances = []
        fraction = []
        dic_compatible = {}
        for i in range(len(self.instrumentID)):
            if self.instrumentID[i] == 0:
                moc = self.moc_TESS
            else:
                telescope, FOV1, FOV2, reg = find_instrument(self.instrumentID[i])
                instrument.append(telescope)
                center = SkyCoord(self.ra[i], self.dec[i], unit="deg", frame="icrs")
                if telescope == 'ZTF':
                    moc = get_ztf_quadrant_moc(
                        self.ra[i], self.dec[i], max_depth=self.moc_order
                    )
                elif telescope == "DeCAM":
                    moc = get_decam_quadrant_moc(
                        self.ra[i], self.dec[i], max_depth=self.moc_order
                    )
                else:
                    if reg == 'circle':
                        region = regions.CircleSkyRegion(center, FOV1 * u.deg)
                    else:
                        region = regions.RectangleSkyRegion(center, FOV1 * u.deg, FOV2 * u.deg)
                        
                    moc = MOC.from_astropy_regions(region,
                                                   max_depth=self.moc_order)  # changed from 18
            mocs.append(moc)

            
        return mocs #, median_distances, fraction, dic_compatible

class S230627c_observations():
    def __init__(self, light_curves, tmin, tmax, skymap):
        self.tmin = tmin
        self.tmax = tmax
        self.light_curves = light_curves
        self.skymap = skymap
        self.T0 = '2023-06-27T01:53:37.834547'
        self.event = 'S230627c'
        self.angle = [0.00, 25.84, 36.87, 45.57, 53.13, 60.00, 66.42, 72.54, 78.46, 84.26, 90.00]
        self.moc_order = 10

        BASE = "https://treasuremap.space/api/v1"
        TARGET = "pointings"
        API_TOKEN = "3AkLCduZwOAaW_1-5uFbMa_j4ULIXNwnFyQZsA"

        json_params = {
            "api_token": API_TOKEN,
            "status": "completed",
            "graceid": self.event,
            "depth_unit": "ab_mag"
        }

        url = "{}/{}".format(BASE, TARGET)

        r = requests.get(url=url, json=json_params)
        t = []
        depth = []
        ra = []
        dec = []
        instrumentID = []
        filters = []

        for i in range(len(json.loads(r.text))):
            pointing = json.loads(r.text)[i]
            # print(pointing)
            filters.append(pointing['band'])
            t.append(pointing['time'])
            depth.append(pointing['depth'])
            instrumentID.append(int(pointing['instrumentid']))
            ra.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[1]))
            dec.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[2]))

        t = Time(t, format='isot', scale='utc')
        t_mjd = t.to_value('mjd')
        T02 = Time(self.T0, format='isot', scale='utc')
        T0_mjd = T02.to_value('mjd')
        depth = np.array(depth)
        t_from_T0 = np.array(t_mjd - T0_mjd)
        ra = np.array(ra)
        dec = np.array(dec)
        instrumentID = np.array(instrumentID)
        filters = np.array(filters)

        idx_not_7dt = np.array([], dtype=int)
        for m in range(len(instrumentID)):
            if instrumentID[m] != 94:
                    idx_not_7dt = np.append(idx_not_7dt, m)

        t_from_T0 = t_from_T0[idx_not_7dt]
        depth = depth[idx_not_7dt]
        instrumentID = instrumentID[idx_not_7dt]
        ra = ra[idx_not_7dt]
        dec = dec[idx_not_7dt]                                                                                                                                     
        filters = filters[idx_not_7dt]
        #s = filters[idx_not_7dt]

        idx_not_goto = np.array([], dtype=int)
        for m in range(len(instrumentID)):
            if instrumentID[m] != 93:
                idx_not_goto = np.append(idx_not_goto, m)
        t_from_T0 = t_from_T0[idx_not_goto]
        depth = depth[idx_not_goto]
        instrumentID = instrumentID[idx_not_goto]
        ra = ra[idx_not_goto]
        dec = dec[idx_not_goto]
        filters = filters[idx_not_goto]
        
        css_idx = np.flatnonzero(instrumentID == 11)
        for i in css_idx:
            depth[i] = depth[i] - 0.125
            filters[i] = 'G'
                    

        df = pd.read_csv('atlas_exposures.csv', index_col=None)
        # instrumentID_ATLAS = df['instrument_ID'].to_numpy()
        times_ATLAS = df['mjd_t0'].to_numpy()
        filter_ATLAS = df['filter'].to_numpy()
        ra_ATLAS = df['raDeg'].to_numpy()
        dec_ATLAS = df['decDeg'].to_numpy()
        mag_ATLAS = df['mag5sig'].to_numpy()
        instrumentID_ATLAS = 6000 * np.ones(len(mag_ATLAS))

        ra = np.append(ra, ra_ATLAS)
        dec = np.append(dec, dec_ATLAS)
        t_from_T0 = np.append(t_from_T0, times_ATLAS)
        depth = np.append(depth, mag_ATLAS)
        instrumentID = np.append(instrumentID, instrumentID_ATLAS)
        filters = np.append(filters, filter_ATLAS)

        df = pd.read_csv('GECKO_obs_S230627C.csv', index_col=None)

        telescope_GECKO = df['instrument'].to_numpy()
        times1_GECKO = df['t-from-T0'].to_numpy()
        r_a_GECKO = df['RA'].to_numpy()
        decl_GECKO = df['DEC'].to_numpy()
        m_GECKO = df['Mag'].to_numpy()
        filter_GECKO = df['filter'].to_numpy()

        t_from_T0 = np.append(t_from_T0, times1_GECKO)
        depth = np.append(depth, m_GECKO)

        ra_G = np.array([])
        dec_G = np.array([])
        for i in range(len(r_a_GECKO)):
            c_G = SkyCoord(r_a_GECKO[i].replace(' ', '') + ' ' + decl_GECKO[i].replace(' ', ''), unit=(u.hourangle, u.deg))
            ra_G = np.append(ra_G, c_G.ra.degree)
            dec_G = np.append(dec_G, c_G.dec.degree)

            if telescope_GECKO[i] == 'LOAO':
                instrumentID = np.append(instrumentID, 2000)
            elif telescope_GECKO[i] == 'KHAO':
                instrumentID = np.append(instrumentID, 7000)
            elif telescope_GECKO[i] == 'CBNUO':
                instrumentID = np.append(instrumentID, 4000)

        ra = np.append(ra, ra_G)
        dec = np.append(dec, dec_G)
        filters = np.append(filters, filter_GECKO)

        df = pd.read_csv('GRANDMA_obs_S230627C.csv', index_col=None)

        telescope_GRANDMA = df['instrument'].to_numpy()
        times1_GRANDMA = df['time'].to_numpy()
        r_a_GRANDMA = df['RA'].to_numpy()
        decl_GRANDMA = df['DEC'].to_numpy()
        m_GRANDMA = df['Mag'].to_numpy()
        filter_GRANDMA = df['filter'].to_numpy()

        t_from_T0 = np.append(t_from_T0, times1_GRANDMA)
        depth = np.append(depth, m_GRANDMA)
        ra = np.append(ra, r_a_GRANDMA)
        dec = np.append(dec, decl_GRANDMA)
        filters = np.append(filters, filter_GRANDMA)

        for telesc in telescope_GRANDMA:
            if telesc == 'Abastumani-T70':
                instrumentID = np.append(instrumentID, 8000)
            elif telesc == 'UBAI-AZT-22':
                instrumentID = np.append(instrumentID, 1000)
            elif telesc == 'UBAI-NT60':
                instrumentID = np.append(instrumentID, 9000)
            elif telesc == 'UBAI-ST60':
                instrumentID = np.append(instrumentID, 10000)
            elif telesc == 'OST-CDK':
                instrumentID = np.append(instrumentID, 12000)
            elif telesc == 'OPD60':
                instrumentID = np.append(instrumentID, 11000)

        df = pd.read_csv('/Users/Virgo/Desktop/GW_GOTO_coverage/230627c_images.csv', index_col=None)
        time_GOTO = df['mjd'].to_numpy()
        filter_GOTO = np.array([])
        for i in range(len(time_GOTO)):
            filter_GOTO = np.append(filter_GOTO, 'L')
        #ra_GOTO = df['corners_ra_dec'].to_numpy()
        #dec_GOTO = df['corners_ra_dec'].to_numpy()
        ra_GOTO = df['centre_ra'].to_numpy()
        dec_GOTO = df['centre_dec'].to_numpy()
        radec_center_GOTO = df['corners_ra_dec'].to_numpy()
        depth_GOTO = df['5sig_lim'].to_numpy()
        instrumentID_GOTO = 93*np.ones(len(depth_GOTO))

        radec_center = np.zeros(len(ra))
        radec_center = np.append(radec_center, radec_center_GOTO)
        ra = np.append(ra, ra_GOTO)
        dec =  np.append(dec, dec_GOTO)
        t_from_T0 = np.append(t_from_T0, time_GOTO-T0_mjd)
        depth = np.append(depth, depth_GOTO)
        instrumentID = np.append(instrumentID, instrumentID_GOTO)
        filters = np.append(filters, filter_GOTO)
        
        idxs = np.intersect1d(np.flatnonzero(t_from_T0 >= self.tmin),
                              np.flatnonzero(t_from_T0 < self.tmax))

        self.ra = ra[idxs]
        self.dec = dec[idxs]
        self.filters = filters[idxs]
        self.instrumentID = instrumentID[idxs]
        self.depth = depth[idxs]
        self.t_from_T0 = t_from_T0[idxs]
        self.radec_center = radec_center[idxs]

    def mocs(self):
        skmp = read_sky_map(self.skymap, moc=True, distances=True)
        instrument = []
        mocs = []
        median_distances = []
        fraction = []
        dic_compatible = {}
        for i in range(len(self.instrumentID)):
            telescope, FOV1, FOV2, reg = find_instrument(self.instrumentID[i])
            instrument.append(telescope)
            center = SkyCoord(self.ra[i], self.dec[i], unit="deg", frame="icrs")
            if telescope == 'ZTF':
                moc = get_ztf_quadrant_moc(
                    self.ra[i], self.dec[i], max_depth=self.moc_order
                )
            elif telescope == "DeCAM":
                moc = get_decam_quadrant_moc(
                    self.ra[i], self.dec[i], max_depth=self.moc_order
                )
            else:

                if reg == 'circle':
                    region = regions.CircleSkyRegion(center, FOV1 * u.deg)
                elif telescope == 'GOTO':
                    couple = self.radec_center[i].replace('(',')').split(')')
                    radec1 = couple[2]
                    radec2 = couple[4]
                    radec3 = couple[6]
                    radec4 = couple[8]
                    
                    radec = [radec1, radec2, radec3, radec4]
                    
                    ra_G = np.array([])
                    dec_G = np.array([])
                    for coord in radec:
                        ra1 =  float(coord.split(',')[0])
                        dec1 = float(coord.split(',')[1])

                        ra_G = np.append(ra_G, ra1)
                        dec_G = np.append(dec_G, dec1)

                    vertices = SkyCoord(ra_G, dec_G, unit='deg', frame='icrs')
                    region = PolygonSkyRegion(vertices=vertices)
                else:
                    region = regions.RectangleSkyRegion(center, FOV1 * u.deg, FOV2 * u.deg)

                moc = MOC.from_astropy_regions(region, max_depth=self.moc_order)  # changed from 18

            mocs.append(moc)

        return mocs


class S230529ay_observations():
    def __init__(self, light_curves, tmin, tmax, skymap):
        self.tmin = tmin
        self.tmax = tmax
        self.light_curves = light_curves
        self.skymap = skymap
        self.T0 = '2023-05-29T18:15:00.747222'
        self.event = 'S230529ay'
        self.angle = [0.00, 25.84, 36.87, 45.57, 53.13, 60.00, 66.42, 72.54, 78.46, 84.26, 90.00]
        self.moc_order = 6
        
        #def find_pointings(self):
        BASE = "https://treasuremap.space/api/v1"
        TARGET = "pointings"
        API_TOKEN = "3AkLCduZwOAaW_1-5uFbMa_j4ULIXNwnFyQZsA"
        
        json_params = {
            "api_token": API_TOKEN,
            "status": "completed",
            "graceid": self.event,
            "depth_unit": "ab_mag"
        }
        
        url = "{}/{}".format(BASE, TARGET)
        
        r = requests.get(url=url, json=json_params)
        t = []
        depth = []
        ra = []
        dec = []
        instrumentID = []
        filters = []
        
        for i in range(len(json.loads(r.text))):
            pointing = json.loads(r.text)[i]
            #print(pointing)
            
            filters.append(pointing['band'])
            t.append(pointing['time'])
            depth.append(pointing['depth'])
            instrumentID.append(int(pointing['instrumentid']))
            ra.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[1]))
            dec.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[2]))
            
        t = Time(t, format='isot', scale='utc')
        t_mjd = t.to_value('mjd')
        T02 = Time(self.T0, format='isot', scale='utc')
        T0_mjd = T02.to_value('mjd')
        depth = np.array(depth)
        t_from_T0 = np.array(t_mjd - T0_mjd)
        ra = np.array(ra)
        dec = np.array(dec)
        instrumentID = np.array(instrumentID)
        filters = np.array(filters)


        idx_not_7dt = np.array([], dtype=int)
        for m in range(len(instrumentID)):
            if instrumentID[m] != 94:
                    idx_not_7dt = np.append(idx_not_7dt, m)

        t_from_T0 = t_from_T0[idx_not_7dt]                                                                                                                     
        depth = depth[idx_not_7dt]                                                                                                                             
        instrumentID = instrumentID[idx_not_7dt]                                                                                                              
        ra = ra[idx_not_7dt]                                                                                                                                    
        dec = dec[idx_not_7dt]                                                                                                                                       
        filters = filters[idx_not_7dt]

        idx_not_goto = np.array([], dtype=int)
        for m in range(len(instrumentID)):
            if instrumentID[m] != 93:
                idx_not_goto = np.append(idx_not_goto, m)

        t_from_T0 = t_from_T0[idx_not_goto]
        depth = depth[idx_not_goto]
        instrumentID = instrumentID[idx_not_goto]
        ra = ra[idx_not_goto]
        dec = dec[idx_not_goto]
        filters = filters[idx_not_goto]
        
        css_idx = np.flatnonzero(instrumentID == 11)
        for i in css_idx:
            depth[i] = depth[i] - 0.125
            filters[i] = 'G'

        df = pd.read_csv('atlas_exposures.csv', index_col=None)
        #instrumentID_ATLAS = df['instrument_ID'].to_numpy()
        times_ATLAS = df['mjd_t0'].to_numpy()
        filter_ATLAS = df['filter'].to_numpy()
        ra_ATLAS = df['raDeg'].to_numpy()
        dec_ATLAS = df['decDeg'].to_numpy()
        mag_ATLAS = df['mag5sig'].to_numpy()
        instrumentID_ATLAS = 6000*np.ones(len(mag_ATLAS))
        
        ra = np.append(ra, ra_ATLAS)
        dec = np.append(dec, dec_ATLAS)
        t_from_T0 = np.append(t_from_T0, times_ATLAS)
        depth = np.append(depth, mag_ATLAS)
        instrumentID = np.append(instrumentID, instrumentID_ATLAS)
        filters = np.append(filters, filter_ATLAS)

        j = np.flatnonzero(filters == 'TESS')
        filters[j] = 'i'

        df = pd.read_csv('230529ay_images.csv', index_col=None)
        time_GOTO = df['mjd'].to_numpy()
        filter_GOTO = np.array([])
        for i in range(len(time_GOTO)):
            filter_GOTO = np.append(filter_GOTO, 'L')
        #ra_GOTO = df['corners_ra_dec'].to_numpy()
        #dec_GOTO = df['corners_ra_dec'].to_numpy()
        ra_GOTO = df['centre_ra'].to_numpy()
        dec_GOTO = df['centre_dec'].to_numpy()
        radec_center_GOTO = df['corners_ra_dec'].to_numpy()
        depth_GOTO = df['5sig_lim'].to_numpy()
        instrumentID_GOTO = 93*np.ones(len(depth_GOTO))

        radec_center = np.zeros(len(ra))
        radec_center = np.append(radec_center, radec_center_GOTO)
        ra = np.append(ra, ra_GOTO)
        dec =  np.append(dec, dec_GOTO)
        t_from_T0 = np.append(t_from_T0, time_GOTO-T0_mjd)
        depth = np.append(depth, depth_GOTO)
        instrumentID = np.append(instrumentID, instrumentID_GOTO)
        filters = np.append(filters, filter_GOTO)
        
        idxs = np.intersect1d(np.flatnonzero(t_from_T0 >= self.tmin),
                                  np.flatnonzero(t_from_T0 < self.tmax))
        
        self.ra = ra[idxs]
        self.dec = dec[idxs]
        self.filters = filters[idxs]
        self.instrumentID = instrumentID[idxs]
        self.depth = depth[idxs]
        self.t_from_T0 = t_from_T0[idxs]
        self.radec_center = radec_center[idxs]

    def mocs(self):
        skmp = read_sky_map(self.skymap, moc=True, distances=True)
        instrument = []
        mocs = []
        median_distances = []
        fraction = []
        dic_compatible = {}
        
        for i in range(len(self.instrumentID)):
            telescope, FOV1, FOV2, reg = find_instrument(self.instrumentID[i])
            instrument.append(telescope)
            #if telescope != 'GOTO':
            center = SkyCoord(self.ra[i], self.dec[i], unit="deg", frame="icrs")
            if telescope == 'ZTF':
                    moc = get_ztf_quadrant_moc(
                        self.ra[i], self.dec[i], max_depth=self.moc_order
                    )
            elif telescope == "DeCAM":
                moc = get_decam_quadrant_moc(
                    self.ra[i], self.dec[i], max_depth=self.moc_order
                )
            else:

                if reg == 'circle':
                    region = regions.CircleSkyRegion(center, FOV1 * u.deg)
                elif telescope == 'GOTO':
                    couple = self.radec_center[i].replace('(',')').split(')')
                    radec1 = couple[2]
                    radec2 = couple[4]
                    radec3 = couple[6]
                    radec4 = couple[8]
                    
                    radec = [radec1, radec2, radec3, radec4]
                    
                    ra_G = np.array([])
                    dec_G = np.array([])
                    
                    for coord in radec:
                        ra1 =  float(coord.split(',')[0])
                        dec1 = float(coord.split(',')[1])
                        
                        ra_G = np.append(ra_G, ra1)
                        dec_G = np.append(dec_G, dec1)
                    
                    vertices = SkyCoord(ra_G, dec_G, unit='deg', frame='icrs')
                    region = PolygonSkyRegion(vertices=vertices)
                else:
                    region = regions.RectangleSkyRegion(center, FOV1 * u.deg, FOV2 * u.deg)

                moc = MOC.from_astropy_regions(region, max_depth=self.moc_order) #changed from 18

            mocs.append(moc)
            
        return mocs

def find_FOV(instrument):
    if instrument == 'ZTF':
        FOV1 = 6.86
        FOV2 = 6.86
        reg = 'square'
    if instrument == 'CSS':
        FOV1 = 2.236
        FOV2 = 2.236
        reg = 'square'
    if instrument == 'ATLAS':
        FOV1 = 5.46
        FOV2 = 5.46
        reg = 'square'
    return FOV1, FOV2, reg


def find_instrument(instrumentID):
    if instrumentID == 38:
        FOV1 = 1.1
        FOV2 = 1.1
        telescope = 'DeCAM'
        reg = 'circle'
    if instrumentID == 44:
        FOV1 = 4.96
        FOV2 = 4.96
        reg = 'square'
        telescope = 'Gattini'
    if instrumentID == 47:
        FOV1 = 6.86
        FOV2 = 6.86
        reg = 'square'
        telescope = 'ZTF'
    if instrumentID == 94:
        FOV2 = 0.9
        FOV1 = 1.36
        reg = 'square'
        telescope = '7DT'
    if instrumentID == 96:
        FOV1 = 0.1516
        FOV2 = 0.1516
        reg = 'square'
        telescope = 'Las Cumbres 2m'
    if instrumentID == 11:
        FOV1 = 2.236
        FOV2 = 2.236
        reg = 'square'
        telescope = 'CSS'
    if instrumentID == 9:
        FOV1 = 0.44166
        FOV2 = 0.44166
        reg = 'square'
        telescope = 'Las Cumbres 1m'
    if instrumentID == 78:
        FOV1 = 1.65
        FOV2 = 1.65
        reg = 'square'
        telescope = 'MeerLICHT'
    if instrumentID == 79:
        FOV1 = 1.65
        FOV2 = 1.65
        reg = 'square'
        telescope = 'BlackGEM'
    if instrumentID == 13:
        FOV1 = 0.2
        FOV2 = 0.2
        reg = 'circle'
        telescope = 'Swift XRT'
    if instrumentID == 12:
        FOV1 = 0.28334
        FOV2 = 0.28334
        reg = 'square'
        telescope = 'Swift UVOT'
    if instrumentID == 91:
        FOV1 = 1.2556
        FOV2 = 0.9416
        reg = 'rectangle'
        telescope = 'WINTER'
    if instrumentID == 93:
        FOV1 = 8.24
        FOV2 = 5.7
        reg = 'rectangle'
        telescope = 'GOTO'
    if instrumentID == 95:
        FOV1 = 2
        FOV2 = 2
        reg = 'square'
        telescope = 'KMTNet'
    if instrumentID == 97:
        FOV1 = 1.25
        FOV2 = 1.25
        reg = 'square'
        telescope = 'PRIME'
    if instrumentID == 99:
        FOV1 = 0.1028
        FOV2 = 0.1306
        reg = 'rectangle'
        telescope = 'Magellan'
    if instrumentID == 68:
        FOV1 = 0.4974
        FOV2 = 0.4952
        reg = 'rectangle'
        telescope = 'Swope'
    if instrumentID == 1000:
        # created by me
        FOV1 = 0.305
        FOV2 = 0.305
        reg = 'square'
        telescope = 'UBAI-AZT-22'
    if instrumentID == 2000:
        # create by me
        FOV1 = 0.467
        FOV2 = 0.467
        reg = 'square'
        telescope = 'LOAO'
    if instrumentID == 3000:
        # create by me
        FOV1 = 0.35
        FOV2 = 0.35
        reg = 'square'
        telescope = 'SAO'
    if instrumentID == 4000:
        FOV1 = 1.19
        FOV2 = 1.19
        reg = 'square'
        telescope = 'CBNUO'

        # to check
    if instrumentID == 5000:
        FOV1 = 2
        FOV2 = 4
        reg = 'rectangle'
        telescope = 'MASTER-Net'
    if instrumentID == 6000:
        FOV1 = 5.46
        FOV2 = 5.46
        reg = 'square'
        telescope = 'ATLAS'
    if instrumentID == 7000:
        FOV1 = 0.3833333
        FOV2 = 0.3833333
        reg = 'square'
        telescope = 'KHAO'
    if instrumentID == 8000:
        FOV1 = 0.5
        FOV2 = 0.5
        reg = 'square'
        telescope = 'Abastumani-T70'
    if instrumentID == 9000:
        FOV1 = 0.18
        FOV2 = 0.18
        reg = 'square'
        telescope = 'UBAI-NT60'
    if instrumentID == 10000:
        FOV1 = 0.11
        FOV2 = 0.11
        reg = 'square'
        telescope = 'UBAI-ST60'
    if instrumentID == 11000:
        FOV1 = 0.095
        FOV2 = 0.095
        reg = 'square'
        telescope = 'OPD60'
    if instrumentID == 12000:
        FOV1 = 0.6
        FOV2 = 0.39
        reg = 'rectangle'
        telescope = 'OST-CDK'
    if instrumentID == 0:
        FOV1 = None
        FOV2 = None
        reg = None
        telescope = 'TESS'
    if instrumentID == 13000:
        FOV1 = 0.3
        FOV2 = 0.3
        reg = 'square'
        telescope = 'Les Makes-T60'
    if instrumentID == 14000:
        FOV1 = 0.0191
        FOV2 = 0.0191
        reg = 'square'
        telescope = 'KAO'
    if instrumentID == 15000:
        FOV1 = 1
        FOV2 = 1
        reg = 'square'
        telescope = 'FRAM-Auger'
    if instrumentID == 16000:
        FOV1 = 1
        FOV2 = 1
        reg = 'square'
        telescope = 'ASTEP'
    if instrumentID == 17000:
        FOV1 = 0.166667
        FOV2 = 0.166667
        reg = 'square'
        telescope = 'TRT'
    if instrumentID == 18000:
        FOV1 = 0.27
        FOV2 = 0.2
        reg = 'rectangle'
        telescope = 'KNC-BBO'
    if instrumentID == 19000:
        FOV1 = 0.4633
        FOV2 = 0.6933
        reg = 'rectangle'
        telescope = 'KNC-T30'

    return telescope, FOV1, FOV2, reg

def get_mag(time, time_inf, time_sup, filter, distance):
    filter_app = filter + 5 * np.log10(distance / 10)
    time_plot, filter_app_plot = remove_nan_inf(time, filter_app)

    splr = splrep(time_plot, filter_app_plot, s=20)
    apparent_mag_new_new_smooth = BSpline(*splr)(time_plot)

    idx1 = np.flatnonzero(time_plot >= time_inf)
    idx2 = np.flatnonzero(time_plot < time_sup)
    idx = np.intersect1d(idx1, idx2)                                                                              

    return apparent_mag_new_new_smooth[idx]

def get_properties(txtfile):
    name = txtfile.split('/')[-1].replace("dat", "")
    keySplit = name.split("_")
    dyn_m = keySplit[1]
    dyn_mej = float(dyn_m.replace("mejdyn", ""))
    dyn_v = keySplit[2]
    wind_mej = float(dyn_v.replace("mejwind", ""))
    phi = float(keySplit[3].replace("phi", ""))
    theta = float(keySplit[4].replace("theta", ""))

    return dyn_mej, wind_mej, phi, theta

def interpolate_one(x,xnew, y):
    f = interpolate.interp1d(x, y)
    ynew = f(xnew)

    return xnew, ynew

def remove_nan_inf(time, filter):
    idx1 = np.flatnonzero(np.isnan(filter))
    filter_new = np.delete(filter, idx1)
    time_new = np.delete(time, idx1)
    idx2 = np.flatnonzero(np.isinf(filter_new))
    filter_new_new = np.delete(filter_new, idx2)
    time_new_new = np.delete(time_new, idx2)

    return time_new_new, filter_new_new

def find_most_probable_obs(skmp, observation_class):
    mocs = observation_class.mocs()
    t_from_T0 = observation_class.t_from_T0
    depth = observation_class.depth
    filters = observation_class.filters
    instrumentID = observation_class.instrumentID

    idx_XRT = np.flatnonzero(filters != 'XRT')

    mocs = np.array(mocs)[idx_XRT]
    t_from_T0 = t_from_T0[idx_XRT]
    depth = depth[idx_XRT]
    filters = filters[idx_XRT]
    instrumentID = instrumentID[idx_XRT]
    
    with fits.open(skmp) as hdul:
        hdul.info()
        data = hdul[1].data
        LEVEL = hdul[1].header["MOCORDER"]

    skymap = read_sky_map(skmp, moc=True, distances=True)

    level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])

    shift = 2 * (LEVEL - level)
    hpx = np.array(np.vstack([ipix << shift, (ipix + 1) << shift]), dtype=np.uint64).T
    nside = ah.level_to_nside(level)
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))

    uniq = skymap["UNIQ"]
    distances = skymap['DISTMU']
    prob = skymap['PROBDENSITY']*pixel_area

    for i in range(len(mocs)):
        if i == 0:
            moc_union = mocs[i]
        else:
            moc_union = moc_union.union(mocs[i])

    mask = moc_union.mask_uniq(uniq, fully_covered_only=True)
    dist_pixels = distances[mask]
    pixels_in_moc = uniq[mask]
    prob_pixels = prob[mask]

    max_prob_idx = np.argmax(prob_pixels)
    dist_max_prob = dist_pixels[max_prob_idx]
    max_prob = prob_pixels[max_prob_idx]

    pixel = pixels_in_moc[max_prob_idx]
    level_pixel, ipix_pixel = ah.uniq_to_level_ipix(pixel)
    nside_pixel = ah.level_to_nside(level_pixel)

    ra, dec = ah.healpix_to_lonlat(ipix_pixel, nside_pixel, order='nested')
    center = SkyCoord(ra.deg, dec.deg, unit="deg", frame="icrs")
    print('most prob loc')
    print('ra', ra.deg)
    print('dec', dec.deg)

    print(moc_union.contains_skycoords(center))
    moc_interest = np.array([])
    mag_interest = np.array([])
    t_interest = np.array([])
    filter_interest = np.array([])
    ID_interest = np.array([])

    for i in range(len(mocs)):
        if mocs[i].contains_skycoords(center):
            moc_interest = np.append(moc_interest, mocs[i])
            mag_interest = np.append(mag_interest, depth[i])
            t_interest = np.append(t_interest, t_from_T0[i])
            filter_interest = np.append(filter_interest, filters[i])
            ID_interest = np.append(ID_interest, instrumentID[i])
            

    idx = np.argmax(mag_interest)

    return ID_interest[idx], t_interest[idx], filter_interest[idx], dist_max_prob, mag_interest[idx], max_prob

def pixel_TESS(moc_TESS, uniq, distances):
    mask = moc_TESS.mask_uniq(uniq)
    dist_pixels = distances[mask]
    pixels_in_moc = uniq[mask]
    teles = 'TESS'

    return teles, pixels_in_moc, dist_pixels
    
def find_field_faster(mocs, t_from_T0, depth, filters, light_curves, skmp, instrumentID):
    nb = len(light_curves)
    max_mag = {}
    fraction_dic = {}
    scenario_dic = {}
    most_constraining = {}
    with fits.open(skmp) as hdul:
        hdul.info()
        data = hdul[1].data
        LEVEL = hdul[1].header["MOCORDER"]

    skymap = read_sky_map(skmp, moc=True, distances=True)
    uniq = skymap["UNIQ"]
    distances = skymap['DISTMU']
    median_distance = np.median(distances)

    if event == 'S230518h':
        idx = np.flatnonzero(instrumentID == 0)
        t_TESS = t_from_T0[idx]
        depth_TESS = depth[idx]

        moc_TESS = np.array(mocs)[idx][0]
        _, pixels_TESS, dist_TESS = pixel_TESS(moc_TESS, uniq, distances)
        for i in range(len(pixels_TESS)):
            dic = compute_LC_for_TESS(pixels_TESS[i], light_curves, dist_TESS[i])

            rule_out_TESS, k_TESS, total_TESS = compute_fraction_TESS(t_TESS, depth_TESS, nb, dic)
            fraction_TESS = float(k_TESS / total_TESS)
            max_mag[pixels_TESS[i]] = 16
            fraction_dic[pixels_TESS[i]] = fraction_TESS
            scenario_dic[pixels_TESS[i]] = rule_out_TESS

            most_constraining[pixels_TESS[i]] = '0-'+str(t_TESS[0])+'-tess-'+str(dist_TESS[i])+'-'+str(16)+'-'+str(fraction_TESS)

        try:
            dic_mag_TESS = dict(map(lambda i,j : (float(i),j) , max_mag.keys(),max_mag.values()))
            dic_frac_TESS = dict(map(lambda i,j : (float(i),j) , fraction_dic.keys(),fraction_dic.values()))
            dic_scenar_TESS = dict(map(lambda i,j : (float(i),j) , scenario_dic.keys(),scenario_dic.values()))
            
            with open(plotDir+"TESS_max_mag_S230518h.json", "w") as outfile:
                json.dump(dic_mag_TESS, outfile)
            with open(plotDir+"TESS_fraction_dic_S230518h.json", "w") as outfile:
                json.dump(dic_frac_TESS, outfile)
            with open(plotDir+"TESS_scenario_dic_S230518h.json", "w") as outfile:
                json.dump(dic_scenar_TESS, outfile)
        except:
            print('no dictionnary')
    
    for i in range(len(mocs)):
        print('-----')
        print(str(i)+' moc over '+str(len(np.flatnonzero(instrumentID != 0))))
        teles, _,_,_ = find_instrument(instrumentID[i])
        print(teles)
        if instrumentID[i] != 0:
            mask = mocs[i].mask_uniq(uniq)
            dist_pixels = distances[mask]
            pixels_in_moc = uniq[mask]
            
            for j in range(len(pixels_in_moc)):
                rule_out, k, total = compute_fraction(t_from_T0[i],
                                    depth[i], light_curves, filters[i], dist_pixels[j])
                fraction = float(k / total)
                if pixels_in_moc[j] not in max_mag.keys():
                    max_mag[pixels_in_moc[j]] = depth[i]
                    fraction_dic[pixels_in_moc[j]] = fraction
                    scenario_dic[pixels_in_moc[j]] = rule_out
                
                    most_constraining[pixels_in_moc[j]] = str(instrumentID[i])+'-'+str(t_from_T0[i])+'-'+str(filters[i])+'-'+str(dist_pixels[j])+'-'+str(depth[i])+'-'+str(fraction)
                elif pixels_in_moc[j] in max_mag.keys():
                    mag = max_mag.pop(pixels_in_moc[j])
                    frac = fraction_dic.pop(pixels_in_moc[j])
                    scenar = scenario_dic.pop(pixels_in_moc[j])
                    
                    if fraction == frac:
                        if depth[i] > mag:
                            constrain = most_constraining.pop(pixels_in_moc[j])
                            most_constraining[pixels_in_moc[j]] = str(instrumentID[i])+'-'+str(t_from_T0[i])+'-'+str(filters[i])+'-'+str(dist_pixels[j])+'-'+str(depth[i])+'-'+str(fraction)
                    
                    elif fraction > frac:
                        constrain = most_constraining.pop(pixels_in_moc[j])
                        most_constraining[pixels_in_moc[j]] = str(instrumentID[i])+'-'+str(t_from_T0[i])+'-'+str(filters[i])+'-'+str(dist_pixels[j])+'-'+str(depth[i])+'-'+str(fraction)
                
                    max_mag[pixels_in_moc[j]] = max(mag, depth[i])
                    fraction_dic[pixels_in_moc[j]] = max(frac, fraction)
                    concat_uniq = np.unique(np.concatenate((np.array(scenar), np.array(rule_out))))
                    scenario_dic[pixels_in_moc[j]] = concat_uniq

    frac_max = []
    depth_max = []
    dist_max = []
    for key, val in most_constraining.items():
        fi = val.split('-')[2]
        if fi == 'XRT' or fi == 'q':
            dist_max.append(0)
        else:
            dist_max.append(float(val.split('-')[-3]))
        frac_max.append(float(val.split('-')[-1]))
        depth_max.append(float(val.split('-')[-2]))
        
    uni = list(most_constraining.keys())
    m = np.max(np.array(frac_max))
    idx = np.flatnonzero(frac_max == m)
    frac_max2 = np.array(frac_max)[idx]
    depth_max2 = np.array(depth_max)[idx]
    dist_max2 = np.array(dist_max)[idx]
    uni2 = np.array(uni)[idx]

    idx2 = np.flatnonzero(depth_max2 == np.max(depth_max2))
    frac_max3 = frac_max2[idx2]
    depth_max3 = depth_max2[idx2]
    dist_max3 = dist_max2[idx2]
    uni3 = uni2[idx2]

    diff = np.abs(dist_max3 - median_distance)
    idx3 = np.flatnonzero(diff == np.min(diff))[0]
    frac_max4 =frac_max3[idx3]
    depth_max4 = depth_max3[idx3]
    dist_max4 = dist_max3[idx3]
    uni4 = uni3[idx3]
    
    most_constraining2 = {uni4: most_constraining[uni4]}
    
    return max_mag, fraction_dic, scenario_dic, most_constraining2


def find_fields(skycoord, mocs, t_from_T0, depth, dist, filters, light_curves):
    mag = np.array([])
    moc_interest = np.array([])
    fraction = []
    rule = np.array([])
    t_interest = np.array([])
    filter_interest = np.array([])
    
    for i in range(len(mocs)):
        if mocs[i].contains_skycoords(skycoord):
            moc_interest = np.append(moc_interest, mocs[i])
            mag = np.append(mag, depth[i])
            t_interest = np.append(t_interest, t_from_T0[i])
            filter_interest = np.append(filter_interest, filters[i])

            rule_out, k, total = compute_fraction(t_from_T0[i],
                        depth[i], light_curves, filters[i], dist)
            fraction.append(float(k / total))
            #print(rule_out)

            rule = np.concatenate((rule, np.array(rule_out)))

    if len(mag) != 0:
        idx = np.argmax(mag)
        color_mag = mag[idx]
        moc = moc_interest[idx]
        rule = np.unique(rule)
        pessimistic = np.max(fraction)
    else:
        moc = None
        color_mag = 0
        pessimistic = -1

    return moc, color_mag, pessimistic, rule

def peak_and_time(time, filter, distance):
    time2 = np.copy(time)
    filter2 = np.copy(filter)
    time_new2, filter_new2 = remove_nan_inf(time2, filter2)
    apparent_mag_new_new = filter_new2 + 5 * np.log10(distance / 10)
    splr = splrep(time_new2, apparent_mag_new_new, s=20)
    smooth_r = BSpline(*splr)(time_new2)
    
    return time_new2, smooth_r


def compute_LC_for_TESS(pixels_TESS, light_curves, dist):
    dic = {}

    for i in range(len(light_curves)):
        file = light_curves[i]
        f = np.loadtxt(file)

        time = f[:, 0]
        dyn_mej, wind_mej, phi, theta = get_properties(file)

        filt = f[:, 13]

        if len(filt) == 0 or math.isinf(dist) == True or dist <= 0:
            array_LC_2D = np.vstack((np.array([]), np.array([])))
        else:
            filter_app = filt + 5 * np.log10(dist * 10 ** 6 / 10)
            times, filter_apps = remove_nan_inf(time, filter_app)

            splr = splrep(times, filter_apps, s=20)
            apparent_mag_new_new_smooth = BSpline(*splr)(times)

            array_LC_2D = np.vstack((times, apparent_mag_new_new_smooth))

        #dic[ str(dyn_mej) + '-' + str(wind_mej) + '-' + str(theta)] =  array_LC_2D
        dic[str(dyn_mej) + '-' + str(wind_mej) + '-' + str(theta)] = array_LC_2D

    return dic

def compute_fraction_TESS(time_obs, depth, nb_config, dictionnary):
    rule_out = []
    k = 0

    for key, values in dictionnary.items():
        dyn_mej = key.split('-')[0]
        wind_mej = key.split('-')[1]
        theta = key.split('-')[2]
        
        times = values[0]
        apparent_mag_new_new_smooth = values[1]
        

        if len(times) != 0:
            spl = splrep(times, apparent_mag_new_new_smooth, k=1)
            x_new = np.array([0, times[0]])
            y_new = splev(x_new, spl, der=0)
            
            new_times = np.concatenate((np.array([x_new[0]]), times))
            new_mag = np.concatenate((np.array([y_new[0]]), apparent_mag_new_new_smooth))
            
            tnew, mag_model = interpolate_one(new_times, time_obs, new_mag)

            idx = np.flatnonzero(depth > mag_model)

            if len(idx) != 0:
                rule_out.append(str(dyn_mej) + '-' + str(wind_mej) + '-' + str(theta))
                k = k + 1
            
    return rule_out, k, nb_config

def compute_fraction(time_obs, depth, light_curves, filter, dist):
    rule_out = []
    k = 0
    for i in range(len(light_curves)):
        file = light_curves[i]
        f = np.loadtxt(file)

        time = f[:, 0]
        dyn_mej, wind_mej, phi, theta = get_properties(file)

        if filter == 'r':
            filt = f[:, 3]
        elif filter == 'i':
            filt = f[:, 4]
        elif filter == 'g':
            filt = f[:, 2]
        elif filter == 'u':
            filt = f[:, 1]
        elif filter == 'z':
            filt = f[:, 5]
        elif filter == 'R' or filter == 'Rc':
            filt = f[:, 13]
        elif filter == 'G':
            filt = f[:,-4]
        elif filter == 'tess':
            filt = f[:,-3]
        elif filter == 'L':
            filt = f[:,-2]
        elif filter=='C' or filter == 'open':
            filt = f[:,-1]
        elif filter == 'I':
            filt = f[:, 14]
        elif filter == 'B':
            filt = f[:, 11]
        elif filter == 'o':
            filt = f[:,-5]
        elif filter == 'c':
            filt = f[:,-6]
        elif filter == 'XRT':
            filt = []
        elif filter == 'q':
            filt = []
        elif filter == 'J':
            filt = f[:,7]

        if len(filt) == 0 or math.isinf(dist) == True or dist <= 0:
            k = 0
        else:
            filter_app = filt + 5 * np.log10(dist * 10 ** 6 / 10)
            times, filter_apps = remove_nan_inf(time, filter_app)
        
            splr = splrep(times, filter_apps, s=20)
            apparent_mag_new_new_smooth = BSpline(*splr)(times)
            if time_obs >= times[0]:
                tnew, mag_model = interpolate_one(times, time_obs, apparent_mag_new_new_smooth)

            elif time_obs < times[0]:
                spl = splrep(times, apparent_mag_new_new_smooth, k=1)
                x_new = np.array([0, times[0]])
                y_new = splev(x_new, spl, der=0)

                new_times = np.concatenate((np.array([x_new[0]]), times))
                new_mag = np.concatenate((np.array([y_new[0]]), apparent_mag_new_new_smooth))
                
                tnew, mag_model = interpolate_one(new_times, time_obs, new_mag)
            
            if depth > mag_model:
                rule_out.append(str(dyn_mej) + '-' + str(wind_mej) + '-' + str(theta))
                k = k + 1
    return rule_out, k, len(light_curves)

def coverage_scenario_ruled_out(observations_class, light_curves, skymap):
    skmp = read_sky_map(skymap, moc=True, distances=True)
    _, _, _, dic_compatible = observations_class.mocs()
    dic_coverage = {}
    for m in range(len(light_curves)):
        file = light_curves[m]
        f = np.loadtxt(file)
        time = f[:, 0]
        dyn_mej, wind_mej, phi, theta = get_properties(file)
        scenario = str(dyn_mej) + '-' + str(wind_mej) + '-' + str(theta)

        coverage = []
        for key, values in dic_compatible.items():
            if scenario not in eval(key):
                coverage.append(values)
        if len(coverage) != 0:
            for i in range(len(coverage)):
                if i == 0:
                    prob_cov = coverage[i]
                else:
                    prob_cov = prob_cov.union(coverage[i])
            vals = MOC.probabilities_in_multiordermap([prob_cov], skmp)
            prob = np.round(np.sum(vals) * 100.0, 3)
        else:
            prob = 0

        dic_coverage[scenario] = prob
    return dic_coverage
        
def plot_skymap(skmp, observations_class, output):
    scenario = []
    dic_scenario = {}

    with fits.open(skmp) as hdul:
        hdul.info()
        data = hdul[1].data
        LEVEL = hdul[1].header["MOCORDER"]
    
    skymap = read_sky_map(skmp, moc=True, distances=True)

    level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])

    shift = 2 * (LEVEL - level)
    hpx = np.array(np.vstack([ipix << shift, (ipix + 1) << shift]), dtype=np.uint64).T
    nside = ah.level_to_nside(level)
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))

    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
    skymap_ra = ra.deg
    skymap_dec = dec.deg
    nside_voulu = 1024
    
    depth = observations_class.depth
    ra_center = observations_class.ra[0]
    dec_center = observations_class.dec[0]
    obs_center = SkyCoord(ra_center, dec_center, unit="deg", frame="icrs")
    max_order = observations_class.moc_order
    mocs = observations_class.mocs()
    instrumentID = observations_class.instrumentID
    filters = observations_class.filters
    t_from_T0 = observations_class.t_from_T0

    plt.rcParams.update({'font.size': 11})
    rect = [0, 0, 1, 1]
    
    fig1 = plt.figure(figsize=(8, 6), dpi=100)
    ax = fig1.add_axes(rect, projection='astro hours mollweide')
    ax.grid()
    
    skymap_raster = rasterize(skymap, order=hp.nside2order(nside_voulu))

    hdu = fits.table_to_hdu(skymap_raster)
    extra_header = [
        ("PIXTYPE", "HEALPIX", "HEALPIX pixelisation"),
        ("ORDERING", "NESTED", "Pixel ordering scheme: RING, NESTED, or NUNIQ"),
        ("COORDSYS", "C", "Ecliptic, Galactic or Celestial (equatorial)"),
        ("INDXSCHM", "EXPLICIT", "Indexing: IMPLICIT or EXPLICIT"),
    ]
    hdu.header.extend(extra_header)

    ax.imshow_hpx(hdu, cmap='cylon')

    fig2 = plt.figure(figsize=(8, 6), dpi=100)
    ax2 = fig2.add_axes(rect, projection='astro hours mollweide')
    ax2.grid()
    ax2.imshow_hpx(hdu, cmap='cylon')


    dic_mag, dic_frac, dic_scenar, most_constraining = find_field_faster(mocs, t_from_T0, depth, filters, light_curves, skmp, instrumentID)

    dic_scenar_float = dict(map(lambda i,j : (float(i),list(j)) , dic_scenar.keys(),dic_scenar.values()))
    
    with open(plotDir+"ipix_vs_scenario_S230627c_"+str(t_min)+"_"+str(t_max)+"_.json", "w") as outfile:
        json.dump(dic_scenar_float, outfile)

    
    dic_scenario = {}
    for key, val in dic_scenar.items():

        lev, match_ipix = ah.uniq_to_level_ipix(key)
        pixel_index = np.flatnonzero(ipix == match_ipix)[0]

        prob = skymap[pixel_index]['PROBDENSITY'] * pixel_area[pixel_index].to_value()
        if len(val) != 0:
            for scenar in val:
                if scenar not in dic_scenario.keys():
                    dic_scenario[scenar] = prob
                else:
                    proba = dic_scenario.pop(scenar)
                    dic_scenario[scenar] = prob + proba
    
    pix_uniq = list(dic_mag.keys())
    levs, ipixs = ah.uniq_to_level_ipix(pix_uniq)
    for i in range(len(pix_uniq)):
        k = pix_uniq[i]
        color_mag = dic_mag[k]
        fraction = float(len(dic_scenar[k])/len(light_curves))
        
        match_ipixs = ipixs[i]
        level_ipixs = levs[i] 
        
        moc_pix = MOC.from_healpix_cells([match_ipixs], [level_ipixs], max_depth=max_order)

        plot_type = ['mag', 'fraction']

        for types in plot_type:
            if types == 'mag':
                boundaries = [18, 18.25, 18.5, 18.75, 19, 19.25, 19.50, 19.75, 20, 20.25, 20.5, 20.75, 21, 21.25, 21.50, 21.75, 22.]
                axes = ax
            elif types == 'fraction':
                boundaries = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
                axes = ax2
            cmap = plt.get_cmap('winter', len(boundaries) - 1)
            norm = BoundaryNorm(boundaries, cmap.N)
    
            sm = ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])

            if types == 'mag':
                color = cmap(norm(color_mag))
                if color_mag != 0:
                    moc_pix.fill(
                        ax=axes,
                        wcs=axes.wcs,
                        alpha=0.8,
                        fill=True,
                        facecolor=color,
                        edgecolor=color,
                        linewidth=1
                    )
            elif types == 'fraction':
                if fraction != -1:
                    color = cmap(norm(fraction))
                    moc_pix.fill(
                        ax=axes,
                        wcs=ax.wcs,
                        alpha=0.8,
                        fill=True,
                        facecolor=color,
                        edgecolor=color,
                        linewidth=1
                    )

            plt.rcParams.update({'font.size': 12})
            if types == 'mag':
                ax1 = fig1.add_axes([0.1, 0.83, 0.8, 0.8])
                ax1.set_axis_off()
                cbar = plt.colorbar(sm, ax=ax1, orientation='horizontal', extend='min', pad=0.1, aspect=50)

                cbar.set_label('Magnitude')
            if types == 'fraction':
                ax3 = fig2.add_axes([0.1, 0.83, 0.8, 0.8])
                ax3.set_axis_off()
                cbar = plt.colorbar(sm, ax=ax3, orientation='horizontal', pad=0.1, aspect=50)

                cbar.set_label('Fraction of scenario ruled out')


    fig1.savefig(output+'mag.png')
    fig2.savefig(output+'fraction.png')

    
    return dic_scenario, most_constraining

def plot_lc_S230518h(observations_class, tmin, tmax, light_curves, constraints):
    
    col = ['navy', 'red', 'chocolate']

    for key, val in constraints.items():
        KeyDict = val.split('-')
    print(KeyDict)
    instru_idx = float(KeyDict[0])
    time_idx = float(KeyDict[1])
    filter_idx = KeyDict[2]
    distance_idx = float(KeyDict[3]) * 10 ** 6
    mag_idx = float(KeyDict[4])
    fraction_idx = float(KeyDict[5])
    telescope_idx, _, _, _ = find_instrument(instru_idx)


    instrument_color = {'Swope':'darkred', 'MeerLICHT':'darkblue', 'MASTER-Net': 'purple', 'KMTNet':'slateblue', 'ATLAS':'orange', 'ZTF':'skyblue'}
    mark = 'v'

    axes2.scatter(time_idx, mag_idx, marker=mark, color=instrument_color[telescope_idx], zorder=1000)
    if tmin == 2:
        axes2.scatter([], [], marker=mark, color='slateblue', label='KMTNet', zorder=1000)
        axes2.scatter([], [], marker=mark, color='orange', label='ATLAS', zorder=1000)

    maximum0 = {}
    minimum0 = {}
    middle_time0 = np.array([])
    max_mag0 = np.array([])
    min_mag0 = np.array([])
    
    maximum45 = {}
    minimum45 = {}
    middle_time45 = np.array([])
    max_mag45 = np.array([])
    min_mag45 = np.array([])
    
    maximum90 = {}
    minimum90 = {}
    middle_time90 = np.array([])
    max_mag90 = np.array([])
    min_mag90 = np.array([])

    for t_range in np.arange(0, 20.2, 0.2):
        values_theta0 = np.array([])
        values_theta45 = np.array([])
        values_theta90 = np.array([])
        for m in range(len(light_curves)):
            file = light_curves[m]
            f = np.loadtxt(file)
            time = f[:, 0]
            dyn_mej, wind_mej, phi, theta = get_properties(file)
            if filter_idx == 'r':
                filt = f[:, 3]
            elif filter_idx == 'i':
                filt = f[:, 4]
            elif filter_idx == 'g':
                filt = f[:, 2]
            elif filter_idx == 'u':
                filt = f[:, 1]
            elif filter_idx == 'z':
                filt = f[:, 5]
            elif filter_idx == 'R' or filter_idx == 'Rc':
                filt = f[:, 13]
            elif filter_idx == 'G':
                filt = f[:,-4]
            elif filter_idx == 'tess':
                filt = f[:,-3]
            elif filter_idx == 'L':
                filt = f[:,-2]
            elif filter_idx =='C' or filter_idx == 'open':
                filt = f[:,-1]
            elif filter_idx == 'I':
                filt = f[:, 14]
            elif filter_idx == 'B':
                filt = f[:, 11]
            elif filter_idx == 'o':
                filt = f[:,-5]
            elif filter_idx == 'c':
                filt = f[:,-6]
                
            if theta == 0:
                filt_val_theta0 = get_mag(time, t_range, t_range+0.2, filt, distance_idx)
                values_theta0 = np.append(values_theta0, filt_val_theta0)

            elif theta == 45.57:
                filt_val_theta45 = get_mag(time, t_range, t_range+0.2, filt, distance_idx)
                values_theta45 = np.append(values_theta45, filt_val_theta45)
            elif theta == 90:
                filt_val_theta90 = get_mag(time, t_range, t_range+0.2, filt, distance_idx)
                values_theta90 = np.append(values_theta90, filt_val_theta90)
        if len(values_theta0) > 0:
            maxim0 = np.max(values_theta0)
            minim0 = np.min(values_theta0)
            maximum0[str((2 * t_range + 0.2) / 2)] = maxim0
            minimum0[str((2 * t_range + 0.2) / 2)] = minim0
        if len(values_theta45) > 0:
            maxim45 = np.max(values_theta45)
            minim45 = np.min(values_theta45)
            maximum45[str((2 * t_range + 0.2) / 2)] = maxim45
            minimum45[str((2 * t_range + 0.2) / 2)] = minim45
        if len(values_theta90) > 0:
            maxim90 = np.max(values_theta90)
            minim90 = np.min(values_theta90)
            maximum90[str((2 * t_range + 0.2) / 2)] = maxim90
            minimum90[str((2 * t_range + 0.2) / 2)] = minim90

    for key, values in maximum0.items():
        max_mag0 = np.append(max_mag0, values)
        middle_time0 = np.append(middle_time0,float(key))
    for key, values in minimum0.items():
        min_mag0 = np.append(min_mag0, values)
        
    for key, values in maximum45.items():
        max_mag45 = np.append(max_mag45, values)
        middle_time45 = np.append(middle_time45,float(key))
    for key, values in minimum45.items():
        min_mag45 = np.append(min_mag45, values)
        
    for key, values in maximum90.items():
        max_mag90 = np.append(max_mag90, values)
        middle_time90 = np.append(middle_time90,float(key))
    for key, values in minimum90.items():
        min_mag90 = np.append(min_mag90, values)

    if tmin == 0:
        a = np.flatnonzero(middle_time0 < 1)
        b = np.flatnonzero(middle_time45 < 1)
        c = np.flatnonzero(middle_time90 < 1)
        axes2.fill_between(middle_time0[a], min_mag0[a], max_mag0[a], color=col[0], alpha=.2)
        axes2.fill_between(middle_time45[b], min_mag45[b], max_mag45[b], color=col[1], alpha=.2)
        axes2.fill_between(middle_time90[c], min_mag90[c], max_mag90[c], color=col[2], alpha=.2)
        
    elif tmin == 1:
        a = np.intersect1d(np.flatnonzero(middle_time0 >= 1), np.flatnonzero(middle_time0 < 2))
        b = np.intersect1d(np.flatnonzero(middle_time45 >= 1), np.flatnonzero(middle_time45 < 2))
        c = np.intersect1d(np.flatnonzero(middle_time90 >= 1), np.flatnonzero(middle_time90 < 2))
        axes2.fill_between(middle_time0[a], min_mag0[a], max_mag0[a], color=col[0], alpha=.2)
        axes2.fill_between(middle_time45[b], min_mag45[b], max_mag45[b], color=col[1], alpha=.2)
        axes2.fill_between(middle_time90[c], min_mag90[c], max_mag90[c], color=col[2], alpha=.2)

    elif tmin == 2:
        a = np.intersect1d(np.flatnonzero(middle_time0 >= 2), np.flatnonzero(middle_time0 < 6))
        b = np.intersect1d(np.flatnonzero(middle_time45 >= 2), np.flatnonzero(middle_time45 < 6))
        c = np.intersect1d(np.flatnonzero(middle_time90 >= 2), np.flatnonzero(middle_time90 < 6))
        axes2.fill_between(middle_time0[a], min_mag0[a], max_mag0[a], color=col[0], alpha=.2, label=r'$\theta=0$')
        axes2.fill_between(middle_time45[b], min_mag45[b], max_mag45[b], color=col[1], alpha=.2, label=r'$\theta=45.57$')
        axes2.fill_between(middle_time90[c], min_mag90[c], max_mag90[c], color=col[2], alpha=.2, label=r'$\theta=90$')


    if plot_type == 'most constraining obs':
        end = 2.95
        end2 = 3
    else:
        end = 5.95
        end2 = 6
        
    if tmin == 0:
        if plot_type == 'most prob obs':
            axes2.text(0.12, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(0.28, 19.2, '(' + str(filter_idx) + '-band)')
        else:
            axes2.text(0.30, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(0.33, 19.2, '(' + str(filter_idx) + '-band)')                                                                                                                                        
    elif tmin == 1:
        if plot_type == 'most prob obs':
            axes2.text(1.12, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(1.23, 19.2, '(' + str(filter_idx) + '-band)')
        else:
            axes2.text(1.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(1.29, 19.2, '(' + str(filter_idx) + '-band)')
    elif tmin == 2:
        if plot_type == 'most prob obs':
            axes2.text(3.5, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(3.7, 19.2, '(' + str(filter_idx) + '-band)')
        else:
            axes2.text(2.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(2.29, 19.2, '(' + str(filter_idx) + '-band)')
            
    if tmin == 2 :
        axes2.plot(1 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
        axes2.plot(2 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
        # ax.set_xscale('log')                                                                                                                                                           
        axes2.set_xlim(0.1, end2)
        axes2.set_ylim(18, 26)
        axes2.set_xlabel(r'T-T$_{0}$ [days]')
        axes2.set_ylabel(r'mag$_\mathrm{app}$')
        
        axes2.annotate('', xy=(0.15, 25), xytext=(0.95, 25), arrowprops=dict(arrowstyle='<->'))
        axes2.annotate('', xy=(1.05, 25), xytext=(1.95, 25), arrowprops=dict(arrowstyle='<->'))
        axes2.annotate('', xy=(2.05, 25), xytext=(end, 25), arrowprops=dict(arrowstyle='<->'))
        #axes2.set_xscale('log')
        
        axes2.invert_yaxis()
        axes2.legend(loc='upper center', ncol=5)



def plot_lc_S240422ed(observations_class, tmin, tmax, light_curves, constraints):
    col = ['navy', 'red', 'chocolate']

    for key, val in constraints.items():
        KeyDict = val.split('-')
    print(KeyDict)
    instru_idx = float(KeyDict[0])
    time_idx = float(KeyDict[1])
    filter_idx = KeyDict[2]
    distance_idx = float(KeyDict[3]) * 10 ** 6
    mag_idx = float(KeyDict[4])
    fraction_idx = float(KeyDict[5])
    telescope_idx, _, _, _ = find_instrument(instru_idx)

    print(telescope_idx)

    if telescope_idx == 'KNC-BOO' or telescope_idx == 'TRT' or telescope_idx == 'Les Makes-T60' or telescope_idx == 'ASTEP' or telescope_idx == 'KAO' or telescope_idx == 'FRAM-Auger':
        telescope_idx = 'GRANDMA'
    instrument_color = {'GRANMDA': 'magenta', 'Swope': 'darkred', 'MeerLICHT': 'darkblue', 'MASTER-Net': 'purple', 'KMTNet': 'slateblue',
                        'ATLAS': 'orange', 'ZTF': 'skyblue', 'DeCAM':'gold'}
    mark = 'D'

    axes2.scatter(time_idx, mag_idx, marker=mark, color=instrument_color[telescope_idx], zorder=1000)
    if tmin == 2:
        axes2.scatter([], [], marker=mark, color='gold', label='DECam', zorder=1000)

    maximum0 = {}
    minimum0 = {}
    middle_time0 = np.array([])
    max_mag0 = np.array([])
    min_mag0 = np.array([])

    maximum45 = {}
    minimum45 = {}
    middle_time45 = np.array([])
    max_mag45 = np.array([])
    min_mag45 = np.array([])

    maximum90 = {}
    minimum90 = {}
    middle_time90 = np.array([])
    max_mag90 = np.array([])
    min_mag90 = np.array([])

    for t_range in np.arange(0, 20.2, 0.2):
        values_theta0 = np.array([])
        values_theta45 = np.array([])
        values_theta90 = np.array([])
        for m in range(len(light_curves)):
            file = light_curves[m]
            f = np.loadtxt(file)
            time = f[:, 0]
            dyn_mej, wind_mej, phi, theta = get_properties(file)
            if filter_idx == 'r':
                filt = f[:, 3]
            elif filter_idx == 'i':
                filt = f[:, 4]
            elif filter_idx == 'g':
                filt = f[:, 2]
            elif filter_idx == 'u':
                filt = f[:, 1]
            elif filter_idx == 'z':
                filt = f[:, 5]
            elif filter_idx == 'R' or filter_idx == 'Rc':
                filt = f[:, 13]
            elif filter_idx == 'G':
                filt = f[:,-4]
            elif filter_idx == 'tess':
                filt = f[:,-3]
            elif filter_idx == 'L':
                filt = f[:,-2]
            elif filter_idx =='C' or filter_idx == 'open':
                filt = f[:,-1]
            elif filter_idx == 'I':
                filt = f[:, 14]
            elif filter_idx == 'B':
                filt = f[:, 11]
            elif filter_idx == 'o':
                filt = f[:,-5]
            elif filter_idx == 'c':
                filt = f[:,-6]
            elif filter_idx == 'J':
                filt = f[:,7]

            if theta == 0:
                filt_val_theta0 = get_mag(time, t_range, t_range + 0.2, filt, distance_idx)
                values_theta0 = np.append(values_theta0, filt_val_theta0)

            elif theta == 45.57:
                filt_val_theta45 = get_mag(time, t_range, t_range + 0.2, filt, distance_idx)
                values_theta45 = np.append(values_theta45, filt_val_theta45)
            elif theta == 90:
                filt_val_theta90 = get_mag(time, t_range, t_range + 0.2, filt, distance_idx)
                values_theta90 = np.append(values_theta90, filt_val_theta90)
        if len(values_theta0) > 0:
            maxim0 = np.max(values_theta0)
            minim0 = np.min(values_theta0)
            maximum0[str((2 * t_range + 0.2) / 2)] = maxim0
            minimum0[str((2 * t_range + 0.2) / 2)] = minim0
        if len(values_theta45) > 0:
            maxim45 = np.max(values_theta45)
            minim45 = np.min(values_theta45)
            maximum45[str((2 * t_range + 0.2) / 2)] = maxim45
            minimum45[str((2 * t_range + 0.2) / 2)] = minim45
        if len(values_theta90) > 0:
            maxim90 = np.max(values_theta90)
            minim90 = np.min(values_theta90)
            maximum90[str((2 * t_range + 0.2) / 2)] = maxim90
            minimum90[str((2 * t_range + 0.2) / 2)] = minim90

    for key, values in maximum0.items():
        max_mag0 = np.append(max_mag0, values)
        middle_time0 = np.append(middle_time0, float(key))
    for key, values in minimum0.items():
        min_mag0 = np.append(min_mag0, values)

    for key, values in maximum45.items():
        max_mag45 = np.append(max_mag45, values)
        middle_time45 = np.append(middle_time45, float(key))
    for key, values in minimum45.items():
        min_mag45 = np.append(min_mag45, values)

    for key, values in maximum90.items():
        max_mag90 = np.append(max_mag90, values)
        middle_time90 = np.append(middle_time90, float(key))
    for key, values in minimum90.items():
        min_mag90 = np.append(min_mag90, values)

    if tmin == 0:
        a = np.flatnonzero(middle_time0 < 1)
        b = np.flatnonzero(middle_time45 < 1)
        c = np.flatnonzero(middle_time90 < 1)
        axes2.fill_between(middle_time0[a], min_mag0[a], max_mag0[a], color=col[0], alpha=.2)
        axes2.fill_between(middle_time45[b], min_mag45[b], max_mag45[b], color=col[1], alpha=.2)
        axes2.fill_between(middle_time90[c], min_mag90[c], max_mag90[c], color=col[2], alpha=.2)

    elif tmin == 1:
        a = np.intersect1d(np.flatnonzero(middle_time0 >= 1), np.flatnonzero(middle_time0 < 2))
        b = np.intersect1d(np.flatnonzero(middle_time45 >= 1), np.flatnonzero(middle_time45 < 2))
        c = np.intersect1d(np.flatnonzero(middle_time90 >= 1), np.flatnonzero(middle_time90 < 2))
        axes2.fill_between(middle_time0[a], min_mag0[a], max_mag0[a], color=col[0], alpha=.2)
        axes2.fill_between(middle_time45[b], min_mag45[b], max_mag45[b], color=col[1], alpha=.2)
        axes2.fill_between(middle_time90[c], min_mag90[c], max_mag90[c], color=col[2], alpha=.2)

    elif tmin == 2:
        a = np.intersect1d(np.flatnonzero(middle_time0 >= 2), np.flatnonzero(middle_time0 < 6))
        b = np.intersect1d(np.flatnonzero(middle_time45 >= 2), np.flatnonzero(middle_time45 < 6))
        c = np.intersect1d(np.flatnonzero(middle_time90 >= 2), np.flatnonzero(middle_time90 < 6))
        axes2.fill_between(middle_time0[a], min_mag0[a], max_mag0[a], color=col[0], alpha=.2, label=r'$\theta=0$')
        axes2.fill_between(middle_time45[b], min_mag45[b], max_mag45[b], color=col[1], alpha=.2,
                           label=r'$\theta=45.57$')
        axes2.fill_between(middle_time90[c], min_mag90[c], max_mag90[c], color=col[2], alpha=.2, label=r'$\theta=90$')

    if plot_type == 'most constraining obs':
        end = 2.95
        end2 = 3
    else:
        end = 5.95
        end2 = 6

    if tmin == 0:
        if plot_type == 'most prob obs':
            axes2.text(0.25, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(0.28, 19.2, '(' + str(filter_idx) + '-band)')
        else:
            axes2.text(0.30, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(0.33, 19.2, '(' + str(filter_idx) + '-band)')
    elif tmin == 1:
        if plot_type == 'most prob obs':
            axes2.text(1.21, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(1.23, 19.2, '(' + str(filter_idx) + '-band)')
        else:
            axes2.text(1.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(1.29, 19.2, '(' + str(filter_idx) + '-band)')
    elif tmin == 2:
        if plot_type == 'most prob obs':
            axes2.text(3.7, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(3.7, 19.2, '(' + str(filter_idx) + '-band)')
        else:
            axes2.text(2.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(2.29, 19.2, '(' + str(filter_idx) + '-band)')

    if tmin == 2:
        axes2.plot(1 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
        axes2.plot(2 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
        # ax.set_xscale('log')
        axes2.set_xlim(0.1, end2)
        axes2.set_ylim(18, 26)
        axes2.set_xlabel(r'T-T$_{0}$ [days]')
        axes2.set_ylabel(r'mag$_\mathrm{app}$')

        axes2.annotate('', xy=(0.15, 25), xytext=(0.95, 25), arrowprops=dict(arrowstyle='<->'))
        axes2.annotate('', xy=(1.05, 25), xytext=(1.95, 25), arrowprops=dict(arrowstyle='<->'))
        axes2.annotate('', xy=(2.05, 25), xytext=(end, 25), arrowprops=dict(arrowstyle='<->'))
        # axes2.set_xscale('log')

        axes2.invert_yaxis()
        axes2.legend(loc='upper center', ncol=4)

    
def plot_lc_S230529ay(observations_class, tmin, tmax, color_dic, light_curves, constraints):
    for key, val in constraints.items():
        KeyDict = val.split('-')
    print(KeyDict)
    instru_idx = float(KeyDict[0])
    time_idx = float(KeyDict[1])
    filter_idx = KeyDict[2]
    distance_idx = float(KeyDict[3])* 10 ** 6
    mag_idx =  float(KeyDict[4])
    fraction_idx = float(KeyDict[5])
    telescope_idx, _, _, _ = find_instrument(instru_idx)

    instrument_color = {'CSS': 'darkgreen', 'ZTF': 'skyblue', 'MASTER-Net': 'purple', 'ATLAS':'orange', 'GOTO':'fuchsia'}

    mark = 'D'
    axes2.scatter(time_idx, mag_idx, marker=mark, color=instrument_color[telescope_idx], zorder=1000)

    if tmin == 2:
        #axes2.scatter([], [], marker=mark, color='darkgreen', label='CSS open filter', zorder=1000)
        axes2.scatter([], [], marker=mark, color='skyblue', label='ZTF', zorder=1000)
        axes2.scatter([], [], marker=mark, color='orange', label='ATLAS', zorder=1000)
        axes2.scatter([], [], marker=mark, color='fuchsia', label='GOTO', zorder=1000)

    for m in range(len(light_curves)):
        file = light_curves[m]
        f = np.loadtxt(file)
        time = f[:, 0]
        dyn_mej, wind_mej, phi, theta = get_properties(file)

        if filter_idx == 'r':
            filt = f[:, 3]
        elif filter_idx == 'i':
            filt = f[:, 4]
        elif filter_idx == 'g':
            filt = f[:, 2]
        elif filter_idx == 'u':
            filt = f[:, 1]
        elif filter_idx == 'z':
            filt = f[:, 5]
        elif filter_idx == 'R' or filter_idx == 'Rc':
            filt = f[:, 13]
        elif filter_idx == 'G':
            filt = f[:,-4]
        elif filter_idx == 'tess':
            filt = f[:,-3]
        elif filter_idx == 'L':
            filt = f[:,-2]
        elif filter_idx =='C' or filter_idx == 'open':
            filt = f[:,-1]
        elif filter_idx == 'I':
            filt = f[:, 14]
        elif filter_idx == 'B':
            filt = f[:, 11]
        elif filter_idx == 'o':
            filt = f[:,-5]
        elif filter_idx == 'c':
            filt = f[:,-6]
                    
        t, mag = peak_and_time(time, filt, distance_idx)
        
        if tmin == 0:
            n = np.flatnonzero(t < 1)
            axes2.text(0.30, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(0.33, 19.2, '('+str(filter_idx)+'-band)')
        elif tmin == 1:
            n = np.intersect1d(np.flatnonzero(t >= 1), np.flatnonzero(t < 2))
            axes2.text(1.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(1.29, 19.2, '('+str(filter_idx)+'-band)')
        elif tmin == 2:
            n = np.intersect1d(np.flatnonzero(t >= 2), np.flatnonzero(t < 6))
            axes2.text(2.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
            axes2.text(2.29, 19.2, '('+str(filter_idx)+'-band)')

        if tmin == 2 :
            axes2.plot(t[n], mag[n], c=color_dic[theta], label=r'$\theta$=' + str(theta) + r'$^\circ$')
            if plot_type == 'most prob obs':
                axes2.plot(1 * np.ones(1000), np.linspace(17.5, 26, 1000), c='black', linewidth=7, linestyle='solid')
                axes2.plot(2 * np.ones(1000), np.linspace(17.5, 26, 1000), c='black', linewidth=7, linestyle='solid')
            else:
                axes2.plot(1 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
                axes2.plot(2 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
            # ax.set_xscale('log')
            axes2.set_xlim(0.1, 3)
            if plot_type == 'most prob obs':
                axes2.set_ylim(17.5, 26)
            else:
                axes2.set_ylim(18, 26)
            axes2.set_xlabel(r'T-T$_{0}$ [days]')
            axes2.set_ylabel(r'mag$_\mathrm{app}$')

            axes2.annotate('', xy=(0.15, 25), xytext=(0.95, 25), arrowprops=dict(arrowstyle='<->'))
            axes2.annotate('', xy=(1.05, 25), xytext=(1.95, 25), arrowprops=dict(arrowstyle='<->'))
            axes2.annotate('', xy=(2.05, 25), xytext=(2.95, 25), arrowprops=dict(arrowstyle='<->'))

            axes2.invert_yaxis()
            axes2.legend(loc='upper center', ncol=4)

        else:
            axes2.plot(t[n], mag[n], c=color_dic[theta])


def plot_lc_S230627c(observations_class, tmin, tmax, color_dic, light_curves, constraints):
    for key, val in constraints.items():
        KeyDict = val.split('-')
    instru_idx = float(KeyDict[0])
    time_idx = float(KeyDict[1])
    filter_idx = KeyDict[2]
    distance_idx = float(KeyDict[3]) * 10 ** 6
    mag_idx = float(KeyDict[4])
    fraction_idx = float(KeyDict[5])
    telescope_idx, _, _, _ = find_instrument(instru_idx)

    if telescope_idx == 'Abastumani-T70' or telescope_idx == 'UBAI-AZT-22' or telescope_idx == 'UBAI-NT60' or telescope_idx == 'UBAI-ST60' or telescope_idx == 'OST-CDK' or telescope_idx == 'OPD60':
        telescope_idx = 'GRANDMA'
    if telescope_idx == 'LOAO' or telescope_idx == 'KHAO' or telescope_idx == 'CBNUO':
        telescope_idx = 'GECKO'

    instrument_color = {'GRANDMA': 'magenta', 'ZTF': 'skyblue', 'MASTER-Net': 'purple', 'ATLAS': 'orange',
                        'GECKO':'slateblue', 'GOTO':'fuchsia'}

    mark = 'D'
    axes2.scatter(time_idx, mag_idx, marker=mark, color=instrument_color[telescope_idx], zorder=1000)

    if tmin == 2:
        axes2.scatter([], [], marker=mark, color='skyblue', label='ZTF', zorder=1000)
        axes2.scatter([], [], marker=mark, color='orange', label='ATLAS', zorder=1000)
        axes2.scatter([], [], marker=mark, color='fuchsia', label='GOTO', zorder=1000)

    for m in range(len(light_curves)):
        file = light_curves[m]
        f = np.loadtxt(file)
        time = f[:, 0]
        dyn_mej, wind_mej, phi, theta = get_properties(file)

        if filter_idx == 'r':
            filt = f[:, 3]
        elif filter_idx == 'i':
            filt = f[:, 4]
        elif filter_idx == 'g':
            filt = f[:, 2]
        elif filter_idx == 'u':
            filt = f[:, 1]
        elif filter_idx == 'z':
            filt = f[:, 5]
        elif filter_idx == 'R' or filter_idx == 'Rc':
            filt = f[:, 13]
        elif filter_idx == 'G':
            filt = f[:,-4]
        elif filter_idx == 'tess':
            filt = f[:,-3]
        elif filter_idx == 'L':
            filt = f[:,-2]
        elif filter_idx =='C' or filter_idx == 'open':
            filt = f[:,-1]
        elif filter_idx == 'I':
            filt = f[:, 14]
        elif filter_idx == 'B':
            filt = f[:, 11]
        elif filter_idx == 'o':
            filt = f[:,-5]
        elif filter_idx == 'c':
            filt = f[:,-6]
            
        t, mag = peak_and_time(time, filt, distance_idx)

        if tmin == 0:
            n = np.flatnonzero(t < 1)
            if plot_type == 'most prob obs':
                axes2.text(0.27, 25.3, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
                axes2.text(0.33, 19.2, '(' + str(filter_idx) + '-band)')
            else:
                axes2.text(0.30, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
                axes2.text(0.33, 19.2, '(' + str(filter_idx) + '-band)')
        elif tmin == 1:
            n = np.intersect1d(np.flatnonzero(t >= 1), np.flatnonzero(t < 2))
            if plot_type == 'most prob obs':
                axes2.text(1.23, 25.3, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
                axes2.text(1.29, 19.2, '(' + str(filter_idx) + '-band)')
            else:
                axes2.text(1.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
                axes2.text(1.29, 19.2, '(' + str(filter_idx) + '-band)')
        elif tmin == 2:
            n = np.intersect1d(np.flatnonzero(t >= 2), np.flatnonzero(t < 6))
            if plot_type == 'most prob obs':
                axes2.text(3, 25.3, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
                axes2.text(3, 19.2, '(' + str(filter_idx) + '-band)')
            else:
                axes2.text(2.26, 24.8, 'D = ' + str(int(round(distance_idx / 10 ** 6, 1))) + ' Mpc')
                axes2.text(2.29, 19.2, '(' + str(filter_idx) + '-band)')

        if tmin == 2:
            axes2.plot(t[n], mag[n], c=color_dic[theta], label=r'$\theta$=' + str(theta) + r'$^\circ$')
            # ax.set_xscale('log')
            if plot_type == 'most prob obs':
                axes2.plot(1 * np.ones(1000), np.linspace(17.5, 26, 1000), c='black', linewidth=7, linestyle='solid')
                axes2.plot(2 * np.ones(1000), np.linspace(17.5, 26, 1000), c='black', linewidth=7, linestyle='solid')
                axes2.set_xlim(0.1, 4.5)
                axes2.set_ylim(17.5, 26)
            else:
                axes2.plot(1 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
                axes2.plot(2 * np.ones(1000), np.linspace(18, 26, 1000), c='black', linewidth=7, linestyle='solid')
                axes2.set_xlim(0.1, 3)
                axes2.set_ylim(18, 26)
            axes2.set_xlabel(r'T-T$_{0}$ [days]')
            axes2.set_ylabel(r'mag$_\mathrm{app}$')

            if plot_type == 'most prob obs':
                axes2.annotate('', xy=(0.15, 25.5), xytext=(0.95, 25.5), arrowprops=dict(arrowstyle='<->'))
                axes2.annotate('', xy=(1.05, 25.5), xytext=(1.95, 25.5), arrowprops=dict(arrowstyle='<->'))
                axes2.annotate('', xy=(2.05, 25.5), xytext=(4.45, 25.5), arrowprops=dict(arrowstyle='<->'))
            else:
                axes2.annotate('', xy=(0.15, 25), xytext=(0.95, 25), arrowprops=dict(arrowstyle='<->'))
                axes2.annotate('', xy=(1.05, 25), xytext=(1.95, 25), arrowprops=dict(arrowstyle='<->'))
                axes2.annotate('', xy=(2.05, 25), xytext=(2.95, 25), arrowprops=dict(arrowstyle='<->'))

            axes2.invert_yaxis()
            if plot_type == 'most prob obs':
                axes2.legend(loc='upper center', ncol=4)
            else:
                axes2.legend(loc='upper center', ncol=4)

        else:
            axes2.plot(t[n], mag[n], c=color_dic[theta])


event = 'S230518h'
plot_LC = True
plot_type = 'most prob obs' # can be 'most constraining obs' or 'most prob obs' for now
plotDir = './S230518h_right_filters_no_MASTER_no_Treasuremap/' #'./final_faster_S230627c/'

if not os.path.exists(plotDir):
    os.mkdir(plotDir)

if event == 'S230529ay':

    light_curves = glob.glob('/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta*_dMpc0.dat')
    
    skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder_S230529ay.fits'
    S230529ay = S230529ay_observations(light_curves, 0, 1, skymap)
    print(np.max(S230529ay.depth))
    S230529ay_1 = S230529ay_observations(light_curves, 1, 2, skymap)
    print(np.max(S230529ay_1.depth))
    S230529ay_2 = S230529ay_observations(light_curves, 2, 6, skymap)
    print(np.max(S230529ay_2.depth))
    
    print('-------S230529ay--------')
    print('between 0 and 1 day')

    t_min = 0
    t_max = 1
    
    dic_scenario_0_1, most_constraining01 = plot_skymap(skymap, S230529ay, plotDir+'test_change_fontsize_skmp_S230529ay_0_1_')
    print('scenario 0 and 1 day', dic_scenario_0_1)
    print('most constraints 0 and 1 day', most_constraining01)

    with open(plotDir+"constraints_0_1_day_S230529ay.json", "w") as outfile:                                                                                    
        json.dump(dic_scenario_0_1, outfile)
    
    print('--------------------------')
    print('between 1 and 2 day')
    t_min = 1
    t_max = 2
    dic_scenario_1_2, most_constraining12 = plot_skymap(skymap, S230529ay_1, plotDir+'test_change_fontsize_skmp_S230529ay_1_2_')
    print('scenario 1 and 2 day', dic_scenario_1_2)
    print('most constraints 1 and 2 day', most_constraining12)

    with open(plotDir+"constraints_1_2_day_S230529ay.json", "w") as outfile:
        json.dump(dic_scenario_1_2, outfile)
    
    print('--------------------------')
    print('between 2 and 6 day')
    t_min = 2
    t_max = 6
    dic_scenario_2_6, most_constraining26 = plot_skymap(skymap, S230529ay_2, plotDir+'test_change_fontsize_skmp_S230529ay_2_6_')
    print('scenario 2 and 6 day', dic_scenario_2_6)
    print('most constraints 2 and 6 day', most_constraining26)

    with open(plotDir+"constraints_2_6_day_S230529ay.json", "w") as outfile:
        json.dump(dic_scenario_2_6, outfile)
    
    if plot_LC == True:
        angle = S230529ay.angle
        cmap = mpl.colormaps['jet_r']
        colors = cmap(np.linspace(0, 1, 11))
        c = [color for i,color in enumerate(colors)]

        color_dic = {}
        for h in range(len(c)):
            color_dic[angle[h]] = c[h]

        if plot_type == 'most prob obs':
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230529ay)
            most_constraining01 = {'most probable pixel' : str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                   + str(depth) + '-' + str(prob)}
            print(most_constraining01)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230529ay_1)
            most_constraining12 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining12)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230529ay_2)
            most_constraining26 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining26)
            name = 'most_prob'
            
            #fig3, axes2 = plt.subplots(1, 1, figsize=(8, 8))
        fig3, axes2 = plt.subplots(1, 1, figsize=(6,6))
        name = 'most_prob'
        plot_lc_S230529ay(S230529ay, 0, 1, color_dic, light_curves, most_constraining01)
        plot_lc_S230529ay(S230529ay_1, 1, 2, color_dic, light_curves, most_constraining12)
        plot_lc_S230529ay(S230529ay_2, 2, 6, color_dic, light_curves, most_constraining26)
        fig3.savefig(plotDir+'test_'+name+'_light_curves_S230529ay_0_6_days_test.png')

        

elif event == 'S230518h':

    all_files = glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_*.dat") 
    light_curves = []
    for LC in all_files:
        dyn_mej, wind_mej, phi, theta = get_properties(LC)
        if dyn_mej <= 0.07 and wind_mej <= 0.03:
            light_curves.append(LC)

    #print(light_curves)
    skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder_S230518h.fits'
    S230518h = S230518h_observations(light_curves, 0, 1, skymap)
    print(np.max(S230518h.depth))
    S230518h_1 = S230518h_observations(light_curves, 1, 2, skymap)
    print(np.max(S230518h_1.depth))
    S230518h_2 = S230518h_observations(light_curves, 2, 6, skymap)
    print(np.max(S230518h_2.depth))

    
    print('-------S230518h--------')
    # show only the top 10
    print('between 0 and 1 day')

    t_min = 0
    t_max = 1
    dic_scenario_0_1, most_constraining01 = plot_skymap(skymap, S230518h, plotDir+'test_skmp_S230518h_0_1_')

    print('most constraining 0 1', most_constraining01)
    
    with open(plotDir+"constraints_0_1_day_S230518h.json", "w") as outfile: 
        json.dump(dic_scenario_0_1, outfile)

    t_min = 1
    t_max = 2
    dic_scenario_1_2, most_constraining12 = plot_skymap(skymap, S230518h_1, plotDir+'test_skmp_S230518h_1_2_')

    print('most constraining 1 2', most_constraining12)
    
    with open(plotDir+"constraints_1_2_day_S230518h.json", "w") as outfile:
        json.dump(dic_scenario_1_2, outfile)

    print('between 2 and 6 day')

    t_min = 2
    t_max = 6
    dic_scenario_2_6, most_constraining26 = plot_skymap(skymap, S230518h_2, plotDir+'test_skmp_S230518h_2_6_')
    print('most constraining 2 6', most_constraining26)
    with open(plotDir+"constraints_2_6_day_S230518h.json", "w") as outfile:
        json.dump(dic_scenario_2_6, outfile)


    all_files = np.array([glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_*theta0.00*.dat"),                                                             
    glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_*theta45.57*.dat"),                                                                                  
    glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_*theta90.00*.dat")]).flatten()                                                                        
    all_files
    light_curves = []
    for LC in all_files:
        dyn_mej, wind_mej, phi, theta = get_properties(LC)
        if dyn_mej <= 0.07 and wind_mej <= 0.03:
            light_curves.append(LC)
    
    if plot_LC == True:
        if plot_type == 'most prob obs':
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230518h)
            most_constraining01 = {'most probable pixel' : str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                   + str(depth) + '-' + str(prob)}
            print(most_constraining01)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230518h_1)
            most_constraining12 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining12)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230518h_2)
            most_constraining26 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining26)

            fig3, axes2 = plt.subplots(1, 1, figsize=(10,6))
        else:    
            fig3, axes2 = plt.subplots(1, 1, figsize=(6,6))

        plt.rcParams.update({'font.size': 12})
        
        plot_lc_S230518h(S230518h, 0, 1, light_curves, most_constraining01)
        plot_lc_S230518h(S230518h_1, 1, 2, light_curves, most_constraining12)
        plot_lc_S230518h(S230518h_2, 2, 6, light_curves, most_constraining26)
        fig3.savefig(plotDir+'test_most_constraining_light_curves_S230518h_0_6_days.pdf')
        fig3.savefig(plotDir+'test_most_constraining_light_curves_S230518h_0_6_days.png')
    
elif event == 'S230627c':
    light_curves = [
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta0.00_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta25.84_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta36.87_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta45.57_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta53.13_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta60.00_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta66.42_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta72.54_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta78.46_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta84.26_dMpc0.dat",
        "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir_more_filter/nph1.0e+06_mejdyn0.010_mejwind0.010_phi30_theta90.00_dMpc0.dat"
    ]
    skymap = '/Users/Virgo/Desktop/Bilby.multiorder_S230627c.fits'
    S230627c = S230627c_observations(light_curves, 0, 1, skymap)
    #print(np.max(S230627c.depth))
    
    S230627c_1 = S230627c_observations(light_curves, 1, 2, skymap)
    #print(np.max(S230627c_1.depth))
    
    S230627c_2 = S230627c_observations(light_curves, 2, 6, skymap)
    #print(np.max(S230627c_2.depth))
    
    print('-------S230627c--------')
    print('between 0 and 1 day')
    t_min = 0
    t_max = 1

    dic_scenario_0_1, most_constraining01 = plot_skymap(skymap, S230627c,
                                                        plotDir + 'test_skmp_S230627c_0_1_')
    print('scenario 0 and 1 day', dic_scenario_0_1)
    print('most constraints 0 and 1 day', most_constraining01)

    with open(plotDir+"constraints_0_1_day_S230627c.json", "w") as outfile:
        json.dump(dic_scenario_0_1, outfile)
    
    print('--------------------------')
    print('between 1 and 2 day')
    t_min = 1
    t_max = 2
    dic_scenario_1_2, most_constraining12 = plot_skymap(skymap, S230627c_1,
                                                        plotDir + 'test_skmp_S230627c_1_2_')
    print('scenario 1 and 2 day', dic_scenario_1_2)
    print('most constraints 1 and 2 day', most_constraining12)

    with open(plotDir+"constraints_1_2_day_S230627c.json", "w") as outfile:
        json.dump(dic_scenario_1_2, outfile)
    
    print('--------------------------')
    print('between 2 and 6 day')
    t_min = 2
    t_max = 6
    dic_scenario_2_6, most_constraining26 = plot_skymap(skymap, S230627c_2,
                                                        plotDir + 'test_skmp_S230627c_2_6_')
    print('scenario 2 and 6 day', dic_scenario_2_6)
    print('most constraints 2 and 6 day', most_constraining26)

    with open(plotDir+"constraints_2_6_day_S230627c.json", "w") as outfile:
        json.dump(dic_scenario_2_6, outfile)

    if plot_LC == True:
        angle = S230627c.angle
        cmap = mpl.colormaps['jet_r']
        colors = cmap(np.linspace(0, 1, 11))
        c = [color for i, color in enumerate(colors)]

        color_dic = {}
        for h in range(len(c)):
            color_dic[angle[h]] = c[h]

        if plot_type == 'most prob obs':
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230627c)
            most_constraining01 = {'most probable pixel' : str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                   + str(depth) + '-' + str(prob)}
            print(most_constraining01)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230627c_1)
            most_constraining12 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining12)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S230627c_2)
            most_constraining26 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining26)

            fig3, axes2 = plt.subplots(1, 1, figsize=(6, 6))
        else:
            fig3, axes2 = plt.subplots(1, 1, figsize=(6, 6))
        plot_lc_S230627c(S230627c, 0, 1, color_dic, light_curves, most_constraining01)
        plot_lc_S230627c(S230627c_1, 1, 2, color_dic, light_curves, most_constraining12)
        plot_lc_S230627c(S230627c_2, 2, 6, color_dic, light_curves, most_constraining26)
        fig3.savefig(plotDir + 'test_light_curves_S230627c_0_6_days_most_probable_loc.pdf')
    
elif event == 'S240422ed':
    light_curves = glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir/nph1.0e+06_*.dat")

    skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder.fits'
    S240422ed = S240422ed_observations(light_curves, 0, 1, skymap)
    S240422ed_1 = S240422ed_observations(light_curves, 1, 2, skymap)
    S240422ed_2 = S240422ed_observations(light_curves, 2, 6, skymap)
    

    print('-------S240422ed--------')
    print('between 0 and 1 day')
    dic_scenario_0_1, most_constraining01 = plot_skymap(skymap, S240422ed, plotDir + 'test_skmp_S240422ed_0_1_')

    print('most constraining 0 1', most_constraining01)

    with open("constraints_0_1_day_S240422ed.json", "w") as outfile:
        json.dump(dic_scenario_0_1, outfile)


    print('between 1 and 2 day')
    dic_scenario_1_2, most_constraining12 = plot_skymap(skymap, S240422ed_1, plotDir + 'test_skmp_S240422ed_1_2_')

    print('most constraining 1 2', most_constraining12)

    with open("constraints_1_2_day_S240422ed.json", "w") as outfile:
        json.dump(dic_scenario_1_2, outfile)

    print('between 2 and 6 day')
    dic_scenario_2_6, most_constraining26 = plot_skymap(skymap, S240422ed_2, plotDir + 'test_skmp_S240422ed_2_6_')
    print('most constraining 2 6', most_constraining26)
    with open("constraints_2_6_day_S240422ed.json", "w") as outfile:
        json.dump(dic_scenario_2_6, outfile)
    

    light_curves = np.array([glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir/nph1.0e+06_*theta0.00*.dat"),
                          glob.glob("/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir/nph1.0e+06_*theta45.57*.dat"),
                          glob.glob(
                              "/Users/Virgo/Desktop/2021_anand_Bu2019nsbh/lcdir/nph1.0e+06_*theta90.00*.dat")]).flatten()

    print('---')
    print('most constraining 0 1', most_constraining01)
    print('most constraining 1 2', most_constraining12)
    print('most constraining 2 6', most_constraining26)
    print('---')

    if plot_LC == True:
        if plot_type == 'most prob obs':
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S240422ed)
            most_constraining01 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining01)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S240422ed_1)
            most_constraining12 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining12)
            ID, t, filt, dist, depth, prob = find_most_probable_obs(skymap, S240422ed_2)
            most_constraining26 = {
                'most probable pixel': str(ID) + '-' + str(t) + '-' + str(filt) + '-' + str(dist) + '-'
                                       + str(depth) + '-' + str(prob)}
            print(most_constraining26)

            fig3, axes2 = plt.subplots(1, 1, figsize=(10, 6))
        else:
            fig3, axes2 = plt.subplots(1, 1, figsize=(6, 6))
        ra, dec = plot_lc_S230422ed(S240422ed, 0, 1, light_curves, most_constraining01)
    
