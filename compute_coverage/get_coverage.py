import h5py
import argparse
from shapely.geometry import Polygon
import matplotlib
from matplotlib.patches import Rectangle
matplotlib.rcParams.update({'font.size': 12})
from scipy.interpolate import make_smoothing_spline
from scipy.interpolate import splrep, BSpline
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interpolate as interp
import os, sys
import glob
from matplotlib import ticker as mticker
import matplotlib
from matplotlib.patches import Rectangle
matplotlib.rcParams.update({'font.size': 10})
from scipy.interpolate import make_smoothing_spline
from scipy.interpolate import splrep, BSpline
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
import json
from scipy.interpolate import make_smoothing_spline
from scipy.interpolate import splrep, BSpline
import glob
from astropy.coordinates import Angle

def find_instrument(instrumentID):
    # you can add rows here in case of another telescope
    if instrumentID == 38:
        FOV1 = 1.1
        FOV2 = 1.1
        telescope = 'DeCAM'
        reg = 'circle'
    if instrumentID == 44:
        FOV1 = 4.96
        FOV2 = 4.96
        reg = 'square'
        telescope ='Gattini'
    if instrumentID == 47:
        FOV1 = 6.86
        FOV2 = 6.86
        reg = 'square'
        telescope = 'ZTF'
    if instrumentID == 94:
        FOV2 = 0.9
        FOV1 = 1.36
        rectangle = 'square'
        telescope = '7DT'
    if instrumentID == 96:
        FOV1 = 0.1516
        FOV2 = 0.1516
        reg = 'square'
        telescope = 'Las Cumbres 2m'
    if instrumentID == 11:
        FOV1=2.236
        FOV2 = 2.236
        reg = 'square'
        telescope = 'CSS'
    if instrumentID == 9:
        FOV1=0.44166
        FOV2 = 0.44166
        reg = 'square'
        telescope = 'Las Cumbres 1m'
    if instrumentID == 78:
        FOV1=1.65
        FOV2 = 1.65
        reg = 'square'
        telescope = 'MeerLICHT'
    if instrumentID == 79:
        FOV1=1.65
        FOV2 = 1.65
        reg = 'square'
        telescope = 'BlackGEM'
    if instrumentID == 13:
        FOV1=0.2
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
        telescope = 'KM3Net'
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
        #ID created by me
        FOV1 = 0.305
        FOV2 = 0.305
        reg = 'square'
        telescope = 'UBAI-AZT-22'
        
    return telescope, FOV1,FOV2, reg

def plot_skymap(skymap, moc2, center):
    
    level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])
    LEVEL = MOC.MAX_ORDER
    shift = 2 * (LEVEL - level)
    hpx = np.array(np.vstack([ipix << shift, (ipix + 1) << shift]), dtype=np.uint64).T
    nside = ah.level_to_nside(level)
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
    skymap_ra= ra.deg
    skymap_dec = dec.deg
    nside_voulu=1024
    
    vals = MOC.probabilities_in_multiordermap(moc2,skymap)
    print(vals)
    
    skymap_raster=rasterize(skymap, order=hp.nside2order(nside_voulu))
    hdu = fits.table_to_hdu(skymap_raster)
    extra_header = [
	("PIXTYPE", "HEALPIX", "HEALPIX pixelisation"),
	("ORDERING", "NESTED", "Pixel ordering scheme: RING, NESTED, or NUNIQ"),
	("COORDSYS", "C", "Ecliptic, Galactic or Celestial (equatorial)"),
	("INDXSCHM", "EXPLICIT", "Indexing: IMPLICIT or EXPLICIT"),
    ]
    hdu.header.extend(extra_header)
    plt.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(8, 6), dpi=100)
    

    ax = plt.axes(projection='astro hours mollweide')
    ax.grid()
    ax_inset = plt.axes(
	[0.1, 0.45, 0.4, 0.4],
	projection='astro zoom',
	center=center[0],
	radius=8*u.deg)
    #'astro hours mollweide'
    ax.imshow_hpx(hdu, cmap="cylon")
    ax_inset.imshow_hpx(hdu, cmap='cylon')

    g=0
    for moc_tile in moc2:
        moc_tile.fill(
	    ax=ax_inset,
	    wcs=ax_inset.wcs,
	    alpha=vals[g],
	    fill=True,
	    color="black",
	    linewidth=1
        )
        moc_tile.border(ax=ax_inset, wcs=ax_inset.wcs, alpha=1, color="black")
        g=g+1

    plt.savefig('skymap_'+event+'_'+'max_mag.png')
    
def get_distance_of_center(filename, ra,dec):
    map1 = Table.read(filename)
    #print(map1)
    map1.sort('PROBDENSITY', reverse=True)
    max_level = 29
    max_nside = ah.level_to_nside(max_level)
    level, ipix = ah.uniq_to_level_ipix(map1['UNIQ'])
    index = ipix * (2 ** (max_level - level)) ** 2
    sorter = np.argsort(index)
    ra = float(ra)*u.deg
    dec = float(dec)*u.deg
    match_ipix = ah.lonlat_to_healpix(ra, dec, max_nside, order='nested')
    i = sorter[np.searchsorted(index, match_ipix, side='right', sorter=sorter) - 1]

    return map1[i]['DISTMU']


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--event", required=True, help="S230518h, S230529ay, S230627c or S240422ed")
parser.add_argument("--filter", required=True)
parser.add_argument("--telescope", required=True)
parser.add_argument("--skmp_plot", default=False, help="True if you want to plot a skymap with "
                                                      "tiles corresponding to observations")
args = parser.parse_args()

event = args.event
telescope = args.telescope
filter = args.filter
skmp_plot = args.skmp_plot

if event == "S240422ed" and telescope == 'DeCAM' and filter == "r" or filter == "z":
    T0 = '2024-04-22T21:35:13.411201'
    skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder.fits'

    observations_Decam = '/Users/Virgo/Desktop/results_decam.csv'
    f = pd.read_csv(observations_Decam, index_col=None)
    
    t = []
    times = list(f['date_obs'])
    for time in times:
        obs_t = time.split('+')[0]
        t.append(obs_t)
        depth = f['depth']
        ra = f['ra'].to_numpy()
        dec = f['dec'].to_numpy()
        band = list(f['ifilter'])
        
    t = Time(t, format='isot', scale='utc')
    t_mjd = t.to_value('mjd')
    T0 = Time(T0, format='isot', scale='utc')
    T0_mjd = T0.to_value('mjd')
    depth = np.array(depth)
    t_from_T0 = np.array(t_mjd - T0_mjd)
    
    indices_rband = []
    for i in range(len(band)):
        if filter == 'r':
            if 'r DECam' in band[i]:
                indices_rband.append(i)
        elif filter == "z":
            if 'z DECam' in band[i]:
                indices_rband.append(i)
                
    t_from_T0 = t_from_T0[indices_rband]
    depth = depth[indices_rband]
    ra = ra[indices_rband]
    dec = dec[indices_rband]
    
    indices_notnan = []
    for i in range(len(depth)):
        if np.isnan(depth[i]) == False:
            indices_notnan.append(i)

    t_from_T0 = t_from_T0[indices_notnan]
    depth = depth[indices_notnan]
    ra = ra[indices_notnan]
    dec = dec[indices_notnan]
    instrumentID = 38*np.ones(len(ra))

else:
        
    BASE = "https://treasuremap.space/api/v1"
    TARGET = "pointings"
    API_TOKEN = "3AkLCduZwOAaW_1-5uFbMa_j4ULIXNwnFyQZsA"
    
    json_params = {
        "api_token":API_TOKEN,
        "band":filter,
        "status":"completed",
        "graceid":event,
        "instruments":[telescope],
        "depth_unit":"ab_mag"
        ""
    }

    url = "{}/{}".format(BASE, TARGET)
    
    r = requests.get(url=url, json=json_params)
    t = []
    depth = []
    ra = []
    dec = []
    instrumentID = []

    if event == 'S230518h':
        T0 = '2023-05-18T12:59:08.151155'
        skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder_S230518h.fits'
    elif event == 'S230529ay':
        T0 = '2023-05-29T18:15:00.747222'
        skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder_S230529ay.fits'
    elif event == 'S230627c':
        T0 = '2023-06-27T01:53:37.834547'
        skymap = '/Users/Virgo/Desktop/Bilby.multiorder_S230627c.fits'
    elif event == 'S240422ed':
        T0 = '2024-04-22T21:35:13.411201'
        skymap = '/Users/Virgo/Desktop/Bilby.offline0.multiorder.fits'
    
    for i in range(len(json.loads(r.text))):
        pointing = json.loads(r.text)[i]
        t.append(pointing['time'])
        depth.append(pointing['depth'])
        instrumentID.append(int(pointing['instrumentid']))
        ra.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[1]))
        dec.append(float(pointing['position'].replace('(', '').replace(')', '').split(' ')[2]))

    t = Time(t, format='isot', scale='utc')
    t_mjd = t.to_value('mjd')
    T0 = Time(T0, format='isot', scale='utc')
    T0_mjd = T0.to_value('mjd')
    depth = np.array(depth)
    t_from_T0 = np.array(t_mjd - T0_mjd)
    r_a = np.array(ra)
    decli = np.array(dec)

skmp = read_sky_map(skymap, moc=True, distances=True)

for j in range(len(t_from_T0)):
    instrument, FOV1, FOV2, reg = find_instrument(instrumentID[j])
    center = SkyCoord(r_a[j], decli[j], unit="deg", frame="icrs")
    if reg == 'circle':
        region = regions.CircleSkyRegion(center, FOV1 * u.deg)
    else:
        region = regions.RectangleSkyRegion(center, FOV1 * u.deg, FOV2 * u.deg)
    if skmp_plot == True:
        moc2.append(MOC.from_astropy_regions(region, max_depth=18))
        center_all.append(center)
    if j == 0:
        mocs = MOC.from_astropy_regions(region, max_depth=18)
    else:
        mocs = mocs.union(MOC.from_astropy_regions(region, max_depth=18))

vals = MOC.probabilities_in_multiordermap([mocs], skmp)
prob = np.round(np.sum(vals) * 100.0, 3)

print('Pourcentage coverage:', prob)

if skmp_plot == True:
    plot_skymap(skmp, moc2, center_all)
