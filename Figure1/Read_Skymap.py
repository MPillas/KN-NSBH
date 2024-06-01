import healpy as hp
import numpy as np
from astropy.table import QTable
from astropy import units as u
import astropy_healpix as ah
import os
import csv
import pandas as pd


# Fonction pour traiter un fichier FITS et calculer l'aire à 90%
def calculate_area_90(filename):
    skymap = QTable.read(filename)
    skymap.sort('PROBDENSITY', reverse=True)
    level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * skymap['PROBDENSITY']
    cumprob = np.cumsum(prob)
    i = cumprob.searchsorted(0.9)
    area_90 = pixel_area[:i].sum()
    return np.round(area_90.to_value(u.deg**2), 2)

# Répertoire contenant les fichiers FITS
directory = "./IGWN-GWTC3and2-1v2-PESkyLocalizations/"


# Lire les données du fichier CSV
csv_file = "./forchristopher.csv"
df = pd.read_csv(csv_file)


# Liste pour stocker les résultats
results = []

# Parcourir les fichiers du répertoire
for filename in os.listdir(directory):
    if (":Mixed." in filename and filename.endswith(".fits")) or ("GW190521_030229_PEDataRelease_cosmo_reweight_C01:SEOBNRv4PHM" in filename) or ('GW190828_063405_PEDataRelease_cosmo_reweight_C01:IMRPhenomXPHM' in filename) or ('GW190926_050336_PEDataRelease_cosmo_reweight_C01:IMRPhenomXPHM' in filename):
        filepath = os.path.join(directory, filename)
        event_name = filename.split('-')[3].split("_PE")[0]  # Exemple de récupération du nom de l'événement
        area_90 = calculate_area_90(filepath)
        
        # Extraire les données correspondantes du CSV
        event_data = df[df['commonName'] == event_name]
        if not event_data.empty:
            mass_1_source = event_data['mass_1_source'].values[0]
            mass_2_source = event_data['mass_2_source'].values[0]
            luminosity_distance = event_data['luminosity_distance'].values[0]
            luminosity_distance_lower = event_data['luminosity_distance_lower'].values[0]
            luminosity_distance_upper = event_data['luminosity_distance_upper'].values[0]
            source_class=""
            if mass_1_source < 2.4:
            	source_class="BNS"
            elif mass_2_source > 2.3 :
            	source_class="BBH"
            else:
            	source_class="NSBH"
            results.append([event_name, source_class,area_90, mass_1_source, mass_2_source, luminosity_distance, luminosity_distance_lower, luminosity_distance_upper])


# Trier les résultats par event_name
results.sort(key=lambda x: x[0])

# Identifier les commonName manquants dans les résultats
missing_common_names = df[~df['commonName'].isin([result[0] for result in results])]['commonName'].tolist()


# Écrire les résultats dans un fichier CSV
output_csv = "resultsGWTC3.csv"
with open(output_csv, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Event Name", "Classification", "Area 90 (deg^2)", "mass_1_source", "mass_2_source", "luminosity_distance", "luminosity_distance_lower", "luminosity_distance_upper"])
    writer.writerows(results)

#print(f"Résultats écrits dans {output_csv}")
print(f"Common Names manquants: {missing_common_names}")

