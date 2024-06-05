import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import matplotlib.pyplot as plt
import numpy as np
import h5py
import statistics

#import em_bright as emb
from computeDiskMass import *
from matplotlib.patches import Rectangle

def mchirp(m1, m2):
    # Define your mchirp function here
    return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)


def compute_hasNSEMB(M_rem,max_mass):
    prediction_ns = np.sum(mass_2 <= max_mass)/len(mass_2)
    prediction_em = np.sum(M_rem > 0)/len(M_rem)
    return prediction_em,prediction_ns

def discretize_values(matrix, categories):
    discretized_matrix = np.zeros_like(matrix)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not np.isnan(matrix[i, j]):
                idx = np.abs(categories - matrix[i, j]).argmin()
                discretized_matrix[i, j] = categories[idx]
            else:
                discretized_matrix[i, j] = np.nan
    return discretized_matrix

# Définition des paramètres
mass_2 = np.arange(1.0, 2.0, 0.04)
mass_1 = np.arange(3.0, 7.0, 0.5)
#spin_1z = [-0.99, -0.9, -0.8, -0.6, -0.3, 0.0, 0.3, 0.6, 0.8, 0.9, 0.99]
spin_1z = [ 0.0, 0.3, 0.6]
spin_2z = 0.0

# Initialisation des matrices
M_rem_Sly_both_max = np.ones((len(mass_1), len(mass_2))) * 0.01
M_rem_Sly_both_min = np.ones((len(mass_1), len(mass_2))) * 0.09
M_rem_Sly_both_median = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]


# Définition des catégories pour discrétiser
categories = np.array([0.0, 0.01,0.02, 0.03, 0.04, 0.05,0.06,0.07,0.08, 0.09, 0.10])

for k in range(len(spin_1z)):
    for i in range(len(mass_2)):
        for j in range(len(mass_1)):
            M_rem_Sly_both_val = computeEjectae(mass_1[j], mass_2[i], spin_1z[k], spin_2z, eosname='SLy', mode="computeDynaEjecta")
            # Clipping des valeurs hors des bornes
            if M_rem_Sly_both_val < 0.01 or M_rem_Sly_both_val > 0.09:
                M_rem_Sly_both_val = 0.0

            # Mise à jour des matrices min et max
            if M_rem_Sly_both_val!=0.0:
                M_rem_Sly_both_max[j, i] = max(M_rem_Sly_both_max[j, i], M_rem_Sly_both_val)
                M_rem_Sly_both_min[j, i] = min(M_rem_Sly_both_min[j, i], M_rem_Sly_both_val)
                #print(M_rem_Sly_both_val)
                M_rem_Sly_both_median[j][i].append(M_rem_Sly_both_val)
                


# Discrétiser les matrices

M_rem_Sly_both_median_final = np.zeros((len(mass_1), len(mass_2)))

mask_matrix_min = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
mask_matrix_max = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
mask_matrix_median = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
for i in range(len(mass_2)):
	for j in range(len(mass_1)):
		if M_rem_Sly_both_max[j, i] == 0.01 or M_rem_Sly_both_max[j, i] == 0.09:
			mask_matrix_max[j, i] = True
		if M_rem_Sly_both_min[j, i] == 0.01 or M_rem_Sly_both_min[j, i] == 0.09:
			mask_matrix_min[j, i] = True
		if M_rem_Sly_both_median[j][i]:  
			M_rem_Sly_both_median_final[j, i] = np.median(M_rem_Sly_both_median[j][i])
		else:
			M_rem_Sly_both_median_final[j, i] = 0.0
		#print(M_rem_Sly_both_median[j, i])
		#M_rem_Sly_both_median[j, i]=np.median(M_rem_Sly_both_median[j][i])
		if M_rem_Sly_both_median_final[j, i] == 0.0:
			mask_matrix_median[j, i] = True
    	
       


M_rem_Sly_both_max_discretized = discretize_values(M_rem_Sly_both_max, categories)
M_rem_Sly_both_min_discretized = discretize_values(M_rem_Sly_both_min, categories)
M_rem_Sly_both_median_discretized = discretize_values(M_rem_Sly_both_median_final, categories)
print(M_rem_Sly_both_median)
print(M_rem_Sly_both_median_discretized)

# Création des graphiques
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

cx1 = ax1.imshow(M_rem_Sly_both_min_discretized, extent=[mass_2.min() - 0.5 * (mass_2[1] - mass_2[0]), mass_2.max() + 0.5 * (mass_2[1] - mass_2[0]), 
                                       mass_1.min() - 0.5 * (mass_1[1] - mass_1[0]), mass_1.max() + 0.5 * (mass_1[1] - mass_1[0])], 
           origin='lower', aspect='auto', cmap='viridis')
           
           
for i in range(len(mass_2)):
    for j in range(len(mass_1)):
        if mask_matrix_min[j, i]:
            rect = Rectangle((mass_2[i] - 0.5 * (mass_2[1] - mass_2[0]), mass_1[j] - 0.5 * (mass_1[1] - mass_1[0])),
                             mass_2[1] - mass_2[0], mass_1[1] - mass_1[0],
                             linewidth=1, edgecolor='none', facecolor='none', hatch='//', color='grey', fill=False)
            ax1.add_patch(rect)
            
cbar = fig.colorbar(cx1, ax=ax1, label='$M_{ejecta}$ ($M_{sun}$)')
ax1.set_xlabel('Secondary mass $(M_{sun})$')
ax1.set_ylabel('Primary mass $(M_{sun})$')
ax1.set_title('Min $M_{ejecta}$ ($M_{sun}$) for $Spin_{1z} = [-0.99,0.99]$')


cx2 = ax2.imshow(M_rem_Sly_both_median_discretized, extent=[mass_2.min() - 0.5 * (mass_2[1] - mass_2[0]), mass_2.max() + 0.5 * (mass_2[1] - mass_2[0]), 
                                       mass_1.min() - 0.5 * (mass_1[1] - mass_1[0]), mass_1.max() + 0.5 * (mass_1[1] - mass_1[0])], 
           origin='lower', aspect='auto', cmap='viridis')
cbar = fig.colorbar(cx2, ax=ax2, label='$M_{ejecta}$ ($M_{sun}$)')

for i in range(len(mass_2)):
    for j in range(len(mass_1)):
        if mask_matrix_median[j, i]:
            rect = Rectangle((mass_2[i] - 0.5 * (mass_2[1] - mass_2[0]), mass_1[j] - 0.5 * (mass_1[1] - mass_1[0])),
                             mass_2[1] - mass_2[0], mass_1[1] - mass_1[0],
                             linewidth=1, edgecolor='none', facecolor='none', hatch='//', color='grey', fill=False)
            ax2.add_patch(rect)

ax2.set_xlabel('Secondary mass $(M_{sun})$')
ax2.set_ylabel('Primary mass $(M_{sun})$')
ax2.set_title('Median $M_{ejecta}$ ($M_{sun}$) $Spin_{1z} = [-0.99,0.99]$')


cx3 = ax3.imshow(M_rem_Sly_both_max_discretized, extent=[mass_2.min() - 0.5 * (mass_2[1] - mass_2[0]), mass_2.max() + 0.5 * (mass_2[1] - mass_2[0]), 
                                       mass_1.min() - 0.5 * (mass_1[1] - mass_1[0]), mass_1.max() + 0.5 * (mass_1[1] - mass_1[0])], 
           origin='lower', aspect='auto', cmap='viridis')
cbar = fig.colorbar(cx3, ax=ax3, label='$M_{ejecta}$ ($M_{sun}$)')

for i in range(len(mass_2)):
    for j in range(len(mass_1)):
        if mask_matrix_max[j, i]:
            rect = Rectangle((mass_2[i] - 0.5 * (mass_2[1] - mass_2[0]), mass_1[j] - 0.5 * (mass_1[1] - mass_1[0])),
                             mass_2[1] - mass_2[0], mass_1[1] - mass_1[0],
                             linewidth=1, edgecolor='none', facecolor='none', hatch='//', color='grey', fill=False)
            ax3.add_patch(rect)

ax3.set_xlabel('Secondary mass $(M_{sun})$')
ax3.set_ylabel('Primary mass $(M_{sun})$')
ax3.set_title('Max $M_{ejecta}$ ($M_{sun}$) $Spin_{1z} = [-0.99,0.99]$')

plt.show()



