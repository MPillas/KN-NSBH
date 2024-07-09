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
    
    
def Calculate_ejecta(mass_2,mass_1,spin_1z,spin_2z,mode,EOS):
	print(mode,EOS)
	M_rem_Sly_both_max = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	M_rem_Sly_both_min = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	M_rem_Sly_both_median = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	
	file1 = open("Formular.txt","w")
	

	for i in range(len(mass_2)):
		for j in range(len(mass_1)):
			M_rem_Sly_both_val,C_ns,Risco = computeEjectae(mass_1[j], mass_2[i], spin_1z, spin_2z, eosname='SLy', mode=mode)
			if M_rem_Sly_both_val < 0.01 :
			    M_rem_Sly_both_val = 0.0
			M_rem_Sly_both_max[j][i].append(M_rem_Sly_both_val)
			M_rem_Sly_both_min[j][i].append(M_rem_Sly_both_val)
			M_rem_Sly_both_median[j][i].append(M_rem_Sly_both_val)
	
	M_rem_Sly_both_median_final = np.zeros((len(mass_1), len(mass_2)))
	M_rem_Sly_both_min_final = np.zeros((len(mass_1), len(mass_2)))
	M_rem_Sly_both_max_final = np.zeros((len(mass_1), len(mass_2)))

	mask_matrix_min = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
	mask_matrix_max = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
	mask_matrix_median = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
	for i in range(len(mass_2)):
		for j in range(len(mass_1)):
			M_rem_Sly_both_median_final[j, i] = np.median(M_rem_Sly_both_median[j][i])
			M_rem_Sly_both_min_final[j, i] = np.min(M_rem_Sly_both_min[j][i])
			M_rem_Sly_both_max_final[j, i] = np.max(M_rem_Sly_both_min[j][i])
			if M_rem_Sly_both_max_final[j, i] == 0.00 :#or M_rem_Sly_both_max[j, i] == 0.09:
				mask_matrix_max[j, i] = True
			if M_rem_Sly_both_min_final[j, i] == 0.00 :#or M_rem_Sly_both_min[j, i] == 0.09:
				mask_matrix_min[j, i] = True
			if M_rem_Sly_both_median_final[j, i] == 0.0:
				mask_matrix_median[j, i] = True
			#L = [str(mass_2[i])+"\t"+str(mass_1[i])+"\t"+spin_1z+"\t"+str(M_rem_Sly_both_median_final[j, i])+"\t"+str(M_rem_Sly_both_median_final[j, i])+"\t"+str(M_rem_Sly_both_min_final[j, i])+"\t"++str(M_rem_Sly_both_max_final[j, i])]
			#file1.writelines(L)
			
		
	#M_rem_Sly_both_max_discretized = discretize_values(M_rem_Sly_both_max_final, categories)
	#M_rem_Sly_both_min_discretized = discretize_values(M_rem_Sly_both_min_final, categories)
	#M_rem_Sly_both_median_discretized = discretize_values(M_rem_Sly_both_median_final, categories)
	
	M_rem_Sly_both_max_discretized = M_rem_Sly_both_max
	M_rem_Sly_both_min_discretized = M_rem_Sly_both_min_final
	M_rem_Sly_both_median_discretized = M_rem_Sly_both_median_final
	print(M_rem_Sly_both_median_discretized)
	max_value=np.max(M_rem_Sly_both_median_discretized)
	return M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,max_value

def plot_ejecta(mass_2,mass_1,spin_1z,spin_2z,Mej,mask_matrix,ax,label_ejecta,max_value):

	cx = ax.imshow(Mej, extent=[mass_2.min() - 0.5 * (mass_2[1] - mass_2[0]), mass_2.max() + 0.5 * (mass_2[1] - mass_2[0]), 
		                                   mass_1.min() - 0.5 * (mass_1[1] - mass_1[0]), mass_1.max() + 0.5 * (mass_1[1] - mass_1[0])], 
		       origin='lower', aspect='auto', cmap='viridis', vmax=max_value)
		       
		       
	for i in range(len(mass_2)):
		for j in range(len(mass_1)):
		    if mask_matrix[j, i]:
		        rect = Rectangle((mass_2[i] - 0.5 * (mass_2[1] - mass_2[0]), mass_1[j] - 0.5 * (mass_1[1] - mass_1[0])),
		                         mass_2[1] - mass_2[0], mass_1[1] - mass_1[0],
		                         linewidth=1, edgecolor='none', facecolor='none', hatch='//', color='grey', fill=True)
		        ax.add_patch(rect)
		        
	cbar = fig.colorbar(cx, ax=ax, label=label_ejecta)
	ax.set_xlabel('Secondary mass $(M_{sun})$')
	ax.set_ylabel('Primary mass $(M_{sun})$')
	if mode=="computeDynaEjecta":
		ax.set_title('Median $M_{ejecta}$ ($M_{sun}$) for $Spin_{1z} =$'+str(np.round(spin_1z,2)))

def mask_mchirp(mass_1,mass_2,mchirp_cut):
	mask_matrix_mchirp = np.zeros((len(mass_1), len(mass_2)), dtype=bool)
	for i in range(len(mass_2)):
		for j in range(len(mass_1)):
			mchirp_value=mchirp(mass_1[j], mass_2[i])
			if np.abs(mchirp_cut-mchirp_value) > 0.1:
				mask_matrix_mchirp[j, i] = True
	return mask_matrix_mchirp
				
def plot_mchirp(mass_2,mass_1,mask_matrix,ax):
	for i in range(len(mass_2)):
		for j in range(len(mass_1)):
		    if mask_matrix[j, i]:
		        rect = Rectangle((mass_2[i] - 0.5 * (mass_2[1] - mass_2[0]), mass_1[j] - 0.5 * (mass_1[1] - mass_1[0])),
		                         mass_2[1] - mass_2[0], mass_1[1] - mass_1[0],
		                         linewidth=1, edgecolor='none', facecolor='none', hatch='//', color='blue', fill=True)
		        ax.add_patch(rect)
		


#Définition des paramètres

EOS='SLy'
#M_rem,_,_ = computeEjectae(2.0, 1.4, 0.8, 0.0, eosname='SLy', mode="both")
#print("Mrem",M_rem)
#stop
mass_2 = np.arange(1.0, 1.6, 0.1)
mass_1 = np.arange(2.1, 8.0, 0.5)
#spin_1z = [-0.99, -0.9, -0.8, -0.6, -0.3, 0.0, 0.3, 0.6, 0.8, 0.9, 0.99]
spin_1z = [-0.3,0.0,0.3,0.8]
spin_2z = 0.0
mchirp_cut=2.1


# Définition des catégories pour discrétiser
#categories = np.array([0.0, 0.01,0.02, 0.03, 0.04, 0.05,0.06,0.07,0.08, 0.09, 0.10])
precision=0.002
categories = np.arange(0.0,0.2,precision)
print(categories)

#fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
fig, axs = plt.subplots(3, 4, figsize=(18, 18))

fig.suptitle('Dynamical Ejecta, Wind Ejecta, and Total, xi 0.3', fontsize=20)

#axs[0, 1].set_title('Dynamical ejecta')

mask_matrix_mchirp=mask_mchirp(mass_1,mass_2,mchirp_cut)

#initial
a,b,c,d,e,f,max_value=Calculate_ejecta(mass_2,mass_1,0.8,0.0,"both",EOS)
print("max_value",max_value)

mode="computeDynaEjecta"

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[0],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[0],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[0, 0],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[0, 0])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[1],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[1],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[0, 1],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[0, 1])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[2],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[2],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[0, 2],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[0, 2])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[3],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[3],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[0, 3],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[0, 3])
      

#axs[1, 1].set_title('Wind Ejecta')

mode="computeDiskMass"
M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[0],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[0],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[1, 0],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[1, 0])


M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[1],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[1],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[1, 1],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[1, 1])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[2],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[2],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[1, 2],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[1, 2])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[3],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[3],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[1, 3],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[1, 3])

#axs[2, 1].set_title('Ejecta Remnant')

mode="both"
M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[0],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[0],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[2, 0],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[2, 0])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[1],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[1],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[2, 1],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[2, 1])

M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[2],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[2],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[2, 2],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[2, 2])


M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[3],spin_2z,mode,EOS)
plot_ejecta(mass_2,mass_1,spin_1z[3],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[2, 3],'$M_{ejecta}$ ($M_{sun}$)',max_value)
#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[2, 3])

#plt.title("Dynamical Ejecta, Ejecta Wind, Total") 





plt.show()


