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


sys.path.append("./data/")


def max_mass_NS(EOS):
	Max_mass={
	"SLy": 2.054,
	"BHF_BBB2": 1.922,
	"KDE0V": 1.96,
	"KDE0V1": 1.969,
	"SKOP": 1.973,
	"H4": 2.031,
	"HQC18": 2.045,
	"SLY2": 2.054,
	"SLY230A": 2.099,
	"SKMP": 2.107,
	"RS": 2.117,
	"SK255": 2.144,
	"SLY9": 2.156,
	"APR4_EPP": 2.159,
	"SKI2": 2.163,
	"SKI4": 2.17,
	"SKI6": 2.19,
	"SK272": 2.232,
	"SKI3": 2.24,
	"SKI5": 2.24,
	"MPA1": 2.469,
	"MS1B_PP": 2.747,
	"MS1_PP": 2.753
	}
	return Max_mass[EOS]

def mchirp(m1, m2):
    # Define your mchirp function here
    return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)


def compute_hasNSEMB(M_rem,max_mass):
    prediction_ns = np.sum(mass_2 <= max_mass)/len(mass_2)
    prediction_em = np.sum(M_rem > 0)/len(M_rem)
    return prediction_em,prediction_ns

# def discretize_values(matrix, categories):
#     discretized_matrix = np.zeros_like(matrix)
#     for i in range(matrix.shape[0]):
#         for j in range(matrix.shape[1]):
#             if not np.isnan(matrix[i, j]):
#                 idx = np.abs(categories - matrix[i, j]).argmin()
#                 discretized_matrix[i, j] = categories[idx]
#             else:
#                 discretized_matrix[i, j] = np.nan
#     return discretized_matrix
    
    
def Calculate_ejecta(mass_2,mass_1,spin_1z,spin_2z,mode,EOS):
	print(mode,EOS)
	M_rem_Sly_both_max = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	M_rem_Sly_both_min = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	M_rem_Sly_both_median = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	Risco_median = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]
	Cns_median = [[[] for _ in range(len(mass_2))] for _ in range(len(mass_1))]

	

	for i in range(len(mass_2)):
		for j in range(len(mass_1)):
			M_rem_Sly_both_val,C_ns,Risco = computeEjectae(mass_1[j], mass_2[i], spin_1z, spin_2z, eosname=EOS, mode=mode)
			if M_rem_Sly_both_val < 0.001 :
			    M_rem_Sly_both_val = 0.0
			if Risco < 1.0:
				M_rem_Sly_both_val = 0.0
				print("Risco < 1.0", Risco)
				stop
			if Risco > 9.0:
				M_rem_Sly_both_val = 0.0
				print("Risco > 9.0", Risco)
			"""
			if C_ns < 0.14:
				M_rem_Sly_both_val = 0.0
				print("C < 0.14", C_ns)
			"""
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
	#print(M_rem_Sly_both_median_discretized)
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
				
# def plot_mchirp(mass_2,mass_1,mask_matrix,ax):
# 	for i in range(len(mass_2)):
# 		for j in range(len(mass_1)):
# 		    if mask_matrix[j, i]:
# 		        rect = Rectangle((mass_2[i] - 0.5 * (mass_2[1] - mass_2[0]), mass_1[j] - 0.5 * (mass_1[1] - mass_1[0])),
# 		                         mass_2[1] - mass_2[0], mass_1[1] - mass_1[0],
# 		                         linewidth=1, edgecolor='none', facecolor='none', hatch='//', color='blue', fill=True)
# 		        ax.add_patch(rect)
		


#Définition des paramètres

EOS='H4'
m2_max=max_mass_NS(EOS)
print("NS bundary",m2_max)
#M_rem = computeEjectae(8.0, m2_max, 0.8, 0.0, eosname=EOS, mode="both")
#print("Mrem max",np.round(M_rem,3))
#stop
mass_1_d = 0.1
mass_2_d = 0.1

#bundary
m2_max_ceil= np.ceil(m2_max * 10) / 10
mass_2 = np.arange(1.2, m2_max_ceil, mass_2_d)
mass_1 = np.arange(m2_max_ceil, 9.0, mass_1_d)
#spin_1z = [-0.99, -0.9, -0.8, -0.6, -0.3, 0.0, 0.3, 0.6, 0.8, 0.9, 0.99]
spin_1z = [-0.3,0.0,0.3,0.8]
spin_2z = 0.0
mchirp_cut=2.1


# Définition des catégories pour discrétiser
#categories = np.array([0.0, 0.01,0.02, 0.03, 0.04, 0.05,0.06,0.07,0.08, 0.09, 0.10])
# precision=0.002
# categories = np.arange(0.0,0.2,precision)
# print(categories)

#fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
fig, axs = plt.subplots(3, 4, figsize=(18, 18))

fig.suptitle('Dynamical Ejecta (top), Wind Ejecta (middle), and Total (bottom),  for EOS :'+str(EOS)+" given xi=0.3", fontsize=20)

#axs[0, 1].set_title('Dynamical ejecta')

mask_matrix_mchirp=mask_mchirp(mass_1,mass_2,mchirp_cut)

#initial
a,b,c,d,e,f,max_value=Calculate_ejecta(mass_2,mass_1,0.8,0.0,"both",EOS)
print("max_value",np.round(max_value,3))

ejecta_summary = pd.DataFrame([[m1, m2, sp] for m1 in mass_1 for m2 in mass_2 for sp in spin_1z], columns=["mass_1", "mass_2", "spin_1z"])

for i1, mode in enumerate(["computeDynaEjecta", "computeDiskMass", "both"]):
	ejecta_summary[mode] = [0]*len(list(ejecta_summary["mass_1"].iloc))
	for i2 in range(4):
		print("Coords:", i1, i2)
		M_rem_Sly_both_max_discretized,M_rem_Sly_both_min_discretized,M_rem_Sly_both_median_discretized,mask_matrix_min,mask_matrix_max,mask_matrix_median,_=Calculate_ejecta(mass_2,mass_1,spin_1z[i2],spin_2z,mode,EOS)
		for y, a in zip(mass_1, M_rem_Sly_both_median_discretized.tolist()):
			for x, e in zip(mass_2, a):
				ejecta_summary.loc[(ejecta_summary["mass_1"]==y) & (ejecta_summary["mass_2"]==x) & (ejecta_summary["spin_1z"]==spin_1z[i2]), mode]=e
				if (ejecta_summary.loc[(ejecta_summary["mass_1"]==y) & (ejecta_summary["mass_2"]==x) & (ejecta_summary["spin_1z"]==spin_1z[i2]), mode].iloc[0]!=e): exit()
		plot_ejecta(mass_2,mass_1,spin_1z[i2],spin_2z,M_rem_Sly_both_median_discretized,mask_matrix_median,axs[i1, i2],f'{i1},{i2}'+'$M_{ejecta}$ ($M_{sun}$)',max_value)
		#plot_mchirp(mass_2,mass_1,mask_matrix_mchirp,axs[i1, i2])

#plt.title("Dynamical Ejecta, Ejecta Wind, Total") 

ejecta_summary["mchirp_min"] = mchirp(ejecta_summary["mass_1"]-mass_1_d/2, ejecta_summary["mass_2"]-mass_2_d/2)
ejecta_summary["mchirp_max"] = mchirp(ejecta_summary["mass_1"]+mass_1_d/2, ejecta_summary["mass_2"]+mass_2_d/2)
ejecta_summary["spin_2z"] = [spin_2z]*len(list(ejecta_summary["mass_1"].iloc))
ejecta_summary["EOS"] = [EOS]*len(list(ejecta_summary["mass_1"].iloc))
ejecta_summary[""] = [0.0]*len(list(ejecta_summary["mass_1"].iloc))

for c in ["computeDynaEjecta", "computeDiskMass", "both", "mchirp_min", "mchirp_max", "mass_1", "mass_2"]:
	ejecta_summary[c] = np.round(ejecta_summary[c], 3)
	print(ejecta_summary[c])
ejecta_summary.to_csv("ejecta_summary_"+EOS+".csv", index=False)


print("Map_ejecta_"+EOS+".png")

plt.show()

plt.savefig("Map_ejecta_"+EOS+".png")
plt.close(fig)


