import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import matplotlib.pyplot as plt
import numpy as np
import h5py


#import em_bright as emb
from computeDiskMass import *

def mchirp(m1, m2):
    # Define your mchirp function here
    return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)


def compute_hasNSEMB(M_rem,max_mass):
    prediction_ns = np.sum(mass_2 <= max_mass)/len(mass_2)
    prediction_em = np.sum(M_rem > 0)/len(M_rem)
    return prediction_em,prediction_ns


def describe_stat(Mrem,EOS,name):
	mean = np.round(np.mean(Mrem),3)
	median = np.round(np.median(Mrem),3)
	std_dev = np.round(np.std(Mrem),3)
	min_value = np.round(np.min(Mrem),3)
	max_value = np.round(np.max(Mrem),3)
	Mrem=pd.Series(Mrem)
	a25th_percentile = np.round(np.percentile(Mrem,25),3)
	#a50th_percentile = np.median(Mrem,50)
	a75th_percentile = np.round(np.percentile(Mrem,25),3)
	print(name,EOS,"mean",mean,"std",std_dev,"median",median,"min",min_value,"max",max_value,"25%",a25th_percentile,"75%",a75th_percentile)

#def treat_superevent(superevent):
# Load data from HDF5 file


#/home/antier/Téléchargements/posterior_samples.h5
#Bilby.posterior_samples.hdf5

with h5py.File("./posterior_samples.h5", 'r') as data:
	print(data)
	table = data['IMRPhenomXP/priors/samples/']
	#mass_1 = 3.9*np.ones(len(table['spin_1z'][:]))#table['mass_1_source'][:]
	mass_1 = table['mass_1_source'][:]
	mass_2 = table['mass_2_source'][:]
	spin_1z = table['spin_1z'][:]
	spin_2z = table['spin_2z'][:]
	indice=np.where((mass_1<4.5) & (mass_2>1.2) & (spin_1z<0.8) & (spin_2z<0.8) & (spin_1z>-0.8) & (spin_2z>-0.8) )[0]
	mass_1=mass_1[indice]
	mass_2=mass_2[indice]
	spin_1z=spin_1z[indice]#0.8*np.ones(len(mass_1))
	spin_2z=spin_2z[indice]#0.0*np.ones(len(mass_1))
	"""
	spin_1z = table['spin_1z'][:]
	spin_1z[spin_1z > 0.8] = 0
	spin_1z[spin_1z < -0.8] = 0
	spin_2z = table['spin_2z'][:]
	spin_2z[spin_2z > 0.8] = 0
	spin_2z[spin_2z < -0.8] = 0
	"""


name="Dyna"
M_rem_Sly_dyn,_,_  = computeEjectae(mass_1, mass_2, spin_1z, spin_2z, eosname='SLy', mode="computeDynaEjecta")
describe_stat(M_rem_Sly_dyn,"SLy",name)
name="Wind"
M_rem_Sly_disk,_,_ = computeEjectae(mass_1, mass_2, spin_1z, spin_2z, eosname='SLy', mode="computeDiskMass")
describe_stat(M_rem_Sly_disk,"SLy",name)
name="Total"
M_rem_Sly_both,_,_  = computeEjectae(mass_1, mass_2, spin_1z, spin_2z, eosname='SLy', mode="both")
describe_stat(M_rem_Sly_both,"SLy",name)


name="Dyna"
M_rem_Sly_dyn_H4,_,_  = computeEjectae(mass_1, mass_2, spin_1z, spin_2z, eosname='H4', mode="computeDynaEjecta")
describe_stat(M_rem_Sly_dyn_H4,"H4",name)
name="Wind"
M_rem_Sly_disk_H4,_,_ = computeEjectae(mass_1, mass_2, spin_1z, spin_2z, eosname='H4', mode="computeDiskMass")
describe_stat(M_rem_Sly_disk_H4,"H4",name)
name="Total"
M_rem_Sly_both_H4,_,_  = computeEjectae(mass_1, mass_2, spin_1z, spin_2z, eosname='H4', mode="both")
describe_stat(M_rem_Sly_both_H4,"H4",name)



# Create a figure with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))

# Subplot 1: Masses and mchirp
ax1.hist(mass_1, alpha=0.2, color='b', bins=10, label='$m_1^{\\mathrm{PE,IMRPhenomXP}}$ ' + str(np.round(np.median(mass_1), 2))+'$\pm$' + str(np.round(np.std(mass_1), 2))+" $M_\odot$")
ax1.hist(mass_2, alpha=0.2, color='r', bins=10, label='$m_2^{\\mathrm{PE,IMRPhenomXP}}$ ' + str(np.round(np.median(mass_2), 2))+'$\pm$' + str(np.round(np.std(mass_2), 2))+" $M_\odot$")
mchirp_source = mchirp(mass_1, mass_2)
ax1.hist(mchirp_source, alpha=0.2, color='g', bins=10,
         label='$m_{chirp}^{\\mathrm{PE}}$ ' + str(np.round(np.median(mchirp_source), 2)) + '$\pm$' + str(np.round(np.std(mchirp_source), 2)))

ax1.legend(fontsize=10)
ax1.set_title("GW230529", fontsize=16)
ax1.set_xlabel("sources masses source IMRPhenomXP PE", fontsize=16)

# Subplot 2: Spin information
ax2.hist(spin_1z, alpha=0.2, color='b', bins=10, label='$\\mathrm{spin}_{1z,IMRPhenomXP}$ ' + str(np.round(np.median(spin_1z), 2))+'$\pm$' + str(np.round(np.std(spin_1z), 2)))
ax2.hist(spin_2z, alpha=0.2, color='r', bins=10, label='$\\mathrm{spin}_{2z,IMRPhenomXP}$ ' + str(np.round(np.median(spin_2z), 2))+'$\pm$' + str(np.round(np.std(spin_2z), 2)))
ax2.legend(fontsize=10)
ax2.set_xlabel("Spin components aligned with the orbital angular momentum", fontsize=14)


# Subplot 3: M_rem
bins_massrem=np.arange(0,0.3,0.02)
ax3.hist(M_rem_Sly_dyn, alpha=1.0, color='green', bins=bins_massrem,  linewidth=2.0, histtype="step",fill=False,label='$M_{dynam,Sly (dynamical)}$:' + str(np.round(np.median(M_rem_Sly_dyn), 2))+'$\pm$' + str(np.round(np.std(M_rem_Sly_dyn), 2)) + " $M_\odot$")
ax3.hist(M_rem_Sly_dyn_H4, alpha=0.5, color='green', bins=bins_massrem, linewidth=2.0, histtype="step",fill=False, label='$M_{dynam,H4 (dynamical)}$:' + str(np.round(np.median(M_rem_Sly_dyn_H4), 2))+'$\pm$' + str(np.round(np.std(M_rem_Sly_dyn_H4), 2)) + " $M_\odot$")
ax3.hist(M_rem_Sly_disk , alpha=1.0, color='purple', bins=bins_massrem, linewidth=2.0, histtype="step",fill=False,label='$M_{wind,SLy (wind)}$: '+ str(np.round(np.median(M_rem_Sly_disk), 3))+'$\pm$' + str(np.round(np.std(M_rem_Sly_disk), 3)) +" $M_\odot$")
ax3.hist(M_rem_Sly_disk_H4 , alpha=0.5, color='purple', bins=bins_massrem, linewidth=2.0, histtype="step",fill=False, label='$M_{wind,H4 (wind)}$: '+ str(np.round(np.median(M_rem_Sly_disk_H4), 3))+'$\pm$' + str(np.round(np.std(M_rem_Sly_disk_H4), 3)) +" $M_\odot$")
#ax3.hist(M_rem_Sly_both, alpha=0.2, color='blue', bins=bins_massrem, label='$Total, SLy$:' + str(np.round(np.median(M_rem_Sly_both), 2))+'$\pm$' + str(np.round(np.std(M_rem_Sly_both), 2)) + " solar mass")
ax3.legend()
ax3.set_xlabel("Ejectae mass (solar mass)", fontsize=16)

# Adjust layout
plt.tight_layout()

plt.savefig("./PE_GW230529.png")
# Show the plots
plt.show()

