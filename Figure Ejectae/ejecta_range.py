import pandas as pd
import numpy as np

# Load the CSV file
EOS="H4"
spin1=0.8

mchirp = [1.0,1.2,1.6,2.0,2.4,2.8,3.2,3.6]
#mchirp=[1.9]
uncertainty=0.2
i=0
for mchirp_value in mchirp:


	csv_file_path = "ejecta_summary_"+EOS+".csv"
	df = pd.read_csv(csv_file_path)

	# Define mchirp range




	df2 = df[((df['spin_1z'] == spin1))]

	max_val_both = np.round(df2[['both']].max(),3)
	min_val_both = np.round(df2[['both']].min(),3)
	median_val_both=np.round(df2[['both']].median(),3)

	max_val_wind = df2[['computeDiskMass']].max()
	min_val_wind = df[['computeDiskMass']].min()

	max_val_dyn = df2[['computeDynaEjecta']].max()
	min_val_dyn = df2[['computeDynaEjecta']].min()


	mchirp_df = df[(df['mchirp_min'] >= (mchirp_value - uncertainty)) & (df['mchirp_max'] <= (mchirp_value + uncertainty))]

	# If no rows are found, use the full table
	if mchirp_df.empty:
		print("no find compatible mchirp")
		continue

	#print(mchirp_df)


	mchirp_df = mchirp_df[((mchirp_df['spin_1z'] == spin1))]

	"""
	max_val_both_mchirp = np.round(mchirp_df[['both']].max(),3)
	min_val_both_mchirp = np.round(mchirp_df[['both']].min(),3)

	max_val_wind_mchirp = mchirp_df[['computeDiskMass']].max()
	min_val_wind_mchirp = mchirp_df[['computeDiskMass']].min()

	max_val_dyn_mchirp = mchirp_df[['computeDynaEjecta']].max()
	min_val_dyn_mchirp = mchirp_df[['computeDynaEjecta']].min()

	#print("Upper limit for total ejecta with MSI PP:")
	#print(mchirp_df[(mchirp_df['both'] == float(max_val_both))] )
	filt_max=mchirp_df[(mchirp_df['both'] == float(max_val_both))]
	filt_min=mchirp_df[(mchirp_df['both'] == float(min_val_both))]

	sentenceA="For "+EOS+" "+", for 1.2 < m2 < max_massNS given by the EOS, and max_massNS < m1 < 9 Msun, for a scenario with none spining BH : total mass ejecta < "+str(float(max_val_both))+" Msun and > "+str(float(min_val_both))+" Msun."
	#print(sentenceA)

	reduce_space=np.round((1-(float(max_val_both_mchirp) - float(min_val_both_mchirp)) / (float(max_val_both) - float(min_val_both)))*100,0)

	sentenceB="Assuming an estimation of mchirp of "+str(mchirp_value)+ " +/- "+str(uncertainty)+ " and for a scenario with none spining BH : total mass ejecta < "+str(float(max_val_both_mchirp))+" Msun and > "+str(float(min_val_both_mchirp))+" Msun."+" This helped to decrease the uncertainty around the ejecta mass by"+str(reduce_space)+" %."
	#print(sentenceB)
	"""
	df2=df2.rename(columns={'both': 'Total Ejecta', 'computeDynaEjecta': 'Dyn','computeDiskMass': 'Wind'})
	mchirp_df=mchirp_df.rename(columns={'both': 'Total Ejecta', 'computeDynaEjecta': 'Dyn','computeDiskMass': 'Wind'})
	if i==0:
		print("EOS: ",EOS,"spin BH",spin1)
		print(df2[['Total Ejecta','Dyn','Wind']].describe())
	print("EOS: ",EOS,"spin BH",spin1,"mchirp ",mchirp_value)
	print(mchirp_df[['Total Ejecta','Dyn','Wind']].describe())
	#print(mchirp_df)
	i=i+1

