import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_event_data(events, label, classification):
    """
    Plots event data with error bars and specific formatting based on label and classification.
    
    Args:
    events (list of dict): List of event dictionaries with keys: 'distance', 'error', 'cr', 'name'.
    label (str): Label for the plot (e.g., 'O4a', 'O4b', 'O4a/b').
    classification (str): Classification of the events (e.g., 'NSBH').
    """
    # Extract data for plotting
    distances = [event["distance"] for event in events]
    errors = [event["error"] for event in events]
    cr_values = [event["cr"] for event in events]
    labels = [event["name"] for event in events]

    # Print statistics
    print(f"{label} Distance mean/std: {np.round(np.mean(distances))} / {np.round(np.std(distances))}")

    # Plotting based on the label
    if label == "O4a":
        plt.errorbar(distances, cr_values, xerr=errors, fmt='^', alpha=0.8, capsize=5, color="k", label=f"{classification} alerts {label}")
    elif label == "O4b":
        plt.errorbar(distances, cr_values, xerr=errors, fmt='^', alpha=1.0, capsize=5, color="k", markeredgecolor="rosybrown", label=f"{classification} alerts {label}")
    elif label == "O4a/b" and classification == "NSBH":
        plt.errorbar(distances, cr_values, xerr=errors, fmt='^', alpha=1.0, capsize=5, color="r", label=f"{classification} alerts {label}")
        # Add specific text annotations
        plt.text(distances[0], cr_values[0] + 8000, "GW230529ay", color='r', fontsize=8, ha='center')
        plt.text(distances[1], cr_values[1] + 20, "S230627c", color='r', fontsize=8, ha='center')
        plt.text(distances[2], cr_values[2] + 50, "S240422ed", color='r', fontsize=8, ha='center', bbox=dict(facecolor='grey', alpha=0.3, pad=10))

def plot_cumulative_data(events, label, cumulative_type):
    """
    Plots cumulative data with error bars.

    Args:
    events (list of dict): List of event dictionaries with keys: 'distance', 'error'.
    label (str): Label for the plot (e.g., 'O4a').
    cumulative_type (str): Type of cumulative data to plot (e.g., 'distances').
    """
    if cumulative_type == "distances":
        sorted_events = sorted(events, key=lambda event: event['distance'])
        distances = [event['distance'] for event in sorted_events]
        errors = [event['error'] for event in sorted_events]
        # Convertir les distances et les erreurs en tableaux NumPy
        distances = np.array(distances)
        print(distances)
        errors = np.array(errors)
        
        print("ratio BBH with median distance below 2000 Mpc", np.where(distances<2000)[0][-1]/len(distances))
        
        
        #print(np.log10(bins))
        if label=="O4a BBH":
            bins = np.logspace(2., 4.,9)
            print(np.log10(bins))
            plt.hist(distances, bins=bins,weights=errors,log=False,color="k",histtype="step",cumulative=True,density=True,label='O4a BBH alerts ' + str(int(np.round(np.median(distances/100.0), 0)*100))+' Mpc $\pm$' + str(int(np.round(np.std(distances/100), 0)*100))+" Mpc")
        if label=="O4b BBH":
            bins = np.logspace(2., 4.,9)
            print(np.log10(bins))
            plt.hist(distances, bins=bins,weights=errors,log=False,color="b",histtype="step",cumulative=True,density=True,label='O4b BBH alerts ' + str(int(np.round(np.median(distances/100.0), 0)*100))+' Mpc $\pm$' + str(int(np.round(np.std(distances/100), 0)*100))+" Mpc")
        
        plt.gca().set_xscale("log")
  
        plt.legend()
        plt.xlabel("Distance (Mpc)")
        plt.ylabel("Cumulative")

    if cumulative_type == "cr":
        sorted_events = sorted(events, key=lambda event: event['cr'])
        creg = [event['cr'] for event in sorted_events]
        creg = np.array(creg)
        print(creg)
        print("ratio alerts with median credible below 500 deg2", np.where(creg<500)[0][-1]/len(creg))

        if label=="O4a BBH":
            bins = np.logspace(1., 4.,10)
            plt.hist(creg, bins=bins,log=False,color="k",histtype="step",cumulative=False,density=True,label='O4a all alerts 90$\%$ cr ' + str(int(np.round(np.median(creg/100.0), 0)*100))+' Mpc $\pm$' + str(int(np.round(np.std(creg/100), 0)*100))+" Mpc")
        if label=="O4b BBH":
            bins = np.logspace(1., 4.,10)
            print(np.log10(bins))
            plt.hist(creg, bins=bins,log=False,color="b",histtype="step",cumulative=False,density=True,label='O4b all alerts 90$\%$ cr ' + str(int(np.round(np.median(creg/100.0), 0)*100))+' Mpc $\pm$' + str(int(np.round(np.std(creg/100), 0)*100))+" Mpc")
        
        plt.gca().set_xscale("log")
  
        plt.legend()
        plt.xlabel("90% credible region area (square degrees)")
        plt.ylabel("Density")
        
        
def Weizmann_LV_O4():
   mpc = [
    100000, 120595, 146487, 192611, 274152, 365724, 428554, 498574, 567658, 
    655718, 720188, 850152, 900717, 1003721, 1071064, 1159627, 1255576, 
    1349638, 1429984, 1570895, 1664378, 1776334, 1937252, 2112749, 2304074, 
    2477085, 2760498, 2903811, 3054626, 3213048, 3503978, 3766558, 4048651, 
    4352047, 4991703, 5288384, 5766706, 6333940, 7369417, 8636418, 9691873
   ]
   mpc=np.array(mpc)*0.001

   cumulative = [
    -0.001, -0.003, -0.001, 0.002, 0.010, 0.020, 0.030, 0.039, 0.055, 0.076, 
    0.096, 0.128, 0.149, 0.181, 0.196, 0.226, 0.262, 0.291, 0.319, 0.367, 
    0.392, 0.430, 0.469, 0.508, 0.543, 0.595, 0.634, 0.661, 0.690, 0.710, 
    0.743, 0.776, 0.802, 0.834, 0.872, 0.888, 0.909, 0.934, 0.955, 0.980, 
    0.987
   ]
   plt.plot(mpc, cumulative, linestyle='-', color='r', label="LV Prospects")
   plt.legend()

# Plot data

# O1 BBH data
bbh_o1_distances = [421, 445, 1046]
bbh_o1_crs = [(177 + 180 + 177) / 3, (1007 + 1006 + 1006) / 3, (1521 + 1540 + 1518) / 3]
bbh_o1_errors = [(421 - 267, 567 - 421), (445 - 254, 585 - 445), (1046 - 598, 1490 - 1046)]
bbh_o1_errors = ([x[0] for x in bbh_o1_errors], [x[1] for x in bbh_o1_errors])
plt.errorbar(bbh_o1_distances, bbh_o1_crs, xerr=bbh_o1_errors, fmt='o', alpha=0.3, capsize=5, color="k")

# O1 BNS data
bns_distances = [40]
bns_crs = [16]
bns_errors = [(40 - 30, 49 - 40)]
bns_errors = ([x[0] for x in bns_errors], [x[1] for x in bns_errors])
plt.errorbar(bns_distances, bns_crs, xerr=bns_errors, fmt='v', capsize=5, color="m")
plt.text(bns_distances[0], bns_crs[0] + 5, "GW170817", color='m', fontsize=8, ha='center')

# O2 BBH data
bbh_o2_distances = [1035, 579, 316, 969, 933, 2732, 1799]
bbh_o2_crs = [38, 87, 383, 342, 907, 1014, 1640]
bbh_o2_errors = [
    (1035 - 659, 1381 - 1035), (579 - 373, 747 - 579), (316 - 206, 403 - 316),
    (969 - 611, 1294 - 969), (933 - 530, 1260 - 933), (2732 - 1414, 3732 - 2732), (1799 - 932, 2634 - 1799)
]
bbh_o2_errors = ([x[0] for x in bbh_o2_errors], [x[1] for x in bbh_o2_errors])
plt.errorbar(bbh_o2_distances, bbh_o2_crs, xerr=bbh_o2_errors, fmt='v', alpha=0.3, capsize=5, color="k")

# O3 BBH data
data = pd.read_csv('resultsGWTC3.csv')
bbh_o3_data = data[data['Classification'] == 'BBH']
plt.errorbar(
    bbh_o3_data['luminosity_distance'], bbh_o3_data['Area 90 (deg^2)'],
    xerr=[-bbh_o3_data['luminosity_distance_lower'], bbh_o3_data['luminosity_distance_upper']],
    fmt='*', alpha=0.3, capsize=5, color="k", label="BBH: O1 (o), O2 (v), O3 (*)"
)

# O3 NSBH data
nsbh_o3_data = data[data['Classification'] == 'NSBH']
plt.errorbar(
    nsbh_o3_data['luminosity_distance'], nsbh_o3_data['Area 90 (deg^2)'],
    xerr=[-nsbh_o3_data['luminosity_distance_lower'], nsbh_o3_data['luminosity_distance_upper']],
    fmt='*', capsize=5, color="r", label="NSBH: O2 (v), O3 (*)"
)


# Annotate specific NSBH O3 events
plt.text(nsbh_o3_data['luminosity_distance'].values[2], nsbh_o3_data['Area 90 (deg^2)'].values[2] + 120, "GW200115", color='r', fontsize=6, ha='center')

# O3 BNS data
bns_o3_data = data[data['Classification'] == 'BNS']
plt.errorbar(
    bns_o3_data['luminosity_distance'], bns_o3_data['Area 90 (deg^2)'],
    xerr=[-bns_o3_data['luminosity_distance_lower'], bns_o3_data['luminosity_distance_upper']],
    fmt='*', capsize=5, color="m", label="BNS: O2 (v), O3 (*)"
)
plt.text(bns_o3_data['luminosity_distance'].values[0], bns_o3_data['Area 90 (deg^2)'].values[0] + 1200, "GW190425", color='m', fontsize=8, ha='center')

events_O4a = [
    {"name": "S230601bf", "distance": 3565, "error": 1260, "cr": 2531},
    {"name": "S230605o", "distance": 1067, "error": 333, "cr": 1077},
    {"name": "S230606d", "distance": 2545, "error": 874, "cr": 1221},
    {"name": "S230608as", "distance": 3447, "error": 1079, "cr": 1694},
    {"name": "S230609u", "distance": 3390, "error": 1125, "cr": 1287},
    {"name": "S230624av", "distance": 2124, "error": 682, "cr": 1024},
    {"name": "S230628ax", "distance": 2047, "error": 585, "cr": 705},
    {"name": "S230630am", "distance": 5336, "error": 2001, "cr": 3965},
    {"name": "S230630bq", "distance": 999, "error": 286, "cr": 1215}, 
    {"name": "S230702an", "distance": 2428, "error": 849, "cr": 2267},
    {"name": "S230704f", "distance": 2759, "error": 992, "cr": 1700},
    {"name": "S230706ah", "distance": 1962, "error": 594, "cr": 1497},
    {"name": "S230702an", "distance": 2428, "error": 849, "cr": 2267},
    {"name": "S230704f", "distance": 2759, "error": 992, "cr": 1700},
    {"name": "S230706ah", "distance": 1962, "error": 594, "cr": 1497},
    {"name": "S230707ai", "distance": 4074, "error": 1485, "cr": 3181},
    {"name": "S230708cf", "distance": 3336, "error": 1076, "cr": 2032},
    {"name": "S230708t", "distance": 3010, "error": 988, "cr": 1227},
    {"name": "S230708z", "distance": 4647, "error": 1696, "cr": 3372},
    {"name": "S230708cf", "distance": 3336, "error": 1076, "cr": 2032},
    {"name": "S230709bi", "distance": 4364, "error": 1585, "cr": 2644},
    {"name": "S230723ac", "distance": 1551, "error": 436, "cr": 1117},
    {"name": "S230726a", "distance": 2132, "error": 714, "cr": 27774},
    {"name": "S230729z", "distance": 1495, "error": 444, "cr": 1945}, #issue
    {"name": "S230731an", "distance": 1001, "error": 242, "cr": 598},
    {"name": "S230802aq", "distance": 576, "error": 246, "cr": 25885},
    {"name": "S230805x", "distance": 3305, "error": 1113, "cr": 2094}, 
    {"name": "S230806ak", "distance": 5423, "error": 1862, "cr": 3715},
    {"name": "S230807f", "distance": 5272, "error": 1900, "cr": 5436},
    {"name": "S230811n", "distance": 1905, "error": 672, "cr": 810},
    {"name": "S230814r", "distance": 3788, "error": 1416, "cr": 3389},
    {"name": "S230814ah", "distance": 330, "error": 105, "cr": 25260}, #issue
    {"name": "S230819ax", "distance": 4216, "error": 1645, "cr": 4044},
    {"name": "S230820bq", "distance": 3600, "error": 1437, "cr": 1373}, #issue
    {"name": "S230822bm", "distance": 5154, "error": 1771, "cr": 3974}, 
    {"name": "S230824r", "distance": 4701, "error": 1563, "cr": 3279},
    {"name": "S230825k", "distance": 5283, "error": 2117, "cr": 3012},
    {"name": "S230831e", "distance": 4900, "error": 2126, "cr": 3803},
    {"name": "S230904n", "distance": 1095, "error": 327, "cr": 2015},
    {"name": "S230911ae", "distance": 1623, "error": 584, "cr": 27759},
    {"name": "S230914ak", "distance": 2676, "error": 827, "cr": 1532},
    {"name": "S230919bj", "distance": 1491, "error": 402, "cr": 708}, #missing
    {"name": "S230920al", "distance": 3139, "error": 1003, "cr": 2180}, #missing
    {"name": "S230922g", "distance": 1491, "error": 443, "cr": 332},
    {"name": "S230922q", "distance": 6653, "error": 2348, "cr": 4658},
    {"name": "S230924an", "distance": 2558, "error": 596, "cr": 835},
    {"name": "S230927l", "distance": 2966, "error": 1041, "cr": 1177},
    {"name": "S230927be", "distance": 1059, "error": 289, "cr": 298},
    {"name": "S230928cb", "distance": 4060, "error": 1553, "cr": 3102},
    {"name": "S230930al", "distance": 4902, "error": 1671, "cr": 3166},
    {"name": "S231001aq", "distance": 4425, "error": 1946, "cr": 3181},
    {"name": "S231005j", "distance": 6417, "error": 2246, "cr": 5480},
    {"name": "S231005ah", "distance": 3707, "error": 1335, "cr": 2497},
    {"name": "S231008ap", "distance": 3531, "error": 1320, "cr": 3102},
    {"name": "S231014r", "distance": 2857, "error": 903, "cr": 1807},
    {"name": "S231020ba", "distance": 1168, "error": 361, "cr": 1339}, #issue
    {"name": "S231020bw", "distance": 2620, "error": 694, "cr": 386},
    {"name": "S231028bg", "distance": 4221, "error": 923, "cr": 1207},
    {"name": "S231029y", "distance": 3292, "error": 1313, "cr": 29973}, #issue
    {"name": "S231102w", "distance": 3493, "error": 1015, "cr": 2343}, #issue
    {"name": "S231104ac", "distance": 1357, "error": 321, "cr": 759}, #issue
    {"name": "S231108u", "distance": 1986, "error": 494, "cr": 949},
    {"name": "S231110g", "distance": 1849, "error": 533, "cr": 636},
    {"name": "S231113bb", "distance": 3260, "error": 1181, "cr": 2172},
    {"name": "S231113bw", "distance": 1186, "error": 376, "cr": 1713},
    {"name": "S231114n", "distance": 1317, "error": 407, "cr": 1267},
    {"name": "S231118ab", "distance": 4353, "error": 1588, "cr": 2898},
    {"name": "S231118an", "distance": 1337, "error": 347, "cr": 1107},
    {"name": "S231118d", "distance": 2109, "error": 585, "cr": 956},
    {"name": "S231119u", "distance": 6597, "error": 2556, "cr": 5212},
    {"name": "S231123cg", "distance": 1148, "error": 338, "cr": 2714},
    {"name": "S231127cg", "distance": 4425, "error": 1718, "cr": 3450},
    {"name": "S231129ac", "distance": 3964, "error": 1513, "cr": 3089},
    {"name": "S231206ca", "distance": 3230, "error": 1141, "cr": 2335},
    {"name": "S231206cc", "distance": 1467, "error": 264, "cr": 342},
    {"name": "S231213ap", "distance": 3861, "error": 1257, "cr": 1469},
    {"name": "S231223j", "distance": 4468, "error": 1602, "cr": 3520},
    {"name": "S231224e", "distance": 863, "error": 213, "cr": 394},
    {"name": "S231226av", "distance": 1218, "error": 171, "cr": 199},
    {"name": "S231231ag", "distance": 1066, "error": 339, "cr": 27061},
    {"name": "S240104bl", "distance": 1978, "error": 618, "cr": 27949},
    {"name": "S240107b", "distance": 6089, "error": 2429, "cr": 4143},
    {"name": "S240109a", "distance": 1594, "error": 567, "cr": 28049},
    {"name": "S240104bl", "distance": 1978, "error": 618, "cr": 27949},
    {"name": "S240107b", "distance": 6089, "error": 2429, "cr": 4143},
    {"name": "S240109a", "distance": 1594, "error": 567, "cr": 28049},
]

plot_event_data(events_O4a,"O4a","BBH")

events_O4b = [   
    #{"name": "S240406aj", "distance": 2449, "error": 692, "cr": 1724},
    {"name": "S240413p", "distance": 526, "error": 101, "cr": 34},
    {"name": "S240421ar", "distance": 7702, "error": 2899, "cr": 2601},
    {"name": "S240426dl", "distance": 5886, "error": 2242, "cr": 3469},
    {"name": "S240426s", "distance": 3452, "error": 1295, "cr": 3050},
    {"name": "S240428dr", "distance": 831, "error": 145, "cr": 186},
    {"name": "S240430ca", "distance": 6212, "error": 2593, "cr": 4061},
    {"name": "S240501an", "distance": 4022, "error": 1460, "cr": 1079},
    {"name": "S240505av", "distance": 4570, "error": 1415, "cr": 1469},
    {"name": "S240507p", "distance": 1328, "error": 370, "cr": 279},
    {"name": "S240511i", "distance": 1906, "error": 404, "cr": 85},
    {"name": "S240512r", "distance": 1082, "error": 266, "cr": 216},
    {"name": "S240513ei", "distance": 2254, "error": 458, "cr": 37},
    {"name": "S240514c", "distance": 4182, "error": 1833, "cr": 30758},
    {"name": "S240514x", "distance": 2594, "error": 587, "cr": 142},
    {"name": "S240515m", "distance": 3559, "error": 976, "cr": 978},
    {"name": "S240520cv", "distance": 1289, "error": 332, "cr": 370},
    {"name": "S240525p", "distance": 4337, "error": 1519, "cr": 1517},
    {"name": "S240527en", "distance": 7238, "error": 2059, "cr": 1779},
    {"name": "S240527fv", "distance": 1119, "error": 188, "cr": 15},
    {"name": "S240530a", "distance": 1229, "error": 393, "cr": 984}
]

plot_event_data(events_O4b,"O4b","BBH")


#NSBHS230529ay
events_O4NSBH = [ 
    {"name": "S230529ay", "distance": 201, "error": 100, "cr": 24534}, #credible region in alert while distance from the paper
    {"name": "S230627c", "distance": 291, "error": 64, "cr": 82},
    {"name": "S240422ed", "distance": 188, "error": 43, "cr": 259},
]


events_O4NSBHO4a = [ 
    {"name": "S230529ay", "distance": 201, "error": 100, "cr": 24534}, #credible region in alert while distance from the paper
    {"name": "S230627c", "distance": 291, "error": 64, "cr": 82},
]

events_O4NSBHO4b = [ 
    {"name": "S240422ed", "distance": 188, "error": 43, "cr": 259},
]

plot_event_data(events_O4NSBH,"O4a/b","NSBH")

plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 100000)
plt.ylim(10, 100000)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.xlabel("Distance (Mpc)")
plt.ylabel("90% credible region area (square degrees)")



plt.grid(True)
plt.clf()


plot_cumulative_data(events_O4a,"O4a BBH","distances")
Weizmann_LV_O4()
plot_cumulative_data(events_O4b,"O4b BBH","distances")
#plt.xscale('log')
plt.legend()
plt.show()


plot_cumulative_data(events_O4a+events_O4NSBHO4a,"O4a BBH","cr")
plot_cumulative_data(events_O4b+events_O4NSBHO4b,"O4b BBH","cr")
#plt.xscale('log')
plt.legend()
plt.show()
