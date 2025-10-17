import pandas as pd
import os

import numpy as np

from AngleBuilding import ConstructBilayer, calcmuD, calcmuA, calcKappa2
from keT import calckeT
from kFRET import calckFRET
from VdwRadius import get_vdw_volume_from_smiles
from BilayerAnalysisInput import Processing, Donor, Acceptor, Conditions, Constants, Solvents

# --- Processing ---
Directory = Processing['Directory']
Angle_Input_CSV_Filename = Processing['AngleInputName']
Data = Processing['BilayerOutput']
RawData = Processing['RawBilayerOutput']
Subdivisions = Processing['Subdivisions']  #Scales O(n^3)

# --- Constants ---
Na = Constants['Na']
c = Constants['c']
h = Constants['h']
hbar = Constants['hbar']
kB = Constants['kB']
echarge = Constants['echarge']
epsnaught = Constants['epsnaught']
eVConv = Constants['eVConv']
preRo6factor = Constants['preRo6factor']
Wavenumber_to_J = Constants['WavenumbertoJ']
Angstrom_to_m = Constants['Angstromtom']

# --- Conditions ---
Temperature = Conditions['Temperature']
Solvent = Conditions['Solvent']
Refractive_Index = Solvents[f'{Solvent}']['RefractiveIndex']
Dielectric = Solvents[f'{Solvent}']['Dielectric']
print(f'{Solvent=}')
print(f'{Refractive_Index=}, {Dielectric=}')
Electronic_Coupling = Conditions['ElectronicCoupling']
Attenuation_Factor = Conditions['AttenuationFactor']
Linking_Ion_Radius = Conditions['LinkingIonRadius']

# --- Donor Properties ---
Donor_Core_SMILES = Donor['SMILESCore']
theta_nmuD = Donor['theta_nmuBL']
theta_muDD = Donor['theta_muBLrBL']
Donor_Length = Donor['length']
Donor_Dipole_Length = Donor['DipoleLength']
Quantum_Yeild = Donor['QuantumYield']
E00 = Donor['E00']
Donor_Oxidation_Potential = Donor['OxidationPotential']
Fluorescent_Lifetime = Donor['FluorescentLifetime']

# --- Acceptor Properties ---
Acceptor_Core_SMILES = Acceptor['SMILESCore']
theta_muAA = Acceptor['theta_muTLrTL']
Acceptor_Length = Acceptor['length']
Acceptor_Dipole_Length = Acceptor['DipoleLength']
Acceptor_Reduction_Potential = Acceptor['ReductionPotential']
# ----------------------


# --- Import Data ---

df1 = pd.read_csv(os.path.join(Directory, Angle_Input_CSV_Filename))


# --- Precompute Necessary Values ---

Donor_Length += Linking_Ion_Radius
Donor_Dipole_Length += Linking_Ion_Radius
Acceptor_Length += Linking_Ion_Radius
Acceptor_Dipole_Length += Linking_Ion_Radius

theta_muBLrBL_rad = np.deg2rad(theta_muDD)
theta_muTLrTL_rad = np.deg2rad(theta_muAA)

Donor_Radius = get_vdw_volume_from_smiles(Donor_Core_SMILES)
Acceptor_Radius = get_vdw_volume_from_smiles(Acceptor_Core_SMILES)

J_integral = Conditions['J_integral']

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm
def calcAngleSet(row):
    r_length_valid = []
    kappa2_valid = []
    kFRET_all = []
    k_eT_all = []
    Reorganization_Energy_all = []
    G_eT_all = []
    Attenuated_Coupling_all = []
    Bottom_FRET_all = []
    Top_FRET_all = []
    Bottom_eT_all = []
    Top_eT_all = []

    theta_beta = row['beta']
    theta_nmuBL = row['theta_nmuBL']
    theta_nmuTL = row['theta_nmuTL']

    theta_beta_rad = np.deg2rad(theta_beta)
    theta_nmuBL_rad = np.deg2rad(theta_nmuBL)
    theta_nmuTL_rad = np.deg2rad(theta_nmuTL)


    muBL = calcmuD(theta_nmuBL_rad)
    muTL = calcmuA(theta_beta_rad, theta_nmuBL_rad, theta_nmuTL_rad)
 
    r_vector_all, r_length_all, Bottom_Vec_all, Top_Vec_all = ConstructBilayer(Donor_Length, Acceptor_Length, Donor_Dipole_Length, Acceptor_Dipole_Length,
                                                  muBL, muTL, theta_muBLrBL_rad, theta_muTLrTL_rad, Subdivisions)
 
    kappa2_all = calcKappa2(r_vector_all, theta_beta_rad, muTL, muBL)
 
    for index in range(len(r_length_all)):
        if r_length_all[index] < (Donor_Radius + Acceptor_Radius):
 
             continue
        else:
            kFRET = calckFRET(r_length_all[index], J_integral, kappa2_all[index], preRo6factor, Quantum_Yeild, Refractive_Index, Fluorescent_Lifetime)

            k_eT, Reorganization_Energy, G_eT, Attenuated_Coupling = calckeT(Donor_Oxidation_Potential, Acceptor_Reduction_Potential, r_length_all[index], Donor_Radius, Acceptor_Radius,
                                                                             E00, Electronic_Coupling, Attenuation_Factor, Refractive_Index, Dielectric, Temperature, echarge,
                                                                             kB, hbar, Na, epsnaught, Angstrom_to_m, Wavenumber_to_J)
            kFRET_all.append(kFRET)
            k_eT_all.append(k_eT)
            Reorganization_Energy_all.append(Reorganization_Energy)
            G_eT_all.append(G_eT)
            Attenuated_Coupling_all.append(Attenuated_Coupling)
            r_length_valid.append(r_length_all[index])
            kappa2_valid.append(kappa2_all[index])
            Bottom_FRET_all.append(Bottom_Vec_all[index])
            Top_FRET_all.append(Top_Vec_all[index])
            Bottom_eT_all.append(Bottom_Vec_all[index])
            Top_eT_all.append(Top_Vec_all[index])

    Bottom_FRET = Bottom_FRET_all[kFRET_all.index(max(kFRET_all))]
    Top_FRET    = Top_FRET_all[kFRET_all.index(max(kFRET_all))]
    Bottom_eT   = Bottom_eT_all[k_eT_all.index(max(k_eT_all))]
    Top_eT      = Top_eT_all[k_eT_all.index(max(k_eT_all))]
    kFRET_avg = np.average(kFRET_all)
    k_eT_avg = np.average(k_eT_all)
    Reorganization_Energy_avg = np.average(Reorganization_Energy_all)
    G_eT_avg = np.average(G_eT_all)
    Attenuated_Coupling_avg = np.average(Attenuated_Coupling_all)
    kappa2_avg = np.average(kappa2_valid)
    r_length_avg = np.average(r_length_valid)

    r_length_valid = [float(item) for item in r_length_valid]
    kappa2_valid = [float(item) for item in kappa2_valid]
    kFRET_all = [float(item) for item in kFRET_all]
    k_eT_all = [float(item) for item in k_eT_all]
    results = [
        kFRET_avg,
        k_eT_avg,
        Reorganization_Energy_avg,
        G_eT_avg,
        Attenuated_Coupling_avg,
        kappa2_avg,
        r_length_valid,
        kappa2_valid,
        kFRET_all,
        k_eT_all,
        r_length_avg,
        Bottom_FRET,
        Top_FRET,
        Bottom_eT,
        Top_eT
    ]

    return results

results = df1.apply(calcAngleSet, axis=1, result_type='expand')
results.columns = [
    'kFRET_avg', 'k_eT_avg', 'Reorganization_Energy_avg', 'G_eT_avg', 'Attenuated_Coupling_avg',
    'kappa2_avg', 'r_length_valid', 'kappa2_valid', 'kFRET_all', 'k_eT_all', 'r_length_avg', 'Bottom_FRET', 'Top_FRET', 'Bottom_eT', 'Top_eT'
]

df1 = pd.concat([df1, results[['kFRET_avg', 'k_eT_avg', 'Reorganization_Energy_avg', 'G_eT_avg', 'Attenuated_Coupling_avg', 'kappa2_avg', 'r_length_avg', 'Bottom_FRET', 'Top_FRET', 'Bottom_eT', 'Top_eT']]], axis=1)
df2 = results[['r_length_valid', 'kappa2_valid', 'kFRET_all', 'k_eT_all']]

print(df1)
print(df2)

df1.to_csv(os.path.join(Directory, Data), index=False)
df2.to_csv(os.path.join(Directory, RawData), index=False)

# Flatten rlen vs kFRET and kET into a single list of points
r_all = []
kFRET_all = []
kET_all = []
kappa2_all = []
for _, row in df2.iterrows():
    r_all.extend(row['r_length_valid'])
    kFRET_all.extend(row['kFRET_all'])
    kET_all.extend(row['k_eT_all'])
    kappa2_all.extend(row['kappa2_valid'])

# Ensure numeric conversion
df1['kFRET_avg'] = pd.to_numeric(df1['kFRET_avg'], errors='coerce')
df1['k_eT_avg'] = pd.to_numeric(df1['k_eT_avg'], errors='coerce')

# Log-transform with filtering to avoid -inf from log(0)
z1_log = np.log10(df1['kFRET_avg'].where(df1['kFRET_avg'] > 0))
z2_log = np.log10(df1['k_eT_avg'].where(df1['k_eT_avg'] > 0))

x = df1['beta']
y = df1['theta_nmuTL']

# --- Plot setup ---
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
import matplotlib as mpl

# Create grid
z1_log = np.log10(df1['kFRET_avg'].where(df1['kFRET_avg'] > 0))
z2_log = np.log10(df1['k_eT_avg'].where(df1['k_eT_avg'] > 0))

x = df1['beta']
y = df1['theta_nmuTL']
# Create grid
xi = np.linspace(min(x), max(x), 500)
yi = np.linspace(min(y), max(y), 500)
xi, yi = np.meshgrid(xi, yi)

# Interpolate z1 and z2 values over the grid
zi1 = griddata((x, y), z1_log, (xi, yi), method='cubic')
zi2 = griddata((x, y), z2_log, (xi, yi), method='cubic')


# Interpolated difference
diff = zi1 - zi2

# Mask invalid values (NaNs)
valid_mask = ~np.isnan(diff)
valid_diff = diff[valid_mask]

# Ensure range includes 0
vmin = min(valid_diff.min(), -.00001)
vmax = max(valid_diff.max(), 0.00001)

# Define colormap normalization centered at 0
divnorm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

plt.rcParams['savefig.dpi'] = 600

fig = plt.figure(figsize=(8, 6))
gs = fig.add_gridspec(2, 2, wspace=0.1, hspace=0.3)

# --- Plot 1: log10(k_FRET_avg) ---
ax1 = fig.add_subplot(gs[0, 0], projection='3d')
surf1 = ax1.plot_surface(xi, yi, zi1, cmap='viridis', edgecolor='none')
ax1.set_title(r'$k_{FRET}$ $(s^{-1})$')
ax1.set_xlabel(r'$\beta(\degree)$')
ax1.set_ylabel(r'$\theta_{n \mu_{TL}}(\degree)$')
cbar1 = fig.colorbar(surf1, ax=ax1, shrink=0.7, pad=0.1)
cbar1.set_label(r'$\log{k_{FRET}}$', fontsize=9)

# --- Plot 2: log10(k_et_avg) ---
ax2 = fig.add_subplot(gs[0, 1], projection='3d')
surf2 = ax2.plot_surface(xi, yi, zi2, cmap='plasma', edgecolor='none')
ax2.set_title(r'$k_{eT}$ $(s^{-1})$',)
ax2.set_xlabel(r'$\beta(\degree)$')
ax2.set_ylabel(r'$\theta_{n \mu_{TL}}(\degree)$')
cbar2 = fig.colorbar(surf2, ax=ax2, shrink=0.7, pad=0.1)
cbar2.set_label(r'$\log{k_{eT}}$', fontsize=9)

# --- Plot 3: both surfaces ---
ax3 = fig.add_subplot(gs[1, 0], projection='3d')
surf3_1 = ax3.plot_surface(xi, yi, zi1, cmap='viridis', alpha=0.7, edgecolor='none')
surf3_2 = ax3.plot_surface(xi, yi, zi2, cmap='plasma', alpha=0.7, edgecolor='none')
ax3.set_title(r'$k_{FRET}$ and $k_{eT}$ $(s^{-1})$ ')
ax3.set_xlabel(r'$\beta(\degree)$')
ax3.set_ylabel(r'$\theta_{n \mu_{TL}}(\degree)$')
ax3.set_zlabel(r'$\log{(k)}$')

# --- Plot 4: 2D topology ---
divnorm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
ax4 = fig.add_subplot(gs[1, 1])
contour = ax4.contourf(xi, yi, diff, levels=100, cmap='seismic', norm=divnorm)
ax4.set_title(r'$\log{k_{FRET}} - \log{k_{eT}}$')
ax4.set_xlabel(r'$\beta(\degree)$')
ax4.set_ylabel(r'$\theta_{n \mu_{TL}}(\degree)$')

linear_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
sm = mpl.cm.ScalarMappable(norm=linear_norm, cmap='seismic')
cbar4 = fig.colorbar(sm, ax=ax4, shrink=0.9, pad=0.05)
# Apply TwoSlopeNorm colors to match plots
cbar4.solids.set_cmap('seismic')
cbar4.solids.set_norm(divnorm)

cbar4.ax.minorticks_on()
cbar4.set_label(r'$\log{k_{FRET}} - \log{k_{eT}}$', fontsize=9)
ax4.set_aspect('equal', adjustable='box')
major_ticks = cbar4.get_ticks()
if len(major_ticks) > 1:
    minor_spacing = (major_ticks[1] - major_ticks[0]) / 2
    cbar4.ax.yaxis.set_minor_locator(mticker.MultipleLocator(minor_spacing))
    cbar4.ax.yaxis.set_minor_formatter(mticker.NullFormatter())
fig.subplots_adjust(top=0.95, bottom=0.08, left=0.05, right=0.95)
plt.show()