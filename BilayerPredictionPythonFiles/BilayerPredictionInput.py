
Processing = {
    'Subdivisions': 10, #number of stuctures to sample for BL and TL recomended range (20-200)
    'Directory': rf'C:\[venv dir]', #change this to your working directory <--------
    'AngleInputName': 'AngleDict.csv', # name of angles list produced by AllAnglesDefinemuBL.py
    'BilayerOutput': 'BilayerData.csv', #name of outputfile with averaged data
    'RawBilayerOutput': 'RawData.csv' #name of outputfile without averaged data
}

Donor = {
    'SMILESCore': "C12=CC=CC=C1C=C3C(C=CC=C3)=C2",   #structural formula in SMILES format
    'theta_nmuBL': 31,    #angle (deg)
    'theta_muBLrBL': 75,    #angle of bonding axis to dipole(deg)
    'length': 13.2,        #length from dipole center to binding group (A)
    'DipoleLength': 6.6,
    'QuantumYield': 0.14,  #Quantum yield [0,1]
    'E00': 2.955,           #E00 of Donor (eV)
    'OxidationPotential': 0.9, #Reduction potential of radical cation to neutral (V)\
    'FluorescentLifetime': 11.1e-9 #fluorescent lifetime in the absence of acceptor (s)
}

Acceptor = { 
    "SMILESCore": "F[B-]([N+]1=CC=CC1=C2)(F)N3C2=CC=C3",      #structural formula in SMILES format
    'theta_muTLrTL': 90,    #angle of bonding axis to dipole (deg)
    'length': 6.1,        #length from dipole center to binding group (A)
    'DipoleLength': 3.7,
    'ReductionPotential': -1.12     #Reduction potential of radical cation to neutral (V)
}

Conditions = {
    'Temperature': 298, #Temperature (K)
    'Solvent': 'Acetonitrile', # edit "Custom" in solvent dictionary if a new one needs to be added
    'ElectronicCoupling': 75, #electronic coupling of the donor and acceptor at contact radius (cm^-1)
    'AttenuationFactor': 0.3, #distance dependant attenuation factor of electronic coupling (A^-1)
    'LinkingIonRadius': 2.0, #radius of linking ion (A)
    'J_integral': 2.52e14 #J itegral as defined in the FRET equation
}

Constants = {
    'Na': 6.022e23, #Avogadro's number (mol^-1)
    'c': 2.998e8, #Speed of Light in a vacuum (m/s)
    'h': 6.626e-34, #Planck's constant (Js)
    'hbar': 1.05457e-34, #reduced Planck's constant (Js)
    'kB': 1.380649e-23, #Boltzmann's Constant (J/K)
    'echarge': 1.602e-19, #elementary unit of charge (C)
    'epsnaught': 8.854e-12, #permativity of free space (C^2/Jm)
    'eVConv': 1240, #conversion factor of eV and nm (eVnm)
    'preRo6factor': 8.875e-5, #factor for calculating Ro6 in A^6 (A^2)
    'WavenumbertoJ': 1.986446e-23, #factor for converting cm^-1 to J (Jcm)
    'Angstromtom': 1e-10 #factor for converting angstrom to m (m/A)
}

Solvents = {
    'Custom':{
        'RefractiveIndex': 1,
        'Dielectric': 1,
    },
    'Acetone':{
        'RefractiveIndex': 1.3587, 
        'Dielectric': 20.7
    },
    'Acetonitrile':{
        'RefractiveIndex': 1.3441, 
        'Dielectric': 37.5
    },
    'Benzene':{
        'RefractiveIndex': 1.5011, 
        'Dielectric': 2.27
    },
    'Chloroform':{
        'RefractiveIndex':1.4458 , 
        'Dielectric': 4.81
    },
    'DCM':{
        'RefractiveIndex': 1.4241, 
        'Dielectric': 8.39
    },
    'Diethyl Ether':{
        'RefractiveIndex': 1.3524, 
        'Dielectric': 4.33
    },
    'DMF':{
        'RefractiveIndex': 1.4384, 
        'Dielectric': 37.8
    },
    'DMSO':{
        'RefractiveIndex': 1.4783, 
        'Dielectric': 46.7
    },
    'Dioxane':{
        'RefractiveIndex': 1.4224, 
        'Dielectric': 2.25
    },
    'EtOH':{
        'RefractiveIndex': 1.3614, 
        'Dielectric': 24.5
    },
    'EtOAc':{
        'RefractiveIndex': 1.3724, 
        'Dielectric': 6.02
    },
    'IPA':{
        'RefractiveIndex': 1.3772, 
        'Dielectric': 17.9
    },
    'MeOH':{
        'RefractiveIndex': 1.3284, 
        'Dielectric': 32.7
    },
    'Nitrobenzene':{
        'RefractiveIndex': 1.5562, 
        'Dielectric': 34.82
    },
    'Nitromethane':{
        'RefractiveIndex': 1.3817, 
        'Dielectric': 35.87
    },
    'Pyridine':{
        'RefractiveIndex': 1.5102, 
        'Dielectric': 12.4
    },
    'THF':{
        'RefractiveIndex': 1.4072, 
        'Dielectric': 7.58
    },
    'Toluene':{
        'RefractiveIndex': 1.4969, 
        'Dielectric': 2.38
    },
    'Water':{
        'RefractiveIndex': 1.333, 
        'Dielectric': 80.1
    }
}
