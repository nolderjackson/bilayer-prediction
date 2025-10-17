import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
#Ro will be in angstrom
#extinction are in M-1cm-1
#conc in mol/L
#wavelength in nm

def calcRo6(Jint, kappa2, preRo6factor, QuantumYield, refractiveindex):
    Ro6 = preRo6factor*Jint*kappa2*QuantumYield/(refractiveindex**4)
    return Ro6

def calcForsterR(Jint, kappa2, preRo6factor, QuantumYield, refractiveindex):
    ForsterR = (calcRo6(Jint, kappa2, preRo6factor, QuantumYield, refractiveindex))**(1/6)
    return ForsterR

def calckFRET(r, Jint, kappa2, preRo6factor, QuantumYield, refractiveindex, FluorescentLifetime):
    Ro6 = calcRo6(Jint, kappa2, preRo6factor, QuantumYield, refractiveindex)
    kFRET = ((Ro6/((r)**(6)))*(1/FluorescentLifetime))
    return kFRET


if __name__ == "__main__":

    Na= 6.022e23 #Avogadro's number (mol^-1)
    c= 2.998e8 #Speed of Light in a vacuum (m/s)
    h= 6.626e-34 #Planck's constant (Js)
    hbar= 1.05457e-34 #reduced Planck's constant (Js)
    kB= 1.380649e-23 #Boltzmann's Constant (J/K)
    echarge= 1.602e-19 #elementary unit of charge (C)
    epsnaught= 8.854e-12 #permativity of free space (C^2/Jm)
    eVConv= 1240 #conversion factor of eV and nm (eVnm)
    preRo6factor= 8.875e-5 #factor for calculating Ro6 in A^6 (A^2)
    WavenumbertoJ= 1.986446e-23 #factor for converting cm^-1 to J (Jcm)
    Angstromtom= 1e-10 #factor for converting angstrom to m (m/A)