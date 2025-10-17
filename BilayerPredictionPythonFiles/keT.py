import numpy as np
from numpy import pi, sqrt

def calcOutReorgEnergy(rD, rA, R, n, dielectric, echarge, eps0, Na, Angtom):
    rD = rD * Angtom
    rA = rA * Angtom
    R = R * Angtom
    #m, m, m, AU, AU, C, C^2/Nm^2, mol^-1

    p1 = (echarge**2)/2 #C^2
    p2 = (1/(n**2) - (1/dielectric))  #AU
    p3 = (1/rD + 1/rA - 2/R) # m^-1
    OutReorgEnergy = Na*(p1 * p2 * p3 / (4*pi*eps0)) #mol^-1 (C^2 * m^-1 * Nm^2/C^2) = J/mol
    return OutReorgEnergy

def calcDeltaGeT(Dpot, Apot, R, E00, dielectric,
                  echarge, eps0, Na, Angtom):
    #V, V, Ang, eV, AU, C, C^2/Nm^2, AU
    R = R * Angtom
    #V, V, m, eV, AU, C, C^2/Nm^2, AU, mol^-1

    emf    = echarge*(Dpot - Apot)                        # J 
    esWork = (echarge**2)/(4*pi*eps0 * R * dielectric)  # J
    E00Energy = E00*echarge #J
    # convert photon energy to mol^-1
    GeT = Na * (emf + esWork - E00Energy)
    return GeT

def calcCouplingAttenuation(coupling, AttFac, rD, rA, R, wavetoJ):
    # cm^-1, Å^-1, Å, Å, Å
    coupling = coupling*wavetoJ
    # J, Å^-1, Å, Å, Å
    """
    Exponential distance attenuation:
      J = J0 * exp[ -AttFac * (R - (rD + rA)) ]
    """
    return coupling * np.exp(-(AttFac/2) * (R-(rD+rA)))

def calckeT(Dpot, Apot, R, rD, rA, E00, coupling, AttFac, n, dielectric, temp,
                echarge, kB, hbar, Na, eps0, Angtom, wavetoJ):
    J      = calcCouplingAttenuation(coupling, AttFac, rD, rA, R, wavetoJ)
    λ       = calcOutReorgEnergy(rD, rA, R, n, dielectric, echarge, eps0, Na, Angtom)
    ΔG_et   = (calcDeltaGeT(Dpot, Apot, R, E00, dielectric, echarge, eps0, Na, Angtom))
    # print(ΔG_et, λ)
    prefac  = (2*pi/hbar) * (1/ sqrt(4*pi*kB*temp*(λ/Na)))
    # print(prefac)
    expfac  = (-((λ + ΔG_et)**2) / (4*λ*kB*temp))/(Na)
    # print(expfac)
    keT     = prefac * J**2 * np.exp(expfac)
    return keT, λ, ΔG_et, (J/wavetoJ)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
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

