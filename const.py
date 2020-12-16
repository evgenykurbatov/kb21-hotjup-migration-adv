# -*- coding: utf-8 -*-

#import scipy.constants as spconst


##
## Fundamental physical constants
##

## Planck's constant [cm^2/s]
h = 6.626e-27
hbar = 1.055e-27
## Gravitational constant [erg cm/g^2]
G = 6.67e-8
## Speed of light [cm/s]
c = 2.99e10
## Electron charge [(erg cm)^(1/2)]
e = 4.803e-10
## Electron mass [g]
m_e = 9.109e-28
## Proton mass [g]
m_p = 1.673e-24
##Neutron mass [g]
m_n = 1.675e-24
## Hydrogen atom mass
m_H = m_p
## Avogagro constant [mol^{-1}]
N_Avogadro = 6.02214076e23
## Daltons [g]
## It is 1/12 of the mass of an unbound neutral atom of carbon-12
## in its nuclear and electronic ground state and at rest.
## It's very close ot 1/N_Avogadro.
## https://en.wikipedia.org/wiki/Dalton_(unit)
m_u = 1.660e-24
## Boltzmann constant [erg/K]
k_B = 1.38e-16


##
## Fundamentsl units conversion constants
##

## Length units
nm = 1e-7        ## 1e-9 m [cm]
Angstrom = 1e-8  ## [cm]
## Energy units
Joile = 1e7      ## kg m^2/s^2 [erg]
Watt = 1e7       ## Joile/s [erg/s]
EV = 1.6022e-12  ## Electron-Volt [erg]
eV = EV
## Force units
Newton = 1e5     ## kg m/s^2 [g cm/s^2]
## Pressure units
Pascal = 10.0    ## J/m^3 [erg/cm^3]
mbar = 1e3       ## 100 Pascal [erg/cm^3]


##
## Derived physical constants
##

## Bohr radius [cm]
## hbar^2 / (m_e e^2)
a_Bohr = 5.2917721067e-9
## Universal gas constant (in MKT) [erg/mol/K]
## k_B N_Avogadro
R_gas = 8.31446262e7
## Universal gas constant [erg/g/mol/K]
## R_gas/m_H/N_Avogadro
RR_gas = 8.2525343e7

## Stefan-Boltzmann constant [erg/cm^2/s/K^4]
## 2 pi^5 k_B^4 / (15 c^2 h^3)
sigma_SB = 5.6704e-5
## Radiation constant (or radiation density constant) [erg/cm^3/K^4]
## 4 sigma_SB/c
a_rad = 7.58e-15

## Elemantary scatter section in bound-free transition of hydrogen atom [cm^2]
sigma_bf_1 = 7.92e-18


##
## Astronomical constants
##

## Solar units
M_sol = 1.99e33
R_sol = 6.96e10
L_sol = 3.827e33
Magn_sol_V   = 4.83
Magn_sol_B   = 5.48
Magn_sol_U   = 5.61
Magn_sol_bol = 4.75

## Planets
M_jup = 9.54e-4 * M_sol
R_jup = 6.99e9
M_earth = 3.00e-6 * M_sol
R_earth = 6.37e8

## Parsec
pc = 3.0857e18
## Astronomical Unit
AU = 1.4960e13
au = AU
## Year
year = 3.1557e7
yr = year

## Energy flux density unit in radio astronomy [(erg/s)/(cm^2*Hz)]
## In SI it's equals to 10e-26 watts per square metre per hertz
## https://en.wikipedia.org/wiki/Jansky
Jansky = 1e-23
