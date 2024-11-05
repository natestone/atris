# AtRIS

AtRIS: Phantoms
============================================= 
Author: Saša Banjac

 + For description and validation, please read the AtRIS validation papers by Banjac et al (2018).
 + Documentation can be found on https://et-wiki.physik.uni-kiel.de/atris/atris
 + Author: Saša Banjac

PHANTOM REQUIREMENTS
====================
 * Except for the first line, every line of the phantom has 650 columns describing the efective ionization of a specific particle
 * The binning is always fixed and cannot be changed. It can be recreated with numpy using: 10**np.linspace(-6,7,651). The unit is MeV. Therefore the range extends from 1eV to 10TeV.
 * The first line of the phantom holds the PDG codes of all particles and nothing else. No header starting "#" sign. Ne ending empty space.
 * PDG codes are separated with a single space. 
 * PDG codes do not repeat
 * PDG order gives the order of the rows corresponding to the particles.
 * Both phantoms must have the same PDG code order. Otherwise it will not work properly.
 

Standard phantoms:
========================
Good for atmospheric aplications at column depths below 10 g/mc^2. Using the following particles:
22 2212 11 1000010020 1000010030 1000020030 1000020040 -11 13 -13 -2212 2112 -2112 211 -211 111 321 -321

For water, I've used .997 g per ccm as desnity, while for TEM it is 1.00 g per ccm

 * icru phantom is a r=15cm sphere of water
   -          composition: h2o
   -               radius: 15 cm
   -                 mass: 
   -    crosssection area:
 * icru_tem phantom is a r=15cm sphere of the icru 4-body soft tissue equivalent material
   -          composition: G4_TISSUE_SOFT_ICRU-4,  composed of 10.1% H, 11.1% C, 2.6% N and 76.2% O 
   -               radius: 15 cm
   -                 mass: 
   -    crosssection area:
 * lbug phantom is the "large bug" phantom with a radius of 1cm and made out of water
   -          composition: h2o
   -               radius: 1 cm
   -                 mass: 
   -    crosssection area:
 * sbug is the "small bug" phantom with r=0.1cm, composed out of water
   -          composition: h2o
   -               radius: 1 mm
   -                 mass: 
   - crosssection area:
 * cell is the "human cell" phantom with r=50 um, made of water
   -          composition: h2o
   -               radius: 0.05 mm
   -                 mass: 
   -    crosssection area:
 * bact is the "bacteria" phantom with r=4 um, made of water
   -          composition: h2o
   -               radius: 0.004 mm
   -                 mass: 
   - crosssection area:
 * dost is the "DOSTEL" detector optimized phantom made of silicon
   -          composition: silicon
   -               radius: 1 mm
   -                 mass: 
   - crosssection area:
 * fred is the "FRED" detector optimized phantom made of silicon
   -          composition: silicon
   -               radius: 1 mm
   -                 mass: 
   -    crosssection area:

Extended phantoms:
==================
These phantoms include response to GCR fully ionized nuclei from Z=3 to Z=28 (most abundant isotope).
 * mslrad_b is a phantom for the MSL/RAD detector B
   -          composition: silicon
   -               radius: 1 mm
   -                 mass: 
   -    crosssection area:
 * icru_b is the icru phantom extended to have the GCR nuclei.
   -          composition: h2o
   -               radius: 15 cm
   -                 mass: 
   -    crosssection area:


Instructions for use
====================
AtRIS calculates dose rate with two phantoms simultaneously. The default phantoms are placed automatically in to
the compilation directory. Furthermore the code knows the mass and crosssection area of the default phantoms. If 
you decide to replace one of the phantoms, you need to also specify the new mass rho*(4./3.)*r**3and the new
 cross section area, which is calculated as pi*r**2. Here are the commands:




List-of-changes
===============

 + DATE:
 + CHAGE:
 
