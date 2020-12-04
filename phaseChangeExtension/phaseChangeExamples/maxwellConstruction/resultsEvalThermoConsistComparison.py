## Python file to communicate with paraview to evaluate the results

import sys

import numpy as np
import os

currentFolder = os.path.dirname(os.path.realpath(__file__))

sys.path.append(currentFolder+'/../EOS python/')
import eos

from matplotlib import pyplot as plt

UseCS = True

if UseCS:
    resultsDatFile = currentFolder+"/maxwellDenstities_CS/a_0.5/predictedDensities.dat"
else:
    resultsDatFile = currentFolder+"/maxwellDenstities_PR/a_0.01/predictedDensities.dat"


trs = np.linspace(0.49,0.85)
if UseCS:
    EOS = eos.CarnahanStarlingInteractionPotential( a=0.5, b=4. )
else:
    EOS = eos.LatticePengRobinson( acentricity = 0.3443, a=0.5/49., b=2./21. )
theory = np.array([ EOS.getEquilibriumDensities(tr)[:2] for tr in trs ])

with open(resultsDatFile) as csvFile:
    lines = csvFile.readlines()
    
tr = []
rho_v = []
rho_l = []
for line in lines:
    tr0, rho1, rho2, deriv = line.strip().split(",")
    tr.append(float(tr0))
    rho_l.append(float(rho1))
    rho_v.append(float(rho2))

##  Density
plt.clf()          
plt.xscale("log")          
plt.plot(rho_v,tr, "b-.", label = "c++")
plt.plot(rho_l,tr, "b-.")
plt.plot(theory[:,0], trs, linestyle="-", linewidth = 1, color="tab:gray", label="Python" )
plt.plot(theory[:,1], trs, linestyle="-", linewidth = 1, color="tab:gray")
plt.grid(True)

if UseCS:
    plt.title("Thermodynamic consistency: Carnahan-Starling EoS C++")
else:
    plt.title("Thermodynamic consistency: Peng-Robinson EoS C++")

plt.xlabel(r"Density $\rho$")
plt.ylabel(r"Temperature ratio [-]")
plt.legend(loc="upper left")

if UseCS:
    plt.savefig(currentFolder+"/maxwell_comparison_CS.jpg")
else:
    plt.savefig(currentFolder+"/maxwell_comparison_PR.jpg")

#https://arxiv.org/ftp/arxiv/papers/1204/1204.4098.pdf


from matplotlib import pyplot as plt


























