## Python file to communicate with paraview to evaluate the results

import sys

import numpy as np
import os

currentFolder = os.path.dirname(os.path.realpath(__file__))

sys.path.append(currentFolder+'/../EOS python/')
import eos

from matplotlib import pyplot as plt

UseCS = False

if UseCS:
    resultsDatFile = currentFolder+"/maxwellDenstities_CS/a_0.5/predictedDensities.dat"
else:
    resultsDatFile = currentFolder+"/maxwellDenstities_PR/a_0.01/predictedDensities.dat"

constantRho = 0.1
trs = np.linspace(0.49,0.85)
if UseCS:
    EOS = eos.CarnahanStarlingInteractionPotential( a=0.5, b=4. )
else:
    EOS = eos.LatticePengRobinson( acentricity = 0.3443, a=0.5/49., b=2./21. )
theory = np.array([ EOS.dPdT(constantRho, tr) for tr in trs] )

with open(resultsDatFile) as csvFile:
    lines = csvFile.readlines()
    
tr = []
deriv = []
for line in lines:
    tr0, rho1, rho2, dP = line.strip().split(",")
    tr.append(float(tr0))
    deriv.append(float(dP))

##  Density
plt.clf()                  
plt.plot(tr, deriv, "b-.", label = r"C++ ($\rho$="+str(constantRho)+")")
plt.plot(trs, theory, linestyle="-", linewidth = 1, color="tab:gray", label=r"Python ($\rho$="+str(constantRho)+")" )
plt.grid(True)

if UseCS:
    plt.title("Pressure derivative: Carnahan-Starling EoS C++")
else:
    plt.title("Pressure derivative: Peng-Robinson EoS C++")

plt.xlabel(r"Temperature ratio [-]")
plt.ylabel(r"dP / dT")
plt.legend(loc="upper right")

if UseCS:
    plt.ylim(0.15,0.154)
    plt.tight_layout()
    plt.savefig(currentFolder+"/pressure_deriv_CS.jpg")
else:
    plt.tight_layout()
    plt.savefig(currentFolder+"/pressure_deriv_PR.jpg")

#https://arxiv.org/ftp/arxiv/papers/1204/1204.4098.pdf


from matplotlib import pyplot as plt


























