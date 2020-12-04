## Python file to communicate with paraview to evaluate the results

import sys
import numpy as np
import os

currentFolder = os.path.dirname(os.path.realpath(__file__))

from matplotlib import pyplot as plt


resultsFolder = currentFolder+"/laplace2D_results/"
diameter = "0.1, 0.2, 0.3, 0.4, 0.5, 0.6".split(", ")

p_v = [ 0.000681442, 0.000638102, 0.000623126, 0.000615578,  0.000611031, 0.000608056]
p_l = [ 0.00695537,  0.00393553,  0.00285009,  0.00230233,   0.00193717,  0.00176017]

length = np.array([ 0.000171755, 0.000330497, 0.000489816, 0.00065016, 0.000838713, 0.000982372])


print( (np.array(p_l) - np.array(p_v))*1e3 )

pDifference = (np.array(p_l) - np.array(p_v) )
# dP = sigma * 1/r
domainSize = 5*125
r = np.array( [ float(d) for d in diameter  ] ) * domainSize/ 2.0
# sigma = dP * r 
surfaceTension = pDifference * r 

# print(r)

convU = 125
dx = 8e-7
convRho = 947.13/8.72518
convM = convRho * dx**3

radius = length/dx / (2*np.pi)
r = radius
# print(radius)

# sigma_wiki = 58.85*1e-3
rs = np.linspace(r[0],r[-1])
# ptheory = sigma_wiki * 2 / rs

z = np.polyfit(1/(r), pDifference, 1)
p = np.poly1d(z)

plt.clf()                
plt.plot( 1/r , pDifference, "b", linestyle="", marker="+", markersize=10, label = "Simulation")
plt.plot([0, 0.04], p([0,0.04]), linestyle="--", linewidth = 1, color="tab:gray", label=r"Linear best fit")

# plt.plot( rs*1e6, ptheory/1e3, linestyle="--", linewidth = 1, color="tab:gray", label="Theory" )
plt.grid(True)

plt.title("Laplace pressure")
plt.xlabel(r"Inverse droplet radius $\frac{1}{R_d}$ [Lattice units]")
plt.ylabel(r"Pressure difference $\Delta p$ [Lattice units]")
plt.legend(loc="upper left")
# plt.xscale('log')
plt.xlim(0)
plt.ylim(0)
plt.tight_layout()
plt.savefig(resultsFolder+"laplace_pressure2.jpg")
        
## Surface tension 
        
avg = np.average(surfaceTension)
print("Avg surface tension (lu): ", avg)
print("Avg surface tension (pu): ", avg*convU**2 *convM/dx/dx )
plt.clf()                
plt.plot( r ,surfaceTension, "b", linestyle="",  marker="x", label = "Simulation")
plt.plot([0, 190 ], [ avg,avg ] , linestyle="--", linewidth = 1, color="tab:gray", label=r"Average")

# plt.plot(1/theory[:,0], trs, linestyle="--", linewidth = 1, color="tab:gray", label="Theory" )
plt.grid(True)

plt.title("Surface tension")
plt.xlabel(r"Droplet radius [Lattice units]")
plt.ylabel(r"Surface tension [Lattice units] ")
plt.legend(loc="upper left")        
plt.xlim(0)
plt.ylim(0,0.3)
plt.savefig(resultsFolder+"laplace_surface_tension2.jpg")

















