## Python file to communicate with paraview to evaluate the results

import sys

import numpy as np
import os

currentFolder = os.path.dirname(os.path.realpath(__file__))

sys.path.append(currentFolder+'/../../EquationOfStatePython/')
import eos

from matplotlib import pyplot as plt



resultsFolder = currentFolder+"/ThermConsistency3D/"



a_tested = "0.01, 0.02, 0.031, 0.041, 0.051, 0.061".split(", ")
a_49     = "0.5, 1.0, 1.5, 2.0, 2.5, 3.0".split(", ")
sigma_tested = ["0.3","0.31","0.32","0.33","0.34","0.35"]

trs = np.linspace(0.6,0.92)
PR = eos.LatticePengRobinson( acentricity = 0.3443, a=1.0/49., b=2./21. )
theory = np.array([ PR.getEquilibriumDensities(tr)[:2] for tr in trs ])

print("Walking through results folder: {}".format(resultsFolder))

for a_current in a_tested:
    for sigma_current in sigma_tested:
        rho_v = []
        rho_l = []
        tr = []
        for root, dirs, files in os.walk(resultsFolder):
            for di in dirs:
                if di==("a_"+a_current) and "sigma_{}/".format(sigma_current) in root:
                    #print("Loading results of:", root)
                    
                    with open(root+"/"+di+"/gnuplotData/data/density.dat") as datafile:
                        for line in datafile:
                            pass
                        last_line = line
                        time, rhoV, rhoL = last_line.split()
                        
                        if int(time) > 5000:
                            print("adding results of ", root, di)
                            
                            tr.append( float(root.split("_")[-1]) )
                            rho_v.append(float(rhoV))
                            rho_l.append(float(rhoL))


## Specific Volume  

        # if a_current == "0.02" and sigma_current == "0.12":
        #     print("Ratio rho_LBM / rho_theory: " )
        #     for i, trc in enumerate(tr):
        #         print("At tr = {}: {:.5f}    | {:.5f}".format(trc, rho_l[i] , np.interp( trc, trs, theory[:,0] )  ) )
        #         print("Ratio for tr = {}: {}".format(trc, rho_l[i]/ np.interp( trc, trs, theory[:,0] )  ) )

        
        plt.clf()                
        plt.plot([1/rho for rho in rho_v],tr, "b-", label = "Simulation")
        plt.plot([1/rho for rho in rho_l],tr, "b-")
        plt.plot(1/theory[:,0], trs, linestyle="--", linewidth = 1, color="tab:gray", label="Theory" )
        plt.plot(1/theory[:,1], trs, linestyle="--", linewidth = 1, color="tab:gray")
        plt.grid(True)
        plt.ylim(0.55,1)
        #plt.xlim(-0.5, 9.5)
        plt.xscale('log')
        plt.title("Thermodynamic consistency: Peng-Robinson EoS ")
        plt.xlabel(r"Simulation specific volume $1/\rho$")
        plt.ylabel(r"Temperature ratio [-]")
        plt.legend(loc="upper right")
        
        ax = plt.gca()
        props = dict( facecolor="white", alpha=0.5)
        plt.text(0.81, 0.82, "$\\sigma$ = {:.2f}\n$a$ = {:.2f}/49".format(float(sigma_current),float(a_49[a_tested.index(a_current)]) ), fontsize=10,
        verticalalignment='top', horizontalalignment='left', transform=ax.transAxes , bbox=props)
        
        plt.savefig(resultsFolder+"maxwell_a_{:.4f}_sigma_{}.jpg".format(float(a_current),sigma_current))
        
##  Density
        
        plt.clf()                
        plt.plot(rho_v,tr, "b-", label = "Simulation")
        plt.plot(rho_l,tr, "b-")
        plt.plot(theory[:,0], trs, linestyle="--", linewidth = 1, color="tab:gray", label="Theory" )
        plt.plot(theory[:,1], trs, linestyle="--", linewidth = 1, color="tab:gray")
        plt.grid(True)
        plt.ylim(0.55,1)
        plt.xlim(-0.5, 9.5)
        #plt.xscale('log')
        plt.title("Thermodynamic consistency: Peng-Robinson EoS")
        plt.xlabel(r"Simulation density $\rho$")
        plt.ylabel(r"Temperature ratio [-]")
        plt.legend(loc="upper right")
        
        ax = plt.gca()
        props = dict( facecolor="white", alpha=0.5)
        plt.text(0.81, 0.82, "$\\sigma$ = {:.2f}\n$a$ = {:.2f}/49".format(float(sigma_current),float(a_49[a_tested.index(a_current)]) ), fontsize=10,
        verticalalignment='top', horizontalalignment='left', transform=ax.transAxes , bbox=props)
        
        plt.savefig(resultsFolder+"maxwell_rho_{:.4f}_sigma_{}.jpg".format(float(a_current),sigma_current))
        
## Ratio
        plt.clf()                
        plt.plot([ rho_l[i]/rho_v[i] for i in range(len(rho_l)) ], tr, "b-", label = "Simulation")
        plt.plot(theory[:,0]/theory[:,1], trs, linestyle="--", linewidth = 1, color="tab:gray", label="Theory" )
        plt.grid(True)
        plt.ylim(0.55,1)
        #plt.xlim(-0.5, 9.5)
        #plt.xscale('log')
        plt.title("Thermodynamic consistency: Peng-Robinson EoS")
        plt.xlabel(r"Density Ratio ")
        plt.ylabel(r"Temperature ratio [-]")
        plt.legend(loc="upper right")
        
        ax = plt.gca()
        props = dict( facecolor="white", alpha=0.5)
        plt.text(0.81, 0.82, "$\\sigma$ = {:.2f}\n$a$ = {:.2f}/49".format(float(sigma_current),float(a_49[a_tested.index(a_current)]) ), fontsize=10,
        verticalalignment='top', horizontalalignment='left', transform=ax.transAxes , bbox=props)
        
        plt.savefig(resultsFolder+"maxwell_ratio_a_{:.4f}_sigma_{}.jpg".format(float(a_current),sigma_current))
            
            
# plt.plot(rhoW, angleW, label=r"Contour = $\rho_w$")
# plt.plot(rhoW, angle4, label=r"Contour = 4")
# plt.plot([0.5,6.5], p([0.5,6.5]), linestyle="--", linewidth = 1, color="tab:gray", label=r"Best fit")
# plt.grid(True)
# plt.ylim(0,180)
# plt.xlim(0,7)
# plt.xlabel(r"Wall density - $\rho_w$")
# plt.ylabel(r"Contact angle - $\theta$ [deg]")
# plt.legend()
# plt.savefig("thesis/wallWettability/wettabilityResults",)
            



#https://arxiv.org/ftp/arxiv/papers/1204/1204.4098.pdf


from matplotlib import pyplot as plt


























