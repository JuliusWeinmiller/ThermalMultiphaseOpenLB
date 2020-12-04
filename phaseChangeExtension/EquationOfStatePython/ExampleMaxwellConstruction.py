#  Python implementation of the Maxwell Construction
#  used in the OpenLB extension
#    
#  Written by Julius Weinmiller for his thesis
#  https://www.openlb.net/forum/users/julius-weinmiller/
#

import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate

# pressure of the LBM Carnahan Starling EOS
def pressure(rho, tr, a=0.5, b=4.0):
    R = 1.0
    c = b*rho/4.
    tc = 0.18727/0.4963*a/b/R
    t = tr * tc
    return rho*R*t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-a*rho*rho

trs = []
rho_v = []
rho_l = []

for tr in np.linspace(0.49,0.998,50):

    P = lambda v: pressure(1/v, tr)
    
    
    guessLeft, guessRight = 4, 9 # guess is dependent on EOS used
    
    VminX = opt.fmin(               P, guessLeft,  disp=False)
    VmaxX = opt.fmin( lambda x: -P(x), guessRight, disp=False)
    
    PressureGuess = P(VmaxX) /2 + 1e-10
    step = PressureGuess/2.
    if PressureGuess< P(VminX):
        PressureGuess = P(VmaxX)/2 + P(VminX)/2
        step = (PressureGuess - P(VminX) ) /2.
    
    diffArea = 1
    accuracy = 1e-8
    # bounds are dependent on EOS used
    lbound, rbound = 1.01, 50000
    
    while abs(diffArea) > accuracy:
        offsetFunction = lambda V: P(V) - PressureGuess
        
        intersectLeft = opt.brentq( offsetFunction, lbound, VminX )
        intersectRight = opt.brentq( offsetFunction, VmaxX, rbound )
        
        diffArea = integrate.quad(offsetFunction, intersectLeft, intersectRight)[0]
        
        if diffArea<0:
            PressureGuess-=step
        else:
            PressureGuess+=step
        step/=2
        
    print("Rho vapour/liquid: {:.6f} / {:.3f}".format(1/intersectRight, 1/intersectLeft))

    trs.append(tr)
    rho_v.append(1/intersectRight)
    rho_l.append(1/intersectLeft)

import matplotlib.pyplot as plt

plt.xscale("log")
plt.plot(np.array(rho_v),trs, label=" Vapour Branch")
plt.plot(np.array(rho_l),trs, label=" Liquid Branch")

plt.ylabel(r"T / T$_c$")
plt.xlabel(r"$\rho$")
#plt.gca().invert_xaxis()
plt.legend()
plt.title("Maxwell construction: Carnahan-Starling")
plt.grid(True)
plt.savefig("./Maxwell Construction")
plt.show()
