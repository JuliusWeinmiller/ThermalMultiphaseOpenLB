#  Python implementation of the Maxwell Construction
#  used in the OpenLB extension
#    
#  Written by Julius Weinmiller for his thesis
#  https://www.openlb.net/forum/users/julius-weinmiller/
#

import numpy as np
import matplotlib.pyplot as plt

import scipy.optimize as opt
import scipy.integrate as integrate

import helperFunc

import sys
import time


####
# Relevant info
#
# https://arxiv.org/ftp/arxiv/papers/1611/1611.00907.pdf
# http://staff.ustc.edu.cn/~huanghb/Huang_Krafczyk_comment.pdf
# https://arxiv.org/ftp/arxiv/papers/1211/1211.6932.pdf
# https://doi.org/10.4208/cicp.OA-2016-0112
# https://doi.org/10.1063/1.2187070
# https://www.sciencedirect.com/science/article/pii/S0898122109001011
# https://www.researchgate.net/publication/323711628_Modeling_thermal_inkjet_and_cell_printing_process_using_modified_pseudopotential_and_thermal_lattice_Boltzmann_methods
# 10.1142/S0129183117501200
# http://verl.npre.illinois.edu/Documents/C-08-02.pdf
# http://dns2.asia.edu.tw/~ysho/YSHO-English/2000%20Engineering/PDF/Ind%20Eng%20Che%20Fun15,%2059.pdf

class InteractionPotential:
    
    def __init__(self):
        
        self._name="BaseClass"
        self.G = None
        self.tc = None
        self.pc = None
        
        self.rhoc = None
        
        
    def pressure(self, rho, tr=None):
        return 1.0
    
    def getName(self):
        return self._name
    
    def psi(self, rho, tr = None):
        "return the pseudopotential, given a density and optionally a temperature ratio" 
        return np.sqrt( 6*( self.pressure(rho,tr) -rho/3. ) / self.G )
    

    def plotPV(self, tr=None):
        self._plotPV(tr)
    
    def _plotPV(self, tr=None, rhoLim = [None, None], ylim = (None,None)):
        
        if rhoLim[0] and rhoLim[1]:
            rhoStart, rhoEnd = rhoLim
        else:
            rhoStart = 1e-3
            rhoEnd = 5e-1
            
        if not ylim[0] and not ylim[1]:
            ylim=(-0.02,0.05)
        
        rhos = np.linspace( rhoStart , rhoEnd ,1000,endpoint=True)
        P = self.pressure(rhos,tr)
        
        plt.xscale("log")
        plt.plot(1/rhos,P, label=self._name+" "+str(tr))
        plt.ylim( ylim )
        plt.ylabel(r"P")
        plt.xlabel(r"v (1/$\rho$)")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        
    def _plotPT(self, rho=None, xLim = [0.5, 1], ylim = (None,None)):
        
        Tr_start, Tr_end = xLim
            
        # if not ylim[0] and not ylim[1]:
        #     ylim=(0,0.2)
        
        Trs = np.linspace( Tr_start, Tr_end, 100, endpoint=True)
        dTr = Trs[1:] - Trs[:-1]
        Ps = np.array([self.pressure(rho,Tr) for Tr in Trs])
        dP = Ps[1:] - Ps[:-1]
        
        dPdT = dP/(dTr*self.tc)
        
        # plt.xscale("lin")
        plt.plot( Trs[1:]-dTr/2, dPdT, label=r"Numerical dP / dT ($\rho$="+str(rho)+")")
        plt.ylim( ylim )
        plt.ylabel(r"dP / dT")
        plt.xlabel(r"T")
        plt.legend()
        plt.grid(True)
        #plt.show()
        
        
    def initNicePrint(self,name="", verbose=True):
        if verbose:
            return lambda *x: print(name+":",*x)
        else:
            return lambda *x: ""
                
        
    def plotPsi(self):
        rhoStart = 0.001
        rhoEnd = 0.5
        rhos = np.linspace( rhoStart , rhoEnd ,100,endpoint=True)
        P = self.psi(rhos)
        
        plt.xscale("log")
        plt.plot(rhos,P,label=self._name)
        plt.legend()
        
    def setTempRatio(self,tr):
        self.t = tr * self.tc
       
    def _superGetEquilibriumDensities(self, tr, lbound, rbound, GuessVc = None , accuracy = 1e-10, iterMax = 200, verbose = False):
        
        qPrint = self.initNicePrint("Calc Eq Densities",verbose)
        
        self.setTempRatio(tr)
        
        # Since we need to integrate over the volume, not density
        P = lambda v: self.pressure(1/v)
        
        if GuessVc:
            Vc = GuessVc
        else:
            Vc = 1/self.rhoc
        
        qPrint("Critical volume: ", Vc)
        ####
        # Pressure Minimum and Maximum
        
        # Get volume at which P is minimum
        VminX = opt.fmin( P, 1 , xtol=1e-12, ftol=1e-12,  disp=False)[0]
        # Get volume at which P is maximum (negative minimum)
        VmaxX = opt.fmin( lambda x: -P(x), Vc, xtol=1e-8, ftol=1e-8, disp=False )[0]
        
        VminX2 = helperFunc.findMin( P, [lbound, Vc] )
        VmaxX2 = helperFunc.findMin(lambda x: -P(x), [Vc, 10000] )
        qPrint("V min/max:   ",VminX,VmaxX)
        qPrint("V min2/max2: ",VminX2,VmaxX2)
        
        VminX = VminX2
        VmaxX = VmaxX2
        #Simulations break if this is true
        assert VminX < VmaxX
        
        ####
        # Initial pressure guess
        # - Has to be higher than 0 to ensure 3 intersects
        # - Has to be between Pmin and Pmax
        # - Should not allow for stepsize to bring it below 0 -> 2 intersects only
        
        PressureGuess = P(VmaxX) /2 + 1e-100
        step = (PressureGuess) /2.
        
        if PressureGuess< P(VminX):
            PressureGuess = P(VmaxX)/2 + P(VminX)/2
            step = (PressureGuess - P(VminX) ) /2.
        
        ####
        # Perform Maxwell Construction
        # - Find equilibrium Pressure such that 
        #       integral between intersects 1&2 is the same as 
        #       integral between intersects 2&3
        
        #t1 = time.time()
        
      
        diffArea = 1
        it = 0
            
        while abs(diffArea) > accuracy:
            
            ##
            # Find intesects:       Pguess = P(V)
            #  - same as find roots: 0 = P(V) - Pguess
            
            offsetFunction = lambda V: P(V) - PressureGuess
            
            # opt.brentq only finds roots 
            # -> with offsetFunction same as intersects
            intersect1 = opt.brentq( offsetFunction, lbound, VminX )
            #intersect2 = opt.brentq( offsetFunction, VminX, VmaxX)
            intersect3 = opt.brentq( offsetFunction, VmaxX, rbound)

            # calculate area enclosed by line Pguess, intersect 1 and interset 2
            #area1 = integrate.quad(offsetFunction, intersect1, intersect2)[0]
            #area2 = -integrate.quad(offsetFunction, intersect2, intersect3)[0]
            area3 = integrate.quad(offsetFunction, intersect1, intersect3)[0]
            
            
            # use difference in area to determine whether Pguess was too low or high
            #diffArea = area1-area2
            #print(area3, diffArea)
            diffArea = area3
            
            if it==0:
                qPrint("Initial pressure is:  ", PressureGuess)
                qPrint("Initial area diff is: ", diffArea)
            
            # Use bisectional search to find P
            # Not the fastest but good enough
            if diffArea<0:
                PressureGuess-=step
            else:
                PressureGuess+=step
            step/=2


            
            # Check for max iterations            
            it+=1
            if it>=iterMax:
                print("Max iterations have been reached.")
                break
            
        # simplification:
        #t2 = time.time()
        
        
        # calcAreaDiff = lambda Poffset: integrate.quad(lambda V: P(V) - Poffset, opt.brentq( offsetFunction, lbound, VminX ), opt.brentq( offsetFunction, VmaxX, rbound))[0]
        # Peq = opt.brentq( calcAreaDiff, max(1e-20, P(VminX)), P(VmaxX) )
        
        # def fnArea(pGuess):
        #     offsetFunction = lambda V: P(V) - pGuess
        #     IntL = opt.brentq( offsetFunction, lbound, VminX )
        #     intR = opt.brentq( offsetFunction, VmaxX, rbound)
        #     return (integrate.quad(offsetFunction, IntL, intR)[0])**2
        
        # step=1e-10
        # #t3 = time.time()
        # pGuess2=helperFunc.findMin(fnArea, [ max(P(VminX),2*step), P(VmaxX) -2*step ], step=step)
        # #t4 = time.time()
        # qPrint("P guess: Bisec, fmin, brentq ", PressureGuess, pGuess2, Peq)
        # #A1 = fnArea(pGuess2)
        # #print("AreaDiff", diffArea, A1)
        # #print("Time:", t2-t1, t4-t3 )
            
        return [1/intersect1, 1/intersect3 , PressureGuess ]
    

class PengRobinsonInteractionPotential(InteractionPotential):
    
    def __init__(self):
        self.G = None
        self.R = None
        self.a = None
        self.b = None
        self.c = None
    
    def setTempRatio(self,tr):
        #print("New TR set")
        self.t = tr * self.tc
        self._calcAlpha(tr)
    
    def _calcAlpha(self, tr):
        _a = 1. + self.c*(1.- np.sqrt(tr))
        self.alpha = _a*_a
    

    def pressure(self, rho, tr=None):
        "Return the pressure of the EOS"
        if not tr==None:
            self.setTempRatio(tr)
            
        part1 = rho * self.R * self.t           /   ( 1- self.b * rho )
        part2 = self.a * self.alpha * rho**2    /   ( 1+ 2*self.b*rho - self.b**2 *rho**2)

        #return (rho * self.R*self.t/(1.-self.b*rho)) - (self.a*self.alpha*rho*rho/(1. + 2.*self.b*rho + self.b*self.b*rho*rho))
    
        return part1 - part2
    
    
    def dPdT(self, rho, tr):
        "return the local pressure derivative wrt temperature at constant density"
        t = tr * self.tc
        
        
        dalpha = +self.c**2/self.tc - (self.c / np.sqrt( t*self.tc) )*(1. + self.c)
    
        # rho* self.R /(1. - self.b*rho) - self.a * rho**2 * dalpha / ( 1. + 2. * self.b*rho - self.b**2 * rho**2)
    
        derivPart1 = rho * self.R             / (1. - self.b * rho) 
        derivPart2 = self.a * rho**2 * dalpha / (1. + 2*self.b*rho - self.b**2 * rho**2)
        return derivPart1 - derivPart2
    
    def dPdTsimple(self, rho, tr):
        "return the local pressure derivative wrt temperature at constant density"
    
        derivPart1 = rho * self.R             / (1. - self.b * rho) 
        return derivPart1
        
    
class LatticePengRobinson(PengRobinsonInteractionPotential):
    
    def __init__(self,G = -1, acentricity = 0.3443, a = 2./49., b = 2./21., tr = 0.9, R=1.0):
        if not "_name" in dir(self): 
            self._name= "Lattice Peng Robinson"        

        self.G = G
        self.R = R
        self.a = a
        self.b = b
        
        self.tc = 0.0778/0.45724*self.a/self.b/self.R
        self.pc = 0.0778 * self.R*self.tc / self.b

        
        #//T pc = 0.0778*_R*tc/_b
        #//T rhoc = pc/0.307/_R/tc
        #// I get 0.30742569
        
        self.acentricity = acentricity
        self.c = (0.37464 + 1.54226*self.acentricity - 0.26992*self.acentricity**2)
        
        self.setTempRatio(tr)
        self.rhoc = self.pc/0.3074/self.R/self.tc
        

    def getEquilibriumDensities(self, tr, accuracy = 1e-14, iterMax = 500, verbose = False):
        return self._superGetEquilibriumDensities(tr, 1e-1, 1e20, accuracy = accuracy, iterMax = iterMax, verbose = verbose)

    def plotPV(self, tr=None):
        self._plotPV( tr, [1e-2,1e1], [-0.1, 0.1] )
        
    
        

class RealPengRobinson(PengRobinsonInteractionPotential):

    def __init__(self, acentricity = 0.344, Tc = 647.4, Pc = 22.10e6, tr = 0.577, R=8314/18):
        
        self.tc = Tc
        self.R = R
        self.acentricity = acentricity
        self.pc = Pc
        self.a = 0.45724 * R**2 * Tc**2 /Pc
        self.b = 0.0778  * R   * Tc    /Pc
        
        self.c = (0.37464 + 1.54226*self.acentricity - 0.26992*self.acentricity*self.acentricity)
        
        self.t = self.tc*tr
        self._calcAlpha(tr)
        self.rhoc = self.pc/0.3074/self.R/self.tc

    
        self._name= "Real Peng Robinson"
        
    def getEquilibriumDensities(self, tr, accuracy = 1e-10, iterMax = 1000, verbose = False):
        return self._superGetEquilibriumDensities(tr, 1.1e-3, 1e20, accuracy = accuracy, iterMax = iterMax, verbose = verbose)

    def plotPV(self, tr=0.9):
        self._plotPV( tr, [1,0.9e3],[-1e7, 1e8] )
        
class RealPengRobinsonAB(LatticePengRobinson):
    
    def __init__(self,G = -1, acentricity = 0.344, a = 599.4, b = 0.01895, R=8314.0/18.0152, tr = 0.9):
        self._name= "Real Peng Robinson AB" 
        super().__init__(G=G,acentricity=acentricity,a=a, b=b, R=R, tr=tr )
        
        
    def getEquilibriumDensities(self, tr, accuracy = 1e-10, iterMax = 1000, verbose = False):
        return self._superGetEquilibriumDensities(tr, 2e-2, 1e20, accuracy = accuracy, iterMax = iterMax, verbose = verbose)

    # def plotPV(self, tr=0.9):
    #     self._plotPV( tr, [1e-2,5e1],[-1e1, 1e5] )
        
class CarnahanStarlingInteractionPotential(InteractionPotential):
    
    def __init__(self,G = -1, a = 1., b = 4., tr = 0.8):
        self._name= "Carnahan Starling"
        self.G = G
        self.a = a
        self.b = b
        
        self.R = 1.0
        
        self.tc = 0.18727/0.4963*self.a/self.b/self.R
        self.t = self.tc * tr
        
        #//T pc = 0.18727*_R*tc/_b
        #//T rhoc = pc/0.35930763/_R/tc
            

    
    def pressure(self, rho, tr=None):
        "Return the pressure of the EOS"
        if tr:
            self.t = tr * self.tc

        c = self.b*rho/4.
        return rho*self.R*self.t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-self.a*rho*rho
    
    def dPdT(self,rho,t):
        c = self.b*rho/4.
        return rho*self.R*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))
    
    def plotPV(self, tr=None):
        self._plotPV( tr, [1e-3,2./3], [-0.01,0.01] )
        
    def getEquilibriumDensities(self, tr, accuracy = 1e-14, iterMax = 500):
        return self._superGetEquilibriumDensities(tr, 1.01, 5e4, accuracy = accuracy, iterMax = iterMax)
    


    
def testConstruction(EOS, name=""):
    trs, rho_v, rho_l = [],[],[]
    for tr in np.linspace(0.49,0.998,50):
        trs.append(tr)
        rhos = EOS.getEquilibriumDensities(tr)[:2]
        rho_v.append(rhos[1])
        rho_l.append(rhos[0])
    
    plt.clf()
    plt.xscale("log")
    plt.plot(np.array(rho_v),trs, label=" Vapour Branch")
    plt.plot(np.array(rho_l),trs, label=" Liquid Branch")

    plt.ylabel(r"T / T$_c$")
    plt.xlabel(r"$\rho$")
    
    plt.legend()
    plt.title("Maxwell construction: "+name+" (Complex)")
    plt.grid(True)
    
def testConstructionPR():
    PR = LatticePengRobinson( acentricity = 0.3443, a=2.0/49., b=2./21. )
    testConstruction(PR, "Peng Robinson")
    plt.savefig("./MaxwellConstructionPRComplexPy")
    
def testConstructionCS():
    CS = CarnahanStarlingInteractionPotential( a=0.5, b=4)
    testConstruction(CS, "Carnahan Starling")
    plt.savefig("./MaxwellConstructionCSComplexPy")

    
def testDeriv(EOS, rho = 0.1, name =""):
    
    plt.clf()
    trs = np.linspace(0.49,0.85,50)
    EOS._plotPT(rho, xLim = [trs[0], trs[-1]])
    
    deriv = np.array( [ EOS.dPdT(rho,tr) for tr in trs] )
    plt.plot(trs,deriv, "-." , label=r"Analytical dP / dT ($\rho$="+str(rho)+")")

    plt.legend()
    plt.title(name+"Pressure derivative wrt Temperature")
    plt.grid(True)
    
def testDerivPR():
    PR = LatticePengRobinson( acentricity = 0.3443, a=2.0/49., b=2./21. )
    testDeriv(PR, name="Peng Robinson: ")
    plt.tight_layout()
    plt.savefig("./DerivativeVerificationPR")
    
def testDerivPRreal():
    PR = RealPengRobinson( )
    testDeriv(PR, 950, name="Peng Robinson Real ")
    print("max rho: ", 1/PR.b)
    plt.tight_layout()
    plt.show()
    # plt.savefig("./DerivativeVerificationPR")
    
def testDerivPRrealAB():
    PR = RealPengRobinsonAB( )
    testDeriv(PR, 0.03, name="Peng Robinson Real ")
    print("max rho: ", 1/PR.b)
    plt.tight_layout()
    plt.show()

def testDerivCS():
    CS = CarnahanStarlingInteractionPotential( a=0.5, b=4)
    testDeriv(CS, name="Carnahan Starling: ")
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.ylim(0.15,0.154)
    plt.tight_layout()
    plt.savefig("./DerivativeVerificationCS")

def enthalphyPRGong(PR, TempRatio):
    #PR = LatticePengRobinson( acentricity = 0.3443, a=0.65/49., b=2./21. , R=1)
    PR.setTempRatio(TempRatio)
    rhol, rhov = PR.getEquilibriumDensities(tempRatio)[:2]

    def enthalpy(rho):
        sqrt2b2 = 2*np.sqrt(2)*PR.b  
        logTop = (2*PR.b**2*rho - 2*PR.b - sqrt2b2)
        logBot = (2*PR.b**2*rho - 2*PR.b + sqrt2b2)
        logCalc = np.log( abs( logTop / logBot ) )
        aTcSqrtAlpha = PR.a*PR.t*PR.c*np.sqrt(PR.alpha)
        h = (  aTcSqrtAlpha /np.sqrt( PR.t * PR.tc ) + PR.a*PR.alpha ) / sqrt2b2 * logCalc
        return h
    
    return enthalpy(rhov) - enthalpy(rhol)
    

if __name__=="__main__":
    
    PR = LatticePengRobinson( acentricity = 0.344, a=2/49., b=2./21. , R=1)
    PRreal = RealPengRobinson( )
    PRrealAB = RealPengRobinsonAB( acentricity = 0.344, a = 599.4, b = 0.01895, R=8314.0/18.0152 )
    # testDerivPRrealAB()
    
    CS = CarnahanStarlingInteractionPotential()
    #print("CS: Tr = 0.7")
    #print("densities CS = ", CS.getEquilibriumDensities(0.8)[:2] )
    tempRatio= 0.7
    print("Temperature = ", tempRatio*647 -273)
    
    # print("densities PR = ", PR.getEquilibriumDensities(tempRatio)[:2] )
    # print("densities PR real = ", PRreal.getEquilibriumDensities(tempRatio)[:2] )
    # print("densities PR real AB = ", PRrealAB.getEquilibriumDensities(tempRatio)[:2] )
    # print("density R PR = ", PR.getEquilibriumDensities(tempRatio)[0]/PR.getEquilibriumDensities(tempRatio)[1] )
    # print("density R PR AB = ", PRreal.getEquilibriumDensities(tempRatio)[0]/PRreal.getEquilibriumDensities(tempRatio)[1] )
    # print("density R PR real AB = ", PRrealAB.getEquilibriumDensities(tempRatio)[0]/PRrealAB.getEquilibriumDensities(tempRatio)[1] )
    #print("densities CS = ", CS.getEquilibriumDensities(tempRatio)[:2] )
    
    RatioCvVL = 1838. / 3800.
    
    cvL = 10
    cvV = cvL * RatioCvVL
    
    # testConstructionCS()
    # testConstructionPR()
    # testDerivPR()
    # testDerivCS()
    
    
    rhol, rhov = PR.getEquilibriumDensities(tempRatio)[:2]
    rholR, rhovR = PRreal.getEquilibriumDensities(tempRatio)[:2]
    rholRab, rhovRab = PRrealAB.getEquilibriumDensities(tempRatio)[:2]
    convRhoR = 958.4/801.34
    
    v_lg = (1/rhov - 1/rhol)
    v_lgR = (1/rhovR - 1/rholR)
    T_eos = tempRatio * PR.tc
    
    Rsp= 8314.0/18.0152
    
    # convRho = 947.13/8.72518
    
    # print("Rho conv:", convRho)
    # print("Rho c ratio:", PRreal.rhoc/PR.rhoc)
    # print("Rho v:", PRreal.rhoc/PR.rhoc*rhov)
    # print("Rho l:", PRreal.rhoc/PR.rhoc*rhol)
    # print("Rho v R:", convRhoR*rhovR)
    
    TempCrit = 647.4
    Tlow = 0.5*TempCrit
    Thigh = 0.7*TempCrit
    dT = Thigh-Tlow
    dTratio = dT/TempCrit
    Nx = 1
    Nt = 7
    convX = 100e-6/Nx
    convT = convX/30/Nt
    convU = convX/convT
    convRho = PRreal.rhoc/PR.rhoc
    convPressure = convRho*convU*convU
    convTemp = dTratio*647
    convCv = convU*convU/(convTemp)
    
    Rdimless = Rsp/(convU**2)*TempCrit
    print("R dimless = ", Rdimless)
    PRfix = LatticePengRobinson( acentricity = 0.344, a=2/49., b=2./21. , R=Rdimless )
    rholRfix, rhovRfix = PRfix.getEquilibriumDensities(tempRatio)[:2]
    print("density PR fix = ", PRfix.getEquilibriumDensities(tempRatio)[:2] )


    L = 2264.705e3
    L_lb = L/(convU**2)
    cv =3800
    cvLB = cv / convCv
    LcvPhys = L/cv
    print("Cv LBM:                      ", cvLB)
    print("L / Cv Phys :                ", LcvPhys)
    print("L / Cv Phys (LB):            ", LcvPhys/convTemp)
    
    latentHeat = T_eos*v_lg*PR.dPdT(rhol,tempRatio)
    # print("L LB:                        ",  latentHeat )
    # print("L / Cv LB:                   ",  latentHeat / cvLB)
    # print("L / Cv LB/Phys Ratio:        ", (latentHeat / cvLB) /(LcvPhys/convTemp) )
    # print("L LB (Phys):                 ",  latentHeat  *convU*convU / 1e3 )
    
    intF2 = lambda rho: rho**-2 * PR.dPdT(rho, tempRatio)
    integral2 = integrate.quad(intF2, rhov, rhol)[0]

    print("Latent heat EoS  scale:      ",tempRatio *integral2 )
    print("Latent heat Phys scale:      ",tempRatio*PR.tc *integral2 * PRreal.pc/PR.pc * PR.rhoc/PRreal.rhoc )
    print("Latent heat lattice scale:   ",tempRatio*PR.tc *integral2 * PRreal.pc/PR.pc /convPressure  )

    
    intFR = lambda rho: PRreal.dPdT(rho, tempRatio)
    intF2R = lambda rho: rho**-2 *PRreal.dPdT(rho, tempRatio)

    # print(intFR(rhovR),intFR(rholR/2) ,intFR(rholR))
    integralR = integrate.quad(intFR, rhovR, rholR)[0]
    integralR2 = integrate.quad(intF2R, rhovR, rholR)[0]

    # print(integrate.quad(intFR, rhovR, rholR))
    # dx = (rholR-rhovR)/10000
    # x = rhovR
    # sumdPdT=0
    # while(x<=rholR):
    #     sumdPdT += dx*intFR(x)
    #     x+=dx
    averageR = integralR / (rholR-rhovR)
    print("avg dPdT real: ",averageR  )
    print("theoretical avg:",L/v_lgR/(tempRatio*647) )
    print("Latent heat using PR: ", tempRatio*647 * integralR2)
    print("L phys            ", L)
    
    integralFunc = lambda rho: rho**-2 *PRreal.dPdT(rho, tempRatio)
    densityLiquidPRreal, densityVapourPRreal = PRreal.getEquilibriumDensities(tempRatio)[:2]
    integralResult = integrate.quad(integralFunc, densityVapourPRreal, densityLiquidPRreal)[0]
     
    print("Latent heat using PR: ", tempRatio*647.4 * integralResult)
    
    
    intF2Rfix = lambda rho: rho**-2 *PRfix.dPdT(rho, tempRatio)
    print("Rho L = ", rholRfix)
    print("Rho V = ", rhovRfix)

    suming=0
    stepwiseResult = integrate.quad(intF2Rfix, 7.5, rholRfix)[0]
    suming+=stepwiseResult
    print("Latent heat using PR: R_l - 7.5 ", tempRatio * stepwiseResult*convU**2 /1e3)
    stepwiseResult = integrate.quad(intF2Rfix, 4, 7.5)[0]
    suming+=stepwiseResult
    print("Latent heat using PR: 7.5 - 4   ", tempRatio * stepwiseResult*convU**2 /1e3)
    stepwiseResult = integrate.quad(intF2Rfix, 1, 4)[0]
    suming+=stepwiseResult
    print("Latent heat using PR: 4 - 1     ", tempRatio * stepwiseResult*convU**2 /1e3)
    stepwiseResult = integrate.quad(intF2Rfix, 0.2, 1)[0]
    suming+=stepwiseResult
    print("Latent heat using PR: 1 - 0.2   ", tempRatio * stepwiseResult*convU**2 /1e3)
    print("Latent heat using PR: sum       ", tempRatio *suming*convU**2 /1e3)
    stepwiseResult = integrate.quad(intF2Rfix, 0.2, rholRfix)[0]
    print("Latent heat using PR: R_l - 0.2 ", tempRatio * stepwiseResult *convU**2 /1e3)
    
    intF2Rfix = lambda rho: rho**-2 *PRfix.dPdT(rho, tempRatio)
    integralRfix2 = integrate.quad(intF2Rfix, rhovRfix, rholRfix)[0]
    
    TLB = 0.5 + (tempRatio*TempCrit - Tlow)/convTemp
    
    print("Nondimensional latent heat using PR:    ", tempRatio * integralRfix2)
    print("Dimensionalized latent heat using PR:   ", tempRatio * integralRfix2 *convU**2 )
    
    print("Gong latent heat PR:             ", enthalphyPRGong( PR, tempRatio)  )
    print("Gong latent heat PR LB:          ", enthalphyPRGong( PR, tempRatio) *PRreal.pc/PR.pc  /(convPressure) )
    print("Gong latent heat PR corrected:   ", enthalphyPRGong( PR, tempRatio) * PRreal.pc/PR.pc *PR.rhoc/PRreal.rhoc )
    print("Gong latent heat PR real crit:   ", enthalphyPRGong( PRreal, tempRatio) )
    print("Gong latent heat PR real ab:     ", enthalphyPRGong( PRrealAB, tempRatio) )
    print("Gong latent heat PR R fix:       ", enthalphyPRGong( PRfix, tempRatio) )
    print("Gong latent heat PR R fix corr:  ", enthalphyPRGong( PRfix, tempRatio) * PRreal.pc/PRfix.pc *PRfix.rhoc/PRreal.rhoc )
    
    print(" Conv Factor: ",  PRfix.tc )
    print(" Conv Factor: ",  PR.tc )
    
    # print("L avg LB:                      ", T_eos/(rhol*rhov)*integral )
    # print("L avg / Cv LB:               ", T_eos/(rhol*rhov)*integral / cvLB)
    # print("L avg / Cv LB/Phys Ratio:    ", T_eos/(rhol*rhov)*integral / cvLB /(LcvPhys/convTemp) )     
    # print("L avg LB (Phys):             ", T_eos/(rhol*rhothatv)*integral *convU*convU / 1e3)
    



###
    # PR = LatticePengRobinson( acentricity = 0.3443, a=3/49., b=2./21. )
    # tempRatio = 0.7
    # rhoL,rhoV = PR.getEquilibriumDensities(tempRatio)[:2] 
    # rhos = np.linspace(rhoV,rhoL,1000)
    # mixture = lambda rho: (rho-rhoV)/(rhoL-rhoV)
    # cv = lambda rho: mixture(rho) * (cvL-cvV) + cvV

    # yfunc = lambda rho, tr: 1- 1/(rho*cv(rho) ) * PR.dPdT(rho,tr)
    
    # plt.plot( rhos, yfunc(rhos, tempRatio) , label="full")
    
    # yfunc = lambda rho, tr: 1- 1/(rho*cv(rho) ) * PR.dPdTsimple(rho,tr)
    
    # plt.plot( rhos, yfunc(rhos, tempRatio) , label="simple")

    # plt.ylim(0)
    # plt.legend()
    # plt.show()
    
    
    #CS._plotPT(rhoDeriv)
    # xs = np.linspace(0.5,1,1000)
    # rhoDeriv = 0.00308
    # plt.plot(xs, (xs+0.5)/(5*rhoDeriv)* CS.dPdT(rhoDeriv,xs), label="Via Deriv CS ({})".format(rhoDeriv))
    # rhoDeriv = 0.0102
    # plt.plot(xs, (xs+0.5)/(5*rhoDeriv)* PR.dPdT(rhoDeriv,xs), label="Via Deriv PR ({})".format(rhoDeriv))
    # plt.legend()
    # plt.show()
    #conversionDensity =   917.0   / PR.getEquilibriumDensities(tempRatio)[0]
    #print("Crit Ratio:", conversionDensity)
    #print("densities PR L = ", conversionDensity * PR.getEquilibriumDensities(tempRatio)[0] )
    #print("densities PR G = ", conversionDensity * PR.getEquilibriumDensities(tempRatio)[1] )
    #print(CS.getEquilibriumDensities(0.5))
    #print(CS.getEquilibriumDensities(0.4))
    #PR.plotPV(0.8)
    # CS.plotPV(0.49)
    # print(CS.getEquilibriumDensities(0.55)[:2] )
    
    