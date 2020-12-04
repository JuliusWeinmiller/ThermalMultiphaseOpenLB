#  Python implementation of the Maxwell Construction
#  used in the OpenLB extension
#    
#  Written by Julius Weinmiller for his thesis
#  https://www.openlb.net/forum/users/julius-weinmiller/
#

# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:19:31 2020

@author: Julius
"""

def findMin( f, bounds, accuracy = 1e-14, step = 1e-4):
    
    Xa = bounds[0]
    Xb = bounds[1]
    dPab = 1
    
    # ensure that Xa is left and Xb is right
    if Xa < Xb: Xa, Xb = Xb, Xa
    
    while( (Xa - Xb ) > accuracy and dPab !=0 ):
        Xab = (Xa+Xb)/2.
        
        # Get derivative magnitude
        # magnitude doesn't matter, only sign does, no need to divide by dx
        #dPab = f(Xab + step) - f(Xab - step)         
        dPab2 = -f(Xab + 2*step) + 8*f(Xab + step) - 8*f(Xab - step) + f(Xab - 2*step) #// magnitude doesn't matter, only sign does, no need to divide by dx
        
        #print(Xa - Xb )
        #print(dPab)

        #// if bigger than 0 -> positive -> right side
        if( dPab2<0  ): Xb = Xab
        #// if smaller than 0 -> negative -> left side
        else:  Xa = Xab
     #// end while
    
    return (Xa+Xb)/2. 



def findRoot( f, bounds, accuracy = 1e-14 ):
    
    Xa = bounds[0]
    Xb = bounds[1]
    

    assert f(Xa)*f(Xb) <= 0, "f(a) and f(b) must have different signs"
    
    # ensure that f(Xa) is positive and f(Xb) is negative
    if f(Xa) < f(Xb): Xa, Xb = Xb, Xa
    
    Xab = (Xa+Xb)/2.
    
    while( (Xa - Xb ) > accuracy and f(Xab)!=0 ):
        # if f(Xab)<0, replace Xb since Xb is negative
        if ( f(Xab)<0 ): Xb = Xab 
        else:            Xa = Xab
        Xab = (Xa+Xb)/2.
            
    return Xab



if __name__=="__main__":
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    f = lambda x: (x+1)**3 + 4*x**2 - 4
    df = lambda x: 3*x**2 + 14*x + 3
    bounds = [-1,3]
    
    x = np.linspace(*bounds)
    plt.plot(x, f(x))
    plt.grid(True)
    
    analytical = 2*np.sqrt(10)/3 - 7/3
    root = findRoot(df,bounds)
    print( "Via Root: ", analytical-root )
    print( "Minimum: 3", analytical-findMin(f, bounds, step=1e-3) )
    print( "Minimum: 4", analytical-findMin(f, bounds, step=1e-4) )
    print( "Minimum: 5", analytical-findMin(f, bounds, step=1e-5) )
    print( "Minimum: 6", analytical-findMin(f, bounds, step=1e-6) )
    print( "Minimum: 7", analytical-findMin(f, bounds, step=1e-7) )
    print( "Minimum: 8", analytical-findMin(f, bounds, step=1e-8) )
    print( "Minimum: 9", analytical-findMin(f, bounds, step=1e-9) )
    print( "Minimum: 10", analytical-findMin(f, bounds, step=1e-10) )
    
