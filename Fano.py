import numpy as np

# Takes as input an energy deposition in keV
# returns a fano adjusted energy in keV
def fano(edep):
    
    fanoFactor = 0.09
    creationEnergy = 4.6
    
    N = 1000*edep/creationEnergy
    
    sigma = np.sqrt(fanoFactor*N)
    
    fanoN = np.random.normal(N,sigma)
    
    fanoEnergy = fanoN*creationEnergy/1000.0
    
    return fanoEnergy