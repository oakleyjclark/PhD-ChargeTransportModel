# Propagation functions


# function list: 

# GetDriftTime - Finds the drift time for the electrons
# GetSigmaZero - Finds original size of the electron charge cloud



import numpy as np
from scipy.integrate import solve_ivp



# Propagation parameters
mobility = 1380 / (100**2) # m^2 / (V.s)


# Function finds the drift time of the electron cloud (assuming it were to
# fully travel to anode).
# Inputs: local Z coordinate (m), abs(bias voltage) (V)
# Outputs: Time (s) of drift
def GetDriftTime(z,bias=1000):
    
    thickness = 0.002 # m    
    distance = thickness - z
    time = (thickness*distance)/(mobility*bias)
    
    return time
    

# Function finds the original sigma of the electron cloud
# Inputs: Energy of deposition (keV)
# Outputs: Sigma (m)
def GetSigmaZero(Edep):
    
    # Polynomial fitted paramaeters based off paper
    # Koch-Merin et al 'Soectroscopic monte-carlo model to sim...'
    a = 1.7582417*10**-10
    b = 2.24175824*10**-8
    c = -8.54148011*10**-22
    
    sigmaFinal = a*Edep**2 + b*Edep + c
    
    return sigmaFinal
    

# Function that is called in GetSigmaDto find the final sigma
def f(t,y,edep,sigmaI):
    
    N = 1000*edep/4.6
    kb = 1.38064852*10**-23 # SI units
    T = 298 # K
    epsilon = 10.9 * 8.854187*10**-12 # SI units
    q = 1.6*10**-19 # C
    
    k1 = 2*mobility*kb*T/q
    k2 = (mobility*N*q)/(12*np.pi**(3/2)*epsilon)
    
    return k1+(k2 / (np.sqrt((sigmaI**2)+y)) )
    

def GetSigmaD(edep,sigmaI,DriftTime):
    
    sol = solve_ivp(f,[0,DriftTime],[0], args=(edep,sigmaI))
    sigmaD = np.sqrt(sol.y[0,-1])
    
    return sigmaD


def GetSigmaF(edep,sigmaI,DriftTime):
    
    sigmaD = GetSigmaD(edep,sigmaI,DriftTime)
    sigmaF = np.sqrt(sigmaD**2 + sigmaI**2)
    
    return sigmaF


