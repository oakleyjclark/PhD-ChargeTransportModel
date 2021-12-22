# Pixel Readout Functions

import numpy as np
import math
from scipy.interpolate import interp1d

from ConvertCoordinates import GetPixelCoordinates




WP = np.load('/Users/oakleyjclark/Desktop/PhD-Work/ElectricFieldModel-V2/WP1D.npy')
z0 = np.linspace(0,0.2,len(WP))
getWP = interp1d(z0,WP)

# Parameters
L = 0.2
E = 1000/L
ue = 940
uh = 114
te = 1.2*10**-6
th = 2.5*10**-6
lamh = uh*th*E
lame = ue*te*E

# Function to get the CCE - accepts z in m 
def CCE(z):
    z *= 100 # convert z to cm
    holes = 1 - np.exp(-(getWP(z))*(L/lamh))
    electrons = 1 - np.exp(-(L-getWP(z)*L)/lame)
    holes *= (lamh/L)
    electrons *= (lame/L)
    return electrons + holes





# This function adds noise to each pixel in a frame that contains charge
# Inputs: Uncorrected frame, ENC (equivalent noise charge (in electrons))
# Outputs: Uncorrected frame with noise
def AddNoise(frame, ENC):
    
    # Get the non zero elements in the frame
    elements = np.nonzero(frame)
    xpix = elements[0]; ypix = elements[1]
    
    # loop over pixels and add noise
    for i in range(len(xpix)):    
        noise = float(np.random.normal(0.0,np.sqrt(ENC),1))
        frame[xpix[i],ypix[i]] += noise
    
    return frame
    

# Function that takes the propagated electrons and converts it to into
# frame format
# Inputs: SigmaFinal (m), x(m), y(m), Edep (keV), scalar
# Outputs: 80x80 Frame with integer electron hits
def DepositionReadout(SigmaFinal,x,y,z,Edep,scalar=5):
    
    N = int(np.round(((Edep*1000)/4.6)/scalar))
    
    xsample = np.random.normal(x,SigmaFinal,N)
    ysample = np.random.normal(y,SigmaFinal,N)
    
    # Initialize frame
    frame = np.zeros((80,80))
    
    # Loop over each carrier
    for i in range(len(xsample)):
        # Find the x and y pixel coords
        xpix, ypix = GetPixelCoordinates(xsample[i],ysample[i])
        
        # Add electron to frame
        try:
            frame[xpix,ypix] += scalar
        except:
            # Charge is lost outside of frame - just ignore these electrons
            pass
    
    # Add noise to the frame (ENC in electrons)
    frame = AddNoise(frame,1000)
    
    return frame






