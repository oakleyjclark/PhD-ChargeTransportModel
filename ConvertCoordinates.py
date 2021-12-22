# This function converts from global to local coordinates. System described
# in ReadMe

def GetPixelCoordinates(x,y):
    
    pitch = 250*10**-6 # m
    xpix = int(x//pitch)
    ypix = int(y//pitch)
    
    return xpix, ypix