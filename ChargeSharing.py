# Correct a 80x80 frame of data  - Oakley Clark 

import numpy as np




import time



def CSD(frame,threshold):
    
    corrected = np.zeros((80,80))
    
    xpix = 80
    ypix = 80
    
    threshold = 1000*threshold/4.6
    
    # set up new slightly bigger frame to avoid problems at edges
    f = np.zeros((xpix+4,ypix+4))
    f[2:xpix+2,2:ypix+2] = frame
    
    # get event list
    [x,y] = np.where(f>threshold)
    
    # iterate throgh events
    for i in range(len(x)):
        
        hitcount = 0
        
        xhit = x[i]; yhit = y[i]
        
        # find neighbours
        neighbours = f[xhit-1:xhit+2,yhit-1:yhit+2]

        # count how many events in neighbours above threshold
        for a in range(3):
            for b in range(3):
                if (neighbours[a,b] > threshold):
                    hitcount += 1
        
        # for CSD, only keep event if the hitcount = 1
        if (hitcount == 1):
            corrected[xhit-2,yhit-2] = f[xhit,yhit] 
    
    return corrected


def CSA(frame,threshold):
    
    corrected = np.zeros((80,80))
    
    xpix = 80
    ypix = 80
    
    threshold = 1000*threshold/4.6
    
    # set up new slightly bigger frame to avoid problems at edges
    f = np.zeros((xpix+4,ypix+4))
    f[2:xpix+2,2:ypix+2] = frame
    
    # get event list
    [x,y] = np.where(f>threshold)
    
    # iterate throgh events
    for i in range(len(x)):
        
        
        neighboursum = 0
        
        xhit = x[i]; yhit = y[i]
        
        # find neighbours
        neighbours = f[xhit-1:xhit+2,yhit-1:yhit+2]
        
        # Loop through neighbours - sum counts above threshold - find max value pixel
        [maxx,maxy] = np.where(neighbours == np.max(neighbours))

        try:
            maxx = int(maxx); maxy = int(maxy)
        except:
            maxx = int(maxx[0]); maxy = int(maxy[0])
            print('Equal maximum signal pixel, attributing to first')
            
        
        for i in range(3):
            for j in range(3):
                if (neighbours[i,j] > threshold):
                    neighboursum += neighbours[i,j]
        
        # put the summed event into the maximum pixel in the corrected frame
        corrected[xhit-2-1+maxx,yhit-2-1+maxy] = neighboursum

    return corrected        
        
        
        
        
        
        