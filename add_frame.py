# This file turns a frame into a hxt style array for a single frame

import numpy as np
import bisect

def add_frame(frame,binEdges):
    
    nChannels = len(binEdges) - 1
    
    arrayOut = np.zeros((80,80,nChannels))
    
    threshold = 50
    
    # find event list
    [x,y] = np.where(frame>threshold)
   
    
    # loop through event list
    for i in range(len(x)):
        
        xhit = x[i]; yhit = y[i]
        
        # for k in range(nChannels):
            
        #     if (frame[xhit,yhit]>binEdges[k] and frame[xhit,yhit]<binEdges[k+1]):
                
        #         arrayOut[xhit,yhit,k] += 1
                
        #         print('my method index: ',k)
        
        insertIndex = bisect.bisect_left(binEdges,frame[xhit,yhit])-1
        arrayOut[xhit,yhit,insertIndex] += 1
        
        
        
        
    return arrayOut