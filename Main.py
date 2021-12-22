# This is the main file
# This file uses the NEWFORMAT CSV file inside DepositionData
# Make sure this file is up to date for each run
#
# It runs through all of the events and does the following:
# 1) Propagates each deposition of the event
# 2) Reads the total event out to a frame
# 3) Performs CSA and CSD charge sharing on each frame
# 4) Bins each event using CSA and CSD correction to make 2 hxt style arrays
# 5) Produces summed spectra for each charge sharing algorithm



import numpy as np
from csv import reader
import matplotlib.pyplot as plt
from progress.bar import IncrementalBar
import os
import time


# import my own functions
from Propagation import GetDriftTime, GetSigmaZero, f, GetSigmaD, GetSigmaF
from PixelReadout import DepositionReadout, CCE
from ChargeSharing import CSD, CSA
from add_frame import add_frame
from Fano import fano



# User options
threshold = 1.5 # keV (for charge sharing)
Emin = 0.0 # keV
Emax = 120.0 # keV
nChannels = 600
EminElec = (Emin*1000.0)/4.6
EmaxElec = (Emax*1000.0)/4.6


# Initialize:
sumFrame = np.zeros((80,80))
oldEventID = -1
rowNumber = 0


# Histogram Initialization
CSDhist = np.zeros((80,80,nChannels))
CSAhist = np.zeros((80,80,nChannels))
summedCSD = np.zeros(nChannels)
summedCSA = np.zeros(nChannels)
binEdges = np.linspace(EminElec,EmaxElec,nChannels+1)
energyEdges = np.linspace(Emin,Emax,nChannels+1)
energyCentres = np.zeros(nChannels)
for i in range(nChannels):
    energyCentres[i] = 0.5*(energyEdges[i]+energyEdges[i+1])


# Input file
inputFile = 'Data.csv'
# Output files
outputPath = 'Output/'



# Progress Bar info and formatting
print('Getting File Length...')
numRows = 0
#for row in open("DepositionData/NEWFORMAT_Analysis_nt_Hits_t0.csv"):
for row in open(inputPath+inputFile):
    numRows += 1
bar = IncrementalBar('Processing Rows...', max=numRows)


count = 0
t0 = time.time()

#with open('DepositionData/NEWFORMAT_Analysis_nt_Hits_t0.csv', 'r') as read_obj:
with open(inputPath+inputFile, 'r') as read_obj:
    # pass the file object to reader() to get the reader object
    csv_reader = reader(read_obj)
    # Iterate over each row in the csv using reader object
    #for row in csv_reader:
    for row in csv_reader:
        
        # Update progress
        bar.next()
        if (count%1000 == 0):
            percent = round(100*(count/numRows),1)
            timePassed = time.time()-t0
            try:
                ETA = (timePassed*100/percent)-timePassed
                ETAmin = ETA//60
                ETAsec = round(ETA%60,1)
            except:
                ETAmin = ''
                ETAsec = ''
            print(count,' / ',numRows,' ... ',percent,'%',' ETA: ',ETAmin,' mins: ',ETAsec,' seconds:')
        count += 1

        # If it is a new event (but not first event):
        if (oldEventID != int(row[0]) and oldEventID!=-1):
            
            # Do charge sharing on complete sum frame of last event
            CSDframe = CSD(sumFrame,threshold)
            CSAframe = CSA(sumFrame,threshold)
            # Get hxt array for a single frame
            hxtCSDframe = add_frame(CSDframe,binEdges)
            hxtCSAframe = add_frame(CSAframe,binEdges)
            
            
            # Add hxt frame arrays to master hxt arrays
            [i,j,k] = np.where(hxtCSDframe != 0)
            for g in range(len(i)):
                i0 = i[g]; j0 = j[g] ; k0 = k[g]
                CSDhist[i0,j0,k0] += hxtCSDframe[i0,j0,k0]
            [i,j,k] = np.where(hxtCSAframe != 0)
            for g in range(len(i)):
                i0 = i[g]; j0 = j[g] ; k0 = k[g]
                CSAhist[i0,j0,k0] += hxtCSAframe[i0,j0,k0]
                
                
                
            # Reinitialize the sum frame to be equal to first line of new event
            x = float(row[1]); y = float(row[2]); z = float(row[3]); E = float(row[4])
            # Apply Hecht equation and fano factor to find E
            #E *= CCE(z)
            E = fano(E)
            driftTime = GetDriftTime(z)
            sigmaZero = GetSigmaZero(E)
            sigmaFinal = GetSigmaF(E,sigmaZero,driftTime)
            frame = DepositionReadout(sigmaFinal, x, y, z, E)
            sumFrame = frame
            #print('frame made')
            oldEventID = int(row[0])
        
        # If it is another line from the same event (or first line of first event)
        else:
            # Update oldEventId
            oldEventID = int(row[0])
            # get parameters
            x = float(row[1]); y = float(row[2]); z = float(row[3]); E = float(row[4])
            # Apply Hecht equation and fano factor to find E
            #E *= CCE(z)
            E = fano(E)
            # Propagate charges
            driftTime = GetDriftTime(z)
            sigmaZero = GetSigmaZero(E)
            sigmaFinal = GetSigmaF(E,sigmaZero,driftTime)
            # Get the frame of that line
            frame = DepositionReadout(sigmaFinal, x, y, z, E)
            #print('frame made')
            # Add to the sum frame
            sumFrame += frame                                                                                 


bar.finish()

# Make summed histograms
for i in range(80):
    for j in range(80):
        summedCSD[:] += CSDhist[i,j,:]
        summedCSA[:] += CSAhist[i,j,:]
        
        
# Make plots
plt.plot(energyCentres,summedCSD)
plt.title('Summed CSD')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.show()
plt.plot(energyCentres,summedCSA)
plt.title('Summed CSA')
plt.xlabel('Energy (keV)')
plt.ylabel('Counts')
plt.show()

 

#  Save the histograms
np.save(outputPath + '/CSD/' + 'DataCSD',summedCSD)
np.save(outputPath + '/CSA/' + 'DataCSA',summedCSA)







