import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

#These are the velocity of the local standard of rest (i.e. velocity of the Sun around the Milky Way) that we need to subtract from the velocities associated with the spectra.
#We measured these when taking our observations
VLSRArray = [-3.687587719, -4.70862069, -2.669926471, 3.273969575, 4.800390126, 1.681310196, 4.323207041, 0, 0.959573699, 2.924186282, 1.172483009, 0.528421053, 1.796052632,
             -1.223047337, 0.313157895, 0.319736842, -0.103947368, -0.572368421, -1.128157895, -1.415789474, 1.333105951, 3.442228739, 1.054107402, 1.129867061, -3.686923077, 2.403936842,
             1.76765896, -1.258421053, 2.418489209, 2.65680758, 5, 5.736945668, -2.067368421, -3.425, -1.364210526, -4.770789474, -1.968421053, 21.09245536, 4.967955224, -1.597631579, 4.766523031,
             5.800178042, -2.419736842, -1.242894737, 3.752485119, 3.490798817, 1.920029412, 2.842088235, 0, 3.416432749]
    
#This function will take the longitude and distance along it to calculate the expected Doppler shifted frequency based on the tangent point method
#The output frequency is a result of geometric properties of the Milky Way and the fact that the the rotational velocity of the Milky Way
#is constant in radius once outside the central bulge (this is in fact part of the evidence for dark matter in the galaxy!)
def TPM(longitude, distance, vlsr):
    R = np.sqrt(distance ** 2 + 8178 ** 2 - 2 * distance * 8178 * np.cos(np.radians(longitude))) #8178pc is roughly the radius of the Sun to the centre of the Milky Way
    Vp = 237600 * np.sin(np.radians(longitude)) * ((8178 / R) - 1) - vlsr * 1000  #237600m/s is the rotational velocity of the Milky Way outside ~0.2kpc
    CVp = (1 - Vp / 299792458) #Doppler shift factor
    ansp = (CVp) * 1420.405751 #Doppler shifted 21-cm emission frequency
    return ansp

#Initialising the arrays we will be putting distances, longitudes, velocities and frequencies in
distancesArray = np.arange(10, 32460, 10) #Goes to ~32kpc because this is roughly 
galacticLongitudesArray = np.arange(0, 250, 5) #We took data for these galactic longitudes
velocitiesArray = np.empty((len(galacticLongitudesArray), len(distancesArray)))
frequenciesArray = np.empty((len(galacticLongitudesArray), len(distancesArray)))


for i in range(len(galacticLongitudesArray)):
    for j in range(len(distancesArray)):
        frequenciesArray[i][j] = TPM(galacticLongitudesArray[i], distancesArray[j], VLSRArray[i])

#This is the spectral data that we are loading for the different galactic longitudes - this data has been cleaned via averaging the frequency bins
#over the different spectra taken (each longitude individually, not averaging over all longitudes too)
frequencies = np.loadtxt('Frequencies.txt', delimiter='\t')
DopplerShiftVelocityRange = []

#Doppler shift calculator based on input Doppler shifted 21-cm emission line
def DopplerShift(frequency):
    return 299792458 * (1420.406 - frequency) / frequency

for i in frequencies:
    DopplerShiftVelocityRange.append(DopplerShift(i) / 1000)
    


#Initialse array for finding velocities associated with peaks
allMaxima = []
allDistances = []
spectralData = pd.read_csv('G0to1000.csv', delimiter=',')
frequencies = np.loadtxt('Frequencies.txt', delimiter='\t')
#Here we are trying to capture the 21-cm emission peaks in the spectral data to obtain the points in the Milky Way at which
#hydrogen is present. This creates the image titled 'PeakMap.png'
for i in range(0, 50):
    row = spectralData.iloc[i]
    maximum, _ = find_peaks(row, height=10, prominence=5) #Useful to adjust for which peaks to pick up
    if len(maximum) == 1:    
        allDistances.append(frequencies[maximum[0]])
    if len(maximum) == 2:
        allDistances.append([frequencies[maximum[0]], frequencies[maximum[1]]])
    if len(maximum) == 3:
        allDistances.append([frequencies[maximum[0]], frequencies[maximum[1]], frequencies[maximum[2]]])
                             

#These arrays are intialized to be appended with the radii and theta that a longitude spans
rDist = []
actualDistances = []
longitudesArr = []


#Thse for loops go through all the maxima that were found to calculate their radii from the centre of the Milky Way
for i in range(0, 50):
    distances = []
    if type(allDistances[i]) != list :
            maxFreq = allDistances[i]
            difArr = np.absolute(frequenciesArray[i] - maxFreq)
            index = difArr.argmin()        
            distance = distancesArray[index]
            actualDistances.append(distance)
            rDist.append(distance)
            longitudesArr.append(np.radians(i * 5))
            
    elif type(allDistances[i]) == list:
        for j in range(len(allDistances[i])): 
            maxFreq = allDistances[i][j]
            difArr = np.absolute(frequenciesArray[i] - maxFreq)
            index = difArr.argmin()        
            distance = distancesArray[index]
            rDist.append(distance)
            longitudesArr.append(np.radians(i * 5))

        

fig = plt.figure()
ax = fig.add_subplot(projection = 'polar')
ax.set_theta_zero_location("N")
ax.scatter(longitudesArr, rDist, marker='.')
plt.show()

frequenciesfromData = []

#This loop is used to add the distances associated (elementwise) with each frequency observed along each longitude
for i in range(0, 50):
    row = spectralData.iloc[i]
    maxIndex = int(frequenciesArray[i].argmax())
    minIndex = int(frequenciesArray[i].argmin())
    usedIndex = minIndex
    if i == 0 or i == 18 or i >= 37:
        usedIndex = 1
    frequenciesArrayRow1 = frequenciesArray[i, 0:usedIndex]
    frequenciesArrayRow2 = frequenciesArray[i, usedIndex:3244]
    if i >= 18 or i == 0:
        frequenciesArrayRow1 = frequenciesArrayRow2
    distancesforData = []
    for j in range(0, 216):
        freq = frequencies[j]
        dif = np.absolute(frequenciesArrayRow2 - freq)
        index = dif.argmin()
        distance = distancesArray[index+usedIndex]
        distancesforData.append(distance)
    frequenciesfromData.append(distancesforData)




#This code is to produce the final map of the Milky Way, found with code in 'MilkyWayMap.png'
rMesh, thetaMesh = np.meshgrid(np.linspace(0, 216, 216), np.radians(galacticLongitudesArray))
fig, ax = plt.subplots(dpi = 300, subplot_kw=dict(projection='polar'))
ax.set_theta_zero_location("N")
ax.contourf(thetaMesh, frequenciesfromData, spectralData, 100, cmap='plasma', norm='linear')
ax.set_rmax(32460)
plt.show()