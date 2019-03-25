import numpy as np
from matplotlib import pyplot as plt
import os
from astropy.io import fits
import sp

path = "data/"
filenamess = os.listdir(path+"light")
filenames = filenamess
# filenamess = np.sort(filenamess)
# filenames = []
# for i in filenamess:#het moest een lijst zijn en ik had geen internet dus...
#     filenames.append(i)
print(filenames)
fluxes = np.ndarray((1, 6))#time, flux, error, fluxref, errorref, testvalue

windowsize = [411, 1431]

def Bias():
    bias = np.ndarray((windowsize[1], windowsize[0]))
    for filename in os.listdir(path+"bias"):
        bias += fits.getdata(path+"bias/"+filename)#sommeren
    bias /= len(os.listdir(path+"bias"))#gemiddelde nemen
    return bias


def Flat(bias):
    flat = np.ndarray((windowsize[1], windowsize[0])) #lege array maken ter grootte van de flat
    for filename in os.listdir(path+"flats"):
        #print(fits.getdata(path+"flats/"+filename)[940:2371, 850:1261].shape)
        flat += (fits.getdata(path+"flats/"+filename)[940:2371, 850:1261]-bias) / fits.open(path+"flats/"+filename)[0].header["exptime"]
        #bij de lege array de gecorrigeerde flat optellen, (flat - bias)/exptime
        #print(np.max(flat))
    flat /= len(os.listdir(path+"flats"))#de som delen door het aantal flats voor een gemiddelde
    flat /= np.mean(flat)#normaliseren
    return flat


def Find(data):
    found = 0
    hmin = 100
    """het idee was om hmin, de minimale threshold, telkens te vergrooten als hij niks vond(ja hij zou moeten verlagen I know) 
    echter vind hij voor elke image nu de goede positie dus ik laat het lekker zo staan. Tevens had ik ook nog voor de 
    zekerheid alleen de twee waardes met de grootste flux over gehouden (de arg1 en arg2), 
    maar hij vind er telkens exact twee dus dat was ook niet nodig.
    """
    while found == 0:
        x, y, fluxtot, sharpness, roundness = sp.find(data, hmin, 5.)
        found = len(x)
        # hmin += 100
        # if hmin >= 300:
        #     print("poging "+str((hmin-100)/100))
    #print(hmin)
    #print(x)
    arg1 = np.argmax(fluxtot)
    arg2 = np.argmax(np.delete(fluxtot, np.argmax(fluxtot)))
    return np.array([x[arg1], x[arg2]]), np.array([y[arg1], y[arg2]]), np.array([fluxtot[arg1], \
    fluxtot[arg2]]), np.array([sharpness[arg1], sharpness[arg2]]), np.array([roundness[arg1], roundness[arg2]])


def Time(timestr):#tijd (UT) van waarnemen is een string van vorm hh:mm:ss.s, dus die wordt hier omgezet in een getal(in dagen)
    time = 0
    time += float(timestr[1:3])/24#uren
    time += float(timestr[4:6])/24/60#minuten
    time += float(timestr[7:11])/24/60/60#seconden
    return time #in dagen

print(Time(" 23:25:00.0"))

bias = Bias()
flat = Flat(bias)
print(bias)

testfile = fits.open(path+"light/"+filenames[11])
testdata = testfile[1].data
print(testfile[0].header["exptime"])
print(np.max(testdata))
# plt.imshow(testdata)
# plt.show()

#asdfhsgioidnfseodj



for filename in filenames:

    filetje = fits.open(path+"light/"+filename, uint=False)#uint moet false omdat anders 10 - 25 = 2**64 - 15 ipv gewoon -15

    exptime = filetje[0].header["exptime"]
    time = Time(filetje[0].header["utobs"]) #tijd omgerekend naar dagen
    #print(time)

    data = filetje[1].data
    data -= bias
    #print(data)
    data /= exptime
    data /= flat



    x, y, fluxtot, sharpness, roundness = Find(data)
    #print((np.sum(bias)/(len(bias[:,0])*len(bias[0,:]))))
    g = 8 #electrons per count?
    flux, fluxerr, sky, skyerr = sp.aper(data, xc=x, yc=y, phpadu=g, apr=[12], skyrad=[12, 15], flux=True, silent=True, setskyval = 0)
    #, setskyval=np.sum(bias)/(len(bias[:,0])*len(bias[0,:])))
    #apr = Vector of up to 12 REAL photometry aperture radii.
    #skyad = skyrad -- Two element vector giving the inner and outer radiiaper to be used for the sky annulus. Ignored if the SETSKYVA keyword is set.


    

    flux = np.ndarray.flatten(flux)#/flux[1]
    fluxerr = np.ndarray.flatten(fluxerr)
    skyerr = np.ndarray.flatten(skyerr)

    print(exptime, flux[0], fluxerr[0], flux[1], fluxerr[1])
    fluxes = np.append(fluxes, np.array([[time, flux[0], fluxerr[0], flux[1], fluxerr[1], int(filename[1:8])]]), axis=0)
    #print(fluxes)
    print(filenames.index(filename))
    if filenames.index(filename) == 100:
        break


fluxes = fluxes[1:, :]
#fluxes[:, 1] /= np.mean(fluxes[:, 1])

plt.scatter(fluxes[:, 0], fluxes[:, 1])
#plt.errorbar(fluxes[:, 0], fluxes[:, 1], fluxes[:, 2])
plt.xlabel("time (days)")
plt.ylabel("flux")
#for i in range(len(fluxes[:, 0])):
    #plt.text(fluxes[i, 0], fluxes[i, 3], fluxes[i, 5])
#plt.ylim(0.9, 1.1)
plt.xlim(0.87, 0.98) #plt.xlim(0.873523148148, 0.975694444444)
plt.show()

