# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:37:51 2020

@author: Ranjan
"""
"""
Fourier Transform of 1009 set of readings
"""

import xlrd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft #, fftfreq
from scipy.signal import savgol_filter
import xlsxwriter
from astropy.io import ascii



#%%

# <<<READING THE DATA>>>


loc= "/Users/pragya/Desktop/Gale_Crater_Rhythmites/1009_thickness.xlsx"

wb=xlrd.open_workbook(loc)
sheet=wb.sheet_by_name('Sheet1')


data=np.zeros((127,1010))

for i in range(sheet.nrows): 
    for j in range(sheet.ncols):
           data[i,j] = ((sheet.cell_value(i, j)))
           


# raw_data=np.delete(raw_data, (22,23), axis=0)
# raw_data=np.delete(raw_data, np.s_[1000:], axis=1)
# data=np.copy(raw_data)
 
#%%           
# data=ascii.read("/Users/pragya/Desktop/Gale_Crater_Rhythmites/1009_individual_periodicities.csv")
# print(data[63][200])

           
#%%          

# <<Converting the number of pixel to length in cm(?)>>
          
#plt.plot(data[:,3])

s1=0.00695 #0-21   [rows 22 & 23 = 0] 
d1=47.6
s2=0.00603 #24-36 
d2=50.19
s3=0.00789 #37-51 
d3=16.259 
s4=0.00915 #52-80 
d4=56.31
s5=0.0089  #81-95 
d5=9.462
s5a=0.0273 #96-108 
s6=0.0184  #109-126  
d6=28.543


lm=np.zeros((127,1010))
lm[0:24,:]=np.sin(d1*np.pi/180)*data[0:24,:]*s1
lm[24:37,:]=np.sin(d2*np.pi/180)*data[24:37,:]*s2
lm[37:52,:]=np.sin(d3*np.pi/180)*data[37:52,:]*s3
lm[52:81,:]=np.sin(d4*np.pi/180)*data[52:81,:]*s4
lm[81:96,:]=np.sin(d5*np.pi/180)*data[81:96,:]*s5
lm[96:109,:]=np.sin(d5*np.pi/180)*data[96:109,:]*s5a
lm[109:127,:]=np.sin(d6*np.pi/180)*data[109:127,:]*s6


lm=np.delete(lm, (22,23), axis=0)
lm=np.delete(lm, np.s_[1000:], axis=1)

#lm=np.round(lm,3)
lm_avg=np.mean(lm, axis=1)


#%%

"""
 <<<<<<< RUN THIS BLOCK >>>>>>>
 The 0th index Fourier coefficient is the DC which is simply the sum of all the values \<nextLIne>\
 in the input signal. i.e np.sqrt(fm[0,i])=lm[*,i].sum().\<nextLIne>\
     
     BUT REMEMBER, THERE ARE SEVERAL CONVENTIONS FOR DEFINING THE DFT AND IDFT, SOME OF THESE
     USE THE 1/N NORMALIZATION WHERE DC BECOMES THE MEAN OF THE INPUT SIGNALS> SCIPY FFT DOES NOT 
     USE THE 1/N NORMALIZATION, SO DC IS EQUAL TO THE SUM OF ALL INPUTS AND NOT THEIR AVERAGE.
     
     
 The square root is because the DC component is also being squared in this code.
"""

freq=np.linspace(0, 3.14, 125) #Frequency values for the FFT
i=0
f=np.zeros(125) #Temporary variable to store the real part of the FFT, i.e., the amplitude. The Fourier coefficients are actually not being saved, plus the DC component is being squared.
fm=np.zeros((125,1000)) #the matrix to store all the 'f' variable values
for i in range (999):
    
    F =fft(lm[:,i])
    f= (np.sqrt(np.power(F.real[0:125],2)+np.power(F.imag[0:125],2))/len(freq))*2 # Calculation of amplitude and scaling it (see comments at the bottom of this block), but remember you are squaring the DC component too.
    fm[:,i]=f[:]
    plt.figure(1)
    plt.plot(freq[1:125], f[1:125]**2) #Power spectra (obtained by squaring the amplitudes) with the DC value removed

"""
 <<< CALCULATION OF AMPLITUDE AND CORRECT SCALING OF THE AMPLITUDES >>>

-> The amplitudes of each Fourier coefficients = abs(F). 
-> The amplitudes have to be scaled as such: abs(F) / (total no of samples in the input signal).
-> But the amplitudes are split equally between the positive and negative frequencies, so each value has to be multiplied by 2.
-> Hence, the correct scaled amplitudes are given by: Scaled_Amp= abs(F)/len(freq)*2
CHECK FREQUENCY FOLDING DUE TO NYQUIST FREQUENCY. WHEELSPOKE EFFECT.

{Note that: np.power(F.real[0:125],2)+np.power(F.imag[0:125],2) is the same as np.square(abs). 
That is: np.power(F.real[0:125],2)+np.power(F.imag[0:125],2) is the power, unscaled.}


Power is equal to Amplitude**2.  Here we have to use the Scaled_amp**2.

However remember that the DC is euqual to the sum of all the values in the input signal. So do not square it.

"""
#%%

"""    
  ORIGINAL PIECE OF CODE <<< DO NOT CHANGE >>>
"""

# freq=np.linspace(0, 3.14, 64) #Frequency values for the FFT
# i=0
# f=np.zeros(64) #Temporary variable to store the real part of the FFT, i.e., the amplitude
# fm=np.zeros((64,1009)) #the matrix to save all the 'f' variable values
# for i in range (1009):
#     F=fft(lm[:,i])
#     f= (np.power(F.real[0:64],2)+np.power(F.imag[0:64],2))
#     fm[:,i]=f[:]
#     plt.figure(1)
#     plt.plot(freq[1:64], f[1:64])


































#favg=np.zeros((64))
#
#favg=np.mean(fm, axis=1)
#
##favg_SG=np.zeros((64))
##fstd=np.std(fm, axis=1)
##favg_SG = savgol_filter(favg, 17, 10)
#
##print(favg)
#plt.figure(1)
##plt.plot(freq[1:64], favg[1:64])
##plt.plot(fstd[1:64])   
#
#"""
#FFT of 9 true readings only
#"""
#f9=np.zeros((127,9))
#f9=lm[:,0:9]
#f9avg=np.mean(f9,1)
#FFT9=fft(f9avg)
##fz9= (np.power(FFT9.real[0:64],2)+np.power(FFT9.imag[0:64],2))
#fz9=abs(FFT9[0:64])
#plt.figure(2)
#plt.plot(freq[1:64], fz9[1:64], color="orange")
#
#"""
#FFT of 1000 synthetic readings only
#"""
#f1000=np.zeros((127,1000))
#f1000=lm[:,9:1009]
#f1000avg=np.mean(f1000,1)
#FFT1000=fft(f1000avg)
##fz1000= (np.power(FFT1000.real[0:64],2)+np.power(FFT1000.imag[0:64],2))
#fz1000=abs(FFT1000[0:64])
#plt.figure(3)
#plt.plot(freq[1:64], fz1000[1:64], color="black")



