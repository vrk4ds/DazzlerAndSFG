#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import sys
import os
import re
import math
import matplotlib.pyplot as plt
import scipy.fftpack
import pylab
import h5py
from scipy.interpolate import interp1d
import csv
import pandas as pd
import pickle
import matplotlib.colors as colors

class Dazzler_Pulse_Shaper():  
    def __init__(self, position, width, hole_position, hole_width, hole_depth,
                 delay, sec_order, third_order, fourth_order, c = 299792458.0):
        '''
        Initializes all dazzler parameters.
        These include central wavelength and width of input field; hole position, width, and depth; and 
        dispersion factors (first, second, third, and fourth order)
        '''
        #use params on the order of wavelength in m
        self.position = position 
        self.width = width
        self.hole_position = hole_position
        self.hole_depth = hole_depth
        self.hole_width = hole_width
        
        self.delay = delay #minimum order of magnitude e-13
        self.sec_order = sec_order #minimum order of magnitude e-26
        self.third_order = third_order #minimum order of magnitude e-39
        self.fourth_order = fourth_order #minimum order of magnitude e-52
        
        self.c = c
    
        self.omega0 = 2*np.pi*self.c/self.position
        self.chi0 = self.width/(2*self.position)
        self.del_omega0 = self.omega0*(self.chi0-self.chi0**3)
        self.omega1 = 2*np.pi*self.c/self.hole_position
        self.chi1 = self.hole_width/(2*self.hole_position)
        self.del_omega1 = self.omega1*(self.chi1-self.chi1**3)/2        
        
    def set_saved_pulse(self, position, width, hole_position, hole_width, hole_depth,
                 delay, sec_order, third_order, fourth_order):
        self.savedPulse=Dazzler_Pulse_Shaper(position, width, hole_position, hole_width, hole_depth,
                 delay, sec_order, third_order, fourth_order)
        
    def calculate_amplitude_transfer_function(self, ang_freq_vector, components_return = False):
        '''
        THis function takes an angular frequency vector and uses the class parameters to 
        calculate the amplitude transfer function. This ends up being a product of 
        e^-([(w-w0)/delw0]^6) * (1-k*e^-([(w-w0)/delw0]^2))
        '''
        #calculate f and g
        f =(np.exp(-((ang_freq_vector-self.omega0)/self.del_omega0)**6))
        g = 1 - self.hole_depth*(np.exp(-((ang_freq_vector-self.omega1)/self.del_omega1)**2))
        
        # return either f,g, and A or just A, where A = f*g
        if (components_return):
            return f, g, f*g
        return f*g
    
    def calculate_phase_transfer_function(self, ang_freq_vector):
        '''
        This function calculates the phase transfer function. 
        '''
        #calculate freq difference
        omega_dif = ang_freq_vector-self.omega0
        #return h as shown in equations in text
        return -(self.delay*omega_dif + self.sec_order/2 * omega_dif**2 + self.third_order/6 * omega_dif**3 + 
                 self.fourth_order/24 * omega_dif**4)
    
    def calculate_full_transfer_function(self,ang_freq_vector, S_saved=0, a_saved=0, components_return=False):
        '''
        This function calculates the full transfer function which is the product of the 
        phase and amplitude transfer functons. 

        Parameters
        ----------
        ang_freq_vector : TYPE
            DESCRIPTION.
        S_saved : TYPE, optional
            DESCRIPTION. The default is 0.
        a_saved : TYPE, optional
            DESCRIPTION. The default is 0.
        components_return : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        f,g,A_dial = None,None, None
        
        
         #calculate A dial
        if (components_return):
            f,g,A_dial = self.calculate_amplitude_transfer_function(ang_freq_vector, components_return=True)
        else:
            A_dial = self.calculate_amplitude_transfer_function(ang_freq_vector, components_return=False)
            
        #calculate phi dial
        phi_dial = self.calculate_phase_transfer_function(ang_freq_vector)
        
        #calculate S
        S_full = A_dial*np.exp(1j*phi_dial) + a_saved*np.exp(1j*phi_dial)*S_saved
        
        
        if (components_return):
            return f,g,A_dial,phi_dial, S_full
        
        return S_full
 
    
    def shape_input_pulse(self, E_field_input, time_vector, sampling_rate, wavelength_vec_limits,
                          S_saved=0, a_saved=0, components_return=False):
        '''
        This function takes the input electric field, time vector, and sampling rate and calls
        on other functions in the class in order to calculate output electric field
        '''
        #take FT of input electric field
        E_field_input_ft = np.fft.fft(E_field_input)
        
        #get frequency vector and angular frequency vector based on sampling rate
        freq_vector = np.fft.fftfreq(n=(E_field_input_ft.size), d=1/sampling_rate) #for our laser, sampling_rate needs to be at least 240/pulsewidth
        ang_freq_vector = 2*np.pi*freq_vector

        #calculate transfer function
        f,g,A_dial,phi_dial, S_full = None, None, None, None, None
        if (components_return):
            f,g,A_dial,phi_dial, S_full = self.calculate_full_transfer_function(ang_freq_vector, S_saved, a_saved, components_return)
            components_dict = {"f":f,"g":g,"A_dial":A_dial,"phi_dial":phi_dial,"S_full":S_full}

        else :
            S_full = self.calculate_full_transfer_function(ang_freq_vector, S_saved, a_saved, components_return)
            components_dict = {"S_full":S_full}
            
        #calculate output frequency domain electric field
        E_field_output_ft = E_field_input_ft*S_full
        
        #E_field_output = E_field_input*np.fft.ifft(S_full)
        
        #calculate IFT of output field
        E_field_output = np.fft.ifft(E_field_output_ft)
        
        
        return E_field_input, E_field_input_ft, E_field_output, E_field_output_ft,time_vector, freq_vector, components_dict
    
    def calculate_spectrum_phase(self, field):
        '''
        THis function separates spectrum and phase, unwrapping phase
        '''
        phase = -1*np.arctan2(np.imag(field),np.real(field))
        phase = np.unwrap(phase)
        spectrum = (np.abs(field))**2
        return spectrum, phase
    
    def constraints_wavelength_vec_limits(self, freq_vec): #use as check if convert_to_wavelength interpolation is failing (likely need to increase sampling rate)
        maxwave=self.c/freq_vec[1]
        minwave=self.c/max(np.abs(freq_vec))
        print("The wavelengths must be between {} and {}".format(minwave, maxwave))
    
    def convert_to_wavelength(self,I_freq, phase_freq, freq_vec,c, wavelength_vec_limits):
        '''
        This function converts the frequency domain signal into a wavelength domain. 
        Defaults to a wavelength vector that goes from 200nm below to 200nm above the central 
        wavelength. This should be adjusted if width of distribution is very large. 

        Parameters
        ----------
        I_freq : TYPE
            DESCRIPTION.
        phase_freq : TYPE
            DESCRIPTION.
        freq_vec : TYPE
            DESCRIPTION.
        c : TYPE
            DESCRIPTION.
        wavelength_vec_range : TYPE, optional
            DESCRIPTION. The default is 200e-9.

        Returns
        -------
        None.

        '''

        wavelength_vector = np.linspace(wavelength_vec_limits[0],wavelength_vec_limits[1],num=len(freq_vec))
        I_freq_interp = interp1d(2*np.pi*freq_vec, I_freq)
        I_wavelength = (2*np.pi*c/(wavelength_vector**2))*I_freq_interp(2*np.pi*c/wavelength_vector)
        phase_freq_interp = interp1d(2*np.pi*freq_vec, phase_freq)
        phase_wavelength = phase_freq_interp(2*np.pi*c/wavelength_vector)

        return wavelength_vector, I_wavelength, phase_wavelength    
    
    def calculate_total_spectrum_phase(self,E_field_input, E_field_input_ft, E_field_output, E_field_output_ft, S_full,
                             time_vector, freq_vector, wavelength_vector_limits):
        # Input 
        intensity_input_td, phase_input_td = self.calculate_spectrum_phase(E_field_input) # Time domain (td)
        spectrum_input_fd, phase_input_fd = self.calculate_spectrum_phase(E_field_input_ft) # Frequency domain (fd)
        wavelength_vector, spectrum_input_wd, phase_input_wd = self.convert_to_wavelength(spectrum_input_fd,
                                                                                     phase_input_fd,freq_vector,
                                                                                     self.c,wavelength_vec_limits=wavelength_vector_limits) #Wavelength domain (wd)
        
        input_functions = {"time_vector": time_vector,
                           "freq_vector":freq_vector,
                            "wavelength_vector":wavelength_vector,
                            "intensity_input_td":intensity_input_td,
                            "phase_input_td":phase_input_td,
                            "spectrum_input_fd":spectrum_input_fd,
                            "phase_input_fd":phase_input_fd,
                            "spectrum_input_wd":spectrum_input_wd,
                            "phase_input_wd":phase_input_wd}
                
        
        # Transfer function 
        S_full_td = np.fft.ifft(S_full)
        intensity_transfer_td, phase_transfer_td = self.calculate_spectrum_phase(S_full_td) # Time domain (td)
        spectrum_transfer_fd, phase_transfer_fd = self.calculate_spectrum_phase(S_full) # Frequency domain (fd)
        wavelength_vector,spectrum_transfer_wd, phase_transfer_wd = self.convert_to_wavelength(spectrum_transfer_fd,phase_transfer_fd,
                                                                        freq_vector, self.c,
                                                                        wavelength_vec_limits=wavelength_vector_limits) #Wavelength domain (wd) 
        
        transfer_functions = {"time_vector": time_vector,
                              "freq_vector": freq_vector,
                              "wavelength_vector": wavelength_vector, 
                              "intensity_transfer_td": intensity_transfer_td,
                              "phase_transfer_td": phase_transfer_td,
                              "spectrum_transfer_fd":spectrum_transfer_fd,
                              "phase_transfer_fd":phase_transfer_fd,
                              "spectrum_transfer_wd":spectrum_transfer_wd,
                              "phase_transfer_wd":phase_transfer_wd}
              
        # Output  
        intensity_output_td, phase_output_td = self.calculate_spectrum_phase(E_field_output) # Time domain (td)
        spectrum_output_fd, phase_output_fd = self.calculate_spectrum_phase(E_field_output_ft) # Frequency domain (fd)
        wavelength_vector, spectrum_output_wd, phase_output_wd = self.convert_to_wavelength(spectrum_output_fd,phase_output_fd,
                                                                        freq_vector, self.c,
                                                                        wavelength_vec_limits=wavelength_vector_limits) #Wavelength domain (wd)
        output_functions = {"time_vector": time_vector,
                            "freq_vector":freq_vector,
                            "wavelength_vector":wavelength_vector,
                            "intensity_output_td":intensity_output_td,
                            "phase_output_td":phase_output_td,
                            "spectrum_output_fd":spectrum_output_fd,
                            "phase_output_fd":phase_output_fd,
                            "spectrum_output_wd":spectrum_output_wd,
                            "phase_output_wd":phase_output_wd}
        
        return input_functions, transfer_functions, output_functions
    

    def outputforSFGCode(self, time_vector, E_field_output, wavelength_vector_limits, sampling_rate, filename, pulsewidth0, central_wavelength=None):        
        
        wavelength_vector = np.linspace(wavelength_vector_limits[0],wavelength_vector_limits[1],num=len(time_vector))
        if central_wavelength==None:
            central_wavelength=(wavelength_vector_limits[1]-wavelength_vector_limits[0])/2
    
        newdict={"E_field": E_field_output,
                "time_vector": time_vector,
                "wavelength_vector": wavelength_vector}
        
        def space_vecs(dictionary, modulus): #cut down time_vector spacing so it is usable for SFG code 
            #sampling_rate is currently a function of pulsewidth
            for key, value in dictionary.items():
                newlist=[]
                for i in range(len(value)):
                    if i%modulus==0:
                        newlist.append(value[i])
                dictionary[key]=newlist
            return dictionary
        
        for i in range(10,31):
            if (sampling_rate*pulsewidth0) % i == 0:
                modulus=i
                break
        else:
            print("Warning: the SFG code prefers a sampling rate between 10 and 30 per pulsewidth. Consider making your sampling rate a multiple in that range.")

        newFrames = space_vecs(newdict, modulus)
        newdict["central_wavelength"]=central_wavelength

        file1 = open(filename,"wb")
        pickle.dump(newFrames, file1)
        file1.close()
     
    def make_spectrogram(self, E_field, time_vector, gate_starting_index):
        def swap_halves(array):
            a=len(array)
            b=a//2
            if a%2==0:
                list2=list(array[0:b])
                list3=list(array[b:a])
            if a%2==1:
                list2=list(array[0:b+1])
                list3=list(array[b+1:a])
            list_swapped_halves=list3+list2
            return list_swapped_halves
        def EofTGateNew(tcut, tau):
            returnval=[]
            i1=tau*sampling_rate #make integer #this is the amount of indices that our tau displaces the time vector
            for i2 in range(len(tcut)):
                istart=-round((t[0]-tcut[0])*sampling_rate) #make integer
                a=int(istart+i2-i1)
                out=output.get("intensity_output_td")[a]
                returnval.append(out)
            return returnval

        taulist=time_vector[gate_starting_index:-gate_starting_index]
        tcut_index=(len(time_vector)/2)-gate_starting_index
        tcut=time_vector[tcut_index:-tcut_index]

        istart=-round((t[0]-tcut[0])*sampling_rate)
        E_field_output_cut=E_field_output[istart:-istart]
        intensityarray=np.zeros((len(tcut),len(taulist)))
        i=0
        for tau in taulist:
            intensity_FROG=(np.fft.fft(E_field_output_cut*(EofTGateNew(tcut, tau))))**2 #may need EofTGate**2 #this gives 1D array of all frequency for given tau
            intensityarray[:,i]=swap_halves(np.abs(intensity_FROG))
            i+=1
        freq_vector = np.fft.fftfreq(n=intensity_FROG.size, d=1/sampling_rate) #for our laser, sampling_rate needs to be at least 240/pulsewidth
        freq_vector = swap_halves(freq_vector)

        #plt.xlim(-3e-12,3e-12)
        #plt.ylim(2.5e14,3.5e14)
        plt.pcolormesh(taulist,freq_vector,intensityarray, norm=colors.SymLogNorm(linthresh=0.003, linscale=0.003,vmin=1e-27, vmax=57000), shading='gouraud')
        plt.colorbar()
        plt.show()

#our laser: 1022nm, 439.3fs pulse width ?


# In[4]:


storedWaveformDemo=Dazzler_Pulse_Shaper(1010,10,1030,10,0,0,0,0,0)
storedWaveformDemo.set_saved_pulse(1050,10,1030,10,0,0,0,0,0)

wavelength_vector=np.arange(1030-200, 1030+200)
freq_vector=3e8/wavelength_vector
a_saved=storedWaveformDemo.savedPulse.calculate_amplitude_transfer_function(2*np.pi*freq_vector)
S_saved=storedWaveformDemo.savedPulse.calculate_full_transfer_function(2*np.pi*freq_vector)

transferfunction=storedWaveformDemo.calculate_full_transfer_function(2*np.pi*freq_vector, S_saved, a_saved, components_return=False)
plt.plot(wavelength_vector, np.abs(transferfunction))


# In[5]:


c=299792458
lambda0=1030e-9
omega0=2*np.pi*c/lambda0
pulsewidth=85e-15
sampling_rate=1000/pulsewidth
t=np.arange(-10e-12,10e-12,1/sampling_rate)
AofT=np.exp(-(t/pulsewidth)**2)
EofT=AofT*np.exp(1j*omega0*t)


# In[24]:


c=299792458
lambda0=1030e-9
omega0=2*np.pi*c/lambda0

pulsewidth=330e-15
sampling_rate=240/pulsewidth

t=1/sampling_rate*np.arange(-600000,600000) #1/sampling_rate*np.arange(-16384, 16384) #for Amy 
print(len(t))
AofT=np.exp(-1.386*(t/pulsewidth)**2)
EofT=AofT*np.exp(1j*omega0*t)
plt.grid()
plt.plot(t, EofT)


# In[58]:


ourLaserTrial=Dazzler_Pulse_Shaper(1030e-9,20,1030e-9,10e-9,0,0, 5e-26,0,0)
wavelength_vector_limits=(lambda0-200e-9,lambda0+200e-9)
E_field_input, EofW, E_field_output, E_field_output_ft,time_vector, freq_vector, components_dict = ourLaserTrial.shape_input_pulse(EofT, t, sampling_rate,wavelength_vector_limits)
S_full = ourLaserTrial.calculate_full_transfer_function(2*np.pi*freq_vector)
input1, transfer, output = ourLaserTrial.calculate_total_spectrum_phase(E_field_input, EofW, E_field_output, E_field_output_ft, S_full, time_vector, freq_vector,wavelength_vector_limits)

#ourLaserTrial.outputforSFGCode(t, E_field_output, wavelength_vector_limits, sampling_rate, "Output_NewPulsewidth_NoPhaseShape", 330e-15, 1030e-9)   


def swap_halves(array):
    a=len(array)
    b=a//2
    if a%2==0:
        list2=list(array[0:b])
        list3=list(array[b:a])
    if a%2==1:
        list2=list(array[0:b+1])
        list3=list(array[b+1:a])
    list_swapped_halves=list3+list2
    return list_swapped_halves

print(freq_vector)

fig, freqplots = plt.subplots(4,figsize=(5,11))
#freqplots[0].plot(t, np.abs(E_field_input))
freqplots[0].set_title("E-Field Input")
freqplots[0].plot(t, input1['phase_input_td'])
freqplots[1].set_title("E-Field Output")
#freqplots[2].set_xlim(2e14,4e14)
freqplots[1].plot(t, output['phase_output_td'])
freqplots[2].set_title("Frequency Domain Output")
freqplots[2].plot(swap_halves(output['freq_vector']), input1['phase_input_fd'])
freqplots[3].plot(swap_halves(output['freq_vector']), output['phase_output_fd'])


def space_vecs(dictionary, modulus): #cut down time_vector spacing so it is usable for SFG code 
            #sampling_rate is currently a function of pulsewidth
            for key, value in dictionary.items():
                newlist=[]
                for i in range(len(value)):
                    if i%modulus==0:
                        newlist.append(value[i])
                dictionary[key]=newlist
            return dictionary
        
        
newdict={"E_field": E_field_output,
        "time_vector": t,
        "wavelength_vector": output['wavelength_vector']}

newdictcondensed=space_vecs(newdict, 2)

print(len(newdictcondensed["E_field"]))
newdictcondensed['central_wavelength'] = 1030e-9

#import pickle

#file1 = open("OutputforAmy_SamplingRate200_noshape.txt","wb")
#pickle.dump(newdictcondensed, file1)
#file1.close()


# In[56]:


t=10/sampling_rate*np.arange(-800,800)
AofT=np.exp(-1.386*(t/pulsewidth)**2)
EofT=AofT*np.exp(1j*omega0*t)
#plt.plot(t,EofT)

phase = -1*np.arctan2(np.imag(EofT),np.real(EofT))
phase = np.unwrap(phase)
spectrum = (np.abs(EofT))**2
E_field_input_ft=np.fft.fft(EofT)
freq_vector = np.fft.fftfreq(n=(E_field_input_ft.size), d=10/sampling_rate) #for our laser, sampling_rate needs to be at least 240/pulsewidth
#plt.xlim(2.5e14,3.5e14)
plt.plot(swap_halves(freq_vector), swap_halves(np.abs(E_field_input_ft)**2))
#np.fft.fft(EofT)


# In[135]:


import matplotlib.colors as colors


AofT=np.exp(-(t/pulsewidth)**2)
EofT=AofT*np.exp(1j*omega0*t)

taulist=np.arange(-3e-12, 3e-12, 1/sampling_rate)
#AofTGate=np.exp(-((t-tau)/pulsewidth)**2)
#EofTGate=AofTGate*np.exp(1j*omega0*(t-tau))
#Esig=EofT*EofTGate**2
#plt.plot(t, Esig)
#plt.plot(t, EofTGate)
#intensity_FROG=(np.fft.fft(E_field_input*EofTGate**2, axis=t))**2


def swap_halves(array):
    a=len(array)
    b=a//2
    if a%2==0:
        list2=list(array[0:b])
        list3=list(array[b:a])
    if a%2==1:
        list2=list(array[0:b+1])
        list3=list(array[b+1:a])
    list_swapped_halves=list3+list2
    return list_swapped_halves

def EofTGate(time, tau):
    EofTGate=np.exp(-((time-tau)/pulsewidth)**2)*np.exp(1j*omega0*(time-tau))
    return EofTGate

intensityarray=np.zeros((len(EofT),len(taulist)))
i=0
for tau in taulist:
    intensity_FROG=(np.fft.fft(E_field_input*(EofTGate(t, tau)**2)))**2 #this gives 1D array of all frequency for given tau
    intensityarray[:,i]=np.abs(swap_halves(intensity_FROG))
    i+=1
freq_vector = np.fft.fftfreq(n=intensity_FROG.size, d=1/sampling_rate) #for our laser, sampling_rate needs to be at least 240/pulsewidth
#print(freq_vector)


sumAlongFreq=[]
for i in range(len(EofT)):
    jsum=0
    for j in range(len(taulist)):
        jsum+=intensityarray[i, j]
        #if intensityarray[i,j]>=0.2:
            #print(i,j,intensityarray[i,j])
    sumAlongFreq.append(jsum)
wavelength_vector = np.linspace(830e-9,1230e-9,num=len(freq_vector))
#I_freq_interp = interp1d(2*np.pi*freq_vector, sumAlongFreq)
#I_wavelength = (2*np.pi*c/(wavelength_vector**2))*I_freq_interp(2*np.pi*c/wavelength_vector)
#plt.plot(wavelength_vector, I_wavelength)



freq_vector=swap_halves(freq_vector)
print(len(freq_vector))
#plt.plot(freq_vector, sumAlongFreq)
#plt.xlim(-1.5e-12, 1.5e-12)
plt.ylim(2.5e14,3.5e14)
plt.pcolormesh(taulist,freq_vector,intensityarray, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=1e-27, vmax=57000), shading='gouraud')
plt.colorbar()
plt.show()


# In[181]:


plt.plot(t, np.abs(EofT))
total0=[0]*len(t)
plt.plot(t,total0)
gate1=1/sampling_rate*np.arange(1000,1000)
gate2=1/sampling_rate*np.arange(800,800)
total1=[0.1]*len(gate1)
total2=[0.2]*len(gate2)
plt.plot([-1000/sampling_rate, 1000/sampling_rate], [0.05,0.05])
plt.plot([-600/sampling_rate, 600/sampling_rate], [0.1,0.1])

#plt.plot(gate1, total1)
#plt.plot(gate2, total2)


# In[146]:


def EofTGateNew(tcut, tau):
    returnval=[]
    i1=tau*sampling_rate #make integer #this is the amount of indices that our tau displaces the time vector
    for i2 in range(len(tcut)):
        istart=-round((t[0]-tcut[0])*sampling_rate) #make integer
        a=int(istart+i2-i1)
        #print(type(a))
        #print(i1,i2,a, istart)
        out=output.get("spectrum_output_fd")[a]
        returnval.append(out)
    returnvalsquared=[number**2 for number in returnval]
    return returnvalsquared

#taulist=1/sampling_rate*np.arange(-1200,1200)
#tcut=1/sampling_rate*np.arange(-400,400)
#print(len(time_vector))

tcut_index=(len(time_vector)/2)-600
print(tcut_index)
    
taulist=time_vector[600:-600]
tcut=time_vector[1000:-1000]
istart=-round((t[0]-tcut[0])*sampling_rate)
E_field_output_cut=E_field_output[istart:-istart]
intensityarray=np.zeros((len(tcut),len(taulist)))
i=0
for tau in taulist:
    intensity_FROG=(E_field_output_cut*(EofTGateNew(tcut, tau)))**2 #may need EofTGate**2 #this gives 1D array of all frequency for given tau
    #intensity_FROG=(np.fft.fft(E_field_output_cut*(EofTGateNew(tcut, tau))))**2 #may need EofTGate**2 #this gives 1D array of all frequency for given tau
    intensityarray[:,i]=swap_halves(np.abs(intensity_FROG))
    i+=1
freq_vector = np.fft.fftfreq(n=intensity_FROG.size, d=1/sampling_rate) #for our laser, sampling_rate needs to be at least 240/pulsewidth
freq_vector = swap_halves(freq_vector)

sumAlongTime=[]
for i in range(len(taulist)):
    isum=0
    for j in range(len(tcut)):
        isum+=intensityarray[j, i]
        #if intensityarray[i,j]>=0.2:
            #print(i,j,intensityarray[i,j])
    sumAlongTime.append(isum)

plt.plot(taulist,sumAlongTime)

#plt.xlim(-3e-12,3e-12)
#plt.ylim(2.5e14,3.5e14)
#plt.pcolormesh(taulist,freq_vector,intensityarray, norm=colors.SymLogNorm(linthresh=0.003, linscale=0.003,vmin=1e-27, vmax=57000), cmap='hsv', shading='gouraud')
#plt.colorbar()
#plt.show()


# In[ ]:





# In[31]:


import matplotlib.animation as animation
#from matplotlib.animation import FuncAnimation

#%matplotlib qt
E_field_outputs=[]
for a in np.arange(5e-26,50e-26,5e-26):
    eachframe=Dazzler_Pulse_Shaper(1030e-9,2000000,1030e-9,10,0,0,a,0,0)
    E_field_input, E_field_input_ft, E_field_output, E_field_output_ft,time_vector, freq_vector, components_dict = eachframe.shape_input_pulse(E_field_input, time_vector, sampling_rate,wavelength_vector_limits)
    E_field_outputs.append(E_field_output)
    
fig, ax = plt.subplots()
line, = ax.plot(time_vector, np.abs(E_field_outputs[0]))


def animate(i):
    line.set_ydata(np.abs(E_field_outputs[i]))  # update the data.
    return line,


ani = animation.FuncAnimation(fig, animate), #interval=1), #blit=True, save_count=50) 


# In[116]:


newFrames=Dazzler_Pulse_Shaper(1030e-9,2000000,1030e-9,10,0,0,30e-26,0,0)
E_field_input, E_field_input_ft, E_field_output, E_field_output_ft,time_vector, freq_vector, components_dict = newFrames.shape_input_pulse(E_field_input, time_vector, sampling_rate,wavelength_vector_limits)
wavelength_vector = np.linspace(1030e-9,1230e-9,num=len(t))

plt.plot(t,np.abs(E_field_output))

newdict={"E_field": E_field_output,
        "time_vector": t,
        "wavelength_vector": output['wavelength_vector']}

newFrames.space_vecs(newdict)
newdict["central_wavelength"]=1030e-9

#import pickle

#file1 = open("Dazzler4KeyOutput_noShape.txt","wb")
#pickle.dump(newdict, file1)
#file1.close()


# In[ ]:




