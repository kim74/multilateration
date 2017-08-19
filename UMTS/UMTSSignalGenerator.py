# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 15:45:24 2015

@author: kimfeld
"""

import numpy as np
import scipy as sc
import scipy.signal as sig
from numpy.random import choice


from matplotlib import pyplot as plt

########physical channel parameters#########################################
slotsPerFrame=15
chipsPerSlot=2560
fc=0#10.0e6#15.0e6
overSamplingFact=2 #for the emulation of an analog signal
chipRate=3.84e6
Tsamp_adc=1/(overSamplingFact*chipRate)
slotPeriod=chipsPerSlot/chipRate
framePeriod=0.01
framePeriod=slotsPerFrame*slotPeriod
############################################################################ 
sscGroupSequence=np.array([1,1,2,8,9,10,15,8,10,16,2,7,15,7,16])

########Signals############################################################# 
physicalChannel1=np.zeros(2560)

############################################################################ 

############################################################################   
def Display(t,s,s_filtered):
    def get_color():
        for item in ['r', 'g', 'b', 'c', 'm', 'y', 'k','r']:
            yield item
    
    fig=plt.figure(1)        
    color = get_color()
    
    ax1=plt.subplot(2,1,1)
    acolor=next(color)
    plt.plot(t,s,color=acolor)#,label="$W={0:.2f}$".format())
    plt.subplot(2,1,2,sharex=ax1)    
    acolor=next(color)
    plt.plot(t,s_filtered,color=acolor)
    
    #plt.legend( loc='left center', bbox_to_anchor=(1.0,0.0,0.2,1.0),mode="expand", ncol=1, shadow=True, borderaxespad=0.,prop={'size':7}) 
      
    plt.title("transmitted UMTS signal")   
    plt.ylabel("amplitude [V]") 
    plt.xlabel("t [s]")

    plt.show() 
    return 0
############################################################################  
############################################################################  
def OverSample(s):
    aux=np.empty((0,))
    for k in range(overSamplingFact):
        aux=np.append(aux,s)
        
    aux=aux.reshape((overSamplingFact,len(s)))
    aux=aux.T
    T=aux.reshape((len(s)*overSamplingFact,))
    return T
############################################################################ 
############################################################################ 
def ReadFromTxtFile(fileName): 
    I=np.empty((0,))
    Q=np.empty((0,))
    f=open(fileName,'r') 
    for line in f:
        ll=line.strip()
        data=ll.split(" ")
        I=np.append(I,int(data[0]))
        Q=np.append(Q,int(data[1]))
    
    f.close()
    return I , Q
############################################################################ 
############################################################################ 
def GeneratePSC_SSC():
    #PSC-generation
    psc=np.empty((0,))    
    a=np.array([1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1])
    aModulate=np.array([1,1,1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1])

    for k in aModulate:
        if k==1:
            psc=np.append(psc, a)
        else:
            psc=np.append(psc,-a)

    psc=(1+1j)*psc 
    #SSC-generation
    b=np.empty((16,))
    b[0:8]=a[0:8]
    b[8:16]=-a[8:16]
    z=np.empty((0,))
    zModulate=np.array([1,1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1])
    
    for k in zModulate:
        if k==1:
            z=np.append(z,b)
        else:
            z=np.append(z,-b)
            
    H=np.empty((1,1))
    H[0,0]=1
    for k in range(1,9):
        row1=np.hstack((H,H))
        row2=np.hstack((H,-H))
        H=np.vstack((row1,row2))    

    ssc=np.empty((16,256))
    for k in range(16):
        m=16*k
        ssc[k,:]=z*H[m,:]
    
    ssc=(1+1j)*ssc
     
    return psc, ssc
############################################################################  
############################################################################
def GenerateScramblingCode(scramblingCodeGroup, codeNumber ):
    x=np.zeros(18)
    x[0]=1
    y=np.ones(18)
    z=np.empty((0,))
    for i in range(2**18-19):
        a=(x[i+7]+x[i])%2
        b=(y[i+10]+y[i+7]+y[i+5]+y[i])%2      
        x=np.append(x,a)
        y=np.append(y,b)
    mod=2**18-1
    for i in range(mod):
        index=(i+codeNumber)%mod
        a=(x[index]+y[i])%2
        if a==0:
            a=-1
        z=np.append(z,a)
        
    S=np.zeros(38400)
    
    for i in range(38400):
        index=(i+2**17)%mod
        S[i]=z[i]+1j*z[index]
        
    return S
############################################################################  
############################################################################ 
def ChannelCombinerPerSlot(psc,ssc,dataChips=0):
    #zero-pad psc and ssc to chipsPerSlot samples; data and PSCH/SSCH do not overlap
    z_pssc=np.zeros(chipsPerSlot-len(psc))
    z_data=np.zeros(len(psc))
    data=np.concatenate((z_data,dataChips))
    PSC=np.concatenate((psc,z_pssc))
    SSC=np.concatenate((ssc,z_pssc))
    return PSC+SSC+data
############################################################################ 
############################################################################
def Channelisation(data,chipLength):
    seqVal=np.array([-2,0,2])
    dataChips=choice(seqVal,chipLength)
    #dataChips=np.zeros(chipLength)
    return dataChips
############################################################################ 
############################################################################
def SlotGenerator(data,PSC,SSC):
    payloadChipsPerSlot=chipsPerSlot-len(PSC)
    dataChips=Channelisation(data,payloadChipsPerSlot)
    chipsInOneSlot=ChannelCombinerPerSlot(PSC,SSC,dataChips)
    return chipsInOneSlot,PSC,SSC
############################################################################  
############################################################################
def FrameGenerator(data):
    PSC,SSC=GeneratePSC_SSC()
    chipsPerFrame=slotsPerFrame*chipsPerSlot
    chipsInOneFrame=np.empty((0,))
    if len(data)<chipsPerFrame:
        data=np.concatenate((data,np.zeros(chipsPerFrame-len(data))))       
    for i in range(slotsPerFrame):
        slotData=data[i*chipsPerFrame:i*chipsPerFrame+chipsPerSlot]
        chipsInSlot,psc,ssc=SlotGenerator(slotData,PSC,SSC[i,:])
        chipsInOneFrame=np.concatenate((chipsInOneFrame,chipsInSlot))
    return chipsInOneFrame,PSC
############################################################################  
############################################################################     
def Modulator(data,t):
    I=np.real(data)
    Q=np.imag(data)
    carrier_I=np.cos(2*np.pi*fc*t)
    carrier_Q=-np.sin(2*np.pi*fc*t)
    tx=I*carrier_I+Q*carrier_Q
    return tx
############################################################################
############################################################################ 
def DelayAndPadChips(data,delayInChips,randomPaddingInChips):
    #shift composite signal by a known delay
    seqVal=np.array([-2,0,2])
    header=choice(seqVal,delayInChips)
    tailer=choice(seqVal,randomPaddingInChips)    
    dataPadded=np.concatenate((header,data,tailer))
    return dataPadded
############################################################################ 
def Demodulator(data, t, freq, phase):

    I=np.cos(2*np.pi*freq*t+phase)*data
    Q=-np.sin(2*np.pi*freq*t+phase)*data
    b,a = sc.signal.butter(4,2.5*chipRate/(2*overSamplingFact*chipRate),'low',analog=False)
    I=sig.lfilter(b,a,I)
    Q=sig.lfilter(b,a,Q)    
    return I,Q
############################################################################
############################################################################ 
    
delayInChips=500  #number of Chips of which the UMTS signal is initially delayed, filled with random sequence
randomPaddingInChips=1000 #random chip sequence added at the end of the UMTS signal
      
dataPayload=np.array([0])
#PSC,SSC=GeneratePSC_SSC()
chips,psc=FrameGenerator(dataPayload)

chipsData=DelayAndPadChips(chips,delayInChips,randomPaddingInChips)
chipsData_oversampled=OverSample(chipsData)
psc_oversampled=OverSample(psc)
h_mf=np.conj(psc_oversampled[::-1]) #create matched filter
#h_mf=np.conj(psc[::-1])#test without appropriate oversampling (such as the psc itself) of the matched filter 

time=np.arange(0,len(chipsData_oversampled)*Tsamp_adc,Tsamp_adc)
#data_modulated=Modulator(chipsData_oversampled,time)

freqReceiver=fc
phaseReceiver=0
#I,Q=Demodulator(data_modulated,time,freqReceiver,phaseReceiver)

s=sig.lfilter(h_mf,1,chipsData_oversampled)
#
block=40000
Display(time[0:block],chipsData_oversampled[0:block],s[0:block])

I,Q=ReadFromTxtFile("C:\\Users\\kimfeld\\Documents\\Projects\\Multilateration\\software\\VC++\\UMTSFrameExtractor\\UMTSFrameExtractor\\result.txt") 


