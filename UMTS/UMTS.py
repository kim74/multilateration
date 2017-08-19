# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 20:09:39 2015

@author: kimfeld
"""

import numpy as np
import scipy as sc
import scipy.signal as sig
from numpy.random import choice
from matplotlib import pyplot as plt
import math
from numpy.fft import fft, ifft

def SaveToTxtFile(filename,a,b):
    f=open(filename,'w')          
    for i in range(len(a)):
        s="{0} {1}\n".format(int(a[i]),int(b[i]))
        f.write(s)
    f.close()

def ReadFromTxtFile(fileName, nbPoints): 
    I=np.zeros(nbPoints)
    Q=np.zeros(nbPoints)
    f=open(fileName,'r') 
    i=0
    for line in f:
        ll=line.strip()
        data=ll.split(" ")
        if i<nbPoints:
            I[i]=int(data[0])
            Q[i]=int(data[1])
            i+=1
        else:
            break
    f.close()
    #'i' is number of samples read from file
    return I , Q, i
#####################################################
#####################################################

class CPhysicalChannel:
    def __init__(self):
        self.slotsPerFrame=15
        self.chipsPerSlot=2560
        self.chipRate=3.84e6
        self.slotPeriod=self.chipsPerSlot/self.chipRate
        self.framePeriod=self.slotsPerFrame*self.slotPeriod
        self.data=0
    
    def Channelisation(self):
        """
        Implements first part of Figure 8 in TS25.213
        returns complex data stream 
        return I+jQ
        """
        return self.data
        
        
class C_CPICH(CPhysicalChannel):
    def __init__(self):
        CPhysicalChannel.__init__(self)
        self.nbSamplesOfZeroTX=256
        seqVal=np.array([-1+0j,1+0j])
        self.data=choice(seqVal,self.slotsPerFrame*self.chipsPerSlot)
        for i in range(self.slotsPerFrame):
            self.data[self.chipsPerSlot*i:self.chipsPerSlot*i+self.nbSamplesOfZeroTX]=np.zeros(self.nbSamplesOfZeroTX,dtype=complex)
        self.data=(1+1j)*self.data
        
    def Channelisation(self):
        self.data=np.ones(self.chipsPerSlot*self.slotsPerFrame,dtype=complex)#Channelisation Code for CPICH is C(ch,256,0); all 1 according to Figure 4 in TS 25.213
        return self.data
        
class C_SCH(CPhysicalChannel):
    def __init__(self):
        CPhysicalChannel.__init__(self)
        self.a=np.array([1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1])

class C_PSC(C_SCH):
    def __init__(self):
        C_SCH.__init__(self)
        self.psc=np.empty((0,),dtype=complex) 
        aModulate=np.array([1,1,1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1])
        
        for k in aModulate:
            if k==1:
                self.psc=np.append(self.psc, self.a)
            else:
                self.psc=np.append(self.psc,-self.a)

        self.psc=(1+1j)*self.psc
        
        self.data=np.zeros(self.slotsPerFrame*self.chipsPerSlot,dtype=complex)
        for i in range(self.slotsPerFrame):
            self.data[self.chipsPerSlot*i:self.chipsPerSlot*i+len(self.psc)]=self.psc
        
class C_SSC(C_SCH):
    def __init__(self,codeGroupSequence):
        C_SCH.__init__(self)
        self.b=np.empty((16,))
        self.b[0:8]=self.a[0:8]
        self.b[8:16]=-self.a[8:16]
        self.z=np.empty((0,))
        zModulate=np.array([1,1,1,-1,1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1])
        
        for k in zModulate:
            if k==1:
                self.z=np.append(self.z,self.b)
            else:
                self.z=np.append(self.z,-self.b)
                
        self.H=np.empty((1,1))
        self.H[0,0]=1
        for k in range(1,9):
            row1=np.hstack((self.H,self.H))
            row2=np.hstack((self.H,-self.H))
            self.H=np.vstack((row1,row2))    
    
        self.ssc=np.empty((16,256),dtype=complex)
        for k in range(16):
            m=16*k
            self.ssc[k,:]=self.z*self.H[m,:]
        
        self.ssc=(1+1j)*self.ssc  

        self.data=np.zeros(self.slotsPerFrame*self.chipsPerSlot,dtype=complex)
        for i in range(self.slotsPerFrame):
            #-1 is necessary since the secondary SCH sequences are defined from 1..16, and not as indices from 0 to 15 
            self.data[self.chipsPerSlot*i:self.chipsPerSlot*i+self.ssc.shape[1]]=self.ssc[codeGroupSequence[i]-1]

class CScramblingCode:
    def __init__(self,codeNumber):
        """
        Implements second part of Figure 8 in TS25.213
        return S
        """
        self.sdl=0
        x=np.zeros(2**18)
        x[0]=1
        y=np.ones(2**18)
        z=np.zeros(2**18)
        for i in np.arange(2**18-19):
            a=(x[i+7]+x[i])%2.0
            b=(y[i+10]+y[i+7]+y[i+5]+y[i])%2.0     
            x[i+18]=a
            y[i+18]=b
        
        #print "{0} , {1}".format(x[84],y[1045])
        #print "completed 1"
        mod=2.0**18-1
        for i in np.arange(mod):
            index=(i+codeNumber)%mod
            a=(x[index]+y[i])%2.0
            if a==0:
                a=1
            else:
                a=-1
            z[i]=a
        #print "completed 2"   
        S=np.zeros(38400,dtype=complex)
        
        for i in range(38400):
            index=(i+2**17)%mod
            S[i]=z[i]+1j*z[index]
        
        #SaveToTxtFile("result_scramblingCode_python.txt",S.real,S.imag)
        self.sdl=S

class CDownlink:
    def __init__(self):
        self.G_psc=1/2.0
        self.G_ssc=1/2.0
        self.G_cpich=1

    
    def PhysicalChannelCombiner(self,psc,ssc,S):
        """
        Implements Figure 9 in TS25.213
        return T
        """
        psc_temp=np.zeros(len(S),dtype=complex)
        ssc_temp=np.zeros(len(S),dtype=complex)
        psc_temp[0:len(psc)]=psc
        ssc_temp[0:len(ssc)]=ssc
        return self.G_psc*psc_temp+self.G_ssc*ssc_temp+S
    
    def Scrambling(self,data,code):
        I=np.array([1,1,1,-1,1,-1,1,1,1,-1])
        Q=np.array([-1,-1,-1,1,1,1,1,-1,-1,-1])
        d=I+Q*1j
        nb=len(code)
        y=np.zeros(nb,dtype=complex)
        nbChipsPerSymbol=len(y)/len(d)
        i=0
        for dd in d:
            y[i*nbChipsPerSymbol:(i+1)*nbChipsPerSymbol]=dd*code[i*nbChipsPerSymbol:(i+1)*nbChipsPerSymbol]
            i=i+1        
        return y

        
class CUMTSCell:
    def __init__(self,overSamplingFact,codeNumber):
        self.fc=0#10.0e6#15.0e6
        self.sscGroupSequence=np.array([1,2,3,1,8,6,5,2,5,8,4,4,6,3,7])
        self.psc=C_PSC();
        self.ssc=C_SSC(self.sscGroupSequence)
        self.cpich=C_CPICH()
        self.Tsamp_adc=1/(overSamplingFact*self.psc.chipRate)        
        self.umtsGenerator=CDownlink()
        self.codeNumber=codeNumber#400
        scrambCode=CScramblingCode(codeNumber)
        self.scramblingCode=scrambCode.sdl
               
    def Frame(self):
        #print scramblingCode[34550:34560]
        s=self.umtsGenerator.Scrambling(self.cpich.Channelisation(),self.scramblingCode)
        self.t=self.umtsGenerator.PhysicalChannelCombiner(self.psc.data,self.ssc.data,s)  
        #print s[34550:34560]
               
        return self.t
        
    def OverSampling(self,s,overSamplingFact):
        aux=np.empty((0,))
        for k in range(overSamplingFact):
            aux=np.append(aux,s)
            
        aux=aux.reshape((overSamplingFact,len(s)))
        aux=aux.T
        T=aux.reshape((len(s)*overSamplingFact,))
        time=np.arange(0,len(T)*self.Tsamp_adc,self.Tsamp_adc)
        return time,T
        
    def GenerateNFrames(self,N):
        frameLength=len(self.t)
        signalNFrames=np.zeros(N*frameLength,dtype=complex)      
        for i in range(N):
            signalNFrames[i*frameLength:(i+1)*frameLength]=self.t
            
        self.t=signalNFrames
        time=np.arange(0,frameLength*N*self.Tsamp_adc,self.Tsamp_adc)
        return time,self.t
        
    def StreamFromFile(self, filename, nbSamples, samplingFrequency):
        I,Q,nbSamplesRead=ReadFromTxtFile(filename, nbSamples)
        self.t=I+Q*1j
        self.Tsamp_adc=1/samplingFrequency
        time=np.arange(0,nbSamplesRead*self.Tsamp_adc,self.Tsamp_adc)
        return time, self.t
 
#########################            
def get_color():
    for item in ['r', 'g', 'b', 'c', 'm', 'y', 'k','r']:
        yield item               
#########################  
def Display(s1,s2,s3):

    
    fig=plt.figure(1)        
    color = get_color()
    
    ax1=plt.subplot(3,1,1)
    acolor=next(color)
    plt.plot(np.arange(len(s1)),s1,color=acolor)
    plt.subplot(3,1,2,sharex=ax1)    
    acolor=next(color)
    plt.plot(np.arange(len(s2)),s2,color=acolor)
    plt.subplot(3,1,3,sharex=ax1)    
    acolor=next(color)
    plt.plot(np.arange(len(s3)),s3,color=acolor)

    plt.show() 
    return 0

def filter_FFT(h,s):
    y=sig.fftconvolve(s,h,mode='same')
    return y    
    
def Energy(time,signal):
    E=0
    delta=time[1]-time[0]
    for i in signal:
        E+=delta*i**2
    return E/(time[-1]-time[0])
    
def Matched_Filter(h,s):
    if len(s)>len(h):
        N=len(s)
    else:
        N=len(h)
    
    H=fft(h,N)
    S=fft(s,N)
    Y=np.conj(H)*S#*(W**(len(h)*k))
    y=ifft(Y)
    return y
    
############
#
#nbSamples=76800
#umtsFreq=3840000
#samplingFreq=7679975.0
#nFrames=5
#codeNumber=400#6064
#
#overSamplingFact=int(round(samplingFreq/umtsFreq)) #for the emulation of an analog signal
#
#umts=CUMTSCell(overSamplingFact,codeNumber)
#signal=umts.Frame()
#
#time,ssc=umts.OverSampling(umts.ssc.ssc[3],overSamplingFact)
#time,psc=umts.OverSampling(umts.psc.psc,overSamplingFact)
#
#time,scramblingCode=umts.OverSampling(umts.scramblingCode,overSamplingFact)
#time,signal=umts.StreamFromFile("C:\\Users\\kimfeld\\Documents\\Projects\\Multilateration\\software\\VC++\\UMTSFrameExtractor\\UMTSFrameExtractor\\data\\TSMW-300k-Samples-2147.6MHz_6.txt",nbSamples,samplingFreq)
#
##time,signal=umts.GenerateNFrames(nFrames)
##SaveToTxtFile("scramblingCode.txt",scramblingCode.real,scramblingCode.imag)
##nb=153600
##signal=signal[nbSamples-nb:]    
##print signal.shape
#nbPoints=len(signal)
#
#S=fft(signal,nbPoints)
#scramblingCodeLength=len(scramblingCode)
#
#scramblingCodePartial=scramblingCode[0:scramblingCodeLength]
#
#H=fft(scramblingCodePartial,nbPoints)
##H=fft(psc,nbPoints)
#Y=np.conj(H)*S
#y_fft=ifft(Y)
#
#
##nbPoints_lin=len(signal)+len(scramblingCode)-1
##S_lin=fft(signal,nbPoints_lin)
##H_lin=fft(scramblingCodePartial,nbPoints_lin )
##Y_lin=np.conj(H_lin)*S_lin
##y_fft_lin=ifft(Y_lin)
##y_fft_lin=np.fft.fftshift(y_fft_lin)
#
##h_mf_psc=np.conj(psc[::-1])
##h_mf_ssc=np.conj(ssc[::-1])
#h_mf_psc=np.conj(psc[::-1])
#h_mf_ssc=np.conj(ssc[::-1])
#
#s_psc=filter_FFT(h_mf_psc,signal)
#s_ssc=filter_FFT(h_mf_ssc,signal)
#
#nbChipsInFrame=overSamplingFact*38400
#max_corr=np.zeros(nbChipsInFrame)
#nbDataPerFrame=150
#nbChipsPerDataBit=nbChipsInFrame/nbDataPerFrame;
#
#start=49614
##nbDataBits=(nbSamples-start)/nbChipsInFrame*nbDataPerFrame
##print "number of data bits : {0}".format(nbDataBits)
##nbCompleteFrames=nbDataBits/nbDataPerFrame
##scrambledSignal=signal[start:start+nbCompleteFrames*overSamplingFact*38400]
#
##scramblingCodeComplete=np.zeros(overSamplingFact*nbCompleteFrames*38400,dtype=complex)
#scramblingCodeComplete=np.zeros(len(signal),dtype=complex)
#
##for i in range(overSamplingFact*nbCompleteFrames*38400):
#for i in range(len(signal)):
#    scramblingCodeComplete[i]=scramblingCode[i%(overSamplingFact*38400)]
#
#scramblingCodeComplete=np.roll(scramblingCodeComplete,start)
#    
#descrambledSignal=scramblingCodeComplete.real*signal
#
#
#despreadedData=np.zeros(len(signal)/512,dtype=complex)
##SaveToTxtFile("descrambledSignal.txt",descrambledSignal.real, descrambledSignal.imag)
#
#for i in np.arange(len(signal)/512):
#    despreadedData[i]=np.sum(descrambledSignal[i*nbChipsPerDataBit:(i+1)*nbChipsPerDataBit])
#
#print despreadedData.shape  
#
#Display(abs(y_fft),abs(s_psc),abs(s_ssc))
#############################################
#fig=plt.figure(2)        
#color = get_color()   
#acolor=next(color)
#ax1=plt.subplot(2,1,1)
#plt.plot(despreadedData.real,color=acolor)
#ax1=plt.subplot(2,1,2)
#acolor=next(color)
#plt.plot(despreadedData.imag,color=acolor)
#############################################
