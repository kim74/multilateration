# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 13:00:25 2015

@author: kimfeld
"""

import numpy as np
from numpy.fft import fft, ifft
from matplotlib import pyplot as plt
import UMTS
import csv


class SCDetector:
    def __init__(self):
        self.samplingFreq=7679975.0
        self.umtsFreq=3840000
        self.overSamplingFact=int(round(self.samplingFreq/self.umtsFreq))
        self.nbChipsPerFrame=38400
        self.nbSamplesPerFrame=self.nbChipsPerFrame*self.overSamplingFact
        self.nbSlotsPerFrame=15
        self.nbChipsPerSlot=self.nbChipsPerFrame/self.nbSlotsPerFrame
        self.nbSamplesPerSlot=self.nbChipsPerSlot*self.overSamplingFact
        self.guardIntervall=self.nbChipsPerSlot/4 #the entire guard intervall spans half of the slot length
    
    def ReadScramblingCodeGroups(self):
        self.sscGroups=np.genfromtxt("ScramblingCodeGroups.txt", delimiter=' ')
        
    def GeneratePSCSSC(self):
        self.a=np.array([1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,-1,-1,1])
        
        #PSC
        self.psc=np.empty((0,),dtype=complex) 
        aModulate=np.array([1,1,1,-1,-1,1,-1,-1,1,1,1,-1,1,-1,1,1])
        
        for k in aModulate:
            if k==1:
                self.psc=np.append(self.psc, self.a)
            else:
                self.psc=np.append(self.psc,-self.a)

        self.psc=(1+1j)*self.psc
        
        #SSC        
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
    
    def OverSampling(self,s):
        aux=np.empty((0,))
        for k in range(self.overSamplingFact):
            aux=np.append(aux,s)
            
        aux=aux.reshape((self.overSamplingFact,len(s)))
        aux=aux.T
        T=aux.reshape((len(s)*self.overSamplingFact,))
        return T    
    
    def PeakFinder(self,data,nbPeaks):
        peaks=np.zeros(nbPeaks)
        for i in range(nbPeaks):
            maxIndex=np.argmax(data)
            peaks[i]=maxIndex
            data[maxIndex-self.guardIntervall:maxIndex+self.guardIntervall]=0
        return peaks
        
    def ConsistencyCheck(self,peaks):
        cleanedPeaks=peaks
        consistencyCounts=np.zeros(len(peaks))
        for i in range(len(peaks)):
            for j in range(len(peaks)):
                diff=abs(peaks[i]-peaks[j])
                a=diff/float(self.nbSamplesPerSlot)
                a_round=round(a)
                if (abs(a-a_round)<0.001):
                    consistencyCounts[i]+=1     
        
        j=0
        #print consistencyCounts
        for v in consistencyCounts:
            if (v<=1):
                cleanedPeaks=np.delete(cleanedPeaks,j)
            else:
                j+=1
                    
        return cleanedPeaks
        
    def PeakDetection(self,data, slope):
        d=np.diff(data)
        #print d
        peak=[]
        peak_detected=0
        for i in range(len(d)):
            
            #if the negative slope is high, peak is detected at the beginning
            if (i==0 and d[i]<-slope):
                peak.append(0)
                peak_detected=0
                
            #if the slope is higher than the limit, peak is detected
            elif (d[i]>slope):
                peak.append(i+1)
                peak_detected=1
                
            #if there is a peak right after a peak, it is detected here;
            #the threshold for discrimination is 0.333*slope
            elif (abs(d[i])<slope/3.0 and peak_detected==1):
                peak.append(i+1)
                peak_detected=1
                
            else:
                peak_detected=0
            
        return peak
        
    def DetectSSCSequences(self, pks, signal, slopeThreshold):
        pksSorted=np.sort(pks)        
        start=pksSorted[0]
        k=range(self.nbSlotsPerFrame)
        v=[[] for x in range(self.nbSlotsPerFrame)]
        z=zip(k,v)
        self.sscSequence={k:v for (k,v) in z}#initialize a dict with 15 keys and [0]
        
        for sscNumber in range(15):#for each SSC-code
            #create SSC-signal with zero-padded over one slot; since only the first 256 elements are non-zero
            sscPadded=np.zeros(self.nbSamplesPerSlot,dtype=complex)
            ssc=self.OverSampling(self.ssc[sscNumber])
            sscPadded[0:len(ssc)]=ssc
            
            signalAligned=signal[start:]#take correlation from synchronized signal; start from the beginning of a slot given by the PSC peak occurence
            sscAligned=self.ExtendSSC(sscPadded,signalAligned)#repeat the above zeropadded signal over the entire acquisition
            
            
            #descramble and integrate over one slot -> result in a
            #'a' contains the correlation value btw ssc-code and signal at the instant of the PSC-peaks
            s=sscAligned*signalAligned
            
            a=np.array([],dtype=complex)
 
            for i in range(len(s)/self.nbSamplesPerSlot-1):     
                s_sum=np.sum(s[i*self.nbSamplesPerSlot:(i+1)*self.nbSamplesPerSlot])
                a=np.append(a,s_sum)
            
            
            a=np.abs(a)
            
            
#            if (sscNumber==4):
#                fig=plt.figure(sscNumber+2)
#                plt.plot(a,'o')
#                plt.show()
          
            b=self.PeakDetection(a,slopeThreshold)
            print "sscNumber: {0} {1}\n".format(sscNumber+1,b)
            #add the peaks found to the sscSequence-dictionary
            i=0
            for bb in b:
                self.sscSequence[bb].append(sscNumber+1)#add 1 to be coherent with the slot numbers starting at 1 in TS25.213

        return self.sscSequence
                    
        
        
    def ExtendSSC(self,ssc,signal):
        signalLength=len(signal)
        nbSSCsInSignal=int(np.ceil(signalLength/float(len(ssc))))
        sscExtended=np.array([])
        for i in range(nbSSCsInSignal):
            sscExtended=np.append(sscExtended,ssc)
        sscExtended=sscExtended[0:signalLength]
        return sscExtended   
        
    def UnpackSequence(self,subSequence):      
        if len(subSequence)==1:
            a=[]
            if subSequence[0]==[]:
                return [[0]] #the missing peak instants are filled with 0
                #so it is "neutral" when computing the cross-correlation
            else:
                for element in subSequence[0]:
                    a.append([element])               
                return a
        else:  
            if subSequence[0]==[]:
                a=[0]
            else:
                a=subSequence[0]
            subSequence.pop(0)
            b=self.UnpackSequence(subSequence)
#            print "a= {0}  b= {1}\n".format(a,b)
            bLength=len(b)           
            c=[x[:] for x in b]
            for i in range(len(a)-1):
                d=[x[:] for x in c]
                b.extend(d)#create new sequences, as many as elements in 'a'    
            
            i=0            
            for aa in a:
                for j in range(bLength):
                    b[j+i*bLength].insert(0,aa)  
                i+=1

            return b
            
    def PotentialSequences(self, subsequence):
        dictlist=[]
        for key,value in subsequence.iteritems():
            dictlist.append(value)
        self.listOfPotentialSSCs=dictlist[0:]
        self.listOfPotentialSSCs=self.UnpackSequence(self.listOfPotentialSSCs)
        return self.listOfPotentialSSCs
    
    def ScramblingCodeGroupIdentification(self):
        matchingMatrix=np.zeros((len(self.listOfPotentialSSCs),len(self.sscGroups)))       
        listOfPotentialSSCs=np.asarray(self.listOfPotentialSSCs)   
        i=0
        #process algorithm for each sequence from the pool of potential sequences
        for sequence in listOfPotentialSSCs:
        #test through the SSC-group table 
            j=0
            for template in self.sscGroups:#
                counts=np.zeros(len(sequence))#matrix containing the number of matches with the sequences aligned at a given position
                for startingPosition in range(len(sequence)):#shift through each starting position and find the number of matches
                    shiftedSequence=np.roll(sequence,startingPosition)
                    for k in range(len(sequence)):
                        if shiftedSequence[k]==template[k]:
                            counts[startingPosition]+=1
                #get the maximum number of matches
                matchingMatrix[i][j]=np.amax(counts)
                j+=1
            i+=1
        maxIndex=np.argmax(matchingMatrix)
        
        print "code group = {0}\n".format(maxIndex%len(self.sscGroups))
        
                
    def SaveToTxtFile(self):
        f=open("ssc_codes_real.txt",'w')       
        for row in range(self.ssc.shape[0]):
            s=""
            for col in range(self.ssc.shape[1]):
                a=np.asscalar(np.real(self.ssc[row,col])) 
                s=s+"{0} ".format(repr(a))
            f.write(s+"\n")
        f.close()  

        f=open("ssc_codes_imag.txt",'w')       
        for row in range(self.ssc.shape[0]):
            s=""
            for col in range(self.ssc.shape[1]):
                a=np.asscalar(np.imag(self.ssc[row,col])) 
                s=s+"{0} ".format(repr(a))
            f.write(s+"\n")
        f.close()

        f=open("psc_codes_real.txt",'w')       
        s=""
        for col in range(len(self.psc)):
            a=np.asscalar(np.real(self.psc[col])) 
            s=s+"{0} ".format(repr(a))
        f.write(s+"\n")
        f.close()      

        f=open("psc_codes_imag.txt",'w')       
        s=""
        for col in range(len(self.psc)):
            a=np.asscalar(np.imag(self.psc[col])) 
            s=s+"{0} ".format(repr(a))
        f.write(s+"\n")
        f.close()               
                
        

nbSamples=76800
umtsFreq=3840000
samplingFreq=7679975.0
nFrames=5
codeNumber=400#6064

overSamplingFact=int(round(samplingFreq/umtsFreq)) #for the emulation of an analog signal

umts=UMTS.CUMTSCell(overSamplingFact,codeNumber)
signal=umts.Frame()

time,ssc=umts.OverSampling(umts.ssc.ssc[3],overSamplingFact)
time,psc=umts.OverSampling(umts.psc.psc,overSamplingFact)

time,scramblingCode=umts.OverSampling(umts.scramblingCode,overSamplingFact)
time,signal=umts.StreamFromFile("C:\\Users\\kimfeld\\Documents\\Projects\\Multilateration\\software\\VC++\\UMTSFrameExtractor\\UMTSFrameExtractor\\exp1\\acquisition_IQ_2147.599900MHz.txt",nbSamples,samplingFreq)

nbPoints=len(signal)

start=73

S=fft(signal,nbPoints)
scramblingCodeLength=len(scramblingCode)

scramblingCodePartial=scramblingCode[0:scramblingCodeLength]

H=fft(scramblingCodePartial,nbPoints)

Y=np.conj(H)*S
y_fft=ifft(Y)

PSC=fft(psc,nbPoints)
S_PSC=np.conj(PSC)*S
s_psc=ifft(S_PSC)

SSC=fft(ssc,nbPoints)
S_SSC=np.conj(SSC)*S
s_ssc=ifft(S_SSC)

#UMTS.Display(abs(y_fft),abs(s_psc),abs(s_ssc))

scDetector=SCDetector()
scDetector.ReadScramblingCodeGroups()
scDetector.GeneratePSCSSC()
pks=scDetector.PeakFinder(abs(s_psc),15)
#pks=np.array([ 40945. , 46065.,  56305.,  35825.,  20465.,  51185.,  71665.,  66544.,  61424.,  5105.,  15344.,  30705.,  10225.,  25585.])
pks_checked=scDetector.ConsistencyCheck(pks)
slopeThreshold=500000
subsequence=scDetector.DetectSSCSequences(pks_checked,signal,slopeThreshold)

#subsequence={0: [12], 1: [3], 2: [0], 3: [5,8,9], 4: [8], 5: [0], 6: [0], 7: [14], 8: [12,9], 9: [9,13], 10: [8], 11: [9], 12: [14], 13: [], 14: [2]}
subsequence={0: [13], 1: [2], 2: [9], 3: [10], 4: [12], 5: [16], 6: [8], 7: [5], 8: [3], 9: [15], 10: [6], 11: [1], 12: [11], 13: [14], 14: [4]}
#print subsequence
sequenceList=scDetector.PotentialSequences(subsequence)
print sequenceList
#print scDetector.sscGroups[8]

sequenceList=scDetector.ScramblingCodeGroupIdentification()
scDetector.SaveToTxtFile()



