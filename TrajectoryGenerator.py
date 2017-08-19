# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 10:03:09 2014

@author: Kili
"""

from matplotlib import pyplot as plt
import numpy as np
import scipy as sc
import scipy.io
from CoordinatesReader import * 
from numpy.random import normal
import re


class NavRecord:
    def __init__(self, Coord, Speed, TOA, FOA, toff, weightTOA=1, weightFOA=1, time=0):
        if (len(Coord)==3):
            self.xi=Coord
            self.v=Speed
            self.toa=TOA #time of arrival in Microseconds us
            self.foa=FOA #frequency of arrival
            self.toff=toff #time offset of the source clock with respect to the receiver clock
            self.weightTOA=weightTOA
            self.weightFOA=weightFOA
            self.t=time
        else:
            print "the variable 'Coord' does not have the right dimension"

class TrajectoryGenerator:
            
    def __init__(self,navigationFile,carrierFrequency,normParamLoc,normParamTime,normParamSpeed):

        #######################################################################
              
        self.x=np.array([])
        self.y=np.array([])
        self.z=np.array([])
        self.vx=np.array([])
        self.vy=np.array([])
        self.vz=np.array([])

        self.speedOfLight=300000000.0
        self.carrierFrequency=carrierFrequency
        self.toa=np.array([])
        self.foa=np.array([])
        self.weightTOA=np.array([])
        self.weightFOA=np.array([])
        self.normParamLoc=normParamLoc
        self.normParamTime=normParamTime
        self.normParamSpeed=normParamSpeed

        #take external data from CoordinatesReader.py
        nbMeasurementsPerSecond=0.5    
        speedOfReceiver=50     


        a=CoordinatesReader(speedOfReceiver,nbMeasurementsPerSecond)
        #a.ReadFile("testflight11_withAltitude.kml",1)
        a.ReadFile(navigationFile,self.CheckFileExtension(navigationFile))
        d=a.GeodeticToECEF()
        (self.position, self.translationVector)=a.CreateMeasurementPoints()#the positions refer to the origin which is the first point of the data set
        #print "Translation Vector : {0}".format(self.translationVector)        
        speeds=a.CreateSpeedVectors()
        self.time=self.position[:,0]
        
        for i in range(len(self.time)):#Positions of all receivers in the local coordinate system, where (0,0,0) is the FIRST receiver point
            self.x=np.append(self.x,self.position[i,1])
            self.y=np.append(self.y,self.position[i,2])
            self.z=np.append(self.z,self.position[i,3])
            self.vx=np.append(self.vx,speeds[i,1])
            self.vy=np.append(self.vy,speeds[i,2])
            self.vz=np.append(self.vz,speeds[i,3]) 
           
    def CheckFileExtension(self,filename):
        reg=r'\.([a-zA-Z]{3})\Z' #captures the extension and the delimiter dot
        w=re.findall(reg,filename)#returns a list of strings
        if (w[0]=="txt"):
            return 0
        elif (w[0]=="kml"):
            return 1
        else:
            return "wrong file type"
                   
    def GetTranslationVector(self):
        return self.translationVector
    
    def GetTrajectory(self):
        return (self.x,self.y,self.z)
        
    def GetTime(self):
        return self.time
    
    def GetSpeedVectors(self):
        return (self.vx,self.vy,self.vz)
        
    def ComputeWeightTOA(self,e,factor,sigma):
        
        e_normalized=abs(e/sigma)
        #print "{0} {1} {2}\n".format(e_normalized,factor,sigma)
        if e_normalized <= factor:
            weight=1.0
        elif e_normalized<=2*factor:
            weight=0.5
        else: 
            weight=0.1
        return weight
        
    def ComputeWeightFOA(self,e,factor,sigma):
        e_normalized=abs(e/sigma)
        if e_normalized <= factor:
            weight=1.0
        elif e_normalized<=2*factor:
            weight=0.5
        else: 
            weight=0.1
        return weight
        
    def ComputeTOAWithRelativeVariation(self,source, fact):       
        if (self.normParamTime[1]==0 and self.normParamLoc[1]==0):
            #print "no time variation, no location variation"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])
                self.toa=np.append(self.toa, np.linalg.norm(source-rx_vec)/self.speedOfLight)
                self.weightTOA=np.append(self.weightTOA,1)
        elif self.normParamTime[1]!=0 and self.normParamLoc[1]==0:
            #print "time variation only"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])
                e_time=np.asscalar(normal(self.normParamTime[0],self.normParamTime[1],1)) 
                weight=self.ComputeWeightTOA(e_time,fact,self.normParamTime[1]) 
                trueTOA=np.linalg.norm(source-rx_vec)/self.speedOfLight
                TOA=trueTOA*(1+e_time)
                self.toa=np.append(self.toa,TOA )
                #print "{0}, {1}, {2}".format(trueTOA,TOA,e_time)
                self.weightTOA=np.append(self.weightTOA,weight)
        elif self.normParamLoc[1]!=0 and self.normParamTime[1]==0:
            #print "location variation only"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])+normal(self.normParamLoc[0],self.normParamLoc[1],3)             
                self.toa=np.append(self.toa, np.linalg.norm(source-rx_vec)/self.speedOfLight) 
                self.weightTOA=np.append(self.weightTOA,1)
        else:
            #print "both time and location variation"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])+normal(self.normParamLoc[0],self.normParamLoc[1],3)
                #add time variation contained in normParamTime
                e_time=np.asscalar(normal(self.normParamTime[0],self.normParamTime[1],1))
                weight=self.ComputeWeightTOA(e_time,fact,self.normParamTime[1])
                self.toa=np.append(self.toa, (np.linalg.norm(source-rx_vec)/self.speedOfLight)*(1+e_time)) 
                self.weightTOA=np.append(self.weightTOA,weight)
        return self.toa, self.weightTOA      
        
    def ComputeTOAWithAbsoluteVariationAndDrift(self,source, fact):       
        if (self.normParamTime[1]==0 and self.normParamLoc[1]==0):
            #print "no time variation, no location variation"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])
                self.toa=np.append(self.toa, np.linalg.norm(source-rx_vec)/self.speedOfLight)
                self.weightTOA=np.append(self.weightTOA,1)
        elif self.normParamTime[1]!=0 and self.normParamLoc[1]==0:
            #print "time variation only"
            print "parameters of normal distribution: {0}  {1}\n".format(self.normParamTime[0],self.normParamTime[1])
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])          
                e_time=np.asscalar(normal(self.normParamTime[0]*self.time[i],self.normParamTime[1],1))   #includes drift component             
                weight=self.ComputeWeightTOA(e_time,fact[0],self.normParamTime[1]) 
                trueTOA=np.linalg.norm(source-rx_vec)/self.speedOfLight                
                TOA=trueTOA+e_time
                #print "true toa: {0}   error toa: {1}\n".format(trueTOA,TOA)
                self.toa=np.append(self.toa,TOA )
                #print "{0}, {1}, {2}".format(trueTOA,TOA,e_time)
                self.weightTOA=np.append(self.weightTOA,weight)
        elif self.normParamLoc[1]!=0 and self.normParamTime[1]==0:
            #print "location variation only"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])+normal(self.normParamLoc[0],self.normParamLoc[1],3)             
                self.toa=np.append(self.toa, np.linalg.norm(source-rx_vec)/self.speedOfLight) 
                self.weightTOA=np.append(self.weightTOA,1)
        else:
            #print "both time and location variation"
            for i in range(len(self.x)):
                rx_vec=np.array([self.x[i],self.y[i],self.z[i]])+normal(self.normParamLoc[0],self.normParamLoc[1],3)
                #add time variation contained in normParamTime
                e_time=np.asscalar(normal(self.normParamTime[0],self.normParamTime[1],1))
                weight=self.ComputeWeightTOA(e_time,fact,self.normParamTime[1])
                self.toa=np.append(self.toa, (np.linalg.norm(source-rx_vec)/self.speedOfLight)*(1+e_time)) 
                self.weightTOA=np.append(self.weightTOA,weight)
        return self.toa, self.weightTOA   
        
    def ComputeFOA(self,source,fact):
        if (self.normParamSpeed[1]==0):
            for i in range(len(self.x)):
                vx_vec=np.array([[self.vx[i]],[self.vy[i]],[self.vz[i]]])#speed vector of receiver i
                rx_vec=np.array([[source[0]-self.x[i]],[source[1]-self.y[i]],[source[2]-self.z[i]]])#vector from receiver i to source
                vRadial=np.dot(np.transpose(rx_vec)/np.linalg.norm(rx_vec),vx_vec) #compute radial speed vector between the source and the receiver
                #using the scalar-product (gives the component from one vector projected on the other vector)   
                freqAtReceiver=(self.speedOfLight+vRadial)/self.speedOfLight*self.carrierFrequency
                self.foa=np.append(self.foa,freqAtReceiver)
                self.weightFOA=np.append(self.weightFOA,1)
        else:
            for i in range(len(self.x)):
                vx_vec=np.array([[self.vx[i]],[self.vy[i]],[self.vz[i]]])#speed vector of receiver i  
                rx_vec=np.array([[source[0]-self.x[i]],[source[1]-self.y[i]],[source[2]-self.z[i]]])#vector from receiver i to source
                e_speed=np.asscalar(normal(self.normParamSpeed[0],self.normParamSpeed[1],1))
                vRadial=(1+e_speed)*np.dot(np.transpose(rx_vec)/np.linalg.norm(rx_vec),vx_vec) #compute radial speed vector between the source and the receiver
                #using the scalar-product (gives the component from one vector projected on the other vector)                     
                freqAtReceiver=(self.speedOfLight+vRadial)/self.speedOfLight*self.carrierFrequency
                weight=self.ComputeWeightFOA(e_speed,fact,self.normParamSpeed[1])
                self.foa=np.append(self.foa,freqAtReceiver)
                self.weightFOA=np.append(self.weightFOA,weight)
            
        return self.foa, self.weightFOA
            
        
