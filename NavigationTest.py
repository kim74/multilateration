# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 13:17:55 2014

@author: kimfeld
"""
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
import scipy.io
import scipy.signal as signal
import TrajectoryGenerator as TG



class NavigationTest:  
    """
    This class defines the test case and generates the test data
    the test data is then fed to the computational core
    """
    def __init__(self,navigationFile,sourceLocationECEF,normParamLoc,normParamTime,normParamSpeed,thresholdNumberOfSigmasOfTOA=10.0,thresholdNumberOfSigmasOfFOA=10.0):
    #thresholdNumberOfSigmasOfTOA=10 and thresholdNumberOfSigmasOfFOA=10.0 as default value to ensure that the weights are practically all at 1
    
        def PrintTestData():
#            print "-------------------------------------------"
#            print "Test Data"
#            print "(Ideal) ECEF-Source Location: {0}".format(sourceLocationECEF)
#            print "Location Normal Params: {0} ---- Time Normal Params: {1} ----- Speed Normal Params: {2}".format(self.normParamLoc,self.normParamTime,self.normParamSpeed)
#            print "Wave Speed: {0} ---- Carrier Frequency: {1}".format(self.c_speed,self.carrierFrequency)
            print "Local Source Coordinate : {0}".format(self.SourceLocationInLocalCoordSystem)        
#            print "Test Data generated"
            
        ##############Simulation Setup
        
        carrierFrequency=2147000000.0
        
        #statistical Properties for location, time and speed (mean, variance)
        self.normParamLoc=normParamLoc 
        self.normParamTime=normParamTime
        self.normParamSpeed=normParamSpeed

        #######Simulation Setup End
        
        
        #the trajectory points refer to a coordinate system that is centered in the last receiver point that was read
        #the translation vector represents therefore the vector from the ECEF's origin (center of the earth) to the last receiver point in the ECEF-based system
        a=TG.TrajectoryGenerator(navigationFile,carrierFrequency,self.normParamLoc,self.normParamTime,self.normParamSpeed)
        x,y,z=a.GetTrajectory() 
        time=a.GetTime()
        vx,vy,vz=a.GetSpeedVectors()
        self.translationVectorECEF=a.GetTranslationVector()
        
        #sourceLocationECEF is the ECEF-based coordinate of the transmitter/source
        #therefore, SourceLocationInLocalCoordSystem is the source location in the coordinate system centered in the first receiver point
        #"local" refers to the coordinate system centered in the last position of the receiver data
        self.SourceLocationInLocalCoordSystem=sourceLocationECEF-self.translationVectorECEF
           
        #thus, toa and foa are computed in the coordinate system centered in the last receiver point
        if normParamTime[0]==0:
            toa, weightTOA=a.ComputeTOAWithRelativeVariation(self.SourceLocationInLocalCoordSystem,thresholdNumberOfSigmasOfTOA)
        else:
            toa, weightTOA=a.ComputeTOAWithAbsoluteVariationAndDrift(self.SourceLocationInLocalCoordSystem,thresholdNumberOfSigmasOfTOA)            
            
        foa, weightFOA=a.ComputeFOA(self.SourceLocationInLocalCoordSystem,thresholdNumberOfSigmasOfFOA)
        self.c_speed=a.speedOfLight   
        self.carrierFrequency=a.carrierFrequency
        
        ##############################################################
        
        PrintTestData()
        
        self.data=[]
        #create the data set from each receiver point
        for i in range(len(x)):
            self.data.append(TG.NavRecord(array([[x[i]],[y[i]],[z[i]]]),array([[vx[i]],[vy[i]],[vz[i]]]),array([toa[i]]),array([foa[i]]),0,weightTOA[i],weightFOA[i],array([time[i]])))
    
    def TranslateFromLocalOriginToFirstPoint(self):
        self.firstPoint=self.data[0].xi
        self.SourceLocationInLocalCoordSystem=self.SourceLocationInLocalCoordSystem-self.firstPoint.reshape((3,))
        for i in range(len(self.data)):
            self.data[i].xi=self.data[i].xi-self.firstPoint

    def TranslateFromFirstPointBackToLocalOrigin(self):
        self.SourceLocationInLocalCoordSystem=self.SourceLocationInLocalCoordSystem+self.firstPoint.reshape((3,))
        for i in range(len(self.data)):
            self.data[i].xi=self.data[i].xi+self.firstPoint       
    
    def GetLocalSourceLocation(self):
        return self.SourceLocationInLocalCoordSystem
        
    def GetTranslationVectorECEFToLocal(self):
        return self.translationVectorECEF
    
    def getReceivers(self):
        return self.data
    
    def printTest(self):
        print "----------TEST DATA----------------"
        for i in self.data:
            print "x : %12.5f , y : %12.5f , z : %12.5f ; vx : %12.5f ; vy: %12.5f ; vz : %12.5f ; tdo : %12.10f , fdo: %14.7f" % (i.xi[0],i.xi[1],i.xi[2],i.v[0],i.v[1],i.v[2],i.toa[0],i.foa[0])      
        print "-----------------------------------"
    
    def WriteTestDataToFile(self,filename):
        #print "x y z vx vy vz tdo fdo"
        f=open(filename,'w')   
        #f.write("{0}\n".format(len(self.data)))
        for d in self.data:
            s="%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f %20.15f %20.10f %20.10f\n" % (d.xi[0],d.xi[1],d.xi[2],d.v[0],d.v[1],d.v[2],d.toa[0],d.foa[0],d.t[0])  
            f.write(s)
        f.close()

    def ReadTestDataFromFile(self,filename): 
        self.data=[]        
        f=open(filename,'r')

        i=0
        for line in f:
            ll=line.strip()
            data=ll.split()
            self.data.append(TG.NavRecord(array([[float(data[0])],[float(data[1])],[float(data[2])]]),array([[float(data[3])],[float(data[4])],[float(data[5])]]),array([float(data[6])]),array([float(data[7])]),0,1,1,array([float(data[8])])))
            i+=1

        f.close()
        #'i' is number of samples read from file

