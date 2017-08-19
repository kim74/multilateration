# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 16:54:43 2014

@author: kimfeld
"""
from xml.dom.minidom import parse
import xml.dom.minidom
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy import arctan, arctan2, cos, sin, sqrt, pi
from math import pow, degrees, radians

class CoordinatesReader:
    
    def __init__(self,speed,measPerSecond):
        self.speed=speed
        self.measurementsPerSecond=measPerSecond
     
    
    def ReadFile(self,filename,switch):
        if (switch==1):
            reg=r'[0-9.]{1,}' #looks for numbers; the expression only contains at least 1 number
            #only number or decimal is allowed
            self.ECEFcalculated=0
            DT=xml.dom.minidom.parse(filename)
            collection=DT.documentElement
            coordinates=collection.getElementsByTagName("coordinates")
            a=coordinates[0].childNodes[0].data       
            w=re.findall(reg,a)#returns a list of strings
            self.dataGeodetic=[]
            for i in range(0,len(w),3):
                self.dataGeodetic.append([float(w[i]),float(w[i+1]),float(w[i+2])])
                #converts the strings into floats    
        else:
            f=open(filename,'r') 
            self.dataGeodetic=[]
            for line in f:
                ll=line.strip()
                data=ll.split("\t")
                self.dataGeodetic.append([float(data[1]),float(data[0]),float(data[2])])#take into same format as .kml files (longitude, latitude, altitude)  
            f.close()
        
    def getWSG84Coordinates(self):
        return np.asarray(self.dataGeodetic)
               
    def GeodeticToECEF(self):
        #WGS84
        self.ECEFcalculated=1
        semiMajorAxis=6378137.0 
        firstEccentricitySquared=6.69437999014e-3
        
        self.dataECEF=[]
        for coord in self.dataGeodetic:
            longitude=2*pi/360.0*coord[0]
            latitude=2*pi/360.0*coord[1]
            altitude=coord[2]
            N=semiMajorAxis/sqrt(1-firstEccentricitySquared*(sin(latitude))**2)
            X=(N+altitude)*cos(latitude)*cos(longitude)
            Y=(N+altitude)*cos(latitude)*sin(longitude)
            Z=(N*(1-firstEccentricitySquared)+altitude)*sin(latitude)
            self.dataECEF.append([X,Y,Z])
        self.dataECEF=np.asarray(self.dataECEF)
        return self.dataECEF
        
    def ECEFToGeodetic(self,x, y, z):
        """Convert ECEF coordinates to geodetic.
        J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates \
        to geodetic coordinates," IEEE Transactions on Aerospace and \
        Electronic Systems, vol. 30, pp. 957-961, 1994."""
        def cbrt(x):
            if x >= 0:
                return pow(x, 1.0/3.0)
            else:
                return -pow(abs(x), 1.0/3.0)
        
        a = 6378137.0#6378.137
        b = 6356752.3142#6356.7523142
        esq = 6.69437999014 * 0.001
        e1sq = 6.73949674228 * 0.001
        f = 1 / 298.257223563   

        r = sqrt(x * x + y * y)
        Esq = a * a - b * b
        F = 54 * b * b * z * z
        G = r * r + (1 - esq) * z * z - esq * Esq
        C = (esq * esq * F * r * r) / (pow(G, 3))
        S = cbrt(1 + C + sqrt(C * C + 2 * C))
        P = F / (3 * pow((S + 1 / S + 1), 2) * G * G)
        Q = sqrt(1 + 2 * esq * esq * P)
        r_0 =  -(P * esq * r) / (1 + Q) + sqrt(0.5 * a * a*(1 + 1.0 / Q) - \
            P * (1 - esq) * z * z / (Q * (1 + Q)) - 0.5 * P * r * r)
        U = sqrt(pow((r - esq * r_0), 2) + z * z)
        V = sqrt(pow((r - esq * r_0), 2) + (1 - esq) * z * z)
        Z_0 = b * b * z / (a * V)
        h = U * (1 - b * b / (a * V))
        lat = arctan((z + e1sq * Z_0) / r)
        lon = arctan2(y, x)
        return degrees(lat), degrees(lon)
      
    def CreateMeasurementPoints(self):
        #returns also the translationVector that is needed to shift the origin from the earth-center to the 
        #first point in the set of receiver points
        if (self.ECEFcalculated!=1):
            print "convert data to ECEF first"
        else:
            totalDistance=0
            t=np.zeros(len(self.dataECEF))
            time=0
            for i in range(len(self.dataECEF)-1):
                distance=np.linalg.norm(self.dataECEF[i+1]-self.dataECEF[i])
                t[i]=time
                totalDistance+=distance                
                time+=distance/float(self.speed)
            t[len(self.dataECEF)-1]=time    
                
            t=np.asarray(t)
            self.dataECEF=np.insert(self.dataECEF,0,t,axis=1)#insert the time vector as the first column
            
            #compute number of Measurement Points as a function of the total length of the path and the number of measurements per second

            measurementTime=t
            l=len(measurementTime)
            measurementTime=measurementTime.reshape((l,1))
            xnew=self.dataECEF[:,1]
            ynew=self.dataECEF[:,2]
            znew=self.dataECEF[:,3]
            xnew=xnew.reshape((l,1))      
            ynew=ynew.reshape((l,1))
            znew=znew.reshape((l,1))
            
            self.dataMeasurementPathECEF=np.hstack((measurementTime,xnew,ynew,znew))
            
            M=len(self.dataMeasurementPathECEF)
            ###################################################################
            #TAKE THE FIRST POINT as the  origin of the local coordinate system
            ###################################################################
            translationVectorToLocalCoordinates=np.array(self.dataMeasurementPathECEF[0,1:4])     
            voffset=np.array([translationVectorToLocalCoordinates])
            for i in range(M):
                self.dataMeasurementPathECEF[i,1:4]=self.dataMeasurementPathECEF[i,1:4]-voffset 
                #print "time {0}".format(self.dataMeasurementPathECEF[i,0])
            
            
            return self.dataMeasurementPathECEF, translationVectorToLocalCoordinates
    
    def CreateSpeedVectors(self):
        self.speedMeasurementPath=np.empty((0,4))
        M=len(self.dataMeasurementPathECEF)
        for i in range(M-1):
            ds=self.dataMeasurementPathECEF[i+1,1:4]-self.dataMeasurementPathECEF[i,1:4]#segment from point i to point i+1
            v=ds/np.linalg.norm(ds)*self.speed #speed vector at point i
            v=np.insert(v,0,[self.dataMeasurementPathECEF[i,0]],axis=0)#insert time into the current speed vector at position 0
            self.speedMeasurementPath=np.append(self.speedMeasurementPath,[v],axis=0)
        
        #add the last speed vector from point M-1 to M as the speed vector at point M
        #to generate as many speed vectors as points
        a=self.speedMeasurementPath[M-2,1:4]
        a=np.insert(a,0,self.dataMeasurementPathECEF[M-1,0],axis=0)
        self.speedMeasurementPath=np.append(self.speedMeasurementPath,[a],axis=0)
        return self.speedMeasurementPath
                   
    def SaveToFile(self,filename,altitude):
        f=open(filename,'w')          
        for i in range(len(self.dataMeasurementPathGeodetic)):
            s="{0},{1},{2} ".format(str(self.dataMeasurementPathGeodetic[i,2]),str(self.dataMeasurementPathGeodetic[i,1]),str(altitude))
            f.write(s)
        f.close()

