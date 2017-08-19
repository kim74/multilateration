# -*- coding: utf-8 -*-
"""
Created on Sun Dec 28 16:50:36 2014

@author: kimfeld
"""


from xml.dom.minidom import parse
import xml.dom.minidom
import re
import numpy as np
from scipy.interpolate import interp1d
from scipy import arctan, arctan2, cos, sin, sqrt, pi


class FlightGenerator:
    def __init__(self,filename):
        reg=r'[0-9.]{1,}' #looks for numbers; the expression only contains at least 1 number
        #only number or decimal is allowed

        self.DT=xml.dom.minidom.parse(filename)
        self.collection=self.DT.documentElement
        self.coordinates=self.collection.getElementsByTagName("coordinates")
        a=self.coordinates[0].childNodes[0].data
        
        w=re.findall(reg,a)#returns a list of strings
        self.dataGeodetic=[]
        for i in range(0,len(w),3):
            self.dataGeodetic.append([float(w[i]),float(w[i+1]),float(w[i+2])])
            #converts the strings into floats
            
    def SaveToTxtFile(self,filename,altitude):
        f=open(filename,'w')          
        for i in range(len(self.dataGeodetic)):
            s="{0},{1},{2}\n".format(str(self.dataGeodetic[i][1]),str(self.dataGeodetic[i][0]),str(altitude))#longitude before latitude
            f.write(s)
        f.close()
        
    def ReadFromTxtFile(self,filename):
        f=open(filename,'r')
        self.coordinates_data=[]
        str_data=""
        for line in f:
            line=line.strip()
            data=line.split("\t")
            l=[]
#            for d in data:
#                print "{0}\n".format(d)
#                str_data+=d+','
#                l.append(float(d))
            l.append(float(data[1]))
            l.append(float(data[0]))
            l.append(float(data[2]))
            str_data+=data[1]+","+data[0]+","+data[2]+","
            
            self.coordinates_data.append(l)
            str_data=str_data.strip(',')
            str_data+=" "
        self.coordinate_string=str_data
        
    def WriteDataToXML(self,filename):
        f=open(filename,'w')
        self.coordinates[0].childNodes[0].data=self.coordinate_string
        self.DT.writexml(f)
        f.close()
        
m=FlightGenerator("test_kurve_ueber_thun_2.kml")
m.SaveToTxtFile("test_kurve_ueber_thun_2.txt",1000)
