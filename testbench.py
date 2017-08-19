# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 13:33:56 2014

@author: kimfeld
"""
from matplotlib import pyplot as plt
from numpy import *
import scipy.io
import scipy.signal as signal
import Multilateration as ML
import NavigationTest as NT
import copy

#Test Data Generation
normParamLoc=array([0,0])
normParamTime=array([0,0e-6])
normParamSpeed=array([0,0])

#sourceLocationECEF=array([4289469.0,653614.0,4659746.0])#see sender.kml (wädenswil)
sourceLocationECEF=array([ 4338849.05439033,   578965.37502619,  4624078.79392981])#for spirale.txt

nbOfReceiverPoints=7
thresholdNumberOfSigmasOfTOA=[10.0]
thresholdNumberOfSigmasOfFOA=[10.0]
testData=NT.NavigationTest("test_kurve_ueber_thun_2.txt",sourceLocationECEF, normParamLoc, normParamTime,normParamSpeed,thresholdNumberOfSigmasOfTOA,thresholdNumberOfSigmasOfFOA)
testData.WriteTestDataToFile("test_for_chris.txt")

rxData=testData.getReceivers()
rxData2=copy.deepcopy(rxData)
rxData3=copy.deepcopy(rxData)
rxData4=copy.deepcopy(rxData)

#Navigation Algorithm, Location finding
c_speed=testData.c_speed
carrierFrequency=testData.carrierFrequency

#b=ML.Multilateration(rxData,c_speed,carrierFrequency)
#resultSchau=b.ComputeSphericalIntersection()

#c=ML.Multilateration(rxData2,c_speed,carrierFrequency)
#resultSmith=c.ComputeSphericalInterpolation()

N=len(rxData3)-1
Q=identity(2*N)#weighting matrix for [Ho2004]
d=ML.Multilateration(rxData3,c_speed,carrierFrequency)
Result=d.ComputeSphericalInterpolation()

#e=ML.Multilateration(rxData4,c_speed,carrierFrequency)
#result=e.SolveNonLinearFDOAEquations(ResultHo)

