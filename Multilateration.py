# -*- coding: utf-8 -*-
"""
Created on Mon Dec 08 14:17:22 2014

@author: Kili
"""

from matplotlib import pyplot as plt
from numpy import *
from numpy import transpose, dot
from numpy.linalg import inv
import scipy as sc
import scipy.io
import scipy.signal as signal
import TrajectoryGenerator as TG
import copy
   
class Multilateration:
    """
    This class implements the computational core, which consists of
    the geolocation position estimation algorithm.
    
    It receives the Test Data consisting of receiver locations and
    computes an estimate of the source position.
    
    The source position is available as a numpy array [x,y,z]
    """
    def __init__(self,rxData,c_speed,carrierFrequency):
        #take the receiver data into a convenient format for
        #the numerical computation
        self.receiver=rxData
        self.c_speed=c_speed  #m/us speed of light  
        self.c=c_speed
        self.carrierFrequency=carrierFrequency
        self.translationVector=TG.NavRecord(array([[0],[0],[0]]),array([[0],[0],[0]]),array([0]),array([0]),0)
        
    def FormatData(self,nbPoints):
        #creates the different matrices needed in the least squares computation
        #and shifts the receiver/source location to a new coordinate system with the last receiver point as origin (0,0,0)
        R=array([])        
        d=array([])
        M= empty((0,3))
        DELTA=[]        
        

#        for i in range(0,nbPoints):
#            print "toa: {0}\n".format(self.receiver[i].toa)            
            
        for i in range(0,nbPoints-1):
            self.receiver[i+1].xi=self.receiver[i+1].xi-self.receiver[0].xi #shift coordinate system origin into last receiver
            r= linalg.norm(self.receiver[i+1].xi)   
            R= append(R,r)
            dd=((self.receiver[i+1].toa-self.receiver[i+1].toff)-(self.receiver[0].toa-self.receiver[0].toff))*self.c_speed
            print "{0}-range: {1}".format(i,dd)
            d= append(d,dd)
            #print "vector: [%10.5f %10.5f %10.5f]  ------ Distance to last point: %10.5f ------ range difference wrt to last receiver: %10.5f" %(self.receiver[i].xi[0],self.receiver[i].xi[1],self.receiver[i].xi[2],R[i],d[i])
            M= append(M, transpose(self.receiver[i+1].xi),axis=0)
        
        R=R.reshape(nbPoints-1,1)
        d=d.reshape(nbPoints-1,1)
        DELTA= multiply(R,R)- multiply(d,d)
  
        return d,M,DELTA,R
        
    def ComputeSphericalInterpolation(self):
        print "------------"
        print "Compute data according to [Smith1987]"
        print "------------"
        #define the weighting matrices W and V as identity matrices
        nbDataPoints=len(self.receiver)
        #print "number of receivers : {0}".format(nbDataPoints)
        V=identity(nbDataPoints-1)
        W=V
        for i in range(nbDataPoints-1):
            W[i,i]=self.receiver[i+1].weightTOA
            print "{0}.point : {1}\n".format(i,self.receiver[i+1].xi)
               
        d,S,DELTA,R=self.FormatData(nbDataPoints)
        self.translationVector.xi=self.receiver[0].xi #store the translation vector of the new coordinate system
        
        #self.translationVector.xi=self.receiver[0].xi
        self.receiver[0].xi=array([0,0,0]) #origin now new at the last receiver point
        #self.receiver[0].xi=array([0,0,0])        
        
        H=dot( linalg.inv(dot( transpose(S),dot(W,S))),dot( transpose(S),W))
        
        Ps=dot(S,H)
        I= identity(nbDataPoints-1)
        Ps_v=I-Ps

#        print d
#        print "\n"
#        print Ps_v
#        print "\n"
        #print V
        #self.print_data(DELTA,M,d,R)
        
        num=dot(d,dot( transpose(d),dot( transpose(Ps_v),dot(V,Ps_v))))
        den=dot( transpose(d),dot( transpose(Ps_v),dot(V,dot(Ps_v,d))))
        
#        print "\n"
#        print num
        
        xs_est=0.5*dot(H,dot(I-num/den,DELTA))
       
        num=dot( transpose(d),dot( transpose(Ps_v),dot(V,Ps_v))) 
        den=2*dot( transpose(d),dot( transpose(Ps_v),dot(V,dot(Ps_v,d))))
 
        Rs=dot(num,DELTA)/den
        xs_est_2=0.5*dot(H,(DELTA-2*Rs*d))
        
        #self.print_data(DELTA,S,d,R)        
        print "--------------------------------"
        #print "the solution is : {0}".format( (xs_est+self.translationVector.xi).T)
        print "the solution is : {0}".format(xs_est.T)

        print "--------------------------------"
        #print "the alternative solution is : {0}".format( transpose(xs_est_2+self.translationVector.xi))
        return  xs_est+self.translationVector.xi       
        
    def print_data(self,DELTA,M,d,R):
        print "Matrix DELTA : {0}".format(DELTA)
        print "Matrix M : {0}".format(M)
        print "Matrix d : {0}".format(d)
        print "Matrix R : {0}".format(R)
            
    def ComputeSphericalIntersection(self):
          
        print "------------"
        print "Compute data according to [Schau1987]"
        print "------------"        
        
        if len(self.receiver)<4:
            print "Not enough data points available."
            return
        if len(self.receiver)>4: #only take the first 4 sensors in case there are more
            print "More data points than necessary. We only take the first four points."
        
#        #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
#        #2) compute R1, R2 and R3, as well as d14, d24, d34
#        #3) create matrices M, DELTA und d for the algorithm

        d,M,DELTA,R=self.FormatData(4)       
        
        self.translationVector.xi=self.receiver[3].xi #store the translation vector of the new coordinate system
        self.receiver[3].xi=array([[0],[0],[0]]) #origin now new at receiver 4
        
        #self.print_data(DELTA,M,d,R)
 
        #4) inverse M, compute a, b and c as the coefficients for the quadratic equation
        M_inv=inv(M)

        W=dot(transpose(M_inv),M_inv)
        a=4-4*dot(dot(transpose(d),W),d) #4-d^T*M^-1)^T*M^-1*d
        b=2*dot(transpose(d),dot(W,DELTA))+2*dot(transpose(DELTA),dot(W,d))
        c=-dot(transpose(DELTA),dot(W,DELTA))       

        #5) compute Rs from the quadratic equation --> 2 possible solutions
        Rs=roots(array([a[0][0],b[0][0],c[0][0]]))
        #print "the solutions are {0}".format(Rs)
            
            #6) compute xs from equation (10)  (0.5*M^-1*(DELTA-2*Rs*d))
            
            #7) store target coordinates for the current 4 points

        x1=0.5*dot(M_inv,DELTA-2*Rs[0]*d)
        x2=0.5*dot(M_inv,DELTA-2*Rs[1]*d)
        print ""
        print "solution 1: x_target = {0}".format( transpose(x1+self.translationVector.xi))
        print "solution 2: x_target = {0}".format( transpose(x2+self.translationVector.xi))
        return (transpose(x1+self.translationVector.xi),transpose(x2+self.translationVector.xi))

    def TDOAandFDOA_Ho2004(self,Q):
        print "------------"
        print "Compute data according to [Ho2004]"
        print "------------"
        M=len(self.receiver)
        #print "number of Receiver points: {0}".format(M) 
        #print "number of receivers : {0}".format(M)
        h1_loc=array([])     
        h1_speed=array([])
        
        self.translationVector.xi=self.receiver[0].xi #store the translation vector of the new coordinate system (=the "last" receiver location)
        self.translationVector.v=self.receiver[0].v
        
        s_M=array([[0],[0],[0]])  #the last receiver (index M-1) is the new origin/reference point   
        s_dotM=self.receiver[0].v      
        s_MmalM=s_M.T.dot(s_M)# compute s^T.S
        s_dotMmalM=s_dotM.T.dot(s_M)
     
        #make vector h1 (from equation (9) in [Ho2004])
        r_iM=array([])
        r_dotiM=array([])
        for i in range(M-1):
            self.receiver[i+1].xi=self.receiver[i+1].xi-self.receiver[0].xi #shift coordinate system origin into last receiver (index M-1)
            
            r_iM= append(r_iM,((self.receiver[i+1].toa-self.receiver[i+1].toff)-(self.receiver[0].toa-self.receiver[0].toff))*self.c_speed) #range from i-th sensor to the last sensor, which is the new origin
            s_i=self.receiver[i+1].xi#i-th receiver location
 
            h_i=r_iM[i]**2-s_i.T.dot(s_i)+s_MmalM      
            h1_loc= append(h1_loc,h_i)   
            r_dotiM= append(r_dotiM,self.c_speed/self.carrierFrequency*(-self.receiver[i+1].foa+self.receiver[0].foa))
            h_i=2*(r_dotiM[i]*r_iM[i]-(self.receiver[i+1].v).T.dot(self.receiver[i+1].xi)+s_dotMmalM)
            h1_speed= append(h1_speed,h_i)       
            print self.receiver[i].xi
        
        print self.receiver[M-1].xi
        self.receiver[0].xi=array([[0],[0],[0]]) #origin now new at the last receiver point 
        
        h1_loc=h1_loc.reshape((M-1,1))
        h1_speed=h1_speed.reshape((M-1,1))
        h1= vstack((h1_loc,h1_speed)) #
        #print "h1: {0}".format(h1)
        ##########################################
        
        #make matrix G1 (from equation (9) in [Ho2004])
        G11= empty((0,3))
        G22= empty((0,3))
        r11= empty((0,1))
        r22= empty((0,1))
        for i in range(M-1):
            G11= append(G11, transpose(self.receiver[i+1].xi-s_M),axis=0)
            r11= append(r11,[[r_iM[i]]],axis=0)                    
        for i in range(M-1):
            G11= append(G11, transpose(self.receiver[i+1].v-s_dotM),axis=0)
            r11= append(r11,[[r_dotiM[i]]],axis=0)
                      
        G22= zeros((M-1,3))
        for i in range(M-1):
            G22= append(G22, transpose(self.receiver[i+1].xi-self.receiver[0].xi),axis=0)
        r22= zeros((M-1,1))
        for i in range(M-1):
            r22= append(r22,[[r_iM[i]]],axis=0)
        
        G=-2* hstack((G11,r11,G22,r22))
        
        ############################################
        W1=identity(2*(M-1))#inv(Q)
        for i in range(M-1):
            W1[i,i]=self.receiver[i+1].weightTOA
            W1[i+(M-1),i+(M-1)]=self.receiver[i+1].weightFOA
                
        
        theta1= linalg.inv(G.T.dot(W1.dot(G))).dot(G.T.dot(W1)).dot(h1)
       
        #make h2 (from equation (17) in [Ho2004])
        theta1_u=theta1[0:3]
        theta1_dotu=theta1[4:7]
        h21=(theta1_u-self.receiver[0].xi)**2 #element-wise squaring
        h22=(theta1_dotu-self.receiver[0].v)*(theta1_u-self.receiver[0].xi)  #element-wise multiply
        h2= vstack((h21,[theta1[3]**2],h22,[theta1[3]*theta1[7]]))
        
        #make G2 (from equation (17) in [Ho2004])
        G21= vstack(( identity(3),array([[1,1,1]]), zeros((3,3)),array([[0,0,0]])))
        G22= vstack(( zeros((3,3)),array([[0,0,0]]), identity(3),array([[1,1,1]])))
        G2= hstack((G21,G22))

        W2= identity(8)
        theta2= inv(G2.T.dot(W2.dot(G2))).dot(G2.T.dot(W2)).dot(h2)
        
        signum= sign(theta1[0:3]-self.translationVector.xi)
        #print "signum: {0}".format(signum)
        signum[2]=signum[2]#sometimes signum[2]=-signum[2] to find the alternative solution
        xs_est=signum*( sqrt( absolute(theta2[0:3,0]).reshape(3,1)))+self.translationVector.xi
        xsdot_est=signum*(theta2[3:6]/sqrt(theta2[0:3])).reshape(3,1) +self.translationVector.v
        
        print "--------------------------------"
        print "the location solution is : {0}".format( xs_est.T)
        #print "the speed solution is :{0}".format(xsdot_est.T)
        print "--------------------------------"        
        return xs_est
    
    def SolveNonLinearFDOAEquations(self,initial_guess):
        
        def equations(s,v,fdoa,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                
                func_0=dot(-(s[0,:]-u)/linalg.norm(s[0,:]-u),v[0,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_0=f0/c*func_0-(fdoa[0]-fdoa[3])
                
                func_1=dot(-(s[1,:]-u)/linalg.norm(s[1,:]-u),v[1,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_1=f0/c*func_1-(fdoa[1]-fdoa[3])
                
                func_2=dot(-(s[2,:]-u)/linalg.norm(s[2,:]-u),v[2,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_2=f0/c*func_2-(fdoa[2]-fdoa[3])
                
                return (func_0,func_1,func_2)
            
            f0=self.carrierFrequency
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1 
                    
#        print "------------"
#        print "Compute data by solving the nonlinear FDOA equations numerically with an initial guess"
#        print "------------"        
        
        nbRx= len(self.receiver)       
        sourceLocations=empty((0,3))       
        

        error=1000
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. We only take the first four points."
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            
            #result=initial_guess
            delta=int(floor(nbRx/4.0))
            for j in range(delta):
            #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                rx=empty((0,3))
                translationVector=self.receiver[j+3*delta].xi #store the translation vector of the new coordinate system              
                for k in range(3):                    
                    rx=append(rx,(self.receiver[j+k*delta].xi-translationVector).T,axis=0)
                rx=append(rx,array([[0,0,0]]),axis=0) #origin now new at receiver 4 of the current subgroup
                                                               
                speeds=empty((0,3))
                fdoas=empty((0,))
                
                for i in range(4):   
                    speeds=append(speeds,self.receiver[j+i*delta].v.T,axis=0)
                    fdoas=append(fdoas,self.receiver[j+i*delta].foa)
                
                result=equations(rx,speeds,fdoas,initial_guess)+translationVector.reshape((3,))
                
                e=linalg.norm(result-initial_guess.flatten())
            
                if e<error:
                    sourceLocations=append(sourceLocations,result.reshape((1,3)),axis=0)
                    #print "{0}. run : {1}".format(j,result)
                
#        print "--------------------------------"
#        print "the average position is: {0}".format(mean(sourceLocationAvg,axis=0))
        #print "the alternative solution is : {0}".format( transpose(xs_est_2+self.translationVector.xi))
        #print "\n\n{0}".format(mean(sourceLocationAvg,axis=0))
        return sourceLocations

    def SolveNonLinearTDOAEquations(self,initial_guess):
        
        def equations(s,toas,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                r1=c*toas[3]               

                func_0=dot(s[0,:],s[0,:])-dot(s[3,:],s[3,:])-2*dot((s[0,:]-s[3,:]),u)-(c*(toas[0]-toas[3]))**2-2*r1*c*(toas[0]-toas[3])            
                
                func_1=dot(s[1,:],s[1,:])-dot(s[3,:],s[3,:])-2*dot((s[1,:]-s[3,:]),u)-(c*(toas[1]-toas[3]))**2-2*r1*c*(toas[1]-toas[3])
                
                func_2=dot(s[2,:],s[2,:])-dot(s[3,:],s[3,:])-2*dot((s[2,:]-s[3,:]),u)-(c*(toas[2]-toas[3]))**2-2*r1*c*(toas[2]-toas[3])
                
                return (func_0,func_1,func_2)
            
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1       

#        print "------------"
#        print "Compute data by solving the nonlinear TDOA equations numerically with an initial guess"
#        print "------------"        
        
        nbRx= len(self.receiver)       
        sourceLocations=empty((0,3))       
        
        threshold=linalg.norm(initial_guess)  
        error=1000
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            
            #result=initial_guess
            delta=int(floor(nbRx/4.0))
            for j in range(delta):
            #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                rx=empty((0,3))
                translationVector=self.receiver[j+3*delta].xi #store the translation vector of the new coordinate system              
                for k in range(3):                    
                    rx=append(rx,(self.receiver[j+k*delta].xi-translationVector).T,axis=0)
                rx=append(rx,array([[0,0,0]]),axis=0) #origin now new at receiver 4 of the current subgroup
                                                               
                toas=empty((0,))
                
                for i in range(4):   
                    toas=append(toas,self.receiver[j+i*delta].toa)
                
                result=equations(rx,toas,initial_guess)+translationVector.reshape((3,))
                
                e=linalg.norm(result-initial_guess.flatten())
            
                if e<error:
                #abs(linalg.norm(result)-threshold)<error:
                    #print "the solution is {1} : {0}".format(result,j)
                    sourceLocations=append(sourceLocations,result.reshape((1,3)),axis=0)
                
#        print "--------------------------------"
#        print "the average position is: {0}".format(mean(sourceLocationAvg,axis=0))
        #print "the alternative solution is : {0}".format( transpose(xs_est_2+self.translationVector.xi))
        return sourceLocations
        
        
    def SolveNonLinearTDOAFDOAEquations(self,initial_guess):
        
        def equationsFDOA(s,v,fdoa,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                
                func_0=dot(-(s[0,:]-u)/linalg.norm(s[0,:]-u),v[0,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_0=f0/c*func_0-(fdoa[0]-fdoa[3])
                
                func_1=dot(-(s[1,:]-u)/linalg.norm(s[1,:]-u),v[1,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_1=f0/c*func_1-(fdoa[1]-fdoa[3])
                
                func_2=dot(-(s[2,:]-u)/linalg.norm(s[2,:]-u),v[2,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_2=f0/c*func_2-(fdoa[2]-fdoa[3])
                
                return (func_0,func_1,func_2)
            
            f0=self.carrierFrequency
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1 
            
        def equationsTDOA(s,toas,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                r1=c*toas[3]               

                func_0=dot(s[0,:],s[0,:])-dot(s[3,:],s[3,:])-2*dot((s[0,:]-s[3,:]),u)-(c*(toas[0]-toas[3]))**2-2*r1*c*(toas[0]-toas[3])            
                
                func_1=dot(s[1,:],s[1,:])-dot(s[3,:],s[3,:])-2*dot((s[1,:]-s[3,:]),u)-(c*(toas[1]-toas[3]))**2-2*r1*c*(toas[1]-toas[3])
                
                func_2=dot(s[2,:],s[2,:])-dot(s[3,:],s[3,:])-2*dot((s[2,:]-s[3,:]),u)-(c*(toas[2]-toas[3]))**2-2*r1*c*(toas[2]-toas[3])
                
                return (func_0,func_1,func_2)
            
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1       

#        print "------------"
#        print "Compute data by solving the nonlinear TDOA equations numerically with an initial guess"
#        print "------------"        
        
        nbRx= len(self.receiver)       
        sourceLocationsTDOA=empty((0,3))
        sourceLocationsFDOA=empty((0,3)) 
        
        threshold=linalg.norm(initial_guess)  
        error=1000
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            
            #result=initial_guess
            delta=int(floor(nbRx/4.0))
            for j in range(delta):
            #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                rx=empty((0,3))
                translationVector=self.receiver[j+3*delta].xi #store the translation vector of the new coordinate system              
                for k in range(3):                    
                    rx=append(rx,(self.receiver[j+k*delta].xi-translationVector).T,axis=0)
                rx=append(rx,array([[0,0,0]]),axis=0) #origin now new at receiver 4 of the current subgroup
                                                               
                toas=empty((0,))
                speeds=empty((0,3))
                fdoas=empty((0,))
                
                for i in range(4):   
                    speeds=append(speeds,self.receiver[j+i*delta].v.T,axis=0)
                    fdoas=append(fdoas,self.receiver[j+i*delta].foa)   
                    toas=append(toas,self.receiver[j+i*delta].toa)
                
                resultTDOA=equationsTDOA(rx,toas,initial_guess)+translationVector.reshape((3,))
                resultFDOA=equationsFDOA(rx,speeds,fdoas,initial_guess)+translationVector.reshape((3,))
                
                efdoa=linalg.norm(resultFDOA-initial_guess.flatten())
                etdoa=linalg.norm(resultTDOA-initial_guess.flatten())
            
                if efdoa<error and etdoa<error :                
                    #abs(linalg.norm(result)-threshold)<error:
                    #print "the solution is {1} : {0}".format(result,j)
                    sourceLocationsTDOA=append(sourceLocationsTDOA,resultTDOA.reshape((1,3)),axis=0)
                    sourceLocationsFDOA=append(sourceLocationsFDOA,resultFDOA.reshape((1,3)),axis=0)
                
#        print "--------------------------------"
#        print "the average position is: {0}".format(mean(sourceLocationAvg,axis=0))
        #print "the alternative solution is : {0}".format( transpose(xs_est_2+self.translationVector.xi))
        return sourceLocationsTDOA,sourceLocationsFDOA
        
    def SolveNonLinearTDOAEquationsWithHistory(self,initial_guess):
        
        def equations(s,toas,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                r1=c*toas[3]               

                func_0=dot(s[0,:],s[0,:])-dot(s[3,:],s[3,:])-2*dot((s[0,:]-s[3,:]),u)-(c*(toas[0]-toas[3]))**2-2*r1*c*(toas[0]-toas[3])            
                
                func_1=dot(s[1,:],s[1,:])-dot(s[3,:],s[3,:])-2*dot((s[1,:]-s[3,:]),u)-(c*(toas[1]-toas[3]))**2-2*r1*c*(toas[1]-toas[3])
                
                func_2=dot(s[2,:],s[2,:])-dot(s[3,:],s[3,:])-2*dot((s[2,:]-s[3,:]),u)-(c*(toas[2]-toas[3]))**2-2*r1*c*(toas[2]-toas[3])
                
                return (func_0,func_1,func_2)
            
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1       

#        print "------------"
#        print "Compute data by solving the nonlinear TDOA equations numerically with an initial guess"
#        print "------------"        
        
        nbRx= len(self.receiver)       
        sourceLocationsWithHistory=empty((0,3))       
        

        error=1000
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            
            for tmpNbPoints in range(4,nbRx,4):          
                delta=int(floor(tmpNbPoints/4.0))
                nbResults=0
                result_avg=zeros((3,)) 
                for j in range(delta):
                #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                    rx=empty((0,3))
                    translationVector=self.receiver[j+3*delta].xi #store the translation vector of the new coordinate system              
                    for k in range(3):                    
                        rx=append(rx,(self.receiver[j+k*delta].xi-translationVector).T,axis=0)
                    rx=append(rx,array([[0,0,0]]),axis=0) #origin now new at receiver 4 of the current subgroup
                                                                   
                    toas=empty((0,))
                    
                    for i in range(4):   
                        toas=append(toas,self.receiver[j+i*delta].toa)
                    
                    result=equations(rx,toas,initial_guess)+translationVector.reshape((3,))
                    
                    e=linalg.norm(result-initial_guess.flatten())
                
                    if e<error:
                        nbResults+=1.0
                        if nbResults==1:
                            result_avg=result
                        else:
                            result_avg=(nbResults-1)/nbResults*result_avg+1/nbResults*result
                
                #print "{0},  {1}\n".format(delta,result_avg)
                if nbResults>0:
                    sourceLocationsWithHistory=append(sourceLocationsWithHistory,result_avg.reshape((1,3)),axis=0)
                
        return sourceLocationsWithHistory
        
    def SolveNonLinearFDOAEquationsWithHistory(self,initial_guess):
        
        def equations(s,v,fdoa,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                
                func_0=dot(-(s[0,:]-u)/linalg.norm(s[0,:]-u),v[0,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_0=f0/c*func_0-(fdoa[0]-fdoa[3])
                
                func_1=dot(-(s[1,:]-u)/linalg.norm(s[1,:]-u),v[1,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_1=f0/c*func_1-(fdoa[1]-fdoa[3])
                
                func_2=dot(-(s[2,:]-u)/linalg.norm(s[2,:]-u),v[2,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_2=f0/c*func_2-(fdoa[2]-fdoa[3])
                
                return (func_0,func_1,func_2)
            
            f0=self.carrierFrequency
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1 
                               
        nbRx= len(self.receiver)       
        sourceLocationsWithHistory=empty((0,3))            

        error=1000
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. We only take the first four points."
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            
            for tmpNbPoints in range(4,nbRx,4):
                delta=int(floor(tmpNbPoints/4.0))
                nbResults=0
                result_avg=zeros((3,)) 
                for j in range(delta):
                #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                    rx=empty((0,3))
                    translationVector=self.receiver[j+3*delta].xi #store the translation vector of the new coordinate system              
                    for k in range(3):                    
                        rx=append(rx,(self.receiver[j+k*delta].xi-translationVector).T,axis=0)
                    rx=append(rx,array([[0,0,0]]),axis=0) #origin now new at receiver 4 of the current subgroup
                                                                   
                    speeds=empty((0,3))
                    fdoas=empty((0,))
                    
                    for i in range(4):   
                        speeds=append(speeds,self.receiver[j+i*delta].v.T,axis=0)
                        fdoas=append(fdoas,self.receiver[j+i*delta].foa)
                    
                    result=equations(rx,speeds,fdoas,initial_guess)+translationVector.reshape((3,))
                    
                    e=linalg.norm(result-initial_guess.flatten())
                
                    if e<error:
                        nbResults+=1.0
                        if nbResults==1:
                            result_avg=result
                        else:
                            result_avg=(nbResults-1)/nbResults*result_avg+1/nbResults*result
                
                if nbResults>0:
                    sourceLocationsWithHistory=append(sourceLocationsWithHistory,result_avg.reshape((1,3)),axis=0)
                    
        return sourceLocationsWithHistory
        
        
    def SolveNonLinearTDOAFDOAEquationsWithHistory(self,initial_guess):
        
        def equationsFDOA(s,v,fdoa,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                
                func_0=dot(-(s[0,:]-u)/linalg.norm(s[0,:]-u),v[0,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_0=f0/c*func_0-(fdoa[0]-fdoa[3])
                
                func_1=dot(-(s[1,:]-u)/linalg.norm(s[1,:]-u),v[1,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_1=f0/c*func_1-(fdoa[1]-fdoa[3])
                
                func_2=dot(-(s[2,:]-u)/linalg.norm(s[2,:]-u),v[2,:])-dot(-(s[3,:]-u)/linalg.norm(s[3,:]-u),v[3,:])
                func_2=f0/c*func_2-(fdoa[2]-fdoa[3])
                
                return (func_0,func_1,func_2)
            
            f0=self.carrierFrequency
            c=self.c_speed  
        
            result1=sc.optimize.fsolve(f,initial)
            return result1 
            
        def equationsTDOA(s,toas,initial):
        #s: 4 x 3 matrix ; position vector of receiver i-th in row
        #v: 4 x 3 matrix ; speed vector of receiver i-th in row
        #the 4th-row is normally (0,0,0) because it is the last receiver that was chosen as the new local origin (0,0,0)
        #fdoa: 4x1 vector; foda
            def f(u):
                r1=c*toas[3]               

                func_0=dot(s[0,:],s[0,:])-dot(s[3,:],s[3,:])-2*dot((s[0,:]-s[3,:]),u)-(c*(toas[0]-toas[3]))**2-2*r1*c*(toas[0]-toas[3])            
                
                func_1=dot(s[1,:],s[1,:])-dot(s[3,:],s[3,:])-2*dot((s[1,:]-s[3,:]),u)-(c*(toas[1]-toas[3]))**2-2*r1*c*(toas[1]-toas[3])
                
                func_2=dot(s[2,:],s[2,:])-dot(s[3,:],s[3,:])-2*dot((s[2,:]-s[3,:]),u)-(c*(toas[2]-toas[3]))**2-2*r1*c*(toas[2]-toas[3])
                
                return (func_0,func_1,func_2)
            
            c=self.c_speed      
        
            result1=sc.optimize.fsolve(f,initial)
            return result1       

#        print "------------"
#        print "Compute data by solving the nonlinear TDOA equations numerically with an initial guess"
#        print "------------"        
        
        nbRx= len(self.receiver)       
        sourceLocationsWithHistory=empty((0,3))
        print "nb of receivers: {0}\n".format(nbRx) 
        error=200000
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            
            #result=initial_guess
            for tmpNbPoints in range(4,nbRx,4):
                delta=int(floor(tmpNbPoints/4.0))
                nbResults=0
                result_avg=zeros((3,)) 
                for j in range(delta):
                #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                    rx=empty((0,3))
                    translationVector=self.receiver[j+3*delta].xi #store the translation vector of the new coordinate system              
                    for k in range(3):                    
                        rx=append(rx,(self.receiver[j+k*delta].xi-translationVector).T,axis=0)
                    rx=append(rx,array([[0,0,0]]),axis=0) #origin now new at receiver 4 of the current subgroup
                                                                   
                    toas=empty((0,))
                    speeds=empty((0,3))
                    fdoas=empty((0,))
                    
                    for i in range(4):   
                        speeds=append(speeds,self.receiver[j+i*delta].v.T,axis=0)
                        fdoas=append(fdoas,self.receiver[j+i*delta].foa)   
                        toas=append(toas,self.receiver[j+i*delta].toa)
                    
                    resultTDOA=equationsTDOA(rx,toas,initial_guess)+translationVector.reshape((3,))
                    resultFDOA=equationsFDOA(rx,speeds,fdoas,initial_guess)+translationVector.reshape((3,))
                    
                    efdoa=linalg.norm(resultFDOA-initial_guess.flatten())
                    etdoa=linalg.norm(resultTDOA-initial_guess.flatten())
                    
                    if efdoa<error and etdoa<error : 
                        nbResults+=1.0
                        result=(resultTDOA+resultFDOA)/2
                        if nbResults==1:                      
                            result_avg=result
                        else:
                            result_avg=(nbResults-1)/nbResults*result_avg+1/nbResults*result
                #print "{0}\n".format(result_avg)
                if nbResults>0:
                    sourceLocationsWithHistory=append(sourceLocationsWithHistory,result_avg.reshape((1,3)),axis=0)
              
        return sourceLocationsWithHistory
    
    def dfdx(self,xj,vj):
        den=asscalar(((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2))
        num=asscalar((xj[0]-self.u_hat[0])*(vj[0]*(xj[0]-self.u_hat[0])+vj[1]*(xj[1]-self.u_hat[1])+vj[2]*(xj[2]-self.u_hat[2])))
        return -self.carrierFrequency/self.c*(num/den**1.5-vj[0]/den**0.5)
    
    def dfdy(self,xj,vj):
        den=asscalar(((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2))
        num=asscalar((xj[1]-self.u_hat[1])*(vj[0]*(xj[0]-self.u_hat[0])+vj[1]*(xj[1]-self.u_hat[1])+vj[2]*(xj[2]-self.u_hat[2])))
        return -self.carrierFrequency/self.c*(num/den**1.5-vj[1]/den**0.5)
    
    def dfdz(self,xj,vj):
        den=asscalar(((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2))
        num=asscalar((xj[2]-self.u_hat[2])*(vj[0]*(xj[0]-self.u_hat[0])+vj[1]*(xj[1]-self.u_hat[1])+vj[2]*(xj[2]-self.u_hat[2])))
        return -self.carrierFrequency/self.c*(num/den**1.5-vj[2]/den**0.5)
    
    def dpdx(self,xj):
        num=-asscalar((xj[0]-self.u_hat[0]))
        den=asscalar(sqrt((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2))
        return num/den
    
    def dpdy(self,xj):
        num=-asscalar((xj[1]-self.u_hat[1]))
        den=asscalar(sqrt((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2))
        return num/den
    
    def dpdz(self,xj):
        num=-asscalar((xj[2]-self.u_hat[2]))
        den=asscalar(sqrt((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2))
        return num/den
    
    def p(self,xj):
        return sqrt((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2)
        
    def f(self,xj,vj):
        den=asscalar(sqrt(((xj[0]-self.u_hat[0])**2+(xj[1]-self.u_hat[1])**2+(xj[2]-self.u_hat[2])**2)))*self.c_speed
        num=asscalar(vj[0]*(xj[0]-self.u_hat[0])+vj[1]*(xj[1]-self.u_hat[1])+vj[2]*(xj[2]-self.u_hat[2]))*self.carrierFrequency

        return self.carrierFrequency-num/den
        
    def Test_f(self):
        print "TEST*******\n"
        for i in range(len(self.receiver)):
            pass
            #print "Doppler frequency  stored {0} :  Doppler frequency computed : {1}\n".format(self.receiver[i].foa,self.f(self.receiver[i].xi.T.flatten(),self.receiver[i].v.T.flatten()))
            #print "source {0}\n".format(self.u_hat-self.receiver[i].xi.T.flatten())
    
    def CreateHn(self,xi,vi,t):
        #print "parent CreateHn method\n"
        n=len(xi)-1
        dp1dx=self.dpdx(xi[0,:])
        dp1dy=self.dpdy(xi[0,:])
        dp1dz=self.dpdz(xi[0,:])
        self.Hn=zeros((2*n,4))
        for i in range(n):
            x=xi[i+1,:]
            #print "{0}.th point : {1}\n".format(i,x.flatten())
            v=vi[i+1,:]
            #print "{0}.th point : {1}\n".format(i,v.flatten())
            ax=dp1dx-self.dpdx(x)    
            ay=dp1dy-self.dpdy(x)
            az=dp1dz-self.dpdz(x)
            bx=self.dfdx(x,v)
            by=self.dfdy(x,v)
            bz=self.dfdz(x,v)
            self.Hn[2*i,:]=[ax,ay,az,(t[0]-t[i+1])*self.c_speed]
            self.Hn[2*i+1,:]=[bx,by,bz,0]
            
    def CreateHn_TDOA(self,xi,vi,t):
        #print "parent CreateHn method\n"
        n=len(xi)-1
        dp1dx=self.dpdx(xi[0,:])
        dp1dy=self.dpdy(xi[0,:])
        dp1dz=self.dpdz(xi[0,:])
        self.Hn=zeros((n,4))
        for i in range(n):
            x=xi[i+1,:]
            ax=dp1dx-self.dpdx(x)    
            ay=dp1dy-self.dpdy(x)
            az=dp1dz-self.dpdz(x)
            #print "c*(t1-ti) = {0}\n".format((t[0]-t[i+1])*self.c_speed)
            self.Hn[i,:]=[ax,ay,az,(t[0]-t[i+1])*self.c_speed]

            
    def CreateHn_FDOA(self,xi,vi):
        #print "parent CreateHn method\n"
        n=len(xi)-1
        self.Hn=zeros((n,3))
        for i in range(n):
            x=xi[i+1,:]
            v=vi[i+1,:]
            bx=self.dfdx(x,v)
            by=self.dfdy(x,v)
            bz=self.dfdz(x,v)
            self.Hn[i,:]=[bx,by,bz]
            
    def CreateHnMinus(self):
        #print "parent CreateHnMinus method\n"
        A=dot(self.Hn.T,self.Hn)
        print "\n A= {0}\n".format(A)
        B=linalg.inv(A)
        self.HnMinus=dot(B,self.Hn.T)  
                                   
    
    def CreateRho(self,xi,vi,tdoas,foas):
        #print "parent CreateCovarianceRho method\n"
        n=len(xi)-1
        self.rho=zeros((2*n,1))
        p1=self.p(xi[0,:])
        for i in range(n):
            self.rho[2*i,0]=(tdoas[0]-tdoas[i+1])*self.c_speed-(p1-self.p(xi[i+1,:]))
            self.rho[2*i+1,0]=foas[i+1]-self.f(xi[i+1,:],vi[i+1,:])
    
    def CreateRho_TDOA(self,xi,vi,tdoas,foas):
        #print "parent CreateCovarianceRho method\n"
        n=len(xi)-1
        self.rho=zeros((n,1))
        p1=self.p(xi[0,:])
        for i in range(n):
            self.rho[i,0]=(tdoas[0]-tdoas[i+1])*self.c_speed-(p1-self.p(xi[i+1,:]))

            
    def CreateRho_FDOA(self,xi,vi,tdoas,foas):
        #print "parent CreateCovarianceRho method\n"
        n=len(xi)-1
        self.rho=zeros((n,1))
        for i in range(n):
            self.rho[i,0]=foas[i+1]-self.f(xi[i+1,:],vi[i+1,:])
            
    
    def TaylorApproxTDOAFDOAEquations(self,initial_guess):
        
        self.u_hat=initial_guess.flatten()
           
        nbRx= len(self.receiver)      
        source=zeros((nbRx-3,4))
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            #translationVector=self.receiver[0].xi #store the translation vector of the new coordinate system              
            for j in range(4,nbRx+1,1):
            #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                xi=empty((0,3))
                vi=empty((0,3))
                tdoas=empty((0,))
                fdoas=empty((0,))
                t=empty((0,))
                for k in range(j):        
                    #xi=append(xi,(self.receiver[k].xi-translationVector).T,axis=0)
                    xi=append(xi,self.receiver[k].xi.T,axis=0)
                    vi=append(vi,self.receiver[k].v.T,axis=0)
                    tdoas=append(tdoas,self.receiver[k].toa,axis=0)
                    fdoas=append(fdoas,self.receiver[k].foa,axis=0)
                    t=append(t,self.receiver[k].t,axis=0)
                
                self.CreateRho(xi,vi,tdoas,fdoas)
                self.CreateHn(xi,vi,t)
                self.CreateHnMinus()
                s=dot(self.HnMinus,self.rho)
                source[j-4,:]=s[0:4,0].T

        return source
        
    def TaylorApproxTDOAEquations(self,initial_guess):
        
        self.u_hat=initial_guess.flatten()
           
        nbRx= len(self.receiver)     
        
        source=zeros((nbRx-3,4))
         
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            #translationVector=self.receiver[0].xi #store the translation vector of the new coordinate system              
            for j in range(4,nbRx+1,1):
            #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                xi=empty((0,3))
                vi=empty((0,3))
                tdoas=empty((0,))
                fdoas=empty((0,))
                t=empty((0,))
                for k in range(j):        
                    #xi=append(xi,(self.receiver[k].xi-translationVector).T,axis=0)
                    xi=append(xi,self.receiver[k].xi.T,axis=0)
                    vi=append(vi,self.receiver[k].v.T,axis=0)
                    tdoas=append(tdoas,self.receiver[k].toa,axis=0)
                    fdoas=append(fdoas,self.receiver[k].foa,axis=0)
                    t=append(t,self.receiver[k].t,axis=0)
                
                self.CreateRho_TDOA(xi,vi,tdoas,fdoas)
                self.CreateHn_TDOA(xi,vi,t)
                self.CreateHnMinus()
                s=dot(self.HnMinus,self.rho)
                source[j-4,:]=s[0:4,0].T

        return source
        
    def TaylorApproxFDOAEquations(self,initial_guess):
        
        self.u_hat=initial_guess.flatten()
           
        nbRx= len(self.receiver)      
        source=zeros((nbRx-3,3))
        
        if nbRx<4:
#            print "Not enough data points available."
            return
        if nbRx>=4: #only take the first 4 sensors in case there are more
#            print "More data points than necessary. The total path is segmented into 4 segments."
#            print "Then, a data point is taken from each of the segment and the navigation solution is computed"
#            print "--------------------------------"
            #split set of receivers in subgroups of 4 receivers
            #translationVector=self.receiver[0].xi #store the translation vector of the new coordinate system              
            for j in range(4,nbRx+1,1):
            #1) format data; take the 4th sensor as origin, adjust the sensor data accordingly
                xi=empty((0,3))
                vi=empty((0,3))
                tdoas=empty((0,))
                fdoas=empty((0,))
                for k in range(j):        
                    #xi=append(xi,(self.receiver[k].xi-translationVector).T,axis=0)
                    xi=append(xi,self.receiver[k].xi.T,axis=0)
                    vi=append(vi,self.receiver[k].v.T,axis=0)
                    tdoas=append(tdoas,self.receiver[k].toa,axis=0)
                    fdoas=append(fdoas,self.receiver[k].foa,axis=0)
                
                self.CreateRho_FDOA(xi,vi,tdoas,fdoas)
                self.CreateHn_FDOA(xi,vi)
                self.CreateHnMinus()
                s=dot(self.HnMinus,self.rho)

                source[j-4,:]=s[0:3,0].T

        return source