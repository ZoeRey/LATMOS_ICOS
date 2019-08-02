#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""
import statistics
import rayleigh as ray
import matplotlib.pyplot as plt
import numpy as np
import outils_OHP as outils
from scipy.optimize import curve_fit
import os.path
import scipy.signal
import matplotlib.cm as cm
from scipy import *

#%%Lecture des fichiers de la journée

annee=["2019"]
dates=["20190704"]


#trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2017/CE376_'+dates[0]+'/mli2_'+dates[0]+'.000000.OHP'
trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2019/CE376_'+dates[0]+'/mli2_'+dates[0]+'.000000.OHP'

#file=np.loadtxt(trajet,skiprows=33)
file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]  

#%%
num1=["0","1","2"]
num2=["0","1","2","3","4","5","6","7","8","9"]
num3=["0","1","2","3","4","5"]
num4=["0","1","2"]
hours=[]
minutes=[]
secondes=[]


#data532paral=[[0]*len(altitude) for k in range(0,1440)]
data532paral=[[0]*1440 for k in range(len(altitude))]
data532paral=np.array(data532paral)

i=0
for a in num1:
    for b in num2:
        for c in num3:
            for d in num2:
                for e in num2:
                    for f in num2:
                        #trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2017/CE376_'+dates[0]+'/mli2_'+dates[0]+'.'+a+b+c+d+e+f+'.OHP'
                        trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2019/CE376_'+dates[0]+'/mli2_'+dates[0]+'.'+a+b+c+d+e+f+'.OHP'
                        if os.path.exists(trajet):
                            print(i)
                            hours.append(a+b)
                            minutes.append(c+d)
                            secondes.append(e+f)
                            #file=np.loadtxt(trajet,skiprows=33)
                            file=np.loadtxt(trajet,skiprows=36)
                            col1=file[:,1]
                            col1reverse=col1[::-1]   
                            data532paral[:,i]=col1reverse
                            i=i+1

#%%Rayleigh                          
LRp=50
LRm=(8*np.pi)/3.

for i in range(len(altitude)):
    altitude[i]=altitude[i]+0.68


source_mol='ISA'
rayleigh532=ray.calc_molecular_profile(source_mol,altitude,532)

#%%time et altitude min correction géométrie

time=[]    
for i in range(len(data532paral[0])):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i])/100)/100.)
    
#%%Moyenne 1h
    
datamoy=[[0]*48 for k in range(len(altitude))]
datamoy=np.array(datamoy)
for i in range(0,48):
    data=data532paral[:,30*i:30*i+30]
    datamean=np.mean(data,axis=1)
    datamoy[:,i]=scipy.signal.medfilt(datamean,25)
    

timemoy=np.arange(0,24,0.5)


plt.plot(data532paral[:,19*60+15],altitude)
plt.plot(datamoy[:,19*2],altitude)
    
#%%plt.figure()

lev=np.linspace(0.,400.,800)
plt.figure()
plt.title("PR2_ICOS")
CS = plt.contourf(timemoy[42:], altitude[654:1292],datamoy[654:1292,42:],levels=lev,cmap='jet', extend='min')
colobar=plt.colorbar(CS)
plt.xlabel('time (UT)')
plt.ylabel('altitude (km)')
plt.show()

#%%


lev=np.linspace(0.,56.,20)
plt.figure()
plt.title("PR2_ICOS")
CS = plt.contourf(timemoy[40:46], altitude[225:1359],datamoy[225:1359,40:46],levels=lev,cmap='jet', extend='min')
colobar=plt.colorbar(CS)
plt.xlabel('time (UT)')
plt.ylabel('altitude (km)')
plt.legend()
plt.show()

#%%inversion pour chaque profil



trajet_geom='/net/nfs/home/rey/Documents/code_Zoe/OHP/fct_geom/fct_geom_OHP_20170901.txt'
file=np.loadtxt(trajet_geom)
fgeom532=file[:,1]

geom532_fit=list(fgeom532[::-1])

for i in range(len(altitude)-len(geom532_fit)):
    geom532_fit.append(1)

geom532_fit[0]=geom532_fit[1]



data532cor=[[0.]*48 for k in range(len(altitude))]
data532cor=np.array(data532cor)

for i in range(len(altitude)):
    for j in range(len(timemoy)):
        data532cor[i,j]=datamoy[i,j]/geom532_fit[i]

#%%
        
cste532=25010
izr=1960

beta_aer=[[]*48 for k in range(0,izr)]
#beta_aer=np.array(beta_aer)
beta_aer=[]

for i in range(0,48):
    data=data532cor[:,i]
    beta532_aer=outils.calc_beta_aer_stable2(izr,altitude,rayleigh532,data,LRp,LRm,cste532)
    beta_aer.append(beta532_aer)
    
    
print(beta_aer)

beta_aer=np.array(beta_aer)
beta_aer2=np.transpose(beta_aer)

#%%
    
lev=np.linspace(0.,0.01,20)
plt.figure()
plt.title("beta_aer_ICOS")
CS = plt.contourf(timemoy,altitude[:izr],beta_aer2,levels=lev,cmap='jet', extend='min')
plt.xlabel('time (UT)')
plt.ylabel('altitude (km)')
colobar=plt.colorbar(CS)
plt.show()

lev=np.linspace(0.,0.0065,300)
plt.figure()
plt.title("beta_aer_ICOS")
CS = plt.contourf(timemoy[42:],altitude[654:1292],beta_aer2[654:1292,42:],levels=lev,cmap='jet', extend='min')
colobar=plt.colorbar(CS)
plt.xlabel('time (UT)')
plt.ylabel('altitude (km)')
plt.show()


#%%

beta532_aer_f=scipy.signal.medfilt(beta532_aer,9)


alpha532_aer,AOD532,trans532_aer=outils.calc_alpha_AOD_trans(beta532_aer_f,LRp,altitude[:izr])

