#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""
import matplotlib.pyplot as plt
import os.path
import scipy.signal
import numpy as np

#%%Lecture des fichiers de la journée

annee=["2017"]
dates=["20170901"]


trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2017/CE376_'+dates[0]+'/mli2_'+dates[0]+'.000000.OHP'
#trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2019/CE376_'+dates[0]+'/mli2_'+dates[0]+'.000000.OHP'

file=np.loadtxt(trajet,skiprows=33)
#file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]  

for i in range(len(altitude)):
    altitude[i]=altitude[i]+0.68

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
                        trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2017/CE376_'+dates[0]+'/mli2_'+dates[0]+'.'+a+b+c+d+e+f+'.OHP'
                        #trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2019/CE376_'+dates[0]+'/mli2_'+dates[0]+'.'+a+b+c+d+e+f+'.OHP'
                        if os.path.exists(trajet):
                            print(i)
                            hours.append(a+b)
                            minutes.append(c+d)
                            secondes.append(e+f)
                            file=np.loadtxt(trajet,skiprows=33)
                            #file=np.loadtxt(trajet,skiprows=36)
                            col1=file[:,1]
                            col1reverse=col1[::-1]   
                            data532paral[:,i]=col1reverse
                            i=i+1


#%%time et altitude min correction géométrie

time=[]    
for i in range(len(data532paral[0])):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i])/100)/100.)
    
#%%Moyenne 1h
    
datamoy=[[0]*96 for k in range(len(altitude))]
datamoy=np.array(datamoy)
for i in range(0,96):
    data=data532paral[:,15*i:15*i+15]
    datamean=np.mean(data,axis=1)
    datamoy[:,i]=datamean
    

timemoy=np.arange(0,24,0.25)


plt.plot(data532paral[:300,8*60+8],altitude[:300])
plt.plot(datamoy[:300,8*4],altitude[:300])

#%%recherche du max 

iymax=[]

for i in range(8*4,18*4+1):
    y=list(datamoy[:300,i])
    iymax.append(y.index(max(y)))

#%%filtre médian
    
time_hcla=timemoy[8*4:18*4+1]
data_minmax=[[] for i in range(len(time_hcla))]
xbis=[[] for i in range(len(time_hcla))]

for i in range(0,(18*4+1)-8*4):
    if i<=32 :
        data_minmax[i]=data_minmax[i]+list(datamoy[iymax[i]-15:150,i+8*4])
        xbis[i]=xbis[i]+list(altitude[iymax[i]-15:150])   
    else :
        data_minmax[i]=data_minmax[i]+list(datamoy[15:100,i+8*4])
        xbis[i]=xbis[i]+list(altitude[15:100])   



"""
plt.plot(data532_medianfilter[80-49],altitude[iymax[80-49]-4:iymax4[80-49]+12])
#plt.plot(data532paral[80][iymax[80-49]-4:iymax4[80-49]+12],altitude[iymax[80-49]-4:iymax4[80-49]+12])
plt.grid()
plt.show()
"""
#%% Approximation polynôme de degré 5 et calcul dérivé 2nd de ce polynome

datafit=[[] for i in range(len(time_hcla))]
datafitderiv2=[[] for i in range(len(time_hcla))]

for i in range(0,(18*4+1)-8*4):
    x=xbis[i]
    y=data_minmax[i]
    coeff_poly=np.polyfit(x,y,5)

    for j in range(len(x)):
        datafit[i].append(coeff_poly[0]*x[j]*x[j]*x[j]*x[j]*x[j]+coeff_poly[1]*x[j]*x[j]*x[j]*x[j]+coeff_poly[2]*x[j]*x[j]*x[j]+coeff_poly[3]*x[j]*x[j]+coeff_poly[4]*x[j]+coeff_poly[5])
        datafitderiv2[i].append(20*coeff_poly[0]*x[j]*x[j]*x[j]+12*coeff_poly[1]*x[j]*x[j]+6*coeff_poly[2]*x[j]+2*coeff_poly[3])
    
    print(i)    

    plt.figure()
    plt.plot(datamoy[:500,i+32],altitude[:500],label="PR2")
    plt.plot(datafit[i],x,color='black',linewidth=2.5,linestyle='--',label="approximation polynomiale")
    #plt.plot(y,x)
    plt.xlabel("signal rétrodiffusé, 532// nm")
    plt.ylabel("altitude (km)")
    plt.legend()
    #plt.grid()
    plt.show()
    
    
#%% Recherche point d'inflexion
#il peut y avoir plusieurs points d'inflexion
pt_inflexion=[[] for i in range(0,(18*4+1)-8*4)]

for i in range(len(time_hcla)):
    x=xbis[i]
    for k in range(len(x)-1):
        if(datafitderiv2[i][k+1]>=0 and datafitderiv2[i][k]<=0) or (datafitderiv2[i][k+1]<=0 and datafitderiv2[i][k]>=0):
            pt_inflexion[i].append(x[k])
            
#%%Vérif s'il y a des endroits sans points d'inflexion
                
for i in range(len(pt_inflexion)):
    if len(pt_inflexion[i])==0:
        pt_inflexion[i]=pt_inflexion[i-1]       


#%% test en partant du début

hcla2=[]

for i in range(0,10):
    hcla2.append(altitude[iymax[i]])



for i in range(10,len(pt_inflexion)):
    diff1=100.
    for j in range(len(pt_inflexion[i])):
        diff=abs(pt_inflexion[i][j]-hcla2[i-1])
        if diff<diff1:
            diff1=diff
            jgood=j
    hcla2.append(pt_inflexion[i][jgood])  

#hcla2[33]=pt_inflexion[33][0]
#hcla2[34]=pt_inflexion[34][0]
hcla_f=scipy.signal.medfilt(hcla2,13)

plt.figure()
plt.plot(time_hcla,hcla_f,label="filtre médian",linewidth=0.2)
plt.plot(time_hcla,hcla2,linewidth=0.90)
plt.xlabel("time (UT)")
plt.ylabel("hcla (km)")
plt.legend()
#plt.grid()
plt.show()


