#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""

import matplotlib.pyplot as plt
import numpy as np
import rayleigh as ray
import outils



#%% Read file

trajet_alt="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/LTAOHP20190704altkm.txt"
#trajet_alt="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/OHPLTA20170901altkm.txt"
file_alt=np.loadtxt(trajet_alt)
altitude=list(file_alt)

trajet_time="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/LTAOHP20190704timeUT.txt"
#trajet_time="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/OHPLTA20190901timeUT.txt"
file_time=np.loadtxt(trajet_time)
time=list(file_time)
time=time[:202] #PROBLEME !!!


trajet_PR2="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/LTAOHP20190704pr2k2d.txt"
#trajet_PR2="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/OHPLTA20170901pr2k2d.txt"
file_PR2=np.loadtxt(trajet_PR2,skiprows=2)

PR2_LTA=[[] for i in range(len(altitude))]

for i in range(len(altitude)):
    #print(i)
    PR2_LTA[i]=PR2_LTA[i]+list(file_PR2[:,i])


#%%
hours=[]
minutes=[]
secondes=[]
for i in range(len(time)):
    hours.append(int(time[i]/10000))
    minutes.append(int((time[i]-hours[i]*10000)/100))
    secondes.append(int(time[i]-(hours[i]*10000+minutes[i]*100)))

time=[]   
for i in range(len(hours)):
    time.append(hours[i]+((100./60.)*minutes[i]+secondes[i]/100)/100.)


#%%

PR2_LTA=np.array(PR2_LTA)

lev=np.linspace(0.,13500.,300)

plt.figure()
CS=plt.contourf(time[22:175],altitude[:65], PR2_LTA[:65,22:175],levels=lev,cmap='jet', extend='min')
plt.title("PR2_LTA")
plt.xlabel("time (UT)")
plt.ylabel("altitude (km)")
colobar=plt.colorbar(CS)
plt.legend()
plt.show()


#%%Rayleigh 

rayleigh532=ray.calc_molecular_profile2("ISA",altitude,532)


#%%test constante

data532moy=np.mean(PR2_LTA,axis=1)


cste532=15500000

plt.figure()
plt.plot(data532moy,altitude,linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m'][:]*rayleigh532['tm'][:],altitude[:],linestyle='--',color='k',label="Rayleigh*K")
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='g',label="Rayleigh")
#plt.yscale('Log')
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("(UA)")
#plt.grid()
plt.legend()
plt.show()

print(cste532)




#%%inversion pour chaque profil


izr=135
LRp=50
LRm=(8*np.pi)/3.


beta_aer_LTA=[]

for i in range(0,202):
    data=PR2_LTA[:,i]
    beta532_aer=outils.calc_beta_aer_stable2(izr,altitude,rayleigh532,data,LRp,LRm,cste532)
    print(beta532_aer)
    beta_aer_LTA.append(beta532_aer)
    
    
print(beta_aer_LTA)

beta_aer_LTA=np.array(beta_aer_LTA)
beta_aer_LTA2=np.transpose(beta_aer_LTA)

#%%
lev=np.linspace(0.,0.0004,300)

plt.figure()
CS=plt.contourf(time[22:175],altitude[:65], beta_aer_LTA2[:65,22:175],levels=lev,cmap='jet', extend='min')
plt.title("beta_aer_LTA")
plt.xlabel("time (UT)")
plt.ylabel("altitude (km)")
colobar=plt.colorbar(CS)
plt.legend()
plt.show()

#%% difference LTA ICOS




