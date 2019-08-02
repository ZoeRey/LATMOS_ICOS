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
import outils
#from scipy.optimize import curve_fit
import scipy.signal 


#%% Read file

trajet="/net/nfs/home/rey/Documents/code_Zoe/LTA/data4Zoe/Case2_2019/aerosol_LTA_OHP_20190704.txt"
file=np.loadtxt(trajet,skiprows=17)

altitude=file[:,0]
PR2=file[:,1]
pressure=file[:,2]
temperature=file[:,3]
beta_aer=file[:,5]
beta_mol=file[:,6]
alpha_aer=file[:,7]
SR=file[:,8]

LRp=50
LRm=(8*np.pi)/3.


alpha_mol=[]
for i in range(len(pressure)):
    pressure[i]=pressure[i]*100
    alpha_mol.append(beta_mol[i]*LRm)
tau_mol=[np.trapz(alpha_mol[:i],x=altitude[:i]) for i in range(len(altitude))]
T_mol=[np.exp(-2*tau_mol[i]) for i in range(len(altitude))]



data532moy_LTA=PR2

plt.figure()
plt.plot(data532moy_LTA,altitude)
plt.xlabel("PR2 (UA)")
plt.ylabel("altitude(km)")
#plt.xscale("log")
plt.legend()
plt.show()


#%% calcul profil moleculaire et rayleigh

rayleigh532=ray.calc_molecular_profile(altitude,532,temperature,pressure)
rayleigh532_2=ray.calc_molecular_profile2("ISA",altitude,532)
rayleigh532_3=[beta_mol[i]*T_mol[i] for i  in range(len(altitude))]


plt.figure()
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,label="S")
plt.plot(rayleigh532_2['beta_m']*rayleigh532_2['tm'],altitude,label="Z")
plt.plot(rayleigh532_3,altitude,label="S_2")
plt.legend()
plt.show()





#%% calcul constante
#Attention : ajuster la zone où on fait la moyenne à la main

id8=128
id9=132

cste532=outils.calc_cste(data532moy_LTA,rayleigh532,id8,id9)
cste532_2=outils.calc_cste(data532moy_LTA,rayleigh532_2,id8,id9)



plt.figure()
plt.title('532//')
plt.plot(data532moy_LTA,altitude,linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh*K")
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='g',label="Rayleigh")
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("PR2 (UA)")
#plt.grid()
plt.legend()
plt.show()



#%%Calcul data corrigé de constante


data532cor=[]

for i in range(len(altitude)):
    data532cor.append(data532moy_LTA[i]/cste532)

#%%Calcul inversion

izr=132

cste532=13800000

beta532_aer_LTA=outils.calc_beta_aer_stable2(izr,altitude,rayleigh532,data532moy_LTA,LRp,LRm,cste532)


plt.figure()
plt.title('532//')
plt.plot(data532moy_LTA,altitude,linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh*K")
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='g',label="Rayleigh")
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("PR2 (UA)")
#plt.grid()
plt.legend()
plt.show()


#%%Calcul alpha et transmission

alpha532_aer_LTA,AOD532_LTA,trans532_aer_LTA=outils.calc_alpha_AOD_trans(beta532_aer_LTA,LRp,altitude[:izr])

#%%plot profil rétrodiffusion

plt.figure()
plt.title("Profil de rétrodiffusion 532 ")
plt.plot(beta532_aer_LTA,altitude[:izr],label="Zoé")
plt.plot(beta_aer[:izr],altitude[:izr],label="Serguey")
plt.grid()
plt.xlabel("beta 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
#plt.xlim(0,0.0009)
plt.legend()
plt.show()

"""
plt.figure()
plt.title("Profil de rétrodiffusion 532 ")
plt.plot(beta532_aer_LTA,altitude[:izr],label="Zoé")

plt.grid()
plt.xlabel("beta 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
#plt.xlim(0,0.0009)
plt.legend()
plt.show()
"""
#%% plot LTA+ICOS
plt.figure()
#plt.title("Profil de rétrodiffusion 532 ")
plt.plot(beta532_aer_LTA[:],altitude[:izr],label='LTA')
plt.plot(beta532_aer_OHP[12:1230],altitude_OHP[12:1230],label='ICOS')
#plt.grid()
plt.xlabel("beta 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
#plt.xlim(0,0.0011)
plt.legend()
plt.show()

plt.figure()
plt.title("Profil de rétrodiffusion 532 ")
plt.plot(beta532_aer_OHP[10:],altitude_OHP[10:1965])
plt.grid()
plt.xlabel("beta 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
#plt.xlim(0,0.0009)
plt.legend()
plt.show()


#%%AOD

#i_pic_begin=1070
#i_pic_end=1167
#i_pic_begin=682
#i_pic_end=820
#i_pic_begin=92
#i_pic_end=106
i_pic_begin=1015
i_pic_end=1100

#AOD_pic=np.trapz(alpha532_aer_OHP[i_pic_begin:i_pic_end],x=altitude_OHP[i_pic_begin:i_pic_end])
#print("AOD_pic=",AOD_pic)


AOD_tot_Z=np.trapz(alpha532_aer_LTA,x=altitude[:izr])
AOD_tot_S=np.trapz(alpha_aer[:izr],x=altitude[:izr])
print("AOD_tot Zoé=",AOD_tot_Z)
print("AOD_tot Sergey=",AOD_tot_S)




#%%Rapport de diffusion

R532=[]
for i in range(len(beta532_aer_LTA)):
    R532.append(1.+beta532_aer_LTA[i]/rayleigh532['beta_m'][i])


plt.figure()
plt.title("20170901 19:44-22:42 UT")
plt.plot(R532[:izr],altitude[:izr],label="LTA")
plt.plot(R532_ICOS,altitude_ICOS[:2000],label="ICOS")
#plt.grid()
plt.xlim(1,8)
plt.ylim(6,20)
plt.xlabel("Rapport de diffusion 532// nm")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

