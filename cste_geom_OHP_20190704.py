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
import outils_OHP
from scipy.optimize import curve_fit
import os.path
import scipy.signal

#%%Lecture des fichiers de la journée

annee=["2019"]
dates=["20190704"]


#trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2017/CE376_'+dates[0]+'/mli2_'+dates[0]+'.000000.OHP'
trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/ICOS_OHP/ascii/2019/CE376_'+dates[0]+'/mli2_'+dates[0]+'.235019.OHP'

#file=np.loadtxt(trajet,skiprows=33)
file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]  

for i in range(len(altitude)):
    altitude[i]=altitude[i]+0.680

#%%
num1=["0","1","2"]
num2=["0","1","2","3","4","5","6","7","8","9"]
num3=["0","1","2","3","4","5"]
num4=["0","1","2"]
hours=[]
minutes=[]
secondes=[]


data532paral=[[0]*len(altitude) for k in range(0,1440)]

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
                            data532paral[i]=col1reverse[:]
                            i=i+1
                            
LRp=50
LRm=(8*np.pi)/3.

#%%time et altitude min correction géométrie

time=[]    
for i in range(len(data532paral)):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i])/100)/100.)


#%%recherche pic altitude min

pic=[]
for a in range(0,22):
    time_moy_begin=a
    time_moy_end=a+1
    data532moy=[0]*len(altitude)
    for i in range(len(altitude)):
        for n in range(time_moy_begin*60,time_moy_begin*60+60,1):
            data532moy[i]=data532moy[i]+data532paral[n][i]
        data532moy[i]=data532moy[i]/60
    data=list(data532moy)
    pic.append(altitude[data.index(max(data[:50]))])
    

print(pic)
   
    

#%%moyenne sur 1h
t1=1236
t2=1439
t3=["20:36-23:59 UT"]    


"""
time_moy_begin=10
time_moy_end=11

for i in range(len(altitude)):
    for n in range(time_moy_begin*60,time_moy_begin*60+60,1):
        data532moy[i]=data532moy[i]+data532paral[n][i]
    data532moy[i]=data532moy[i]/60
    
data532moy=scipy.signal.medfilt(data532moy,13)

plt.figure()
plt.title("Vérification moyenne, 532//")
#plt.plot(data532paral[time_moy_begin*60+30],altitude)
plt.plot(data532moy,altitude,label="moy")
plt.legend()
plt.show()
"""
"""
data=list(data532moy)
plt.plot(data532paral[9*60+20][:100],altitude[:100])
print("pic à ",altitude[data.index(max(data[:100]))])
"""


for i in range(len(altitude)):
    for n in range(t1,t2,1):
        data532moy[i]=data532moy[i]+data532paral[n][i]
    data532moy[i]=data532moy[i]/(t2-t1)

data532moy=scipy.signal.medfilt(data532moy,25)

plt.figure()
plt.title('ICOS, 532nm, 20190704, '+t3[0])
#plt.plot(data532paral[1200],altitude)
plt.plot(data532moy,altitude)
plt.grid()
plt.ylabel("altitude (km)")
plt.xlabel("PR2 (UA)")
plt.legend()
plt.show()




#%%looking for blind zone et max du pic

time_moy_begin2=[]
time_moy_begin2.append(str(time_moy_begin))
time_moy_end2=[]
time_moy_end2.append(str(time_moy_end))


plt.figure()
plt.title('532// '+time_moy_begin2[0]+'h-'+time_moy_end2[0]+'h UT')
plt.plot(data532moy[:100],altitude[:100])
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.grid()
plt.show()

plt.figure()
plt.title('ICOS, 532nm '+time_moy_begin2[0]+'h-'+time_moy_end2[0]+'h UT')
plt.plot(data532moy[:20],altitude[:20])
plt.ylabel("altitude (km)")
plt.xlabel("PR2 (UA)")
#plt.grid()
#plt.xscale('log')
plt.show()

iav=7



#%% limite de bruit étendue
z0=0.62

data1=[]
for i in range(len(altitude)):
    data1.append(data532moy[i]/(altitude[i]-z0)**2)
 
bruit532=[]

ectype1=np.std(data1[1300:])

for i in range(len(altitude)):
    bruit532.append(2.*ectype1*(altitude[i]-z0)*(altitude[i]-z0))


#%%read fonction geom

trajet_geom='/net/nfs/home/rey/Documents/code_Zoe/OHP/fct_geom/fct_geom_OHP_20170901.txt'
file=np.loadtxt(trajet_geom)
fgeom532=file[:,1]

geom532_fit=list(fgeom532[::-1])
geom532_fit[0]=geom532_fit[1]

for i in range(len(altitude)-len(geom532_fit)):
    geom532_fit.append(1)

plt.figure()
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(geom532_fit[:119],altitude[:119],label="532nm")
plt.ylim(0,2)
#plt.grid()
plt.legend()
plt.show()


#%% calcul profil moleculaire et rayleigh

source_mol='ISA'
rayleigh532=ray.calc_molecular_profile(source_mol,altitude,532)


#%% calcul constante

idx1=np.where(altitude==20.12)
id1=idx1[0][0]
idx2=np.where(altitude==21.62)
id2=idx2[0][0]


#cste532=outils_OHP.calc_cste(data532moy,rayleigh532,id1,id2)
cste532=25010

plt.figure()
plt.title('532// '+t3[0])
plt.plot(data532moy[:],altitude[:],linewidth=1.,label="PR2")
#plt.plot(bruit532,altitude,label="limite bruit étendue")
plt.plot(cste532*rayleigh532['beta_m'][:]*rayleigh532['tm'][:],altitude[:],linestyle='--',color='k',label="Rayleigh*K")
plt.ylabel("altitude (km)")
#plt.grid()
plt.xlabel("PR2 voie 1 (UA)")
plt.legend()
plt.show()



print(cste532)


#%%Calcul data corrigé de constante

data532cor=[]
for i in range(len(altitude)):
    data532cor.append(data532moy[i]/cste532)


plt.figure()
plt.title('ICOS 532// '+t3[0])
plt.plot(data532cor,altitude,linewidth=1.,label="PR2/K")
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()
 

#%%Verification alignement

print(altitude[850:900])

idx1=np.where(altitude==13.37)
id1=idx1[0][0]
idx2=np.where(altitude==14.09)
id2=idx2[0][0]

id1=1160
id2=1180

cste532=np.log(np.mean(data532moy[id1:id2]))-np.log(np.mean(rayleigh532['beta_m'][id1:id2]*rayleigh532['tm'][id1:id2]))
cste532=np.exp(cste532)
print(cste532)

cste532=228180


plt.figure()
plt.title('ICOS, 532nm, '+t3[0])
plt.plot(data532moy,altitude,linewidth=1.,label="PR2")
plt.plot(cste532*rayleigh532['beta_m'][iav:]*rayleigh532['tm'][iav:],altitude[iav:],linestyle='--',color='k',label="Rayleigh*K")
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='g',label="Rayleigh")
#plt.yscale('Log')
plt.xscale('log')
plt.ylabel("altitude (km)")
plt.xlabel("(UA)")
plt.grid()
plt.legend()
plt.show()

print(cste532)


#%%Calcul data corrigé de constante et géométrie

data532cor=[]
for i in range(len(altitude)):
    data532cor.append(data532moy[i]/cste532/geom532_fit[i])


plt.figure()
plt.title('ICOS 532// '+t3[0])
plt.plot(data532cor[iav:],altitude[iav:],linewidth=1.,label="PR2/K")
plt.plot(rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh")
plt.ylabel("altitude (km)")
#plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()


#%%plot avant après

data532cor2=[]
for i in range(len(altitude)):
    data532cor2.append(data532moy[i]/geom532_fit[i])

"""   
plt.figure()
plt.title('532// '+time_moy_begin2[0]+'h-'+time_moy_end2[0]+'h UT')
plt.subplot(1,2,1)

plt.subplot(121)
plt.plot(data532moy[:300],altitude[:300],linewidth=1.,label="PR2")
plt.xscale('log')
plt.ylim(0,4)
plt.ylabel("altitude (km)")
plt.xlabel("PR2 (UA)")


plt.subplot(122)
plt.plot(data532cor2[iav:300],altitude[iav:300],linewidth=1.)
plt.xlabel("PR2 corrigé de la géométrie (UA)")
#plt.grid()
plt.ylim(0,4)
plt.xscale('log')
plt.legend()
plt.show()
"""

plt.figure()
plt.plot(data532cor2[iav:100],altitude[iav:100],linewidth=1.,label="PR2 corrigé de la géométrie")
plt.plot(data532moy[:100],altitude[:100],linewidth=1.,label="PR2")
plt.xlabel("(UA)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

#%% test pour géométrie


approx=[]
x=altitude[80:90]
y=data532moy[80:90]
coeff_poly=np.polyfit(x,y,1)
for i in range(0,80):
    approx.append(coeff_poly[0]*altitude[i]+coeff_poly[1])
"""
approx=[]
for i in range(0,70):
    approx.append(max(data532moy[:200]))
"""
plt.figure()
plt.title('532// '+time_moy_begin2[0]+'h-'+time_moy_end2[0]+'h UT')
plt.plot(data532moy[:100],altitude[:100],linewidth=1.,label="PR2")
plt.plot(approx,altitude[:80],label="apro")
plt.ylabel("altitude (km)")
plt.xlabel("voie 532")
plt.grid()
plt.legend()
plt.show()
 
#%%
c=[]
for i in range(0,60):
    #print(altitude[i])
    print(data532moy[i]/approx[i])
    c.append(data532moy[i]/approx[i])

datacor2=[]
for i in range(len(c)):
    datacor2.append(data532moy[i]/c[i])
    
for i in range(len(data532moy)-len(datacor2)):
    #print(i)
    datacor2.append(data532moy[i+60])

plt.figure()
plt.title('Vérification superpostition courbes')
plt.plot(datacor2[7:400],altitude[7:400],linewidth=1.,label="PR2 corrigé")
plt.ylabel("altitude (km)")
plt.xlabel("voie 1")
#plt.xscale('log')
plt.grid()
plt.legend()
plt.show()


c=c[::-1]
for i in range(len(c)):
    print(c[i])
"""  
for i in range(55,-1,-1):
    print(altitude[i])
"""
    

#%% data corrigé de géométrie seulement

data532cor2=[]
for i in range(len(altitude)):
    data532cor2.append(data532moy[i]/geom532_fit[i])


plt.figure()
plt.title('ICOS 532// '+t3[0])
plt.plot(data532cor2[15:],altitude[15:],linewidth=1.,label="PR2/K")
plt.ylabel("altitude (km)")
plt.grid()
plt.xlabel("voie1")
#plt.ylim(0,10)
#plt.xlim(-2.5,2.5)
plt.legend()
plt.show()

#%%Calcul inversion

izr=1965
#cste532=255500

beta532_aer=outils_OHP.calc_beta_aer_stable2(izr,altitude,rayleigh532,data532cor2,LRp,LRm,cste532)

beta532_aer_f=scipy.signal.medfilt(beta532_aer,9)


alpha532_aer,AOD532,trans532_aer=outils_OHP.calc_alpha_AOD_trans(beta532_aer_f,LRp,altitude[:izr])


#%%plot beta aer

plt.figure()
plt.title("ICOS 532// 20190704 "+t3[0])
plt.plot(beta532_aer[500:1230],altitude[500:1230])
plt.grid()
plt.xlabel("beta aer 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

print("à comparer avec 5.96*10**-5  ",max(beta532_aer[1000:]))

#%%valable seulement pour le 01092017 le soir

beta532_aer_f=list(beta532_aer_f)
i_pic=beta532_aer_f.index(max(beta532_aer_f[200:]))

i_pic_begin=1070
i_pic_end=1167

plt.figure()
#plt.title("ICOS 532// "+t3[0])
plt.plot(beta532_aer_f[i_pic-80:i_pic+80],altitude[i_pic-80:i_pic+80])
plt.grid()
plt.xlabel("beta aer 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

AOD_pic=np.trapz(alpha532_aer[i_pic_begin:i_pic_end],x=altitude[i_pic_begin:i_pic_end])
print("AOD_pic=",AOD_pic)

AOD_tot=np.trapz(alpha532_aer[iav:1400],x=altitude[iav:1400])
print("AOD_tot=",AOD_tot)


#%%plot alpha,T,AOD


plt.figure()
plt.title("ICOS 532// "+t3[0])
plt.plot(alpha532_aer[iav:izr],altitude[iav:izr])
plt.grid()
plt.xlabel("alpha aer 532// (km⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()


plt.figure()
plt.title("Transmission")
plt.plot(trans532_aer[iav:izr],altitude[iav:izr])
plt.grid()
plt.xlabel("transmission 532// aer")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

plt.figure()
plt.title("Epaisseur optique")
plt.plot(AOD532[iav:izr],altitude[iav:izr])
plt.grid()
plt.xlabel("AOD 532//")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()


#%% plot rapport de diffusion

altitude_ICOS=altitude

R532_ICOS=[]
for i in range(len(beta532_aer_f)):
    R532_ICOS.append(1.+beta532_aer_f[i]/rayleigh532['beta_m'][i])
    
plt.figure()
plt.plot(R532_ICOS,altitude[:izr])
plt.xlim(1,8)
plt.show()


