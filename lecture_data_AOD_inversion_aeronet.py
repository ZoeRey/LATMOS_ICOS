import matplotlib.pyplot as plt
import numpy as np

#%%


trajet='/net/nfs/home/rey/Documents/data_OHP/photomètre/20170901_20170901_OHP_OBSERVATOIRE AEROSOL INVERSION/data_AOD_aeronet_inversion.txt'


file=np.loadtxt(trajet,skiprows=1)
time=file[:,0]
AOD_ET_440=file[:,1]
AOD_ET_675=file[:,2]
AOD_EF_440=file[:,3]
AOD_EF_675=file[:,4]

hours=[]
minutes=[]
secondes=[]
for i in range(len(time)):
    hours.append(int(time[i]/10000))
    minutes.append(int((time[i]-hours[i]*10000)/100))
    secondes.append(int(time[i]-(hours[i]*10000+minutes[i]*100)))


time=[]   
for i in range(len(AOD_ET_440)):
    time.append(hours[i]+((100./60.)*minutes[i]+secondes[i]/100)/100.)
    
#%%

plt.figure()
plt.plot(time,AOD_ET_440,label="440nm")
plt.plot(time,AOD_ET_675,label="675nm")
plt.title("AOD Extinction-Total Aeronet Inversion, 01092017")
plt.xlabel("time")
plt.ylabel("AOD")
plt.legend()
plt.show()
