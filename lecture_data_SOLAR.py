import matplotlib.pyplot as plt
import numpy as np


#%%

trajet='/net/nfs/home/rey/Documents/data_OHP/photomètre/20170901_20170901_OHP_OBSERVATOIRE AOD LUNAR/data_AOD_LUNAR.txt'
file=np.loadtxt(trajet)
date=file[:,0]
time=file[:,1]
AOD_500=file[:,2]

#%%
hours=[]
minutes=[]
secondes=[]
for i in range(11,104):
    hours.append(int(time[i]/10000))
    minutes.append(int((time[i]-hours[i-11]*10000)/100))
    secondes.append(int(time[i]-(hours[i-11]*10000+minutes[i-11]*100)))

time=[]   
for i in range(len(hours)):
    time.append(hours[i]+((100./60.)*minutes[i]+secondes[i]/100)/100.)

AOD_500_LUNAR=AOD_500[11:104]
time_LUNAR=time
    
#%%

plt.figure()
plt.plot(time,AOD_500)
plt.title("AOD 500nm, LUNAR, 01092017")
plt.xlabel("time")
plt.ylabel("AOD 500")
plt.legend()
plt.show()

#%%

trajet='/net/nfs/home/rey/Documents/data_OHP/photomètre/20170901_20170901_OHP_OBSERVATOIRE AOD SOLAR/data_AOD_500_Solar.txt'
file=np.loadtxt(trajet)
time=file[:,0]
AOD_500=file[:,1]


hours=[]
minutes=[]
secondes=[]
for i in range(len(time)):
    hours.append(int(time[i]/10000))
    minutes.append(int((time[i]-hours[i]*10000)/100))
    secondes.append(int(time[i]-(hours[i]*10000+minutes[i]*100)))


time=[]   
for i in range(len(AOD_500)):
    time.append(hours[i]+((100./60.)*minutes[i]+secondes[i]/100)/100.)

AOD_500_SOLAR=AOD_500
time_SOLAR=time
   

plt.figure()
plt.plot(time,AOD_500)
plt.title("AOD 500nm, SOLAR, 01092017")
plt.xlabel("time")
plt.ylabel("AOD 500")
plt.legend()
plt.show()

#%%plot solar lunar

plt.figure()
plt.plot(time_SOLAR,AOD_500_SOLAR,label="SOLAR")
plt.plot(time_LUNAR,AOD_500_LUNAR,label="LUNAR")
plt.title("AOD 500nm, SOLAR+LUNAR, 01092017")
plt.xlabel("time (UT)")
plt.ylabel("AOD 500nm")
plt.legend()
plt.show()


