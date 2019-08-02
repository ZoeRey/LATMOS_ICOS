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
    
#%%

plt.figure()
plt.plot(time,AOD_500_LUNAR)
plt.title("AOD 500nm, LUNAR, 01092017")
plt.xlabel("time")
plt.ylabel("AOD 500")
plt.legend()
plt.show()

#%%Comparaison AOD LTA 20170901 19:44 22:42

AOD_LTA=np.mean(AOD_500[15:68])
print("AOD_LTA=",AOD_LTA)
