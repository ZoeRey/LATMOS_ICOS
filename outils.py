#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 09:55:35 2019

@author: rey
"""
import numpy as np
import statistics
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt


#%%Calcul constante première itération

def calc_cste(datamoy,rayleigh,indice3,indice4):
    cste=statistics.mean(datamoy[indice3:indice4]/(rayleigh['beta_m'][indice3:indice4]*rayleigh['tm'][indice3:indice4]))   
    if cste<0:
        cste=-1*cste
    print(cste)
    return cste


#%%
def calc_beta_aer_stable2(izr,altitude,rayleigh,datamoy,LRp,LRm,cste):
    
    a=[]
    b=[]
    beta_p=[]
    for i in range(0,izr):
        a.append(np.trapz(rayleigh['beta_m'][i:izr],x=altitude[i:izr]))
        b.append(datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))
    c=[np.trapz(b[i:izr],x=altitude[i:izr]) for i in range(0,izr)]
    for i in range(0,izr):
        #beta_p.append(b[i]/(datamoy_zr/rayleigh['beta_m'][izr]+2.*LRp*c[i])-rayleigh['beta_m'][i])
        beta_p.append(b[i]/(cste*rayleigh['tm'][izr]+2.*LRp*c[i])-rayleigh['beta_m'][i])
    return beta_p
    
    """
    altitude2=altitude[:izr]
    altitude3=altitude2[::-1]
    ray=rayleigh['beta_m'][:izr]
    ray2=ray[::-1]
    datamoy2=datamoy[:izr]
    datamoy3=datamoy2[::-1]

    a=[]
    b=[]
    beta_p=[]
    for i in range(len(ray)):
        a.append(np.trapz(ray2[:i],x=altitude3[:i]))
        b.append(datamoy3[i]*np.exp(-2.*(LRp-LRm)*a[i]))
    #print(a)
    c=[np.trapz(b[:i],x=altitude3[:i]) for i in range(len(ray))]
    for i in range(len(ray)):
        #beta_p.append(b[i]/(datamoy_zr/rayleigh['beta_m'][izr]-2.*LRp*c[i])-ray2[i])
        beta_p.append(b[i]/(cste*rayleigh['tm'][izr]-2.*LRp*c[i])-ray2[i])
    beta_p=beta_p[::-1]
    return beta_p
    """
#%%calcul alpha et transmission
    
def calc_alpha_AOD_trans(beta,LR,altitude):
    alpha=[]
    trans=[]
    for i in range(len(beta)):
        alpha.append(beta[i]*LR)
    AOD=[np.trapz(alpha[:i],x=altitude[:i]) for i in range(len(altitude))]
    for j in range(len(AOD)):
        trans.append(np.exp(-1*AOD[j]))
        
    return alpha,AOD,trans
    


















