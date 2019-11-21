# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 09:16:26 2019

@author: pfm
"""

import numpy as np
from matplotlib import pyplot as plt
from math import pi



def get_elipse_dataset(n=10, f0=1, fe=100, scale_noise=1, form_variability=True):
    ds=[]
    Ti=[]
    np.random.seed(42)

    u=1.     #x-position of the center
    v=0.5    #y-position of the center
    a0=2.     #radius on the x-axis
    b0=1.    #radius on the y-axis
    phi0=0
    L=int(2*fe/f0)
    
    t0 = np.linspace(0, 2/f0, L)
    dt=2/f0/L
    t=[0]
    T=0
    a=a0
    b=b0
    for i in range(L-1):
        inc=dt#*(np.random.rand()/8+.85)
        T+=inc
        t.append(T)
    t=np.array(t)
    for i in range(n):
        if form_variability:
            a=a0*(1+np.random.uniform(-.5,.5))
            b=b0*(1+np.random.uniform(-.5,.5))
        f=f0*(1+np.random.uniform(-.1,.1))
        phi=phi0+np.random.uniform(-.25,.25)
        L=int(2*fe/f)
        t = np.linspace(0, 2/f, L)
        dt=2/f0/L
        T=0

        noise=(np.random.rand(2,L)-.5)*scale_noise
        noise=(np.random.normal(size=(2,L)))*scale_noise

        ti=2*pi*f*t+phi
        x=u+a*np.cos(ti)+noise[0,:]
        y=v+b*np.sin(ti)+noise[1,:]
        plt.plot(x,y,'c')
        s=np.array([x,y])
        ds.append(s.T)
        Ti.append(ti)
    L=int(2*fe/f0)
    t = np.linspace(0, 2/f0, L)
    x=u+a0*np.cos(2*pi*f0*t+phi0)
    y=v+b0*np.sin(2*pi*f0*t+phi0)
    C=np.array([x,y])
    T0=2*pi*f0*t+phi0
    #plt.grid(color='lightgray',linestyle='--')
    #plt.show()
    return ds, np.array(Ti), C.T, T0
    
def get_elipse_dataset_df(n=10, f0=1., delta_f0=.1, fe=100, scale_noise=1, form_variability=True):
    ds=[]
    T=[]
    u=1.     #x-position of the center
    v=0.5    #y-position of the center
    a0=2.     #radius on the x-axis
    b0=1.    #radius on the y-axis
    a=a0
    b=b0
    for i in range(n):
        if form_variability:
            a=a0*(1+np.random.uniform(-.5,.5))
            b=b0*(1+np.random.uniform(-.5,.5))
        x=[]
        y=[]
        t=0
        f=f0
        while f*t<=2.:
            T.append(2*pi*f*t)
            x.append(u+a*np.cos(2*pi*f*t)+(np.random.normal())*scale_noise)
            y.append(v+b*np.sin(2*pi*f*t)+(np.random.normal())*scale_noise)
            t+=1/fe
            f+=delta_f0
            if f <.25:
                f=.25
    x=np.array(x)
    y=np.array(y)
    plt.plot(x,y)
    s=np.array([x,y])
    ds.append(s.T)
    #plt.grid(color='lightgray',linestyle='--')
    #plt.show()
    return ds, np.array(T)


def get_elipse_dataset_fm(n=10, f0=1., delta_f0=.01, fe=100, scale_noise=1, form_variability=True):
    ds=[]
    T=[]
    u=1.     #x-position of the center
    v=0.5    #y-position of the center
    a0=2.     #radius on the x-axis
    b0=1.    #radius on the y-axis
    a=a0
    b=b0
    for i in range(n):
        if form_variability:
            a=a0*(1+np.random.uniform(-.5,.5))
            b=b0*(1+np.random.uniform(-.5,.5))
        x=[]
        y=[]
        t=0
        f=f0
        phi=0
        while phi<=4*pi+.1:
            T.append(phi)
            phi+=2*pi*f/fe
            x.append(u+a*np.cos(phi)+(np.random.normal())*scale_noise)
            y.append(v+b*np.sin(phi)+(np.random.normal())*scale_noise)
            t+=1/fe
            f=f0+np.sin(2*pi*1.7*f0*t)*delta_f0
            #print(f, end=' ')
    x=np.array(x)
    y=np.array(y)
    plt.plot(x,y)
    s=np.array([x,y])
    ds.append(s.T)
    #plt.grid(color='lightgray',linestyle='--')
    #plt.show()
    return ds, np.array(T)

