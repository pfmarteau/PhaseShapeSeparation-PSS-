import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
from ekats2 import PyeKats2
import ellipse as el
import sys

EKATS = PyeKats2()

#
# Given set S of time series, returns the index within S
# of the time series with the mean length
def getMeanLengthId(S):
    L=[]
    for n in range(len(S)):
        L.append(len(S[n]))
    lid=np.argsort(L)
    mean=np.mean(L)
    mn=1e10
    nmin=0
    for n in range(len(S)):
        if np.abs(L[n]-mean)<mn:
            mn=np.abs(L[n]-mean)
            nmin=n
    mdl=lid[int(np.ceil(len(S)/2))]
    print(nmin)
    return nmin

#
# get the inertia for estimate ts (sum of the pairwise kernel evaluations between ts and all the elements of S)
def get_kdtw_inertia(ts, S, nu):
    inertia = 0
    for i in range(len(S)):
        inertia = inertia + EKATS.kdtw(ts, S[i], nu)
    return inertia

#
# get the estimated xi function given the estimation of the centroid C
# s is a time series that is aligned on C to provide xi
def alignCentroid(s, C, nu=1):
    Y=EKATS.iTEKA_stdev(C, [s], nu)
    X=np.array(list(Y[0]))
    E=np.array(list(Y[2]))
    errv=E[:,0]/len(C)
    xi=X[:,0]
    return xi, errv

#
# get the iTEKA centroid estimate for the set of time series S
# nu is the smoothing parameter, 
# c0 is the current estimate (e.g. the medoid of S)    
# npass is the max number of itertions 
# OUTPUTS: Cp: the centroid estimate temporaly re-interpolated
#         CCp: the centroid estimate (without temporal interpolation)
#         xip: the mean xi associated to CCp
#         XIp: the estimated xi for each element of S  
def get_iTEKACentroid(S, c0, nu=1, npass=100):
    ii = 0
    inertiap = 1e10
    Y=EKATS.iTEKA_stdev(c0, S, nu)
    XI=Y[1]
    XIp=XI
    X=np.array(list(Y[0]))
    Tstd=X[:,len(X[0])-1]
    Tstdp=Tstd
    xi=X[:,0]
    xip=xi
    X=X[:,0:len(X[0])-1]
    C = EKATS.interpolate(X)
    Cp=C
    dim=len(S[0][0])
    print('DIM=',dim)
    CC=X[:,1:dim+1]
    C0=C[:,0:dim]
    CCp=CC
    inertia=Y[3][0][0]
    n=0
    print(n, "inertia: ", inertia)
    while (not np.isnan(inertia)) and (ii < npass) and (inertia < inertiap):
        inertiap = inertia
        Cp = C
        CCp=CC
        xip=xi
        XIp=XI
        Y=EKATS.iTEKA_stdev(Cp[:,0:dim], S, nu)
        XI=Y[1]
        X=np.array(list(Y[0]))
        E=np.array(list(Y[2]))
        errv=E[:,0]/len(Cp)

        xiTrue=X[:,0]
        X=X[:,0:len(X[0])-1]
        C = EKATS.interpolate(X)

        C0=C[:,0:dim]
        CC=X[:,1:dim+1]

        inertia=Y[3][0][0]

        n+=1
        if not np.isnan(inertia):
           print(n, "inertia: ", inertia)
        ii = ii + 1
    return np.copy(Cp), np.copy(CCp), xip, inertiap, XIp


# Generate the data set
Fe=400
scale_noise=1.5
# 
S, xi, Ctrue, Ttrue = el.get_elipse_dataset(n=20, f0=1, fe=Fe, scale_noise=scale_noise)
print(len(S), "ts loaded")
plt.plot(Ctrue[:,0],Ctrue[:,1], '--r', linewidth=4)


# set the smoothing parameter nu
nu=.1
# use the time series that have the mean length as the initial centroid estimate
idc=getMeanLengthId(S)
# get the iTEKA centroid and the xi (temporal) functions 
# (one for each time series within S)
C, CC, T, inertia, TTp = get_iTEKACentroid(S, S[idc], nu=nu)
plt.plot(C[:,0],C[:,1], 'k', linewidth=4)
plt.savefig("images/dataset.pdf")

# Plot the iTEKA centroid estimate and the true centroid
plt.figure(2)
C=C[:,0:2]
plt.plot(C[:,0],C[:,1], 'k', linewidth=4)
plt.plot(Ctrue[:,0],Ctrue[:,1], '--r', linewidth=4)
plt.savefig("images/centroids.pdf")
plt.show()


# Plot the True xi functions
plt.figure(100)
for i in range(len(xi)):
    plt.plot(xi[i][:], 'b', linewidth=1)
plt.plot(Ttrue, 'r', linewidth=1)
plt.title(r"Fonctions $h_i$ réelles")
plt.grid()
plt.savefig("images/true_hi_ex1.pdf")
plt.show()


# Plot the estimated xi functions
L=len(TTp[0])
TT=np.zeros((L, len(TTp)))
for n in range(L):
    for j in range(len(TTp)):
        TT[n][j]=TTp[j][n]
# xi estimated from the outputs of iTEKA
plt.figure(101)
for i in range(L):
    plt.plot(TT[i][:]*max(Ttrue)/len(TT[i]), 'b', linewidth=1)
plt.plot(Ttrue, 'r', linewidth=1)
plt.title(r"Fonctions $\hat{h}_i$ estimées")
plt.grid()
plt.savefig("images/estimated_hi_ex1.0.pdf")
# xi restimated through the alignment with the estimated centroid
plt.figure(102)
for n in range(len(S)):
    T00, err00 = alignCentroid(S[n], CC, nu=nu)
    plt.plot(T00[:]*max(Ttrue)/len(T00),'b')
plt.plot(Ttrue, 'r', linewidth=1)
plt.title(r"Fonctions $\hat{h}_i$ estimées")
plt.grid()
plt.savefig("images/estimated_hi_ex1.1.pdf")
plt.show()

FV=True

# generate xiTrue (f0=1)
S00, xiTrue, ck00, Tk00 = el.get_elipse_dataset(n=1, f0=1, fe=Fe, scale_noise=scale_noise, form_variability=FV)

# generate a noisy ellipse with a xi function twice slower than the expected xi (f0=.5)
S0, xi0, ck0, Tk0 = el.get_elipse_dataset(n=1, f0=.5, fe=Fe, scale_noise=scale_noise, form_variability=FV)
T0,err0=alignCentroid(S0[0], CC, nu=nu/2)
# generate a noisy ellipse with a xi function twice faster than the expected xi (f0=2)
S1, xi1, ck1,Tk1 = el.get_elipse_dataset(n=1, f0=2, fe=Fe, scale_noise=scale_noise, form_variability=FV)
T1,err1 = alignCentroid(S1[0], CC, nu=nu)
# generate a noisy ellipse with a xi function that accelerates compared to the expected xi (delta_f0>0)
S2, xi2 = el.get_elipse_dataset_df(n=1, f0=1, delta_f0=.0005, fe=Fe, scale_noise=scale_noise, form_variability=FV)
T2,err2 = alignCentroid(S2[0], CC, nu=nu)
# generate a noisy ellipse with a xi function that decelerates compared to the expected xi (delta_f0<0)
S3, xi3 =el.get_elipse_dataset_df(n=1, f0=1, delta_f0=-.0002, fe=Fe, scale_noise=scale_noise, form_variability=FV)
T3,err3 = alignCentroid(S3[0], CC, nu=nu)
# generate a noisy ellipse with a xi function that accelerate and decelerates (frequency modulation with f varying sinusoidally in [.1;1.9]) compared to the expected xi (delta_f0>0)
S4, xi4 = el.get_elipse_dataset_fm(n=1, f0=1, delta_f0=.9, fe=Fe, scale_noise=scale_noise, form_variability=FV)
T4,err4 = alignCentroid(S4[0], CC, nu=nu)

# plot the true xi functions
plt.figure(6)
y=list(range(len(xiTrue[0])))
plt.plot(xiTrue[0]/np.max(xiTrue[0])*len(Ttrue),y, 'r', label='mean xi')
y=list(range(len(xi0[0])))
plt.plot(xi0[0]/np.max(xi0[0])*len(Ttrue),y,'k', label='1/2 slower')
y=list(range(len(xi1[0])))
plt.plot(xi1[0]/np.max(xi1[0])*len(Ttrue),y,'c', label='2 faster')
y=list(range(len(xi2)))
plt.plot(xi2/np.max(xi2)*len(Ttrue),y,'--b', label='acceleration')
y=list(range(len(xi3)))
plt.plot(xi3/np.max(xi3)*len(Ttrue),y,'--b', label='deceleration')
y=list(range(len(xi4)))
plt.plot(xi4/np.max(xi4)*len(Ttrue),y,'-.g', label='freq. modulation')
plt.grid()
plt.legend()
plt.title(r"True $\xi_i$")
plt.savefig("images/true_xi_ex2.pdf")
# plot the estimated xi functions
plt.figure(7)
plt.plot(T,'r',label='mean xi')
plt.plot(T0,'k', label='1/2 slower')
plt.plot(T1,'c', label='2 faster')
plt.plot(T2,'--b', label='acceleration')
plt.plot(T3,'--b', label='deceleration')
plt.plot(T4,'-.g', label='freq. modulation')
plt.grid()
plt.legend()
plt.title(r"Estimated $\xi_i$")
plt.savefig("images/estimated_xi_ex2.pdf")
plt.show()


