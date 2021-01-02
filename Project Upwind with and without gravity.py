# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:53:17 2020

@author: Tristan
"""

# -*- coding: utf-8 -*-



import numpy as np

import matplotlib.pyplot as plt



a,b=0,100

alpha=0.1472   #in m3/s value found on suez website and used on real installation

N=100 #[0,1]

h=(b-a)/N

X=np.linspace(a,b,N)

T=50 #boucles de temps      T*k = temps reel

k=h/(3*alpha) #lié à la CFL

r=k/h

S=np.zeros((len(X),T)) #init saturation vector

S[0,:]=1 #boundary conditions





def fwater(s):

    return s**2



def foil(s):

    return ((1-s)**2)/4



def f(s):

    return alpha*fwater(s)/(fwater(s)+foil(s))



for n in range(0,T-1):

    for i in range(1,len(S)-1):

       S[i,n+1] = S[i,n] - r*(f(S[i,n])-f(S[i-1,n]))





K=41*(10**(-12)) #permeabilite

g=9.81 #accélération de la pensateur, m/s²

S1=np.zeros((len(X),T)) #init saturation vector

S1[0,:]=1 #boundary conditions

rho_w,rho_o = 1000,900  # in kg/m^3

phi=0.35

beta = (rho_w - rho_o)*(K/phi)*g





def G(a,b):

    if -alpha+beta*fwater(a)<=0:

        G=(fwater(a)*(alpha+beta*foil(a)))/(fwater(a)+foil(a))

    else:

        G=(fwater(a)*(alpha+beta*foil(b)))/(fwater(a)+foil(b))

    return(G)





for n in range(0,T-1):

    for i in range(1,len(S)-1):

       a=S1[i,n]

       b=(k/h)

       c=G(S1[i,n],S1[i+1,n])-G(S1[i-1,n],S1[i,n])

       S1[i,n+1] = a - (b*c)



plt.clf()

plt.plot(X,S[:,-1])

plt.plot(X,S1[:,-1])

plt.plot

plt.xlabel('x(m)')

plt.ylabel('s_water(x)')

plt.legend(["without gravity",'with gravity'])

plt.title('Upwind-scheme')

plt.show()