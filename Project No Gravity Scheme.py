#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ Project Ma411 ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤


# Common parameters :
#
# Variables :
import numpy as np
import matplotlib.pyplot as plt
a,b=0,150 #Beginning and end of the space model
alpha=0.1472#Initial conditions
N=100 #points
h=(b-a)/N
X=np.linspace(a,b,N) #Generate N points between a and b
T=50 #boucles de temps (secondes ?)
k=h/(3*alpha) #lié à la CFL
r=k/h # dx/dt
S1=np.zeros((len(X),T)) #init saturation vector for upwind
S2=S1 # for LW
S3=S1 # for Godunov
S1[0,:]=1 #boundary conditions
S2[0,:]=S1[0,:]
S3[0,:]=S1[0,:]
U=np.zeros(N) #for Godunov
U[0]=1 #Initial conditions also apply !

#Functions :
def fw(s):
    return s**2

def fo(s):
    return ((1-s)**2)/4

def f(s):
    return alpha*fw(s)/(fw(s)+fo(s))


# Loop to store values of saturation :

#Upwind :
for n in range(0,T-1):
    for i in range(1,len(S1)-1):
       S1[i,n+1] = S1[i,n] - r*(f(S1[i,n])-f(S1[i-1,n]))

#Lax-Friederich:
for n in range(0,T-1):
    for i in range(1,len(S2)-1):
       S2[i,n+1] = (0.5*(S2[i+1,n]+S2[i-1,n]))-(0.5*r*(f(S2[i+1,n])-f(S2[i-1,n])))

#Godunov:
for n in range(0,T-1):
    for i in range(1,len(S3)-1):
        U[i]=S3[i,n]-r*(f(S3[i,n])-f(S3[i-1,n]))
        S3[i,n+1] = 0.5*(S3[i,n]+U[i]-r*(f(U[i])-f(U[i-1])))


# Plotting the curve :
plt.clf()
plt.plot(X,S1[:,-1],label='Upwind scheme')
plt.plot(X,S2[:,-1],label='Lax-Friederich scheme')
plt.plot(X,S3[:,-1],label='Godunov scheme')
plt.xlabel('length (m) on x-direction')
plt.ylabel('water saturation')
plt.title('Comparaison of numerical scheme for petroleum recovery without gravity')
plt.legend()
plt.show()