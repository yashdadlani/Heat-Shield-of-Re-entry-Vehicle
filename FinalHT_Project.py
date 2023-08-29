import numpy as np
import matplotlib.pyplot as plt
s = 5.67*(10**(-8))  #boltzmann constant
e = 0.9              #emissivity
K = float(input("Enter thermal conductivity: "))                       #1               
rho = float(input("Enter density: "))            #1600
cp = float(input("Enter specific heat: "))                       #500
Li = 0
Lf = float(input("Enter desired thickness in m : "))    #0.1
dx =  float(input("Enter step size: "))                                                  #0.01
m = int((Lf-Li)/dx) + 1 # no of grid points: discretization
dt =   float(input("Enter time step: "))                                                  #1
n = int(720/dt) + 1

alpha = K/(rho*cp)
r = (alpha*dt) / (dx**2)
if r <= 0.5 :
 x = np.linspace(Li,Lf, m)
 time = np.linspace(0,720,n)
 heat_flux = np.empty(n)
 for a in range(0,201):
    heat_flux[a] = 4000*time[a] + 200000
 for b in range(201,551):
    heat_flux[b] = 1000000
 for c in range(551,651):
    heat_flux[c] = (496*(10**4) -7200*time[c])
 for d in range(651,721):
    heat_flux[d] = 2.8*(10**5)

 T = np.zeros((m,n))
 T[:,0] = 300 #initial condition
 for j in range(0,n-1):
     for i in range(1,m-1):
         T[i, j+1] = (r*T[i+1,j]) + ((1-(2*r))*T[i,j]) + r*(T[i-1,j])
     T[m-1,j+1] = T[m-2,j+1] #boundary condition
     T[0,j+1] = T[0,j] + (dt/(rho*cp))*((heat_flux[j]) + (K*(T[1,j] - T[0,j]))  -((e*s)*((T[0,j])**4)))
 max_Temp  = np.amax(T)
 pos_time = np.where(T==max_Temp)
 print(max_Temp)
 print(pos_time)
 plt.plot(T[m-1,:])
 plt.grid()
 plt.show()
 print(T[m-1,:])
else:
   print("Enter smaller time step once you re-run the code ")
   

