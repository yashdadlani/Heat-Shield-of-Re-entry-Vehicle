import numpy as np
import matplotlib.pyplot as plt
s = 5.67*(10**(-8))  #boltzmann constant
e = 0.9 #emissivity
Li = 0
Lf = 0.1
dxp = [0.0025,0.004,0.005,0.01,0.0125,0.025,0.05]
plots = []

for l in range(0,7):
 dx = dxp[l]
 m = int((Lf-Li)/dx) + 1                    # no of grid points: discretization
 dt = 1
 n = int(720/dt) + 1
 K = 1                        # Reinforced Carbon Carbon
 rho = 1600
 cp = 500
 alpha = K/(rho*cp)
 r = (alpha*dt) / (dx**2)
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
#  plots.append(T[:,100])
#  plots.append(T[:,200])
#  plots.append(T[:,300])
#  plots.append(T[:,450])
#  plots.append(T[:,550])
#  plots.append(T[:,650])
#  plots.append(T[:,720])
 plots.append(T[m-1,:])
#  max_Temp  = np.amax(T)
# pos_time = np.where(T==max_Temp)
#  print(max_Temp)
# print(pos_time)


#  plt.plot(time,T[0,:])
#  plt.xlabel('Time(in s)')
#  plt.ylabel('Temp (in K)')
#  plt.title('Grid Test')

#  plt.grid()
#  plt.show()
A, = plt.plot(time,plots[0],color = 'red')
B, = plt.plot(time,plots[1], color = 'green')
C, = plt.plot(time,plots[2], color = 'orange')
D, = plt.plot(time,plots[3], color = 'violet')
E, = plt.plot(time,plots[4], color = 'purple')
F, = plt.plot(time,plots[5], color = 'cyan')
G, = plt.plot(time,plots[6], color = 'magenta')

plt.xlabel('time(in s)')
plt.ylabel('Temp (in K)')
plt.legend((A,B,C,D,E,F,G),['40 div','25 div','20 div ', '10 div','8 div', '4 div', '2 div'])
plt.title('Grid Independence: for different number of elements')
plt.grid()
plt.show()
