import numpy as np
import matplotlib.pyplot as plt

#Variables
v1 = 20
v2 = 24
T = 20
L = 200
dt = 0.001
dx = 1
Nt = int(round(T/dt))
t = np.linspace(0,T, Nt+1)
c = 5
r = 7
a = v1*dt/dx
b = v2*dt/dx

Nx = int(round(L/dx))
x = np.linspace(0,L,Nx+1)
dx = x[1] - x[0]
dt = t[1] - t[0]

u1 = np.zeros(Nx+1)
u2 = np.zeros(Nx+1)
u1_n = np.zeros(Nx+1)
u2_n = np.zeros(Nx+1)

U1 = np.zeros((Nt,Nx+1))
U2 = np.zeros((Nt,Nx+1))

u1_n[0] = 10
u1_n[1] = 10
u1_n[Nx] = 0


#Solve equations
for n in range(0,Nt):
    for i in range(1,Nx):
        u1[i] = u1_n[i] - a*(u1_n[i] - u1_n[i-1]) - dt*(c*u1_n[i] - r*u2_n[i])
        u2[i] = u2_n[i] + b*(u2_n[i+1] - u2_n[i]) + dt*(c*u1_n[i] - r*u2_n[i])
        U1[n,:] = u1 
        U2[n,:] = u2 
        #Boundary conditions
        u1[0] = u1_n[0] - dt*(c*u1_n[0] - r*u2_n[0])
        u2[Nx] = u2_n[Nx] + dt*(c*u1_n[Nx] - r*u2_n[Nx]) 
   
    u1_n, u1  = u1, u1_n 
    u2_n, u2  = u2, u2_n 

#Plot solutions
    
for nn in range(0, Nt, Nt//10):
   print(U1[nn,:])
   plt.plot(x, U1[nn,:])      
   plt.ylabel('Density of growing microtubules (ρ_g)', fontsize=14)
   plt.ylim(0,0.5)
   plt.xlabel('length (l)', fontsize = 14)
   
   


for n in range(0,Nt):
    print('u1', U1[n,:])
    print('u2', U2[n,:])

      
   


    


                       


