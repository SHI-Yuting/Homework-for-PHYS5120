import numpy as np
import matplotlib.pyplot as plt
import math
#define the parameter
g=9.81
l=1.2

# Input initial angle in degrees and convert to radians
deg=input("Please enter the initial angle in degrees: ")
theta_0=math.radians(float(deg))
omega_0=0

#time step
h=float(input("Please enter the time step in seconds: "))
t_max=40

#initial
time=np.arange(0,t_max,h)
theta=np.zeros(len(time))
omega=np.zeros(len(time))
theta[0]=theta_0
omega[0]=0.0

# Define functions for acceleration, energy, kinetic energy, and potential energy
def accelaration(x):
    return -g/l*np.sin(x)
def energy(x,v):
    return 0.5*l*l*v*v+g*l*(1-np.cos(x))
def Kinetic(v):
    return 0.5*l*l*v*v
def V(x):
    return g*l*(1-np.cos(x))

#velocity Verlet
for n in range(1,len(time)):
    theta[n]=theta[n-1]+h*omega[n-1]+0.5*accelaration(theta[n-1])*h*h
    omega[n]=omega[n-1]+0.5*(accelaration(theta[n])+accelaration(theta[n-1]))*h

#calculate energy
energy_cal=energy(theta,omega)
Kinetic_num=Kinetic(omega)
V_num=V(theta)

# Calculate analytical solutions for comparison
theta_linear=theta_0*np.cos(np.sqrt(g/l)*time)
omega_linear=-np.sqrt(g/l)*theta_0*np.sin(np.sqrt(g/l)*time)
energy_theory=np.full(len(time),energy(theta[0],omega[0]),dtype=float)
Kinetic_theory=np.full(len(time),Kinetic(omega_linear),dtype=float)
V_theory=np.full(len(time),V(theta_linear),dtype=float)
appli=np.full(len(time),theta_0)

# Plotting the results
plt.xlabel("time(s)",fontsize=20)
plt.ylabel("angle(rad)",fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.plot(time,appli,linewidth=5,linestyle='--',color='green',label='amplitude')
plt.plot(time,-appli,linewidth=5,linestyle='--',color='green')
plt.plot(time,theta,linewidth=5,label=r'numerical $\theta$')
plt.plot(time,theta_linear,linewidth=5,label=r"analytical $\theta$")
plt.legend(loc=1,fontsize=20)
plt.show()

plt.xlabel("time(s)",fontsize=20)
plt.ylabel("Angular velocity(rad/s)",fontsize=20)
plt.plot(time,omega,linewidth=5,label=r'numerical $\omega$')
plt.plot(time,omega_linear,linewidth=5,label=r"analytical $\omega$")
plt.legend(loc=1,fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

plt.xlabel("time(s)",fontsize=20)
plt.ylabel(r"numerical $\theta$-analytical $\theta$ (rad)",fontsize=20)
plt.plot(time,theta-theta_linear,linewidth=5)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

plt.xlabel("time(s)",fontsize=20)
plt.ylabel(r"numerical $\omega$-analytical $\omega$ (rad/s)",fontsize=20)
plt.plot(time,omega-omega_linear,linewidth=5)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.show()

plt.ylabel("energy per kg(J/kg)",fontsize=20)
plt.xlabel("time(s)",fontsize=20)
plt.plot(time,energy_cal,linewidth=5,label="numerical energy")
plt.plot(time,energy_theory,linewidth=5,label="initial energy")
plt.legend(loc=2,fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
#plt.ylim(17.47,17.49)
plt.show()

plt.ylabel("kinetic energy per kg(J/kg)",fontsize=20)
plt.xlabel("time(s)",fontsize=20)
plt.plot(time,energy_theory,linewidth=5,label="initial energy")
plt.plot(time,Kinetic_num,linewidth=5,label="Numerical kinetic energy")
#plt.plot(time,Kinetic_theory,linewidth=5,label="Analytical kinetic energy")
plt.legend(loc=2,fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()

plt.ylabel("potential energy per kg(J/kg)",fontsize=20)
plt.xlabel("time(s)",fontsize=20)
plt.plot(time,energy_theory,linewidth=5,label="initial energy")
plt.plot(time,V_num,linewidth=5,label="Numerical potential energy")
#plt.plot(time,V_theory,linewidth=5,label="Analytical potential energy")
plt.legend(loc=2,fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()

plt.ylabel("energy per kg(J/kg)",fontsize=20)
plt.xlabel("time(s)",fontsize=20)
plt.plot(time,V_num,linewidth=5,label="Numerical potential energy")
plt.plot(time,Kinetic_num,linewidth=5,label="Numerical kinetic energy")
plt.plot(time,energy_theory,linewidth=5,label="initial energy")
plt.legend(loc=2,fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.show()
# Print final values of theta and omega
#print(theta)
#print(omega)

