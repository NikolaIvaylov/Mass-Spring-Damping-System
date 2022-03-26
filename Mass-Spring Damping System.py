import math
from numpy import linspace,sin,pi,int16
import numpy as np
from pylab import plot,show,axis
from matplotlib import pyplot as plt

#Variables you can play with:
#Time:
t_step = 0.01 #s
n = 200 #number of steps
#Mass-Spring System:
u0 = 0 #initial displacement (m)
v0 = 0 #initial velocity (m/s)
f0 = -3 #initial force F0 (N)
m = 3 #mass (kg)
b = 0.5 #damping constant (kg/s)
k = 1 #spring constant (N/m)
    
#Other variables (that you cannot play with):
#Time:
t = np.linspace(0, (n*t_step), n) #generate x-axis values
#Mass-Spring System:
displacements = [] #list for all displacements in time
velocities = [] #list for all velocities in time
displacements.append(u0) #set initial displacement
velocities.append(v0) #set initial velocity

#function for acceleration y'' (m/s^2):
def func(t, u, v):
        return ((f0 - b*v - k*u) / m) * t
    
#Runge-Kutta Method:
def rungeKutta(t0, u, v, t_step, tn):
    #Count number of iterations using step size:
    n = (int)((tn - t0)/t_step)
    
    #Iterate for n number of iterations:
    for i in range(1, n + 1):
        #Apply Runge-Kutta Formulas to find next values of k1-k4:
        k1 = (1 / 2) * (t_step**2) * func(t0, u, v)
        k2 = (1 / 2) * (t_step**2) * func(t0+(1/2)*t_step, u+(1/2)*t_step*v+(1/4)*k1, v+k1/t_step)
        k3 = (1 / 2) * (t_step**2) * func(t0+(1/2)*t_step, u+(1/2)*t_step*v+(1/4)*k2, v+k2/t_step)
        k4 = (1 / 2) * (t_step**2) * func(t0+t_step, u+t_step*v+k3, v+(2*k3)/t_step)
        
        #Apply Formulas to find next values of P and Q:
        P = (1 / 3) * (k1 + k2 + k3)
        Q = (1 / 3) * (k1 + 2*k2 + 2*k3 + k4)
        
        t0 = t0 + t_step #Update next value of t
        u = u + t_step*v + P #Update the value of u
        v = v + Q/t_step #Update the value of v
        
    return u, v  #return the new displacement and velocity


#Runge-Kutta Loop:
for time in range(1, len(t)):
    previous_displacement = displacements[time-1] #get y(n-1)
    previous_velocity = velocities[time-1] #get y'(n-1)
    
    #Use the rungeKutta() function:
    current_displacement, current_velocity = rungeKutta(0, previous_displacement, previous_velocity, t_step, t[time])
    
    displacements.append((current_displacement-f0)) #save the value of the current displacement
    velocities.append(current_velocity) #save the value of the current velocity

#Displaying the graph:
plt.figure(dpi=300) 
plt.title(f"Runge-Kutta Method - Slides (step={t_step})")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.plot(t, displacements, color="blue")
plt.plot(t, np.zeros(len(t)) , color="gray")
plt.show()
