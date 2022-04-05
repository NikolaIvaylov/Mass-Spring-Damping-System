import math
from numpy import linspace,sin,pi,int16
import numpy as np
from pylab import plot,show,axis
from matplotlib import pyplot as plt

#Variables you can play with:
#Time:
t_step = 0.01 #s
n = 8000 #number of steps
#Mass-Spring System:
u0 = 0 #initial displacement (m)
v0 = 0 #initial velocity (m/s)
f0 = 3 #initial force F0 (N)
m = 3 #mass (kg)
b = 0.5 #damping constant (kg/s)
k = 1 #spring constant (N/m)

#generate a force pulse:
force = np.zeros(n) #list for force pulse values
for i in range(0, 100):
    a = np.pi / 100 * (i - 200)
    force[i] = f0*np.sin(a)
    
t = np.linspace(0, (n*t_step), n) #generate x-axis values

def rungeKutta(t, u0, v0, m, b, k, force):
    """
    t: list of all values of time (x-axis)
    u0: the displacement at t[0]
    v0: the velocity at t[0]
    m: the mass
    b: damping constant
    k: spring constant
    force: list of the values of the force pulse
    
    returns all displacements u, and velocities v in time
    """

    u = np.zeros(t.shape) #list for all displacements in time
    v = np.zeros(t.shape) #list for all velocities in time
    u[0] = u0 #set initial displacement
    v[0] = v0 #set initial velocity
    
    #function for acceleration y'' (m/s^2):
    def func(u, v, force):
        return (force - b * v - k * u) / m

    for i in range(t.size - 1):
        #force at half a time step:
        half_force = (force[i + 1] - force[i]) / 2 + force[i]

        #apply Runge-Kutta formulas to find u1-u4, v1-v4, and a1-a4:
        u1 = u[i]
        v1 = v[i]
        a1 = func(u1, v1, force[i])
        u2 = u[i] + v1 * t_step / 2
        v2 = v[i] + a1 * t_step / 2
        a2 = func(u2, v2, half_force)
        u3 = u[i] + v2 * t_step / 2
        v3 = v[i] + a2 * t_step / 2
        a3 = func(u3, v3, half_force)
        u4 = u[i] + v3 * t_step
        v4 = v[i] + a3 * t_step
        a4 = func(u4, v4, force[i + 1])
        u[i + 1] = u[i] + (t_step / 6) * (v1 + 2 * v2 + 2 * v3 + v4)
        v[i + 1] = v[i] + (t_step / 6) * (a1 + 2 * a2 + 2 * a3 + a4)

    return u, v

#Use Runge-Kutta Method:
u, v = rungeKutta(t, u0, v0, m, b, k, force)

#Displaying the graph:
plt.figure(dpi=300) 
plt.title(f"Runge-Kutta Method - Modified (step={t_step})")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.plot(t, u, color="blue")
plt.plot(t, force/(f0*40), color="red")
plt.plot(t, np.zeros(len(t)) , color="gray")
plt.show()
