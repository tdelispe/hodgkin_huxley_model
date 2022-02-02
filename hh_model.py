#Hodkin-Huxley model
#equation are from "Friedman, Avner & Borisyuk, Alla & Ermentrout, Bard & Terman, David. (2005). Introduction to Neurons. 10.1007/978-3-540-31544-5_1".

import numpy as np

import matplotlib.pyplot as plt


# import of the function that solves the ODE system

from scipy.integrate import odeint
 


#calculation of parameters from equations (4), page 15

def parameters (V):

     an = 0.01 * (- V + 10.0) / (np.exp ((- V + 10.0) / 10.0) -1.0)

     bn = 0.125 * np.exp (-V / 80.0)

     am = 0.1 * (- V + 25.0) / (np.exp ((- V + 25.0) / 10.0) -1.0)

     bm = 4.0 * np.exp (-V / 18.0)

     ah = 0.07 * np.exp (-V / 20.0)

     bh = 1.0 / (np.exp ((- V + 30.0) / 10.0) + 1.0)
     

     return an, bn, am, bm, ah, bh


# Cm represents the capacity of the axon membrane per unit area (mF/cm^2)

# gk, gNa, gL are the conductivities of potassium, sodium and leak channel respectively, per unit area (mS/cm^2)

# Ek, ENa, El are the Nerst equilibrium voltages for each channel in (mV)

# The default values have been taken from "Weinberg, Seth. (2013). High-Frequency Stimulation of Excitable Cells and Networks.

# The default values have been taken from "Cellular function given parametric variation in the Hodgkin and Huxley model of excitability, Hillel Ori, Eve Marder, Shimon Marom,
#Proceedings of the National Academy of Sciences Aug 2018, 115 (35) E8211-E8218; DOI: 10.1073/pnas.1808552115"


def values (Cm = 1.0, gk_star = 36.0, gNa_star = 120.0, gL = 0.3, Ek = -77.0, ENa = 50.0, El = -54.0):

     values = (Cm, gk_star, gNa_star, gL, Ek , ENa , El )

     return values



    
# the hodgkin_huxley function implements the system of Ordinary Differental Equation to be solved (equations 1 & 3)

# y0 is an array with the initial values of action potential V and n, m, h

# t is an array with the timesteps of the simulation



def hodgkin_huxley(y, t, *values):

        Cm, gk_star, gNa_star, gL, Ek , ENa , El = values


        V, n, m, h = y
        

        #we build a system dy/dt using the equations (1) and (3)
    
        #dydt array contains four rows, for dV/dt, dn/dt, dm/dt, dh/dt
        
        dydt = np.zeros((4,))


        #Calculation of conductivities from equations (2)
        
        #gk, gNa are (V,t) dependent, we do not calculate gL because it is constant from the assumption  
        
        gk, gNa = n**4.0 * gk_star, m**3.0 * h * gNa_star

        
        #equation (1)
        #we take dV/dt from equation (1) and we put it in the first row of the dydt array

        
        dVdt =  ( Im(t) - gk * (V-Ek) - gNa * (V-ENa) - gL * (V-El) ) /Cm
        

        dydt[0] = dVdt
        

        # we calculate dn/dt, dm/dt, dh/dt from equations (3) and we put them in the dydt array
        # we call the function "parameters" to compute an, am, ah, bn, bm, bh

        an, bn, am, bm, ah, bh = parameters (V)  

        dndt, dmdt, dhdt =   an * (1-n) - bn * n,   am * (1-m) - bm * m,    ah * (1-h) - bh * h

        dydt[1], dydt[2], dydt[3] = dndt, dmdt, dhdt


        return dydt




#set the current that is being injected into the axon as a function of time
#here we stimulate the membrane with a current pulse (100uA/cm^2) between 1-2ms

def Im(t):

    if 2.0 < t < 3.0:
        Im = 100.0
    else:
         Im = 0
    return Im
  

    
  

# set values of axon membranes capacity, conductivities and Nerst equilibrium voltages for the channels
# values of membrane capacity, channel conductivities and Nerst voltages can be changed from here
# e.g. to set axon membrane capacity equal to 2 mF/cm^2 and Nerst equilibrium voltage for potassium channel equal to -70 mV, type: ch_values(Cm=2, Ek=-70)

##gK = 35
####
##### Average sodium channel conductance per unit area (mS/cm^2)
##gNa = 40
##
### Average leak channel conductance per unit area (mS/cm^2)
##gL = 0.3
##
### Membrane capacitance per unit area (uF/cm^2)
##Cm = 1.0
##
### Potassium potential (mV)
##VK = -77
##
### Sodium potential (mV)
##VNa = -65
##
### Leak potential (mV)
##Vl = 10.613

ch_values = values()
       
#set the initial values for resting potential Vrest and n, m, h

Vrest = -55

an, bn, am, bm, ah, bh = parameters (Vrest)


#steady state solution dn/dt=0...

n0= an / (an + bn)

m0= am / (am + bm)

h0= ah / (ah + bh)



#store initial values in an array


y0 = np.array([ Vrest, n0, m0, h0])


  

#set the time of simulation (start, end, step) in ms

t=np.linspace(0, 30.0, 1000)


#solve ODE system


solution=odeint(hodgkin_huxley,y0,t, args=ch_values)


               



#plot the results


# calculate currents of potassium and sodium channels


Ik=(solution[:,1]**4 * ch_values[1] * (solution[:,0]- ch_values[4]))

INa=(solution[:,2] **3* ch_values[2] *solution[:,3]* (solution[:,0]-ch_values[5]))


I = [Im(i) for i in t]

# plot stimulation

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(t, I)
ax.set_xlabel('time (ms)')
ax.set_ylabel(r'Current density (uA/$cm^2$)')
ax.set_title('stimulation current')
plt.grid()

#plot action potential vs time

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(t, solution[:, 0])
ax.set_xlabel('time (ms)')
ax.set_ylabel('voltage (mV)')
ax.set_title('action potential')
plt.grid()

#plot trajectories with limit cycles, action potential vs n & action potential vs m

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(solution[:, 0], solution[:, 1], '--', label='action potential vs gating variable n')
ax.plot(solution[:, 0], solution[:, 2], '--', label='action potential vs gating variable m')
ax.plot(solution[:, 0], solution[:, 3], '--', label='action potential vs gating variable h')
ax.set_title('limit cycles')
ax.legend()
plt.grid()

# plot channel currents vs time
fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(t, Ik, label='Ik')
ax.plot(t, INa, label='INa')
ax.set_title('channel currents')
ax.legend()
plt.grid()

# plot gating variables vs time

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(t, solution[:,1], label='n')
ax.plot(t, solution[:,2], label='m')
ax.set_title('gating variables')
ax.legend()
plt.grid()


plt.show()


