#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 18:10:33 2018

@author: leonardo
"""
import math as m
import numpy as np
import matplotlib.pyplot as plt
import constants
#-------------------------------------------------------------------------
#       Program Overview
#-------------------------------------------------------------------------
#       The purpose of this program is to calculate the altitude and
#       velocity of a rocket at a given time t, taking into account
#       air resistance, changing air density, changing gravity,
#       the thrust of the rocket, and the changing mass of the rocket.
#       It would be impossible to determine a closed solution for
#       this problem, so using numerical methods is the only possible
#       way to find the solution. In this program, the Runge-Kutta-Fehlberg
#       method is used for integration (4th and 5th order). 
#-------------------------------------------------------------------------
#       Variables
#-------------------------------------------------------------------------
#       Earthradius [m] Radius of Earth
#       Initgravity [m/s²] Gravity at surface
#       Initdensity [kg/m³] Air density at surface
#       R [Joules/kg*K] Gas constant for air
#       Temp [K] Mean temperature of atmosphere
#       Endtime [s] When to stop plotting
#       Rocketmass [kg] Mass of rocket without fuel
#       Fuelmass [kg] Mass of fuel
#       Impulse [s] Impulse of engine
#                   300 - Kerosene/Oxygen
#                   360 - Hydrogen/Oxygen
#                   490 - Hydrogen/Fluorine
#       Burnrate [kg/s] Amount of fuel burned per second
#       Endthrust [s] Total time of thrust
#       Velocity [m/s] Velocity of material ejected from rocket nozzle
#       Initthrust [N] Thrust of engine (constant)
#       Cd [dimensionless] Coefficient of Drag
#       SurfaceArea [m²] Frontal surface area
#       Initheight [m] Initial altitude of rocket
#       Initvelocity [m/s] Initial velocity of rocket
#       time [s] Current time
#       dt [s] Step size
#       Vold,Vnew [m/s] Velocity
#       Hold, Hnew [m] Altitude
#       Used in the Runge-Kutta-Fehlberg (45) integration method
#       function accel computes the overall acceleration of the rocket at
#               a given time, altitude, and velocity.
#
#-------------------------------------------------------------------------   
#-------------------------------------------------------------------------
#       Initialization (Constant Physics)
#-------------------------------------------------------------------------

Earthradius = constants.R_earth        # [m]
Initgravity =  constants.g             # [m/s²] initial gravity at the earths surface
Initdensity = constants.rho_air        # [kg/m³] Air density sea level
R = constants.R_air                    # [J/kgK] Air constant
Temp = constants.Tamb                  # [K] Tempature at time 0. Aprox. 20 Celsius 
M = constants.M                        # kg/mole
R_univ = constants.R_univ              # Joules/Kelvin*Mole 
L = constants.L                        # Kelvin/meter  "Tempature lapse rate"
gamma = constants.gamma_air            # Air coefficient of specific heats 
Po = constants.Pamb                    # Pascal  Pressure at time 0 

Rocketmass = 500.0
Fuelmass = 300.0
Impulse = 300.0
Burnrate = 5.0                         # Mass flow of propellant    
Endthrust = Fuelmass/Burnrate
Velocity = Initgravity*Impulse         # Effective exhaust velocity 
Initthrust = Burnrate*Velocity         # Thrust 
SurfaceArea = 0.1
NozzleArea = 0.025
NozzlePressure = 0.9*Po
Initheight = 0.0
Initvelocity = 0.0
Endtime = 1e1000


# (by setting Endtime = 1.e9, we can watch the rocket crash)
#-------------------------------------------------------------------------
#       Functions
#-------------------------------------------------------------------------

def accel(h,v,t):
# This function computes the total acceleration of the 
# rocket at a given time, height, and velocity.
    Acceleration = (Thrust(t) - Drag(v,h))/Mass(t) - Gravity(h)#*m.sin(Angle(t,v,h))

    return Acceleration

#-------------------------------------------------------------------------

def Thrust(t):   
# This function computes the thrust at a given time. Note that
# a more complex thrust curve could be put here rather than just
# a constant value.
    Thrust = Initthrust
    if (t >= Endthrust):
        Thrust = 0

    return Thrust

#-------------------------------------------------------------------------

def Mass(t):   
# This function computes the mass of the rocket at a given 
# time t.    
    Mass = Rocketmass
    if (t < Endthrust):
        
        Mass = Rocketmass + Fuelmass - Burnrate*t
        
    return Mass
        
#-------------------------------------------------------------------------

def Gravity(h):
# This function computes the acceleration of gravity at an
# altitude h, given that the gravity at the surface is
# Initgravity.
    Gravity = Initgravity*(Earthradius/(Earthradius + h))**2

    return Gravity

#-------------------------------------------------------------------------

def Drag(v,h):
# Function Drag computes the air drag on the rocket for a certain
# velocity and altitude. Note that drag acts in the direction 
# opposite of the given velocity v. 
    Cd = Cdrag(v,h)
    Drag = Cd*(0.5)*Density(h)*(v**2)*SurfaceArea

    return Drag
        
#-------------------------------------------------------------------------

def Density(h):       
# Function Density finds the density at a given altitude h.
# Taken from the Barometric formula. Accurate up to 80,000m above sea level
# The altitude model used below is based on the standard atmospheric model used in modern meteorology.
# It takes into accout the different regression rates and properties of the thermoclines.    
    pb = [1.2250, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064]
    hb = [0.0, 11000, 20000, 32000, 47000, 51000, 71000]
    Tb = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
    Lb = [-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002]

    if h < 11000:
        i = 0
        Density = pb[i]*(Tb[i]/(Tb[i] + Lb[i]*(h - hb[i])))**(1 + (Initgravity*M)/(R_univ*Lb[i]))
        return Density
    elif h > 11000 and h < 20000:
        i = 1
        Density = pb[i]*m.exp((-1*(Initgravity)*M*(h - hb[i]))/(R_univ*Tb[i]))
        return Density
    elif h > 20000 and h < 32000:
        i = 2
        Density = pb[i]*(Tb[i]/(Tb[i] + Lb[i]*(h - hb[i])))**(1 + (Initgravity*M)/(R_univ*Lb[i]))
        return Density
    elif h > 32000 and h < 47000:
        i = 3
        Density = pb[i]*(Tb[i]/(Tb[i] + Lb[i]*(h - hb[i])))**(1 + (Initgravity*M)/(R_univ*Lb[i]))
        return Density
    elif h > 47000 and h < 51000:
        i = 4
        Density = pb[i]*m.exp((-1 * (Initgravity)*M*(h - hb[i]))/(R_univ*Tb[i]))
        return Density
    elif h > 51000 and h < 71000:
        i = 5
        Density = pb[i]*(Tb[i]/(Tb[i] + Lb[i]*(h - hb[i])))**(1 + (Initgravity*M)/(R_univ*Lb[i]))
        return Density
    elif h > 71000 and h < 86000:
        i = 6
        Density = pb[i]*(Tb[i]/(Tb[i] + Lb[i]*(h - hb[i])))**(1 + (Initgravity*M)/(R_univ*Lb[i]))
        return Density
    else:
        return 0.0

    return Density

#-------------------------------------------------------------------------

def Temperature(h):
# Calculates air temperature [Kelvin] at altitude [m]
# from equations at 
# http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    if h <= 11000:
        # Troposphere
        T = Temp - L*h
    elif h <= 25000:
        # Lower Stratosphere
        T = 216.69
    else:
        T = 141.94 + 0.00299*h
        
    return T

#-------------------------------------------------------------------------

def Pressure(h):
# Calculates air pressure [Pa] at altitude [m]"
# from equations at 
# http://www.grc.nasa.gov/WWW/K-12/airplane/atmosmet.html
    rho = Density(h)    
    T = Temperature(h)
    p = rho*R*T

    return p 

#------------------------------------------------------------------------- 

def Speed_sound(h):
# Calculates seeped sound [m/s] at altitude [m]"
    T = Temperature(h)
    speed_sound = np.sqrt(gamma*R*T)
        
    return speed_sound   

#--------------------------------------------------------------------------

def Mach_number(v,c,h):
# Calculates Mach Number
    c = Speed_sound(h)
    Mach = v/c
    
    return Mach

#-------------------------------------------------------------------------

def to_radians(degree):
    
    return degree*np.pi/180

#-------------------------------------------------------------------------

def Cdrag (v,h):
# Calculates drag coefficient [-] at altitude [m] and Mach number[-]"
# Drag function for V2
# derived from Sutton, "Rocket Propulsion Elements", 7th ed, p108
# probably not that relevant to other body types
# use nose cone formula
    theta = to_radians(15)  
    c = Speed_sound(h) 
    mach = Mach_number(v,c,h)
    
    if mach > 5:
        cd = 2*(m.sin(theta)**2)  
            
    elif mach > 1.8 and mach <= 5:
        cd = -0.03125*mach + 0.30625

    elif mach > 1.2 and mach <= 1.8:
        cd = -0.25*mach + 0.7

    elif mach > 0.8 and mach <= 1.2:
        cd = 0.625*mach - 0.35

    elif mach <= 0.8:
        cd = 0.15
        
    return cd

#-------------------------------------------------------------------------

def main():
#-------------------------------------------------------------------------
#       Main Program
#-------------------------------------------------------------------------
##### Coefficients used to compute the independent variable argument of f [time]
    a2  =   2.500000000000000e-01  #  1/4
    a3  =   3.750000000000000e-01  #  3/8
    a4  =   9.230769230769231e-01  #  12/13
    a5  =   1.000000000000000e+00  #  1
    a6  =   5.000000000000000e-01  #  1/2

##### Coefficients used to compute the dependent variable argument of f [velocity and acceleration]
    b21 =   2.500000000000000e-01  #  1/4
    b31 =   9.375000000000000e-02  #  3/32
    b32 =   2.812500000000000e-01  #  9/32
    b41 =   8.793809740555303e-01  #  1932/2197
    b42 =  -3.277196176604461e+00  # -7200/2197
    b43 =   3.320892125625853e+00  #  7296/2197
    b51 =   2.032407407407407e+00  #  439/216
    b52 =  -8.000000000000000e+00  # -8
    b53 =   7.173489278752436e+00  #  3680/513
    b54 =  -2.058966861598441e-01  # -845/4104
    b61 =  -2.962962962962963e-01  # -8/27
    b62 =   2.000000000000000e+00  #  2
    b63 =  -1.381676413255361e+00  # -3544/2565
    b64 =   4.529727095516569e-01  #  1859/4104
    b65 =  -2.750000000000000e-01  # -11/40

##### Coefficients used to compute local truncation error estimate.  These
    # come from subtracting a 4th order RK estimate from a 5th order RK
    # estimate.
    r1  =   2.777777777777778e-03  #  1/360
    r3  =  -2.994152046783626e-02  # -128/4275
    r4  =  -2.919989367357789e-02  # -2197/75240
    r5  =   2.000000000000000e-02  #  1/50
    r6  =   3.636363636363636e-02  #  2/55

##### # Coefficients used to compute 4th order RK estimate
    # c1  =   1.157407407407407e-01  #  25/216
    # c3  =   5.489278752436647e-01  #  1408/2565
    # c4  =   5.353313840155945e-01  #  2197/4104
    # c5  =  -2.000000000000000e-01  # -1/5
    
##### Coefficients used to compute 5th order RK estimate
    d1  =  1.18518518519000000e-1  #  16/135
    d3  =  5.18986354776000000e-1  #  6656/12825
    d4  =  5.06131490342000000e-1  #  28561/56430
    d5  = -1.80000000000000000e-1  # -9/50
    d6  =  3.63636363636363636e-2  #  2/55

    
    # Initialize arrays that will be returned
    Hold = Initheight
    Vold = Initvelocity
    time = 0.0
    dtmax = 0.01; dtmin = 1e-15
    dt = (dtmax - dtmin)/2.0
    tol = dt/(10**6)
    
    # Array's
    HArray = np.array([Hold])
    VArray = np.array([Vold])
    tArray = np.array([time])
    
    while time <= Endtime:
        # Compute values needed to compute truncation error estimate and
        # the 4th and 5th order RK estimate.
        #---------------------------------------------------------------------------------------------------------------
        # Runge-Kutta-Fehlberg (4th and 5th order) Integration Method
        #---------------------------------------------------------------------------------------------------------------
        k12 = Vold
        k15 = accel(Hold, k12, time)
        k22 = Vold + dt*b21*k15
        k25 = accel(Hold + b21*k12, k22, time + a2*dt)
        k32 = Vold + dt*(b31*k15 + b32*k25)
        k35 = accel(Hold + b31*k12 + b32*k22, k32, time + a3*dt)
        k42 = Vold + dt*(b41*k15 + b42*k25 + b43*k35)
        k45 = accel(Hold + b41*k12 + b42*k22 + b43*k32, k42, time + a4*dt)
        k52 = Vold + dt*(b51*k15 + b52*k25 + b53*k35 + b54*k45)
        k55 = accel(Hold + b51*k12 + b52*k22 + b53*k32 + b54*k42, k52, time + a5*dt)
        k62 = Vold + dt*(b61*k15 + b62*k25 + b63*k35 + b64*k45 + b65*k55)
        k65 = accel(Hold + b61*k12 + b62*k22 + b63*k32 + b64*k42 + b65*k52, k62, time + a6*dt)
        #---------------------------------------------------------------------------------------------------------------        
        #Hnew = Hold + (c1*k12 + c3*k32 + c4*k42 + c5*k52)  # RK4
        H_new = Hold + dt*(d1*k12 + d3*k32 + d4*k42 + d5*k52 + d6*k62)  # RK5
        
        #Vnew = Vold + (c1*k15 + c3*k35 + c4*k45 + c5*k55)  # RK4
        V_new = Vold + dt*(d1*k15 + d3*k35 + d4*k45 + d5*k55 + d6*k65)  # RK5
        #---------------------------------------------------------------------------------------------------------------
        
        r_H = (r1*k12 + r3*k32 + r4*k42 + r5*k52 + r6*k62)*dt 
        r_V = (r1*k15 + r3*k35 + r4*k45 + r5*k55 + r6*k65)*dt
        #r_H = (H_new - Hnew)/dt
        #r_V = (V_new - Vnew)/dt
        r = np.sqrt(r_H**2 + r_V**2)

        # Compute the estimate of the local truncation error.  If it's small enough then we accept this step and save the 4th order estimate.
        #if len(np.shape(r)) > 0:
        #    r = max(r)
        
        if r <= tol:
            time += dt
            Hold = H_new  # RK5
            Vold = V_new  # RK5
            HArray = np.append(HArray, Hold)
            VArray = np.append(VArray, Vold)
            tArray = np.append(tArray, time)
            
        # Now compute next step size, and make sure that it is not too big or too small.
        dt = dt*min(max((tol/(2.0*r))**0.25, 0.1), 4.0)
        
        if dt > dtmax:
            dt = dtmax
            
        elif dt < dtmin:
            print "Error: stepsize should be smaller than %e." %(dtmin)
            break
        
        if (H_new < 0.0):
            print 'The rocket has crashed.'
            break
       
        print'\nFlight Data: Time: %.4f'%(time), '[secs]', 'Altitude: %.4f'%(H_new/10**3), '[km]', 'Velocity: %.4f'%(V_new), '[m/s]'   
        print'Erro rel.:%+12.4e'%(r), 'Step:%+12.4e'%(dt)

    # endwhile
    # Plot's
    plt.fig1 = plt.figure(1, figsize = (20, 10))
    plt.plot(tArray, HArray/10**3)
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (km)')
    plt.title('Trajectory')
    plt.grid()
        
    plt.fig1 = plt.figure(2, figsize = (20, 10))
    plt.plot(tArray, VArray)
    plt.xlabel('Time (s)', fontsize = 15)
    plt.ylabel('Velocity (m/s)', fontsize = 15)
    plt.title('Velocity Y', fontsize = 20)
    plt.grid()


    return 

if __name__ == "__main__":
    main()