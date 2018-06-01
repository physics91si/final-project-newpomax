import rocket
import matplotlib.pyplot as plt
import numpy as np

def __init__rocket(rocket_params):
    '''Initializes Rocket object based on input parameters (see Rocket class for explication of parameters).'''
    Rocket=rocket.Rocket(*rocket_params)
    return Rocket

def run_sim(Rocket,dt=0.01,g=9.80665,plot=True,MAX_TIME=1000,parachute_deploy_time=0.5):
    '''Runs simulation based off Rocket object and standard Earth environment. Simulation uses time steps of size dt and a gravitational acceleration g in m/s^2. \nIf plot is True, the rocket's altitude, velocity, and acceleration will be plotted. \nMAX_TIME is a catch in case of infinite loop error or too-long simulation; if rocket is expected to be in flight more than 1000 sec, this parameter should be modified. \nParachute_deploy_time is the time it takes for a parachute to deploy; in the case of extremely large parachutes, this parameter prevents sudden variable overflow at high speeds (to remove effects, make this parameter arbitrarily small).\nIf MAX ACCEL EXCEEDED is printed, then the parameters created an overflow in acceleration.'''
    v=0
    h=0
    time=0
    time_since_apogee = -1
    current_time_step=0
    apogee = False
    chute_deployed = False   
    deploy_at_alt = False
    deploy_at_time =False
    MAX_ACCEL_EXCEEDED = False
    
    if Rocket.deploy_time>0:
        deploy_at_time=True
    elif Rocket.deploy_altitude>0:
        deploy_at_alt=True
        
    time_list = [time]
    altitude_list= [h]
    velocity_list=[v]
    acceleration_list=[0] 
    parachute_multiplier =0.0
    maximums = [0,0,0,0] #height,speed,accel,accel_up
    max_times = [0,0,0,0] #correlated with above
    motor_depletion_rate = 0.4*Rocket.motor_mass/Rocket.motor_curve[-1][0]
    current_motor_mass = Rocket.motor_mass
    
#---------START OF ITERATION LOOP---------------------------------------------------------------------    
    while (h>0 or apogee==False) and time<MAX_TIME: #MAX_TIME is here just a catch for infinite loops in case of error
        Temp = temp(h)
        Pres = pressure(Temp)
        dens = density(Temp,Pres)
        R = reynolds(Rocket,h,v,Temp,dens)
        Motor_force,current_time_step = F_motor(Rocket,time,current_time_step)
        F_tot = 0
        max_accel = abs(v/dt)
        
        if current_motor_mass>Rocket.motor_mass*0.6:
            current_motor_mass = Rocket.motor_mass-(motor_depletion_rate*time)
            if current_motor_mass<Rocket.motor_mass*0.6:
                current_motor_mass=Rocket.motor_mass*0.6
                
        total_mass = Rocket.mass + current_motor_mass
        
        if parachute_multiplier<1 and chute_deployed==True:
            parachute_multiplier += dt/parachute_deploy_time
            
        if (deploy_at_alt==True and apogee==True and h<=Rocket.deploy_altitude) or (deploy_at_time==True and time>=Rocket.deploy_time):
            chute_deployed=True
            F_tot =  -1*total_mass*g + parachute_multiplier*F_chute(Rocket,v,Pres)
            if F_tot>max_accel*total_mass:
                F_tot=max_accel*total_mass
                if MAX_ACCEL_EXCEEDED==False:
                    print "Max. Acceleration was exceeded: Data below may not be accurate. Please check input parameters (parachute/motor may be too large/powerful, or deployment time may be prior to apogee)."
                    MAX_ACCEL_EXCEEDED=True
        elif Motor_force<=total_mass*g and current_time_step<len(Rocket.motor_curve)/2:
            F_tot = 0
        else:
            if v>=0:
                F_tot=-1*F_d(Rocket,h,v,Pres,R) + Motor_force - total_mass*g
            else:
                F_tot = F_d(Rocket,h,v,Pres,R) + Motor_force - total_mass*g
        
        acceleration = F_tot/total_mass
        h,v = update_all(h,v,acceleration,dt)
        maximums,max_times= update_maxs(maximums,max_times,time,apogee,[h,v,acceleration])
        time+=dt
        
        if v<0:
            apogee = True
        time_list.append(time)
        altitude_list.append(h)
        velocity_list.append(v)
        acceleration_list.append(acceleration)
#---------END OF ITERATION LOOP------------------------------------------------------------------   

    if MAX_ACCEL_EXCEEDED==True:
        print "Max. Acceleration was exceeded: Data below may not be accurate. Please check input parameters (parachute/motor may be too large/powerful, or deployment time may be prior to apogee)"
    print 'Apogee =', maximums[0],'m at time =',max_times[0],'sec'
    print 'Max. speed = ', maximums[1],'m/s at time =',max_times[1],'sec'
    print 'Max. accel = ', maximums[2],'m/s^2 at time =',max_times[2],'sec'
    print 'Max. motor accel = ', maximums[3],'m/s^2 at time =',max_times[3],'sec'
    print 'Final height: ',altitude_list[-1], 'm'
    print 'Final vel: ',velocity_list[-1],'m/s'
    print 'Final accel: ', acceleration_list[-1],'m/s^2'
    if plot==True:
        plt.plot(time_list,altitude_list,c='blue',label='Altitude ($m$)')
        plt.plot(time_list,velocity_list,c='green',label ='Velocity ($\\frac{m}{s}$)')
        plt.plot(time_list,acceleration_list,c='red',label = 'Acceleration($\\frac{m}{s^2}$)')
        plt.legend()
        plt.xlabel('Time (sec)')
        plt.title('Rocket Flight Data')

def reynolds(Rocket,altitude,velocity,Temp,dens):
    '''Compute the Reynolds number, R, based off of Rocket characteristics, velocity (m/s), air temperature (C), and air density(kg/m^3).'''
    kinvisc = (0.001792/dens)*np.exp(-1.94-4.80*(273.16/(Temp+273.16))+6.74*np.power((273.16/(Temp+273.16)),2))
    return velocity*Rocket.length/kinvisc

def pressure(T):
    '''Compute the pressure in kPa at a given temperature (C).'''
    if(T>-273):
        return 101.29 * np.power(((T+273.1)/288.08),5.256)
    
    else:
        return 101.29
    
def temp(h):
    '''Compute the air temperature in C at a given altitude, h (m)'''
    return 15.04 - 0.00649*h

def density(T,p):
    '''Compute the air density in kg/m^3 from at a given air temperature (C) and pressure (kPa).'''
    return p/(.2869*(T+273.1))

def F_d(Rocket,altitude,velocity,Pres,R):
    '''Return the drag force on the rocket in N based off of the Rocket profile, altitude (m), velocity (m/s), pressure (kPa) and the Reynolds number, R.'''
    return 2.15*Rocket.C_d(altitude,velocity,R)*np.pi/4.*np.power(Rocket.diameter,2)*Rocket.num_fins*Rocket.fin_height*Rocket.fin_thickness*Pres*np.power(velocity,2)

def F_motor(Rocket,time,current_time_step):
    '''Calculate the force of the motor based off the characteristics of the Rocket as well as the current time and size of time step.'''
    if time>Rocket.motor_curve[len(Rocket.motor_curve)-1][0] or current_time_step>=len(Rocket.motor_curve):
        return 0.0,current_time_step
    
    if time>=Rocket.motor_curve[current_time_step+1][0] and current_time_step<len(Rocket.motor_curve)-1:
        current_time_step= current_time_step+1
        
    if time == Rocket.motor_curve[current_time_step][0]:
        return Rocket.motor_curve[current_time_step][1],current_time_step
    
    else:
        slope = (Rocket.motor_curve[current_time_step][1]-Rocket.motor_curve[current_time_step+1][1])/(Rocket.motor_curve[current_time_step][0]-Rocket.motor_curve[current_time_step+1][0])
        return (time-Rocket.motor_curve[current_time_step][0])*slope + Rocket.motor_curve[current_time_step][1],current_time_step
    
def F_chute(Rocket,velocity,Pres):
    '''Calculate the drag force of the parachute based off the Rocket characteristics, velocity(m/s), and air pressure (kPa).'''
    return 0.005*Rocket.chute_Cd*np.pi/4.*Pres*np.power(Rocket.chute_diameter,2)*np.power(velocity,2)

def update_all(h,v,a,dt):
    '''For the simulation, updating the altitude and velocity based on current acceleration and the size of the time step.'''
    pos_vel = [h,v]
    k1 = np.dot(dt,[v,a])
    k2 = np.dot((dt/2),[v,a]+k1)
    k3 = np.dot((dt/2),([v,a]+k2))
    k4 = np.dot((dt),([v,a]+k3))
    pos_vel = pos_vel + np.dot((1./4.),(k1+2*k2+2*k3+k4))
    return pos_vel[0],pos_vel[1]

def update_maxs(current_maxs,max_times,time,apogee,new_vals):
    '''For the simulation, updating the maximum values of altitude, velocity, and acceleration for reporting at the end of the simulation.'''
    for i in range(0,len(new_vals)):
        if abs(new_vals[i])>current_maxs[i]:
            current_maxs[i]=new_vals[i]
            max_times[i]=time
            
    if new_vals[2]>current_maxs[3] and apogee == False:
        current_maxs[3]=new_vals[2]
        max_times[3]=time
        
    return current_maxs,max_times