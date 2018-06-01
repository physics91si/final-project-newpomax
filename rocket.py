import motor
import numpy as np

class Rocket:
    '''Calculates and stores information about a standard rocket, with rectangular/trapezoidal fins, tangent ogive nosecone, single deploy chute, and little to no surface protrusion. All units are SI (standard m,kg,s).'''
    def r_crit(self):
        '''Determines the critical Reynold's number for the rocket.'''
        return 51*np.power((self.surface_finish/self.length),-1.039)
    
    def Mach(self,velocity):
        '''Determines the rocket's current Mach number.'''
        return velocity/331.2
    
    def C_f(self,altitude,R):
        '''Determines the rocket's Reynolds-based coefficient based off the current Reynold's number, R, and the critical Reynold's number.'''
        R_crit = self.r_c
        if(R<10**4):
            return 1.48*10**(-2)
        elif(R<R_crit):
            return 1/((1.5*np.log(R)-5.6)**2)
        else:
            return 0.032*(self.surface_finish/self.length)**(0.2)
        
    def C_f_c(self,altitude,velocity,R):
        '''Determines the final velocity/Reynolds coefficient based off the current Mach number and the Reynolds-based coefficient.'''
        M = self.Mach(velocity)
        if(M<1):
            return self.C_f(altitude,R)*(1-0.1*M**2)
        else:
            return self.C_f(altitude,R)/((1+.15*M**2)**0.58)
        
    def mean_chord(self):
        '''Calculates the mean chord of the fins.'''
        return (self.inner_chord+self.outer_chord)/2.
    
    def A_fins_wet(self):
        '''Calculates the wet area of the fins (the area in contact with air).'''
        face_area = .5*(self.outer_chord + self.inner_chord)*self.fin_height
        upper_edge_l = self.fin_height/np.cos(self.fin_angle*np.pi/180.)
        lower_edge_l =np.sqrt( (self.inner_chord-(self.fin_height*np.tan(self.fin_angle*np.pi/180.)+self.outer_chord))**2 + self.fin_height**2)
        return self.num_fins*(2*face_area + self.fin_thickness*(upper_edge_l + self.outer_chord + lower_edge_l))
    
    def A_body_we(self):
        '''Calculates the wet area of the rocket body.'''
        tube_area = self.length *np.pi*self.diameter
        ogive_rad = (self.length**2 + self.diameter**2/4)/(self.diameter)
        cone_area = 2*np.pi*ogive_rad*((self.diameter/2-ogive_rad)*np.arcsin(self.nose_length/ogive_rad)+self.nose_length)
        return tube_area+cone_area
    
    def C_d_friction(self,altitude,velocity,R):
        '''Calculates the total coefficient of drag from air-skin friction along the body and fins.'''
        A_body = self.A_body_w
        A_fins = self.A_fins_w 
        t = self.fin_thickness
        c = self.m_chord 
        Cfc = self.C_f_c(altitude,velocity,R)
        return Cfc * ((1+(1./2/((self.length+self.nose_length)/self.diameter)))*A_body + (1+(2*t/c))*A_fins)/(np.pi * self.diameter**2 / 4.)
    
    def C_d_base(self,altitude,velocity):
        '''Calculates the based drag coefficient.'''
        M=self.Mach(velocity)
        if(M<1):
            return (.12+.13*M**2)
        else:
            return (.25/M)
        
    def C_d_fins(self,altitude,velocity):
        '''Calculates the base drag for the fins.'''
        C_d_b = self.C_d_base(altitude,velocity)
        C_d_f = 0.0
        M=self.Mach(velocity)
        if(M<0.9):
            C_d_f = ((1-M**2)**(-1*0.417)) -1.0
        elif(M<1):
            C_d_f = 1-1.785*(M-0.9)
        else:
            C_d_f = 1.214 - (.502/M**2) + (1.095/M**4)
        C_d_tot = C_d_b + C_d_f
        return self.num_fins*self.fin_thickness*self.fin_height*C_d_tot/(np.pi*self.diameter**2/4.)
    
    def C_d(self,altitude,velocity,R):
        '''Calculates the total coefficient of drag for the rocket.'''
        return self.C_d_friction(altitude,velocity,R) + self.C_d_base(altitude,velocity) + self.C_d_fins(altitude,velocity)
    
    def __init__(self,length, diameter,mass,fin_thickness,inner_chord,outer_chord,fin_height,nose_length,motor_name,motor_mass,num_fins=3,surface_finish=60*10**(-6),chute_Cd=1,chute_diameter=1,fin_angle=0,deploy_altitude=-1,deploy_time=-1):
        '''Defines a rocket with the following features: length, length of body tube in m; diameter, diameter of body tube in m; mass, mass of rocket without motor components in kg; fin_thickness, the thickness of the fins in m; inner_chord, the length of the side of the fin closest to the rocket in m; outer_chord, the length of the side of the fin furthest from the rocket in m; fin_height, the distance between the edge of the rocket and the furthest extent of the fin from the rocket in m; nose_length, the length of the nosecone (not including any shoulder) in m; motor_name, the name of the motor in the form 'L2500', with no tags or additional letters included after the average force number; motor_mass, the mass of the motor including its casing or additional hardware in kg; num_fins, the number of fins; surface_finish, the average 'roughness' of the surface finish in m, default is that of finely-sprayed paint; chute_Cd, the coefficient of drag for the parachute; chute_diameter, the diameter of the chute in m; fin_angle, the angle in degrees of the leading edge of the fin (for a rectangle, 0 deg, for more oblique fins, the angle increases); deploy_altitude, the altitude at which the parachute is to be deployed, if applicable; deploy_time, the time after ignition at which the parachute is to be deployed, if there is no deploy altitude (i.e. deploy_altitude=-1).'''
        self.motor_name=motor_name
        self.motor_mass=motor_mass
        self.motor_curve=motor.get_data(self.motor_name)
        self.length=length
        self.diameter=diameter
        self.mass = mass
        self.surface_finish=surface_finish
        self.fin_thickness=fin_thickness
        self.num_fins=num_fins
        self.inner_chord=inner_chord
        self.outer_chord=outer_chord
        self.fin_height=fin_height
        self.fin_angle=fin_angle
        self.nose_length=nose_length
        self.A_body_w=self.A_body_we()
        self.A_fins_w=self.A_fins_wet()
        self.r_c=self.r_crit()
        self.m_chord=self.mean_chord()
        self.deploy_time=deploy_time
        self.deploy_altitude=deploy_altitude
        self.chute_Cd=chute_Cd
        self.chute_diameter=chute_diameter
        
