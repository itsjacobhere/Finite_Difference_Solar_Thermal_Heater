import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation as FA



class Simulation:
    
    
    def __init__(self, solar_radiation = 1000,
                 mass_flow_rate = 1.0,
                 initial_temp = 20, 
                 outside_temp = 20,
                 nodes = 200,
                 days_to_sim = 3):
    
        # General assumptions
        self.Q = solar_radiation # W/m^2 - solar radiation
        self.T0 = initial_temp + 273.15 # (Kelvin) intial condition - Temperature 
        self.T_outside = outside_temp + 273.15
    
        # water characterstics
        self.Cp = 4180 # Joules/kg/kelvin - heat capacity
        self.rho = 1000 # (kg/m^3) - typical density of water
        
        self.n = nodes # number of nodes - granularity or resolution of temperature gradient in pipe system and tank
        self.t_final = int(days_to_sim*86400) # seconds - length of simulation - days*seconds_in_day
        
        # determined by power of pump
        self.m_dot = mass_flow_rate
    
    
    def Panel(self, 
                    length = 1.5,
                    height = 1,
                    efficiency = 0.9):
        
        self.length = length
        self.height = height
        self.eff = efficiency
    
    
    def Pipe(self,
                   length = 10,
                   radius = 0.1):
        
        self.pipe_len = length
        self.pipe_radius = radius
    
    
    def Tank(self,
             height = 2,
             radius = 0.2,
             thickness = 0.1,
             thermal_conductivity = 0.75):
        
        self.tank_height = height
        self.tank_radius = radius
        self.tank_thickness = thickness
        self.K = thermal_conductivity
    
    
    def Run(self):
        
        # set design parameters
        
        # Panel characteristics
        length = self.length # meters (m)
        height = self.height # meters (m)
        Panel_eff = self.eff # 90% of solar radiation converted to heat
        Panel_area = length * height # m^2
        Panel_power = self.Q * Panel_area * Panel_eff # W - total captured power of panel
        
        # Pipe characteristics
        pipe_len = self.pipe_len # (m) - length of heated portion of pipe
        pipe_radius = self.pipe_radius # (m) - radius
        pipe_SA = (2 * np.pi * pipe_radius) * pipe_len # SA of heated portion of pipe
        
        # Tank characteristics
        tank_height = self.tank_height # (m)
        tank_radius = self.tank_radius # (m)
        tank_thickness = self.tank_thickness # (m) 
        tank_SA = 2*np.pi*tank_radius*tank_height + 2*np.pi*tank_radius**2
        K = self.K # tank material thermal conductivity
        
        
        m_dot = self.m_dot #(kg/s) -  mass flow rate - dependent on pump power 
        q_pipe = Panel_power / pipe_SA # (W/m^2) - heat flux applied to pipe from solar panel
        
        n = self.n
        T0 = self.T0
        T_outside = self.T_outside
        Cp = self.Cp
        rho = self.rho
        
        # finite difference parameters
        dt = 1 # seconds - time step
        dx = pipe_len / n # (m) - width of each node in pipe
        dx_tank = tank_height / n
        self.dx = dx
        self.dx_tank = dx_tank
        
        # set initial condiitons
        T = np.ones(n) * T0 # initialize temperature array of each node in pipe
        T_tank = np.ones(n) * T0 # init temp throughout tank
        
        dT_dt = np.zeros(n) # init change in temp of each node as zeros vector
        dT_dt_tank = np.zeros(n) # init temp deriv at each node in tank
        
        
        t = np.arange(0, self.t_final, dt) # define time vector
        self.t = t
        
        t1 = time.time() # start timer
        
        # empty lists for data storage and later visualization
        self.PipeAtTime = []
        self.TankAtTime = []
        self.PipeAvg = []
        
        
        # run simulation
        for timestep in range(0, len(t)):
            
            
            # copy to dereference from updating array, convert to Celsius
            T_copy = np.copy(T - 273.15)
            T_tank_copy = np.copy(T_tank - 273.15)
            
            self.PipeAtTime.append(T_copy) # store pipe system temp profile at time t
            self.TankAtTime.append(T_tank_copy) # stores tank temperature profile
            self.PipeAvg.append(np.mean(T) - 273.15) # store average temp for simple plot
            
            # solar thermal heated pipe system, flow from storage tank  
            if timestep == 0: # use initial conditions at start of simulation
                q_tank = - ((K*tank_SA*(T0 - T_outside)) / tank_thickness)  # dynamic heat transfer of tank to environment
                dT_dt[0] = (m_dot*Cp*(T0 - T[0]) + q_pipe*2*np.pi*pipe_radius*dx) \
                                                / (rho*Cp*dx*np.pi*pipe_radius**2)
            else:   # use flow from storage tank to calculate first element temp change      
                q_tank = - ((K*tank_SA*(np.mean(T_tank) - T_outside)) / tank_thickness) # heat loss from storage tank               
                dT_dt[0] = (m_dot*Cp*(T_tank[n-1] - T[0]) + q_pipe*2*np.pi*pipe_radius*dx) \
                                                    / (rho*Cp*dx*np.pi*pipe_radius**2)
                                                  
            # compute change in temp for all other elements in pipe
            dT_dt[1:n] = (m_dot*Cp*(T[0:n-1] - T[1:n]) + q_pipe*2*np.pi*pipe_radius*dx) \
                                                / (rho*Cp*dx*np.pi*pipe_radius**2)
            
            
            T += dT_dt*dt # update pipe temperature gradient with time
            
            # change in temp of storage tank system, input is from heated pipe system
            dT_dt_tank[0] = (m_dot*Cp*(T[n-1] - T_tank[0]) + q_tank*2*np.pi*tank_radius*dx_tank) \
                                                        / (rho*Cp*dx_tank*np.pi*tank_radius**2)
            dT_dt_tank[1:n] = (m_dot*Cp*(T_tank[0:n-1] - T_tank[1:n]) + q_tank*2*np.pi*tank_radius*dx_tank) \
                                                        / (rho*Cp*dx_tank*np.pi*tank_radius**2)
            
            T_tank += dT_dt_tank*dt # update temp grad of storage tank
            
            
        print('Simulation Run Time: {:.2f} seconds'.format( time.time() - t1))
        
        return self.PipeAtTime, self.TankAtTime
    
    
    def Show_SystemChange(self):
        
        fig = plt.figure(3)
        plt.plot(self.t/(24*3600), np.array(self.PipeAvg), linewidth = 2 )
        plt.xlabel('Days')
        plt.ylabel('Average System Temperature (C)')
        plt.title('Change in Avg. System Temperature')
        
        return fig
        
    
    
    def Show_LineProfile(self, PipeAtTime, speed = 1):
        ## animated temperature profile of heated pipe system
        
        dx = self.dx
        
        fig1 = plt.figure(1)
        ax1 = plt.axes(xlim = (0, self.pipe_len)) 
        plt.xlabel('Length (m)')
        plt.ylabel('Temperature (C)')
        line, = ax1.plot([], [], lw = 2)
        
        def init():
            line.set_data([], [])
            
            return line,
        
        def animate1(i):
            x = np.linspace(dx/2, self.pipe_len - dx/2, self.n)
            y = PipeAtTime[1::speed][i]
            
            plt.title('Time: {:.2f} minutes'.format(speed*i/60))
            ax1.set_ylim(min(y),max(y)) # dynamic y axis
            
            line.set_data(x, y)
            return line,
        
        anim1 = FA(fig1, animate1, init_func = init, frames = self.t_final, interval = 1/1000, blit = False)
        
        plt.show()
        
        return anim1
    
    
    
    def Show_GradientProfile(self, TankAtTime, speed = 1):
    
        ## plot animated temperature profile of storage tank
        
        
        expanded = []
        for i in self.TankAtTime[0]:
            expanded.append(np.ones(self.n)*i)
        expandedArray = np.array(expanded)
        
        fig2 = plt.figure()
        im = plt.imshow(expandedArray, interpolation = 'bilinear', cmap = 'viridis')
        plt.xlabel('Diameter/Width')
        plt.ylabel('Height/Length')
        
        def init():
            im.set_data(expandedArray)
            cbar = plt.colorbar()
            cbar.set_label('Temperature (C)')
            return [im]
        
            
        
        def animate2(i):
            
            expanded = []
            for k in TankAtTime[1::speed][i]:
                expanded.append(np.ones(self.n)*k)
            expandedArray = np.array(expanded)
            
            # dynamic colorbar range
            vmin = min(TankAtTime[1::speed][i])
            vmax = max(TankAtTime[1::speed][i])
            
            plt.title('Time: {:.2f} minutes'.format(speed*i/60))
            im.set_array(expandedArray)
            im.set_clim(vmin, vmax)
            
            return [im]
            
        anim2 = FA(fig2, animate2, init_func = init, frames = self.t_final, interval = 1/100, blit = False)
                
        plt.show()
        
        return anim2
    
    
    