from Solar_Thermal import Simulation 


# Run this file to visualize simulation


    
## Further Improvements:
    
    # add convective heat transfer to model - more significant w/ slow flow rate
    # add heat loss of water traveling through non-heated pipes to storage tank
    # update solar radiation to be a sine function of time of day
        # this could also be done for outside temperature (ie heat loss from storage tank)


# warning: some design parameter combinations may lead to unstable solution

# instantiate simulation with general scenario parameters
sim = Simulation(solar_radiation = 800,
                 mass_flow_rate = 0.75,
                 initial_temp = 15, 
                 outside_temp = 10,
                 nodes = 200,
                 days_to_sim = 3)

# set panel charactersistics
sim.Panel(length = 3.0,
        height = 5.0,
        efficiency = 0.9)

# set pipe design
sim.Pipe(length = 10, radius = 0.1)

# set storage tank design
sim.Tank(height = 2.0,
         radius = 0.40,
         thickness = 0.05,
         thermal_conductivity = 0.4)

# run simulation and store temperature data
PipeData, TankData = sim.Run()

# show summary plot
a = sim.Show_SystemChange()

# show animated plots
# warning: GradientProfile animation is slow

a11 = sim.Show_LineProfile(PipeData, speed = 10) # 10x animation speed
#a12 = sim.Show_GradientProfile(PipeData)


a21 = sim.Show_GradientProfile(TankData, speed = 10) 
#a22 = sim.Show_LineProfile(TankData)
