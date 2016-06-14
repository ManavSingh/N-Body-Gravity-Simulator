import particle
import MyODEs
import math
import csv
import matplotlib.pyplot as plotter
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#### Using units of Solar mass, AU, 1/2pi years, and 29.78 Km/s velocity units ####

G = 1.0
particles = []

def grav_accel(system):                 #Updates the acceleration of each body in "particles" with the gravitational force
    d = [0.0]*3                         #Array to be used for calculations of distances
    particles = system                  #Local variable to hold all the particles in the system
    clear_accels(particles)
    for j in range(0, len(particles) - 1):                                            #Loop through each particle "Force on j..."
        for i in range(j+1, len(particles)):                                          #Loop through each particle "... from i"
            for k1 in range(0, 3):                                                    #Loop through spatial coordinates to fill d[k]
                d[k1] = particles[i].position[1][k1] - particles[j].position[1][k1]   #Calculate position difference in k direction
            r = (d[0]**2 + d[1]**2 + d[2]**2)**(3.0/2)                                #r = (dx^2 + dy^2 + dz^2)^(3/2)
            for k2 in range(0, 3):                                                    #Loop through spatial coordinates to do the physics
                d[k2] /= r                                                            #Self-explanatory
                particles[j].acceleration[k2] += (G*particles[i].mass) * d[k2]        #Acceleration change on j due to i
                particles[i].acceleration[k2] -= (G*particles[j].mass) * d[k2]        #Acceleration change on i due to j
            #end k for    
        #end i for
    #end j for
#end grav_accel

def nBodyRunner(system, tFinal, dt, outputFrequency):
    particles = system               #Local variable to store the system
    tSteps = int((float(tFinal)/dt) + .5)   #Determine the total number of steps needed (Rounded down).
    data = open('output.txt', 'w')
    data.write(str(int(tSteps/outputFrequency)) + ',' + str(len(particles)) + ',' + str(dt) + ',' + str(outputFrequency) + ',' +'\n') #First line of the output file contains the total number of steps and the number of particles
    grav_accel(particles)            #Compute gravitational acceleration on the particles before we move them.
    for p in particles:              #Loop through particles to do Euler step
        data.write(p.state() + '\n') #Write the intial states to the file.
        for k in range(0, 3):        #Loop through dimensions
            p.position[1][k] = MyODEs.euler(p.position[0][k], p.velocity[0][k], dt)
            p.velocity[1][k] = MyODEs.euler(p.velocity[0][k], p.acceleration[k], dt)
    for p in particles: data.write(p.state() + '\n')  #Write the Euler steps to the file.
    for tn in range(int(tSteps-1)):   #Leapfrog Time loop
        grav_accel(particles)         #Accelerate the particles
        for p in particles:           #For each particle
            for k in range(0, 3):     #For each dimension
                p.position[2][k] = MyODEs.leapfrog(p.position[0][k], p.velocity[1][k], dt)
                p.position[0][k] = p.position[1][k]
                p.position[1][k] = p.position[2][k]
                p.velocity[2][k] = MyODEs.leapfrog(p.velocity[0][k], p.acceleration[k], dt)
                p.velocity[0][k] = p.velocity[1][k]
                p.velocity[1][k] = p.velocity[2][k]
            if (tn % outputFrequency == 0): data.write(p.state() + '\n') #Write the particles' states to the file.
    data.close
    print "Complete\n"
    
#end nBodyRunner

def clear_accels(system):
    particles = system
    for p in particles:
        for k in range(3):
            p.acceleration[k] = 0.0
#end clear_accels

def net_momentum(system):
    momentumNet = [0.0 for i in range(3)]
    for p in system:
        for k in range(3):
            momentumNet[k] += p.momentum(k)
    return momentumNet
    
def init_system(inputFile):#Loads in initial condtions for n bodies and converts units.
    del particles[:]
    velConvert = 1731/29.78             #Conversion factor to go from AU/day to 29.78 Km/s
    with open(inputFile, 'r') as initialConditions:
        reader = csv.reader(initialConditions, delimiter=',')
        for row in reader:
            #Each row is a list of strings. Each row is the intial conditions for a particle.
            newBody = particle.Particle([float(row[0]), float(row[1]), float(row[2])],         #Assign particle position
                [float(row[3])*velConvert, float(row[4])*velConvert, float(row[5])*velConvert],#Assign particle velocity
                [0.0, 0.0, 0.0],        #Assign particle acceleration
                float(row[6]), row[7])  #Assign particle mass and name
            particles.append(newBody)
        momentumNet = [0.0 for k in range(3)]
        netMass = 0.0
    for p in particles:
        for k in range(3):
            momentumNet[k] += p.momentum(k)
            netMass += p.mass
    for p in particles:
        for k in range(3):
            p.velocity[0][k] -= momentumNet[k]/netMass
            p.velocity[1][k] = p.velocity[0][k]
#end init_system

def plot_trajectories(file = 'output.txt'):             #Plots trajectories from an output CSV file called output.dat 
    with open(file, 'r') as trajectories:               #Open the output file
        reader = csv.reader(trajectories, delimiter=',')        #Initialize a CSV reader
        firstLine = reader.next() #firstLine stores the first line of the CSV file which contains info about the format of the output
        tMax = int(float(firstLine[0]))  #Read in the total number of time steps in the file
        N = int(firstLine[1])       #Read in the total number of particles
        dt = float(firstLine[2])
        outputFrequency = int(firstLine[3])
        particlePathsx = [[0.0 for i in range(tMax+1)] for j in range(N)]
        particlePathsy = [[0.0 for i in range(tMax+1)] for j in range(N)]
        particlePathsz = [[0.0 for i in range(tMax+1)] for j in range(N)]
        particleNames = ['' for i in range(N)]
        timeSteps = [dt * t * outputFrequency / (2 * math.pi) for t in range(tMax+1)]
        kEnergies = [0.0 for t in range(tMax+1)]
        angMomenta = [0.0 for t in range(tMax+1)]
        symbols = ['k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k']
        markers = ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
        for i in range(tMax+1):    #For each time step...
            for j in range(N):          #For each particle...
                try: row = reader.next()     #Grab the line from the file
                except StopIteration:
                    print "EOF"
                    break
                particlePathsx[j][i] = float(row[1])   #Grab x position of particle
                particlePathsy[j][i] = float(row[2])   #Grab y position of particle
                particlePathsz[j][i] = float(row[3])   #Grab z position of particle
                kEnergies[i] += float(row[10])
                angMomenta[i] += float(row[11])
                if (i == 0): particleNames[j] = row[0]
        mpl.rcParams['legend.fontsize'] = 10
#        fig1 = plotter.figure(1) #Figure for plot of 3D motion
        fig2 = plotter.figure(2) #Figure for plot of energy
        p2 = fig2.add_subplot(1, 1, 1)
        deltaE = kEnergies[0]
        deltaL = angMomenta[0]
        for t in range(tMax + 1): 
            kEnergies[t] -= deltaE* (2*math.pi)**2
            angMomenta[t] -= deltaL*2*math.pi
        p2.plot(timeSteps, kEnergies, 'k.')
        plotter.xlabel(r"Time (Years)", fontsize = 16)
        plotter.ylabel(r"$\Delta$ Kinetic Energy (Solar Masses AU$^2$ / Years$^2$)", fontsize = 16)
        
        fig3 = plotter.figure(3) #Figure for plot of energy
        p3 = fig3.add_subplot(1, 1, 1)
        p3.plot(timeSteps, angMomenta, 'k.')
        plotter.xlabel(r"Time (Years)", fontsize = 16)
        plotter.ylabel(r"$\Delta$ Angular Momentum (Solar Masses AU$^2$ / Years)", fontsize = 16)
        plotter.show()  
        #fig1.clf()
        #p = Axes3D(fig1)
        #for j in range(N):
        #    p.plot3D(particlePathsx[j], particlePathsy[j], particlePathsz[j], zdir = 'z', c = symbols[j], marker=markers[j], label = particleNames[j])
        #p.set_xlabel("x")
        #p.set_ylabel("y")
        #p.set_zlabel("z")
        #zoom_factory(p, base_scale = 1.5)
        plotter.show()
        fig4 = plotter.figure(4)
        p4 = fig4.add_subplot(1, 1, 1)
        for j in range(N):
            p4.plot(particlePathsx[j], particlePathsy[j], c = 'k', marker = markers[j], label = particleNames[j])
        #plotter.legend() ### LEGEND ###
        plotter.xlabel("x-position relative to barycenter (AU)", fontsize = 16)
        plotter.ylabel("y-position relative to barycenter (AU)", fontsize = 16)
        #plotter.legend(loc = 2)

def zoom_factory(ax,base_scale = 2.):   #Not written by me. From user tcaswell on stackoverflow.
    def zoom_fun(event):
        # get the current x and y limits
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
        cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location
        if event.button == 'up':
            # deal with zoom in
            scale_factor = 1/base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            print event.button
        # set new limits
        ax.set_xlim([xdata - cur_xrange*scale_factor,
                     xdata + cur_xrange*scale_factor])
        ax.set_ylim([ydata - cur_yrange*scale_factor,
                     ydata + cur_yrange*scale_factor])
        plotter.draw() # force re-draw

    fig = ax.get_figure() # get the figure of interest
    # attach the call back
    fig.canvas.mpl_connect('scroll_event',zoom_fun)

    #return the function
    return zoom_fun

def run_system(particles, tFinal, dt, outputFrequency, inputFile):
    init_system(inputFile)
    nBodyRunner(particles, tFinal, dt, outputFrequency)
    plot_trajectories()