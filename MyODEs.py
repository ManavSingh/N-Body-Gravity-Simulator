## ODE subroutines. Very prototype. Test test test.

def euler(q, f, dt): #Takes q[i-1], f[i-1], and dt. Computes q[i] with Euler Method
    return q + dt*f

#def eulerPC(q, f, dt, step):#Asks for same arguments as Euler method, but also does a corrector step.
##    qPrediction = euler(q, f, dt, step)    #Needed if f is a function of q and t. Not needed here since f is assumed to be an array.
#    qCorrected = q[step] + .5*dt*(f[step] + f[step+1])
#    return qCorrected

def leapfrog(q, f, dt): #Takes q[i-1], f[i], and dt. Computes q[i+1] using 2nd order Leapfrog
    return q + 2*dt*f