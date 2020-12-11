import math

h = 0.000001
t0 = 0
tmax = 0.5
I0 = 0
q0 = 0
R = 10
L = 1
C = 1*math.pow(10,-5)
beta = 5
omega0 = 316.23
T0 = 1.99*math.pow(10,-2)
omega = 316.0
T = 1.99*math.pow(10,-2)
Emax = 10
Rcr = 632.46

def ODEfunc(t,u):
    y,z = u
    return [ z, -2*beta*z-math.pow(omega0,2)*y+Emax*math.cos(omega*t) ]

def EulerStep(f,t,u,h):
    f1 = f(t,u)
    k1 = [h*f1[0],h*f1[1]]
    f2 = f(t+h, [k1[0]+u[0],k1[1]+u[1]])
    k2 = [h*f2[0],h*f2[1]]
    k1k2 = [0.5*(k1[0]+k2[0]),0.5*(k1[1]+k2[1])]
    return [u[0]+k1k2[0],u[1]+k1k2[1]]

def Euler(f,a,b,x0,h):
    t = a
    x = x0[:]
    fileObj = open("result.txt", "w")
    while t < b+0.5*h:
        x = EulerStep(f,t,x,h)
        fileObj.write("{0} {1} {2} {3} {4}".format(x[0], x[1], x[1]/R, Emax*math.cos(omega*t),'\n'))
        t=t+h
    fileObj.close()

Euler(ODEfunc,t0,tmax,[0.,0.],h)
