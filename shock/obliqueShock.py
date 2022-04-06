
import normalShock as ns
import numpy

def findAngle(g, beta, tRef, r, p, u):

    error = 1
    while(abs(error)>1e-9):

        s0 = oblique(g, r, p, u, beta)
        s1 = oblique(g, r, p, u, beta+0.1)
        
        t0 = s0.theta
        t1 = s1.theta

        error = - (t0-tRef)*0.1/(t1 - t0)
        beta += error
        print(beta)

    return beta

        


class oblique():

    def __init__(self, gas, r1, p1, u1, beta):

        self.gas = gas
        self.r1 = r1
        self.p1 = p1
        self.u1 = u1
        self.beta = beta

        u1t = self.u1*numpy.cos(self.beta)
        u1n = self.u1*numpy.sin(self.beta)

        s = ns.shock(gas, r1 , p1, u1n)

        s.calc()

        self.r2 = s.r2
        self.p2 = s.p2
        u2n = s.u2
        u2t = u1t

        u2 = numpy.sqrt(u2n**2 + u2t**2)

        self.a1 = s.a1
        self.a2 = s.a2

        self.T1 = s.T1
        self.T2 = s.T2

        self.m1 = u1/self.a1
        self.m2 = u2/self.a2

        self.h1 = s.h1
        self.h2 = s.h2

        bt = numpy.arctan(u2n/u2t)

        self.theta = self.beta - bt
        



if __name__=="__main__":

    g = ns.GAS('ideal')    

    # checagem para gas ideal

    s = oblique(g, 0.1 , 1e5, 15000, numpy.pi/4)

    print('mach', s.m1) 
    print('theta', s.theta*180/numpy.pi)
    print('mach2', s.m2, 2.23656895)
    print('pressureRatio', s.p2/s.p1, 93.5833333)
    print('densityRatio', s.r2/s.r1, 5.64853556)
    print('temperatureRatio',s.T2/s.T1, 16.5677160)
    print('entalpyRatio', s.h2/s.h1)

    g = ns.GAS('real')    
    
    r = 0.1
    p = 1e5
    u = 11000

    beta = findAngle(g, numpy.pi/4, 20*numpy.pi/180, r, p, u)
    s = oblique(g, r, p, u, beta)

    print('mach', s.m1) 
    print('theta', s.theta*180/numpy.pi)
    print('mach2', s.m2)
    print('pressure2', s.p2)
    print('density2', s.r2)
    print('temperature2',s.T2)

    r = 0.0180119
    p = 1171.87
    u = 30*301.803

    theta = 20*numpy.pi/180

    g = ns.GAS('ideal')

    beta = findAngle(g, 2*theta, theta, r, p, u)
    sI = oblique(g, r, p, u, beta)

    g = ns.GAS('real')

    beta = findAngle(g, 2*theta, theta, r, p, u)
    sR = oblique(g, r, p, u, beta)

    print('mach1', sI.m1, sR.m1) 
    print('theta', sI.theta*180/numpy.pi, sR.theta*180/numpy.pi)
    print('beta', sI.beta*180/numpy.pi, sR.beta*180/numpy.pi)
    print('mach2', sI.m2, sR.m2)
    print('pressure2', sI.p2, sR.p2)
    print('density2', sI.r2, sR.r2)
    print('temperature2', sI.T2, sR.T2)

