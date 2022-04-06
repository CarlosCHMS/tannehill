import tannehill as tann
import numpy

class GAS():

    def __init__(self, type):

        self.real = False
        if type == 'real':
            self.real = True

    def pressure(self, r, e):

        if self.real:
            a = tann.pressure(r, e)
        else:
            a = rho*e*(1.4-1)

        return a

    def sound(self, r, e):

        if self.real:
            a = tann.sound(r, e)
        else:
            a = numpy.sqrt(e*(1.4-1)*1.4)

        return a

    def temperature(self, r, e):

        if self.real:
            a = tann.temperature(r, e)
        else:
            a = e*(1.4-1)/287.058

        return a

    def entalpy(self, r, p):

        if self.real:
            a = tann.entalpy(r, p)
        else:            
            a = (p/r)*(1.4/(1.4-1))

        return a

    def gamma_c(self, r, p):

        if self.real:
            a = tann.gamma_c(r, p)
        else:            
            a = 1.4

        return a

class shock():

    def __init__(self, gas, r1, p1, u1):

        self.r1 = r1
        self.p1 = p1
        self.u1 = u1

        g = 1.4

        self.gas = gas
        self.m = r1*u1
        self.P = r1*u1*u1 + p1
        self.h0 = gas.entalpy(r1, p1) + 0.5*u1*u1

    def func(self, r, p):

        g = self.gas.gamma_c(r, p)

        a = 0.5*(g+1)/(g-1)
        b = -g*self.P/(self.m*(g-1))
        c = self.h0

        u2 = (-b - numpy.sqrt(b**2 - 4*a*c))/(2*a)

        r2 = self.m/u2
        p2 = self.P - self.m*u2

        return u2, r2, p2
        
    def calc(self):
    
        r2 = self.r1
        p2 = self.p1

        error = 1
        u2old, r2, p2 = self.func(r2, p2)
        while abs(error)>1e-3:  
            u2, r2, p2 = self.func(r2, p2)
            error = u2-u2old
            u2old = u2

        self.r2 = r2
        self.p2 = p2
        self.u2 = u2
    
        self._post()

        return None

    def _post(self):

        self.h1 = self.gas.entalpy(self.r1, self.p1)
        self.h2 = self.gas.entalpy(self.r2, self.p2)

        self.e1 = self.h1 - self.p1/self.r1
        self.e2 = self.h2 - self.p2/self.r2

        self.a1 = self.gas.sound(self.r1, self.e1)
        self.a2 = self.gas.sound(self.r2, self.e2)

        self.T1 = self.gas.temperature(self.r1, self.e1)
        self.T2 = self.gas.temperature(self.r2, self.e2)

        self.m1 = self.u1/self.a1
        self.m2 = self.u2/self.a2

        

if __name__=="__main__":

    g = GAS('real')    

    # problema anderson hypersonics pag 607

    # altitude 170000ft = 51.816 km
    # velocidade 36000ft/s = 10972.8 m/s

    s = shock(g, 0.000783552 , 60.3609, 10972.8)

    s.calc()

    print(s.m1, s.m2)
    print(s.p2/s.p1, 1387)
    print(s.r2/s.r1, 15.19)
    print(s.h2/s.h1, 212.8)
    print(s.T2/s.T1, 41.64)


