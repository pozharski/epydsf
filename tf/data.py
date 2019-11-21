from scipy import exp

class curve(object):
    def __init__(self, p):
        self.setit()
        self.validate(p)
        self.p = p
    def setit(self):
        raise Exception("Undefined method setit")
    def validate(self, p):
        raise Exception("Undefined method validate")
    def values(self, x):
        raise Exception("Undefined method values")

class parcurve(curve):
    def validate(self, p):
        if len(p) != self.dimN:
            raise ValueError('Baseline requries %d parameters, %d supplied' % (self.dimN,len(p)))

class baseline(parcurve):
    def setit(self):
        self.dimN = 2
    def values(self, x):
        return self.p[0]+self.p[1]*x

class melt(parcurve):
    def setit(self):
        self.dimN = 2
    def values(self, x):
        tm, dt = self.p
        k = exp(4*tm**2*(1/tm-1/x)/dt)
        return k/(1+k)
        
