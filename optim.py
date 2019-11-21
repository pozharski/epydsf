from scipy import exp, mean, log, ones, polyfit
from scipy.optimize import fmin
from matplotlib.pyplot import plot
import sys, os
class optimizer:
    def __init__(self, x, y, w=None, p=None):
        self.x = x
        self.y = y
        if w is None:
            self.w = ones(self.x.shape)
        else:
            self.w = w
        if p is None:
            self.guess_params()
        else:
            self.p = p
    def guess_params(self):
        pass
    def f(self,p):
        pass
    def lsq(self,p):
        return sum(self.w*((self.f(p)-self.y)**2))
    def fmin(self):
        fsock = open(os.devnull,'w')
        fbck = sys.stdout
        sys.stdout = fsock
        self.p = fmin(self.lsq, self.p)
        fsock.close()
        sys.stdout = fbck
        return self.p
    def plot(self, dstyle='ro', fstyle='g--'):
        plot(self.xplt(),self.y,dstyle,self.xplt(),self.f(self.p),fstyle)
    def plotnorm(self, dstyle='ro', fstyle='g--'):
        plot(self.xplt(),self.fnorm(self.y),dstyle,self.xplt(),self.fnorm(self.f(self.p)),fstyle)
    def report(self):
        return '%.f '*len(self.p) % self.p
    def getp(self):
        return self.p
    def xplt(self):
        return self.x
    def gety(self):
        return self.f(self.p)
    def getnormy(self, y=None):
        if y is None:
            return self.fnorm(self.f(self.p))
        else:
            return self.fnorm(y)
    def fnorm(self, y, p=None):
        pass
    def setw(self, w=None):
        if w is None:
            self.w = ones(self.x.shape)
        else:
            self.w = w

class linmelt(optimizer):
    def __init__(self, x, y, w=None, p=None, tm=330):
        if p is not None:
            p[0] = p[0]+273.15
        self.tm_guess=tm
        optimizer.__init__(self, x+273.15, y, w, p) 
    def f(self,p):
        tm, dt, a, b, c, d = p
        k = exp(4*tm**2*(1/tm-1/self.x)/dt)
        return (a+b*self.x+(c+d*self.x)*k)/(1+k)
    def fnorm(self, y, p=None):
        if p is not None:
            a, b, c, d = p[-4:]
        else:
            a, b, c, d = self.p[-4:]
        return (y-c-d*self.x)/(a-c+(b-d)*self.x)
    def guess_params(self):
        ind = self.w.astype(bool)
        self.p = [self.tm_guess, 10] + polyfit(self.x[ind][:10],self.y[ind][:10],1)[::-1].tolist() + polyfit(self.x[ind][-10:],self.y[ind][-10:],1)[::-1].tolist()
    def report(self):
        return 'Tm = %.1f deltaT = %.1f' % (self.p[0]-273.15, self.p[1])
    def getp(self):
        return [self.p[0]-273.15] + self.p[1:].tolist()
    def xplt(self):
        return self.x-273.15
    def getG(self, t):
        tm, dt = self.p[:2]
        tk = t+273.15
        return 4*tm**2*(1/tm-1/tk)/dt
    def tm(self):
        return self.p[0]-273.15
    def dt(self):
        return self.p[1]
    def w2delta(self, kfac=4):
        wt = ones(self.x.shape)
        dt = max(self.p[1],10)
        x1 = self.p[0]-kfac*dt
        x2 = self.p[0]+kfac*dt
        wt = wt*(self.x>=x1).astype(float)
        wt = wt*(self.x<=x2).astype(float)
        self.setw(wt)

class linmeltcp(linmelt):
    def f(self,p):
        tm, dt, kappa, a, b, c, d = p
        k = exp(4*tm*(kappa*log(self.x/tm)+(1-kappa)*(1-tm/self.x))/dt)
        return (a+b*self.x+(c+d*self.x)*k)/(1+k)
    def guess_params(self):
        linmelt.guess_params(self)
        self.p.insert(2,0)
    def getG(self, t):
        tm, dt, kappa = self.p[:3]
        tk = t+273.15
        return 4*tm*(kappa*log(tk/tm)+(1-kappa)*(1-tm/tk))/dt
    def report(self):
        return 'Tm = %.1f deltaT = %.1f dCp = %.1f' % (self.p[0]-273.15, self.p[1], 4*1.98e-3*self.p[2]*(273.15+self.p[0])/self.p[1])
    
class linmelt_quadb(linmelt):
    def f(self,p):
        tm, dt, a1, b1, c1, a2, b2, c2 = p
        k = exp(4*tm**2*(1/tm-1/self.x)/dt)
        return (self.bline(a1,b1,c1)+self.bline(a2,b2,c2)*k)/(1+k)
    def bline(self, a, b, c):
        return a+b*self.x+c*self.x**2
    def fnorm(self, y, p=None):
        if p is not None:
            a1, b1, c1, a2, b2, c2 = p[-6:]
        else:
            a1, b1, c1, a2, b2, c2 = self.p[-6:]
        y1 = self.bline(a1,b1,c1)
        y2 = self.bline(a2,b2,c2)
        return (y-y2)/(y1-y2)
    def guess_params(self):
        linmelt.guess_params(self)
        self.p.insert(4,0)
        self.p.append(0)
    
class linmeltcp_quadb(linmelt_quadb):
    def f(self,p):
        tm, dt, a1, b1, c1, a2, b2, c2 = p
        k = exp(4*tm*(kappa*log(self.x/tm)+(1-kappa)*(1-tm/self.x))/dt)
        return (self.bline(a1,b1,c1)+self.bline(a2,b2,c2)*k)/(1+k)
    def guess_params(self):
        linmelt.guess_params(self)
        self.p.insert(2,0)
        self.p.insert(4,0)
        self.p.append(0)
    def getG(self, t):
        tm, dt, kappa = self.p[:3]
        tk = t+273.15
        return 4*tm*(kappa*log(tk/tm)+(1-kappa)*(1-tm/tk))/dt
    def report(self):
        return 'Tm = %.1f deltaT = %.1f dCp = %.1f' % (self.p[0]-273.15, self.p[1], 4*1.98e-3*self.p[2]*(273.15+self.p[0])/self.p[1])
