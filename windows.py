import sys
from matplotlib.figure import Figure
from matplotlib.pyplot import grid, close, title
from optim import linmelt
from scipy import array, argmin
class DSFWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
    def onkeypress(self, event):
        if event.key == 'q':
            close(self)
            sys.exit()
        if event.key == 'right':
            self.well = self.wells[min(self.wells.index(self.well)+1, len(self.wells)-1)]
        if event.key == 'left':
            self.well = self.wells[max(self.wells.index(self.well)-1, 0)]
        if event.key == 'up':
            self.well = self.wells[max(self.wells.index(self.well)-24, 0)]
        if event.key == 'down':
            self.well = self.wells[min(self.wells.index(self.well)+24, len(self.wells)-1)]
        if event.key == 'F':
            if self.wtfit[self.well] is None:
                self.wtfit[self.well] = linmelt(self.wtf[self.well][0], self.wtf[self.well][1], tm=self.tm_guess)
                p = self.wtfit[self.well].fmin()
            self.wtfit[self.well].plot(dstyle='r.')
            titline ='Well #%d'% (self.well)
            if self.well_info:
                titline += ' ('+self.well_info.get_well_value(self.well, 'info')+')'
            titline += ' %s' % (self.wtfit[self.well].report())
            title(titline)
            event.canvas.draw()
        if event.key == 'W':
            if self.wtfit[self.well] is None:
                self.wtfit[self.well] = linmelt(self.wtf[self.well][0], self.wtf[self.well][1], tm=self.tm_guess)
                p = self.wtfit[self.well].fmin()
            self.wtfit[self.well].w2delta(kfac=self.kfac)
            p = self.wtfit[self.well].fmin()
            self.wtfit[self.well].plot(dstyle='r.')
            title('Well #%d: %s' % (self.well,self.wtfit[self.well].report()))
            event.canvas.draw()
        if event.key in ['right','left','up','down']:
            self.plot_well()
    def set_data(self, data, well_info, tm_guess=320.0, wtfit=None, kfac=4.0):
        self.wtf = data
        self.wells = sorted(self.wtf)
        self.set_well(self.wells[0])
        if wtfit is None:
            self.wtfit = dict([(key,None) for key in list(data.keys())])
        else:
            self.wtfit = wtfit
        self.well_info = well_info
        self.tm_guess = tm_guess
        self.kfac=kfac
    def set_well(self, well):
        if well in self.wells:
            self.well = well
    def plot_well(self, well=None):
        if well is None:
            well = self.well
        else:
            self.well = well
        self.add_subplot(111)
        self.axes[0].clear()
        self.pt = self.axes[0].plot(self.wtf[well][0], self.wtf[well][1])
        titline = 'Well #%d' % self.well
        if self.well_info:
            titline += ': '+self.well_info.get_well_value(self.well, 'info')
        if self.wtfit[self.well] is not None:
            titline += ' %s' % (self.wtfit[self.well].report())
            self.wtfit[self.well].plot(dstyle='r.')
        title(titline)
        grid(True)
        self.canvas.draw()
        
class PlateResWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
    def set_data(self, data):
        self.data = data
    def attach_curves(self, figCurve):
        self.figCurve = figCurve
    def plot(self, what, logflag):
        if what in self.data:
            self.add_subplot(111)
            self.axes[0].clear()
            x = array(self.data[what]['x'])
            y = array(self.data[what]['y'])
            fmts = array(self.data[what]['format'])
            fs = set(fmts)
            self.pt = {}
            for f in fs:
                ind = fmts==f
                xx = x[ind]
                yy = y[ind]
                if logflag:
                    self.pt[f] = [self.axes[0].semilogx(xx,yy,f)[0], array(self.data['wells'])[ind]]
                else:
                    self.pt[f] = [self.axes[0].plot(xx,yy,f)[0], array(self.data['wells'])[ind]]
            if 'my' in self.data[what] and 'sy' in self.data[what]:
                self.pmt = self.axes[0].axhspan(self.data[what]['my']-self.data[what]['sy'],self.data[what]['my']+self.data[what]['sy'],facecolor='0.5', alpha=0.5)
            grid(True)
    def set_axlabels(self, xlabel=None, ylabel=None):
        if self.axes:
            if xlabel is not None:
                self.axes[0].set_xlabel(xlabel)
            if ylabel is not None:
                self.axes[0].set_ylabel(ylabel)
    def onkeypress(self, event):
        if event.key == 'q':
            close(self)
            sys.exit()
    def onmouseclick(self, event):
        if event.dblclick:
            if event.xdata is not None:
                if 'figCurve' in dir(self):
                    d2, wells = list(zip(*[(((array(pt.get_data()).T-array([event.xdata,event.ydata]))**2).sum(1).tolist(),wells.tolist()) for pt,wells in list(self.pt.values())]))
                    d2, wells = sum(d2,[]), sum(wells,[])
                    self.figCurve.plot_well(wells[argmin(d2)])
