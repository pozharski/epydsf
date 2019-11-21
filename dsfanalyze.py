#!/usr/bin/env python
headerhelp = \
'''
DSFANALYZE analyzes DSF data
'''
from argparse import ArgumentParser, RawDescriptionHelpFormatter
parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
parser.add_argument('datafile',
                    help='Raw data file name.  Standard format is a text file with a header row and individual columns with tempretaure and sample readings.')
parser.add_argument('--cols',
                    help='Comma-separated list of columns (numbers, starting with first non-temperature column) to read from the file.  By default, all columns are imported.')
parser.add_argument('--ncyc', type=int, default=2,
                    help='Number of optimization cycles, defaults to 2.')
parser.add_argument('--t1', type=float,
                    help='Low temperature cutoff.')
parser.add_argument('--t2', type=float,
                    help='High temperature cutoff.')
parser.add_argument('--show-plots', action='store_true',
                    help='Show infdividual plots')
parser.add_argument('--title', help='Title of the summary figure')
parser.add_argument('--subtract-baseline', action='store_true',
                    help='Subtract baseline when plotting summary results')
parser.add_argument('--skiprows', type=int, default=1,
                    help='Number of rows to skip from the top of the data file')
parser.add_argument('--attempt-cp', action='store_true',
                    help='Try fitting to the equation that takes heat capacity into account')
parser.add_argument('--quadratic-baseline', action='store_true',
                    help='Fit baseline signals to quadratic polynomial')
parser.add_argument('--excludeT', 
                    help='Exclude temperature range (comma-separated, multiple ranges can be separated by underscore)')
parser.add_argument('--plot-partial', action='store_true',
                    help='Only plot the part of the curve restricted by t1/t2')
parser.add_argument('--tm-guess', default=57.0, type=float,
                    help='Initial guess for melting temperature, defaults to 57 C')
parser.add_argument('--plot-colors',
                    help='Specify colors for each curve.  Comma-separated list of the same length as cols parameter')
args = parser.parse_args()

from scipy import loadtxt, array, mean, logical_or, ones, arange
from matplotlib.pyplot import show, grid, figure, title, plot, cla, xlabel, ylabel
if args.quadratic_baseline:
    from optim import linmelt_quadb as linmelt
    if args.attempt_cp:
        from optim import linmeltcp_quadb as linmeltcp
else:
    from optim import linmelt
    if args.attempt_cp:
        from optim import linmeltcp

args.tm_guess += 273.15

rawx = loadtxt(args.datafile, skiprows=args.skiprows)
T = rawx[:,0]
params, fc, fo, g25s, g37s = [], [], [], [], []
if args.attempt_cp:
    paramscp, fc_cp, g25s_cp, g37s_cp = [], [], [], []
if args.cols is None:
    cols = arange(1,rawx.shape[1])
    f = rawx[:,1:]
else:
    cols = array(args.cols.split(',')).astype(int)
    f = rawx[:,cols]
wt = ones(T.shape)
if args.t1 is not None:
    wt = wt*(T>=args.t1).astype(float)
if args.t2 is not None:
    wt = wt*(T<=args.t2).astype(float)
wt1 = wt
if args.excludeT is not None:
    for xrng in args.excludeT.split('_'):
        xt1, xt2 = list(map(float,xrng.split(',')))
        wt = wt*(logical_or(T<=xt1, T>=xt2).astype(float))
for (col,ff) in enumerate(f.T):
    fopt = linmelt(T, ff, wt)
    for i in range(args.ncyc):
        fopt.fmin()
    params.append(fopt.getp())
    print("Column %d: " % cols[col] + fopt.report())
    if args.show_plots:
        cla()
        grid(True)
        fopt.plot(dstyle='r.')
        show()
    g25s.append(fopt.getG(25))
    g37s.append(fopt.getG(37))
    if args.attempt_cp:
        pcp = fopt.getp()
        pcp.insert(2, 0)
        fcp = linmeltcp(T, ff, wt, pcp)
        for i in range(args.ncyc):
            fcp.fmin()
        paramscp.append(fcp.getp())
        if args.show_plots:
            fcp.plot(fstyle='b--', dstyle='r.')
            title('Column %d' % cols[col])
            show()
        print("Column %d: " % cols[col] + fcp.report())
        g25s_cp.append(fcp.getG(25))
        g37s_cp.append(fcp.getG(37))
    if args.subtract_baseline:
        fc.append(fopt.getnormy())
        if args.attempt_cp:
            fc_cp.append(fcp.getnormy())
        fo.append(fopt.getnormy(ff))
    else:
        fc.append(fopt.gety())
        if args.attempt_cp:
            fc_cp.append(fcp.gety())
        fo.append(ff)

pm = array(params).mean(0)
ps = array(params).std(0, ddof=1)
h = 4*1.98e-3*(273.15+array(params)[:,0])**2/array(params)[:,1]
print('Averaged Tm = %.1f +- %.1f' % (pm[0], ps[0]))
print('Averaged delta T = %.1f +- %.1f' % (pm[1], ps[1]))
print('Averaged delta H = %.1f +- %.1f' % (h.mean(), h.std(ddof=1)))
g25s = array(g25s)
g37s = array(g37s)
print('Averaged delta G37 = %.1f +- %.1f' % (g37s.mean(), g37s.std(ddof=1)))
print('Averaged delta G25 = %.1f +- %.1f' % (g25s.mean(), g25s.std(ddof=1)))
print('Recommended temperature range: %d - %d' % (pm[0]-min(3*pm[1],20), pm[0]+min(3*pm[1],20)))

if args.attempt_cp:
    cp_pm = array(paramscp).mean(0)
    cp_ps = array(paramscp).std(0, ddof=1)
    cp_h = 4*1.98e-3*(273.15+array(paramscp)[:,0])**2/array(paramscp)[:,1]
    dcp = array(paramscp)[:,2]*cp_h/(273.15+array(paramscp)[:,0])
    print('Model adjusted for Cp')
    print('Averaged       Tm = %.1f +- %.1f' % (cp_pm[0], cp_ps[0]))
    print('Averaged delta T  = %.1f +- %.1f' % (cp_pm[1], cp_ps[1]))
    print('Averaged delta H  = %.1f +- %.1f' % (cp_h.mean(), cp_h.std(ddof=1)))
    print('Averaged delta Cp = %.1f +- %.1f' % (dcp.mean(), dcp.std(ddof=1)))
    g25s_cp = array(g25s_cp)
    g37s_cp = array(g37s_cp)
    print('Averaged delta G37 = %.1f +- %.1f' % (g37s_cp.mean(), g37s_cp.std(ddof=1)))
    print('Averaged delta G25 = %.1f +- %.1f' % (g25s_cp.mean(), g25s_cp.std(ddof=1)))

cla()
if args.plot_partial:
    ind = wt1.astype(bool)
else:
    ind = ones(T.shape).astype(bool)
if args.plot_colors is not None:
    tt = T[ind]
    fcolors = args.plot_colors.split(',')
    for (col,ff) in enumerate(array(fo)[:,ind]):
        plot(tt, ff, fcolors[col]+'.')
    for (col,ff) in enumerate(array(fc)[:,ind]):
        plot(tt, ff, fcolors[col]+'--')
    if args.attempt_cp:
        for (col,ff) in enumerate(array(fc_cp)[:,ind]):
            plot(tt, ff, fcolors[col]+'-.')
else:
    plot(T[ind],array(fo)[:,ind].T,'.')
    plot(T[ind],array(fc)[:,ind].T,'g--')
    if args.attempt_cp:
        plot(T[ind],array(fc_cp)[:,ind].T,'b--')
grid(True)
title(args.title)
xlabel('Temperature, C')
if args.subtract_baseline:
    ylabel('Fraction of folded protein')
else:
    ylabel('Emission, RFU')
show()
