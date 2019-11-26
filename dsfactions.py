# -*- coding: utf-8 -*-
from windows import DSFWindow, PlateResWindow
from epyparser import viia_parser, exparser
from matplotlib.pyplot import figure, show
from optim import linmelt
from scipy import array, sqrt
import csv

def info(args):
    if args.csv_wells is not None:
        well_info = exparser(args.csv_wells,args.csvexp,args.csvregexp)
        expnames = well_info.get_experiments()
        print('Listed experiments:')
        print('\n'.join(['%30s (%3d wells)' % (k,v) for k,v in expnames.items()]))

def plot(args):
    dataset = viia_parser(args.input_file)
    fig = figure(FigureClass=DSFWindow)
    well_info = None
    if args.csv_wells is not None:
        well_info = exparser(args.csv_wells,args.csvexp,args.csvregexp)
        args.wells = well_info.get_wells()
    fig.set_data(dataset.get_all_readings(args.wells), well_info, tm_guess=args.tm_guess, kfac=args.basespan)
    fig.plot_well()
    fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
    show()
   
def fit(args):
    dataset = viia_parser(args.input_file)
    wtf = dataset.get_all_readings(args.wells)
    if args.csv_output:
        fout = open(args.output_file,'w')
        fout.write('Well,Tm,deltaT\n')
    if args.wells is not None:
        wells = list(map(int, args.wells.split(','))) 
    elif args.csv_wells is not None:
        wells = exparser(args.csv_wells,args.csvexp,args.csvregexp).get_wells()
    else:
        wells = sorted(wtf)
    for well in wells:
        wtfit = linmelt(wtf[well][0], wtf[well][1], tm=args.tm_guess)
        p = wtfit.fmin()
        if args.csv_output:
            fout.write('%d,%f,%f\n' % (well,wtfit.tm(),wtfit.dt()))
        print('Well #%d: %s' % (well,wtfit.report()))
    if args.csv_output:
        fout.close()

def plate(args):
    dataset = viia_parser(args.input_file)
    if args.cwells is not None or (args.csvcontrol is not None and args.csv_wells is not None):
        if args.cwells is not None:
            cwells = list(map(int, args.cwells.split(','))) 
        else:
            cwells = exparser(args.csv_wells,args.csvcontrol).get_wells()
        cwtf = dataset.get_all_readings(cwells)
        c_parms = []
        for well in cwells:
            wtfit = linmelt(cwtf[well][0], cwtf[well][1], tm=args.tm_guess)
            p = wtfit.fmin()
            for i in range(args.nbro):
                wtfit.w2delta(kfac=args.basespan)
                p = wtfit.fmin()
            print('Control well #%03d: %s' % (well,wtfit.report()))
            c_parms.append([wtfit.tm(),wtfit.dt()])
        c_parms = array(c_parms)
        mtm, mdt = c_parms.mean(0)
        stm, sdt = c_parms.std(0,ddof=1)
        print('--------\nAverages:')
        print('          Tm = %.1f +- %.1f' % (mtm, stm))
        print('      deltaT = %.1f +- %.1f' % (mdt, sdt))
        contrflag = True
    else:
        contrflag = False
    print('--------\nResults:')
    wtf = dataset.get_all_readings(args.wells)
    if args.wells is not None:
        wells = list(map(int, args.wells.split(','))) 
    elif args.csv_wells is not None:
        well_info = exparser(args.csv_wells,args.csvexp,args.csvregexp)
        wells = well_info.get_wells()
    else:
        wells = sorted(wtf)
    wtfits = {}
    a_parms = []
    for well in wells:
        if contrflag:
            wtfit = linmelt(wtf[well][0], wtf[well][1], tm=273.15+mtm)
        else:
            wtfit = linmelt(wtf[well][0], wtf[well][1], tm=args.tm_guess)
        p = wtfit.fmin()
        for i in range(args.nbro):
            wtfit.w2delta(kfac=args.basespan)
            p = wtfit.fmin()
        outline = 'Well #%03d: %s' % (well,wtfit.report())
        if contrflag:
            outline += ' ZTm=%8.1f' % ((wtfit.tm()-mtm)/stm/sqrt(1+len(cwells)))
        if args.csv_wells is not None:
            outline += ' : ' + well_info.get_well_value(well, 'info')
            well_info.set_well_value(well, 'tm', wtfit.tm())
            well_info.set_well_value(well, 'dt', wtfit.dt())
        print(outline)
        wtfits[well] = wtfit
        a_parms.append([wtfit.tm(),wtfit.dt()])
    a_parms = array(a_parms)
    mtm, mdt = a_parms.mean(0)
    stm, sdt = a_parms.std(0,ddof=1)
    print('--------\nAverages:')
    print('          Tm = %.1f +- %.1f' % (mtm, stm))
    print('      deltaT = %.1f +- %.1f' % (mdt, sdt))
    if args.csv_wells is not None:
        x,tm,dt,fmt,wellnum = list(zip(*[(float(v.get('x')),v.get('tm'),v.get('dt'),v.get('format','ro'),k) for k,v in well_info.iteritems()]))
        if args.output_file is not None:
            with open(args.output_file,'w') as fout:
                if args.csv_output:
                    fout.write('Well,Tm,deltaT\n')
                    for xx,yy,zz in zip(*(x,tm,dt)):
                        fout.write('%f,%f,%f\n' % (xx,yy,zz))
                else:
                    for xx,yy in zip(*(x,tm)):
                        fout.write('%f %f\n' % (xx,yy))
        fig = figure(FigureClass=PlateResWindow)
        if contrflag:
            fig.set_data({'tm':{'x':x,'y':tm,'format':fmt,'my':mtm,'sy':stm}, 'wells':wellnum})
        else:
            fig.set_data({'tm':{'x':x,'y':tm,'format':fmt}, 'wells':wellnum})
        fig.plot('tm', args.logplot)
        if args.ylabel is None:
            args.ylabel = "Tm, °C"
        fig.set_axlabels(args.xlabel, args.ylabel)
        fig.canvas.mpl_connect('key_press_event', fig.onkeypress)
        fig.canvas.mpl_connect('button_press_event', fig.onmouseclick)
        figCurve = figure(FigureClass=DSFWindow)
        figCurve.set_data(dataset.get_all_readings(wells), well_info, wtfit=wtfits, kfac=args.basespan)
        figCurve.canvas.mpl_connect('key_press_event', figCurve.onkeypress)
        fig.attach_curves(figCurve)
        show()
