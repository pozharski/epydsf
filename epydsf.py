#! /usr/bin/env python

import os, dsfactions
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def main():
    headerhelp = \
'''
Analyze DSF data.

Details on actions:

plot        Plot DSF data.  
fit         Fit DSF data to simple melting curve.
plate       Analyze data from a plate experiment.
info        Print information from well info file
'''
    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter,
                        description=headerhelp)
    parser.add_argument('-a', '--action', action='append',
                    choices = [	'plot',
                                'fit',
                                'plate',
                                'info',
                                ],
                    default = [],
                    metavar = '', help='Action to perform')
    parser.add_argument('-i', '--input-file',help='Input data file.')
    parser.add_argument('-o', '--output-file',help='Output data file.')
    parser.add_argument('--csv-output', action='store_true', help='Produce CSV output written into output file')
    parser.add_argument('--wells', help='Comma separated list of wells to process.')
    parser.add_argument('--csv-wells', help='CSV-file with well information.')
    parser.add_argument('--csvexp', help='Name of experiment(s) to load from CSV file, comma-separated.')
    parser.add_argument('--csvregexp', help='Name of experiment(s) to load from CSV file as regular expression') 
    parser.add_argument('--csvcontrol', help='Name of experiment to load as control wells.')
    parser.add_argument('--cwells', help='Comma separated list of control wells.')
    parser.add_argument('--nbro', type=int, default=3, help='Number of cycles of baseline range optimization.')
    parser.add_argument('--tm-guess', type=float, default=50.0, help='Initial estimate of Tm, degrees Celsius')
    parser.add_argument('--basespan', type=float, default=4.0, help='Extent of baseline range, as dT factor')
    parser.add_argument('--xlabel', type=lambda s: unicode(s, 'utf8'), help='X axis label.')
    parser.add_argument('--ylabel', type=lambda s: unicode(s, 'utf8'), help='Y axis label.')
    parser.add_argument('--logplot', action='store_true', help='Logarithmic scale plots')
    parser.add_argument('--dry-run', action='store_true', help='Dry run.')

    args = parser.parse_args()
    args.tm_guess += 273.15

    for action in args.action:
        dsfactions.__getattribute__(action)(args)

if __name__ == "__main__":
    main()
