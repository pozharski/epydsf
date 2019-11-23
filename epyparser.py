import re
from csv import DictReader
from scipy import array
from collections import Counter

class parser(object):
    def __init__(self, fname):
        self.fname = fname
        self.parse()
    def parse(self):
        self.blocks = {}
    def _getlines(self):
        with open(self.fname) as fin:
            lines = [x.strip() for x in fin.readlines()]
        return lines

class viia_parser(parser):
    def parse(self):
        lines = self._getlines()
        ptrn_head = re.compile('^\*')
        self.blocks = {'Header':[x[1:] for x in lines if ptrn_head.match(x)]}
        ptrn_blck = re.compile('^\[(.*)\]$')
        blck_rows = [x for x in range(len(lines)) if ptrn_blck.match(lines[x])]+[len(lines)]
        for blockey, i, j in [(lines[blck_rows[i]][1:-1],blck_rows[i]+1,blck_rows[i+1]) for i in range(len(blck_rows)-1)]:
            reader = DictReader(lines[i:j], dialect='excel-tab')
            block = dict([(x,[]) for x in reader.fieldnames])
            for row in reader:
                for key, value in row.items():
                    block[key].append(value)
            self.blocks[blockey] = self.typefix(blockey,block)
    def typefix(self, bname, block):
        if bname == 'Raw Data':
            for key in ['Well','Cycle']:
                block[key] = self.keytypefix(block.get(key), 'int')
            for key in ['Well Position']:
                block[key] = self.keytypefix(block.get(key), 'str')
            for exw in ['1','2','3','4','5','6']:
                for emw  in ['1','2','3','4','5','6']:
                    key = 'x'+exw+'-m'+emw
                    block[key] = self.keytypefix(block.get(key), 'float')
        if bname == 'Melt Curve Raw Data':
            for key in ['Well','Reading']:
                block[key] = self.keytypefix(block.get(key), 'int')
            for key in ['Well Position']:
                block[key] = self.keytypefix(block.get(key), 'str')
            for key in ['Temperature','Fluorescence','Derivative']:
                block[key] = self.keytypefix(block.get(key), 'float')
        if bname == 'Multicomponent Data':
            for key in ['Well','Cycle']:
                block[key] = self.keytypefix(block.get(key), 'int')
            for key in ['ROX','SYPRORANGE']:
                block[key] = self.keytypefix(block.get(key), 'float')
        return block
    def keytypefix(self, values, mthd):
        if values:
            if mthd == 'int':
                return [int(x) for x in values]
            if mthd == 'str':
                return [x.strip() for x in values]
            if mthd == 'float':
                return [float('nan' if x is None else x.translate({44: None})) for x in values]
        return values
    def get_block(self, bname):
        return self.blocks[bname]
    def get_block_column(self, bname, colname):
        return self.blocks[bname][colname]
    def get_well_readings(self, well_number):
        t = array(self.get_block_column('Melt Curve Raw Data','Temperature'))[array(self.get_block_column('Melt Curve Raw Data','Well'))==well_number]
        f = array(self.get_block_column('Multicomponent Data','SYPRORANGE'))[array(self.get_block_column('Multicomponent Data','Well'))==well_number]
        return t,f
        #ind = array(self.get_block_column('Melt Curve Raw Data','Well'))==well_number
        #return array(self.get_block_column('Melt Curve Raw Data','Temperature'))[ind], array(self.get_block_column('Melt Curve Raw Data','Fluorescence'))[ind]
    def get_wells(self):
        return sorted(set(self.get_block_column('Melt Curve Raw Data','Well')))
    def get_all_readings(self, wells=None):
        if wells is None:
            wells = self.get_wells()
        elif type(wells) in [str,str]:
            wells = [int(x) for x in wells.split(',')]
        elif type(wells) is not list:
            raise TypeError('Wells could be provided only as a pre-formed list or comma-separated string.')
        wtf = {}
        for well in wells:
            wtf[well] = self.get_well_readings(well)
        return wtf

class exparser(object):
    def __init__(self, fname, expnames=None, regexpnames=None):
        self.fname = fname
        self.parse(expnames, regexpnames)
    def parse(self, expnames=None, regexpnames=None):
        self.well_info = {}
        with open(self.fname) as fcsv:
            for row in DictReader(fcsv):
                row = dict([(k.lower().strip(),v) for k,v in list(row.items())])
                try:
                    self.well_info[int(row['well'])] = row
                except:
                    print(row)
                    raise
        if expnames is not None:
            if type(expnames) is not list:
                expnames = expnames.split(',')
            self.filter('experiment', expnames)
        if regexpnames is not None:
            self.regfilter('experiment', regexpnames)
    def filter(self, key, ptrns):
        self.well_info = dict([(k,v) for k,v in self.well_info.items() if v[key] in ptrns])
    def regfilter(self, key, ptrns):
        ptrn = re.compile(ptrns)
        self.well_info = dict([(k,v) for k,v in self.well_info.items() if ptrn.match(v[key])])
    def get_wells(self):
        return list(self.well_info.keys())
    def get_well_value(self, well, key):
        return self.well_info[well].get(key)
    def set_well_value(self, well, key, value):
        self.well_info[well][key] = value
    def iteritems(self):
        return iter(self.well_info.items())
    def get_experiments(self):
        return Counter([v.get('experiment') for k,v in self.well_info.items()])
