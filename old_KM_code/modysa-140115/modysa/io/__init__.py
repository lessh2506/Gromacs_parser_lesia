import re, glob
import numpy as N
import modysa.statistics as STAT

from gromacs import xtcio as X

class DataFile:

    def __init__(self, filename):
        data = [line.strip() for line in open(filename,'r').readlines() \
                    if line.strip()]
        self.labels = data[0].split()
        self.data = []
        for d in data[1:]:
            numb = map(float,d.split())
            if len(numb) != len(self.labels):
                raise IOError, \
                    ' Inconsistent data in the following line:\n%s' % d
            self.data.append(numb)
        self.data = N.array(self.data)

    def __call__(self, cbeg=1, cend=None, factor=1, rbeg=0, rend=None):
        # [beg, end)
        assert rend is not None and rbeg < rend or rbeg < len(self.data) \
            and rbeg >= 0
        labels = self.labels[cbeg:cend]
        print '     columns: ', '%12s'*len(labels) % tuple(labels)
        data = N.add.reduce(self.data[rbeg:rend, cbeg:cend], -1) * factor
        print ' mean/stddev: %12g +/- %-12g' % \
            (STAT.mean(data), STAT.standardDeviation(data))

    def filter(self, cRE):
        # find labels matching cRE, return indices
        pass

class XVGFile:

    comment_repat = re.compile('^#.*$', re.M)
    header_repat = re.compile('^@.*$', re.M)
    data_repat = re.compile('^\s*\d+\.?\d+.+$', re.M)
    number_repat = re.compile('-?\d+\.?\d+')

    def __init__(self, filename):
        data = open(filename,'r').read()
        self.comment = re.findall(self.comment_repat, data)
        self.header = re.findall(self.header_repat, data)
        self.data = re.findall(self.data_repat, data)

    def getArray(self, skip=None):
        return N.array([map(float, re.findall(self.number_repat, line)) \
                            for line in self.data])

    def getComment(self):
        return '\n'.join(self.comment)

    def getHeader(self):
        return '\n'.join(self.header)

    def __len__(self):
        return len(self.data)

    def validateData(self, step):
        # check for presence of duplicated or skipped time step entries
        pass



class DataBlock:

    def __init__(textfile,recordspec,recordnum):
        """ recordspec: fortran like record format description,
        F for real/double/float, I for integer/long, S for string are
        allowed, e.g. '2F5, 4I6, A10, F2' """

        import caux # this is reminder that caux module needs to be double chkd
        self.data = None # here we keep data read from file
        FIELDPATTERN = '^\s*(\d*)([FIA])(\d+)$'
        if recordnum == 0:
            gj = textfile.readline()
            self.data = []
        check = [re.search(FIELDPATTERN,x) for x in recordspec.split(',')]
        if check.count(None):
            raise ' Bad recordspec: %s' % recordspec
        fields = []
        for m,k,s in [x.groups() for x in check]: # multi,kind,size
            if m: m = string.atoi(m)
            else: m = 1
            fields += [(k,string.atoi(s)) for i in range(m)]
        data = caux.dataBlock(textfile.file,tuple(fields),recordnum)
        if type(data) == type([]): self.data = data
        raise NotImplementedError, 'DataBlockError: %s' % data



if __name__ == '__main__':

    pass
    # t = Trajectory('pope-2.part000[1-3].xtc')
    # print t.first, t.last
    # fcoords = t.getConfiguration()
    # fbox = t.getBox()
    # ftime = t.getTime()
    # t.limit(10000,30000) # we limit reading trajectory from 10ns to 30ns
    #     filename = "/home/murzyn/Projects/Papers/2010-dubai/data/lipAH-EpotTBox.xvg"
    #     xvg = XVGFile(filename)
    #     data = xvg.getArray()
    #     # mdutils
    #     import statistics as S 
    #     print S.mean(data[:,1]), S.standardDeviation(data[:,1])
    #     # data = DataFile('data/hb_pmb-la.dat')
    #     # data(cbeg=2, cend=5, factor=1/6., rbeg=79)
    #     # hb_atoms = [(1,2), (2, 5), (5, 8), (8, 11), (11, 14), (14, 17), (17, 19), (19, 21), (21,24), (24, 27), (27, 30), (30, 31)]
    #     # for b, e in hb_atoms: data(cbeg=b, cend=e, factor=1/6., rbeg=79)
    #     # c_atoms = [(1, 9)
