
from modysa.io import trajectory as T, index as NDX
from modysa.core import geometry as G
from modysa.core import configuration as C
from modysa import statistics as S
from modysa import misc as M

import numpy as N
import os.path as OP

# FIXME: in ndx file add [ molecules ] ??? to select which dihedrals 
# (from the full set defined elsewhere) should be calculated

# FIXME: getConformations mode='path': tagging with numbers 0-7 (eight defined
# conformations, 
# automatic sorting of conformations (a la histogram/searchsorted) for
# arbitrary length of path (now two-dihedrals sequence is hard-coded, i.e.
# gauche-gauche, etc. in a new system: values: -g => 1 vs g- => 

# http://research.stowers-institute.org/efg/R/Statistics/MixturesOfDistributions/index.htm
# http://www.pymix.org/pymix/index.php?n=PyMix.Tutorial#parestim

# from PyMix import mixture
# data = N.loadtxt('phi__400-1199.dat')
# ds = mixture.DataSet()
# a1 = mixture.NormalDistribution(-40,10.)
# a2 = mixture.NormalDistribution(0,8.)
# m = mixture.MixtureModel(2,[.5,.5],[a1,a2])
# xx = [] # density to sample 
# for i,j in data: xx.extend([i]*int(j))
# ds.fromList(xx)
# m.modelInitialization(ds)
# m.EM(ds, 100,.1)
# print m
# clust = m.classify(ds)
# c1 = N.take(xx,N.nonzero(clust))
# print N.mean(c1), N.std(c1)
# c2 = N.take(xx,N.nonzero(clust == 0))
# print N.mean(c2), N.std(c2)


  
EMPTY_TAG = '-'

class StatePath:

    def __init__(self, path=None):
        self.set(path)

    def __update(self, other):
        path = [x is EMPTY_TAG and y or x for x,y in zip(self.path,other)]
        return ''.join(path)
        
    def __add__(self, other):
        return StatePath(self.__update(other))

    def __str__(self):
        return self.path

    def set(self, path):
        if path is None: path = ''
        else: assert type(path) is type('')
        self.path = path

    def update(self, other):
        self.path = self.__update(other)

class DihedralAnalysis:

    named_conformations = {\
        'trans': lambda a: N.logical_or(a < -150, a > 150),
        'gauche_minus': lambda a: N.logical_and(a > -90, a < -30),
        'gauche_plus': lambda a: N.logical_and(a > 30, a < 90),
        'gauche': lambda a: N.logical_and(N.abs(a) > 30, N.abs(a) < 90)
        }

    named_tags = {'trans': 't', 'gauche_minus': 'm',
                  'gauche_plus': 'p', 'gauche': 'g'}

    def getConformations(self, angle, 
                         limits=None, named=None, 
                         mode=None, tag='x'):
        """ return list of values in a given conformation,
        you can either give the name of a particular conformation
        (see named_conformations) or limit dihedral values
        defining a conformation; in the latter case you can set limits with
        either two-element-lists (angle < low or angle > high) or
        two-element-tuple (low < angle < high)
        """
        if type(angle) in [list, tuple]: angle = N.array(angle)
        assert (limits is not None and named is None) or \
            (limits is None and named is not None), \
            'Pick either _limits_ or _named_'
        modes = ['which', 'values', 'path']
        if mode is None: mode = modes[0]
        assert mode in modes
        if named:
            assert self.named_conformations.has_key(named)
            func = self.named_conformations[named]
        else:
            assert len(limits) == 2
            low, high = limits
            if type(limits) == list:
                func = lambda a: N.logical_or(a < low, a > high)
            elif type(limits) == tuple:
                func = lambda a: N.logical_and(a > low, a < high)
            else: func = lambda a: N.ones(len(a), N.int)
        idx = N.flatnonzero(func(angle))
        result = None
        if mode is modes[0]:   result = idx.tolist()
        elif mode is modes[1]: result = N.take(angle, idx).tolist()
        elif mode is modes[2]: 
            result = N.chararray(len(angle),buffer=EMPTY_TAG*len(angle))
            if named: tag = self.named_tags[named]
            N.put(result, idx, tag)
            result = ''.join(result)
        return result

    def mean(self, angle):
        # imaginary
        pass

    def standardDeviation(self, angle):
        pass

    # N.column_stack/vstack



def c5c6(block_steps=5e5,first=200000,last=299999,filename='out.dat',
         named_conformations=None, named_tags=None):

    fn_xtc = '/home/T/lipAH-1a/xtc/lipAH-1a.part*.xtc'
    fn_ndx = OP.join(OP.dirname(__file__), 'data/dihedral.ndx')

    ndx = NDX.IndexFile()
    ndx.read(fn_ndx)

    tors_names = ["alpha:o6'.c6'.c5'.o5'", "beta:o6'.c6'.c5'.c4'"]
    dvec = {}
    for tname in tors_names:
        name, angle = tname.split(':')
        tors = ndx.groups[tname]
        atoms = zip(*tors)
        # renumber indices from 1..n to 0..n-1
        # FIXME: an elegant function needed to renumber 
        #        indices in an array of lists
        for i in range(len(atoms)):
            atoms[i] = [x-1 for x in atoms[i]]
        #
        dvec[name] = [G.DistanceVector(zip(*atoms[a:a+2])) for a in range(3)]

    traj = T.Trajectory(fn_xtc, dt=1, allowed_time_diff=0.1)
    traj.limit(first, last)
    step1st = traj.getStep()

    conf = None
    dihan = DihedralAnalysis()
    if named_conformations:
        dihan.named_conformations.update(named_conformations)
        dihan.named_tags.update(named_tags)

    state_path = StatePath()

    assert not OP.lexists(filename)
    f = open(filename,'w')
    f.write('#  %10s %5s %5s %5s\n' % ('Time [ps]', 'gg', 'gt', 'tg'))
    f.write('# '+'-'*35 + '\n')
    next = True
    result = []
    while next:
        xyz = traj.getConfiguration()
        cell = traj.getCell()
    
        if conf is None: conf = C.Configuration(xyz, cell)
        else: conf.update(xyz, cell)

        conformations = []

        # alpha
        tname = 'alpha'
        dv = dvec[tname]
        [d(conf) for d in dv]
        angles = G.dihedral(dv[0].array, dv[1].array, 
                            dv[2].array)/M.Units.deg
        g = dihan.getConformations(angles, named='agauche', mode='path')
        t = dihan.getConformations(angles, named='atrans', mode='path')
        state_path.set(g)
        state_path.update(t)
        conformations.append(state_path.path)

        # beta
        tname = 'beta'
        dv = dvec[tname]
        [d(conf) for d in dv]
        angles = G.dihedral(dv[0].array, dv[1].array, 
                            dv[2].array)/M.Units.deg
        g = dihan.getConformations(angles, named='bgauche', mode='path')
        t = dihan.getConformations(angles, named='btrans', mode='path')
        state_path.set(g)
        state_path.update(t)
        conformations.append(state_path.path)
        
        # for tname in ['alpha', 'beta']:
        #     dv = dvec[tname]
        #     [d(conf) for d in dv]
        #     angles = G.dihedral(dv[0].array, dv[1].array, 
        #                         dv[2].array)/M.Units.deg
        #     g = dihan.getConformations(angles, named='gauche', mode='path')
        #     t = dihan.getConformations(angles, named='trans', mode='path')
        #     state_path.set(g)
        #     state_path.update(t)
        #     conformations.append(state_path.path)

        union = [''.join(x) for x in zip(*conformations)]
        length = len(union)/100.
        counts = [union.count('gg')/length, 
                  union.count('gt')/length, 
                  union.count('tg')/length]
        result.append([traj.getTime()] + counts)
        step = traj.getStep()
        if len(result) > 1 and (step - step1st) % block_steps == 0:
            data = N.array(result)
            f.write(' %12d %5d %5d %5d\n' % (round(N.mean(data[:,0]),0),
                     round(N.mean(data[:,1]),0), 
                     round(N.mean(data[:,2]),0), 
                     round(N.mean(data[:,3]),0)))
            f.flush()
            result = []
        next = traj.readFrame()

    data = N.array(result)
    f.write(' %12d %5d %5d %5d\n' % (round(N.mean(data[:,0]),0),
             round(N.mean(data[:,1]),0), 
             round(N.mean(data[:,2]),0), 
             round(N.mean(data[:,3]),0)))
    f.close()
    traj.close()


def tors_prof(tors_names=None,block_steps=5e5,first=200000,last=299999,
              filename='out.dat'):

    fn_xtc = '/home/T/lipAH-1a/xtc/lipAH-1a.part*.xtc'
    fn_ndx = OP.join(OP.dirname(__file__), 'data/dihedral.ndx')

    ndx = NDX.IndexFile()
    ndx.read(fn_ndx)

    if tors_names is None:
        #tors_names = ["phi:h1'.c1'.o6.c6", "psi:c1'.o6.c6.c5", "omega:o6.c6.c5.h5"]
        tors_names = ndx.groups.keys()
    TNAMES = [':' in x and x.split(':')[0] or x for x in tors_names]

    for tname in TNAMES:
        out_fnm = '%s_%s' % (tname,filename)
        assert not OP.lexists(out_fnm), 'rename/remove %s' % out_fnm


    dvec = {}
    for tname in tors_names:
        name, angle = tname.split(':')
        tors = ndx.groups[tname]
        atoms = zip(*tors)
        # renumber indices from 1..n to 0..n-1
        # FIXME: an elegant function needed to renumber 
        #        indices in an array of lists
        for i in range(len(atoms)):
            atoms[i] = [x-1 for x in atoms[i]]
        #
        dvec[name] = [G.DistanceVector(zip(*atoms[a:a+2])) for a in range(3)]

    traj = T.Trajectory(fn_xtc)
    traj.limit(first, last)
    step1st = traj.getStep()

    conf = None
    dihan = DihedralAnalysis()

    next = True
    # init storage of results
    result = {}
    for tname in TNAMES:
        result[tname] = S.Histogram(N.arange(-180,180.,1.))

    while next:
        xyz = traj.getConfiguration()
        cell = traj.getCell()
    
        if conf is None: conf = C.Configuration(xyz, cell)
        else: conf.update(xyz, cell)

        for tname in TNAMES:
            dv = dvec[tname]
            [d(conf) for d in dv]
            angles = G.dihedral(dv[0].array, dv[1].array, 
                                dv[2].array)/M.Units.deg
            result[tname].addData(angles)

        next = traj.readFrame()

    for tname in TNAMES:
        out_fnm = '%s_%s' % (tname,filename)
        bins = result[tname].bins.tolist() + [180.]
        data = result[tname](1.0).tolist()
        data.append(data[0])
        N.savetxt(out_fnm, zip(bins,data), fmt='%10.1f %12.5e')
    traj.close()



if __name__ == '__main__':

    named_conformations = {\
        'agauche': lambda a: N.logical_or( \
            N.logical_and( a > 33.56, a < 98.60), 
            N.logical_and( a > -92.41, a < -21.79)),
        'bgauche': lambda a: N.logical_or( \
            N.logical_and( a > -86.59, a < -25.21), 
            N.logical_and( a > 40.42,  a < 98.26)),
        'atrans': lambda a: N.logical_or( a < -140.89, a > 160.37),
        'btrans': lambda a: N.logical_or( a < -142.97, a > 147.03)
        }
    named_tags = {'agauche': 'g', 'atrans': 't',
                  'bgauche': 'g', 'btrans': 't'}

    #c5c6(first=500000,last=599999,filename='gg-gt-tg_500-599.dat',
    #    named_conformations=named_conformations, named_tags=named_tags)
    #c5c6(first=400000,last=499999,filename='gg-gt-tg_400-499.dat',
    #     named_conformations=named_conformations, named_tags=named_tags)

    tors_prof(first=400000,last=1199999,filename='_400-1199.dat')
    #tors_names = ["alpha:o6'.c6'.c5'.o5'", "beta:o6'.c6'.c5'.c4'"]
    #tors_prof(tors_names=tors_names,first=400000,last=1199999,filename='400-1199.dat')

# d = N.loadtxt(filename)
