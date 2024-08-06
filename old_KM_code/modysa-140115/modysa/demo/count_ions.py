import os, sys, itertools
import modysa.io.index as NDX
import modysa.io.trajectory as T
import numpy as N
from modysa.core import pbcgrid as P
import modysa.misc as M

class AtomContacts:

    def __init__(self, ndxfn, groupA_ids, groupB_ids, verbose=1):

        self.verbose = verbose
        self.ndx = NDX.IndexFile()
        self.ndx.read(ndxfn)
        self.ndx.flat('atoms','moieties')
        self.ndx.range('pairs')
        self.atoms = [i-1 for i in self.ndx.groups['atoms']]

        partners = {}
        for i,j in self.ndx.groups['pairs']:
            if partners.has_key(i): partners[i].append(j)
            else: partners[i] = [j]
            if i == j: continue # intra-moiety pairs stored ONCE
            if partners.has_key(j): partners[j].append(i)
            else: partners[j] = [i]

        stored = []
        result = {}
        for gidA in groupA_ids:
            result[gidA] = sorted(set(partners[gidA]).\
                           intersection(groupB_ids).difference(stored))
            stored.append(gidA)

        natoms = [0]
        natoms.extend(N.add.accumulate(self.ndx.groups['moieties']))
        self.pairs = {}
        for k,v in result.items():
            self.pairs[(natoms[k-1],natoms[k])] = \
                [(natoms[i[0]-1],natoms[i[-1]]) for i in M.slicer(v)]        

    def count(self, trajfn, outputfn, r_cutoff=0.31, 
              first=0, last=None, step=1,
              allowed_time_diff=0.05, register=False):
        # FIXME: step
        traj = T.Trajectory(trajfn, allowed_time_diff=allowed_time_diff) # 50 fs
        traj.limit(first,last)

        grid = P.gridsearch()

        if outputfn is None: f = sys.stdout
        else: f = open(outputfn,'w')

        if last is None: last = traj.parts[-1][-1]
        else: last = min(last,traj.parts[-1][-1])
        now = traj.getTime()
        if self.verbose:
            print 'Reading configurations from %d to %d' % (first,last)
        while now <= last:
            configuration = traj.getConfiguration()
            cell = traj.getCell()
            frame = configuration[self.atoms]
            grid.buildgrid(cell, frame, r_cutoff)
            bounded = []
            atom_pairs = {}
            
            for atom_range, intervals in self.pairs.items():
                for ai in range(*atom_range):
                    neighbors = grid.neighbors(ai)
                    if neighbors is None: continue
                    n, s, d =  neighbors
                    if register:
                        atom_pairs[self.atoms[ai]] = N.take(self.atoms,n)
                    bounded.append(len(n))
                    
            occurence = [(i[0], len(list(i[1]))) \
                             for i in itertools.groupby(sorted(bounded))]
            if not register:
                f.write('%12.1f %12d %s\n' % (now,len(bounded),
                              ','.join(['%d:%d' % (i,j) for i,j in occurence])))
            else:
                label = '_'.join(['%d:%s' % \
                                      (i+1,','.join(map(str,
                                                        [k+1 for k in j]))) \
                                      for i,j in atom_pairs.items()])
                f.write('%12.1f %12d %s\n' % (now,len(bounded),label))
            f.flush()
            if self.verbose and now%100 == 0:
                print now, len(bounded)
            traj.readFrame()
            now = traj.getTime()
        traj.close()
        f.close()

                    


def direct_contacts(groupA_ids, groupB_ids,
         r_cutoff=0.31,first=0,last=None,outputfn=None,verbose=1):

    # --- parse IndexFile first ---
    ndxfn = os.path.join(os.path.dirname(__file__), 'data/lipidA-ions.ndx')
    ndx = NDX.IndexFile()
    ndx.read(ndxfn)
    ndx.flat('atoms','moieties')
    ndx.range('pairs')
    atoms = [i-1 for i in ndx.groups['atoms']]
    partners = {}
    for i,j in ndx.groups['pairs']:
        if partners.has_key(i): partners[i].append(j)
        else: partners[i] = [j]
        if i == j: continue # intra-moiety pairs stored ONCE
        if partners.has_key(j): partners[j].append(i)
        else: partners[j] = [i]
    stored = []
    result = {}
    for gidA in groupA_ids:
        result[gidA] = sorted(set(partners[gidA]).\
                                  intersection(groupB_ids).difference(stored))
        stored.append(gidA)
    natoms = [0]
    natoms.extend(N.add.accumulate(ndx.groups['moieties']))
    pairs = {}
    for k,v in result.items():
        pairs[(natoms[k-1],natoms[k])] = \
            [(natoms[i[0]-1],natoms[i[-1]]) for i in M.slicer(v)]
    # --- parsing IndexFile done ---


    traj_fn = '/datac1/murzyn/lipAH-2a/xtc/lipAH-2a.part00*.xtc'
    #traj_fn = '/data/T/lipAH-2/lipAH-2.part004[01].xtc'
    traj = T.Trajectory(traj_fn,allowed_time_diff=0.05) # 50 fs
    # lipAH-1a: (27, 462147.031250, 462146.000000)
    traj.limit(first,last)

    grid = P.gridsearch()

    if outputfn is None: f = sys.stdout
    else: f = open(outputfn,'w')

    if last is None: last = traj.parts[-1][-1]
    else: last = max(last,traj.parts[-1][-1])
    now = traj.getTime()
    print 'Reading configurations from %d to %d' % (first,last)
    while now <= last:
        configuration = traj.getConfiguration()
        cell = traj.getCell()
        frame = configuration[atoms]
        grid.buildgrid(cell, frame, r_cutoff)
        bounded = []
        for atom_range, intervals in pairs.items():
            for ai in range(*atom_range):
                neighbors = grid.neighbors(ai)
                if neighbors is None: continue
                n, s, d =  neighbors
                bounded.append(len(n))
        occurence = [(i[0], len(list(i[1]))) \
                         for i in itertools.groupby(sorted(bounded))]
        f.write('%12.1f %12d %s\n' % (now,len(bounded),
                          ','.join(['%d:%d' % (i,j) for i,j in occurence])))
        f.flush()
        if verbose and now%100 == 0:
            print now, len(bounded)
        traj.readFrame()
        now = traj.getTime()
    traj.close()
    f.close()




if __name__ == '__main__':
    ions = [73,74,75]
    lipidAs = range(1,73)
    #direct_contacts([ions[2]], lipidAs, outputfn='data/count_2a-ca2+_1.dat')
    ndxfn = os.path.join(os.path.dirname(__file__), 'data/lipidA-ions.ndx')
    ion_which = ions[1]
    ac = AtomContacts(ndxfn, [ion_which], lipidAs)
    trajfn = '/data/T/lipAH-2/lipAH-2.part0041.xtc'
    trajfn = '/datac1/murzyn/lipAH-2a/xtc/lipAH-2a.part0043.xtc'
    outputfn='data/test2-mg2+.dat'
    #
    # natoms = [0]
    # natoms.extend(N.add.accumulate(ac.ndx.groups['moieties']))
    # MgIons = range(natoms[ion_which-1],natoms[ion_which])
    # MgIons = (N.array(ac.atoms)[MgIons]+1).tolist()
    #
    ac.count(trajfn, outputfn, register=True)

