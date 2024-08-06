from modysa.analyse import interactions as I
from modysa.analyse.interactions import RadialDistributionFunction as RDF
from modysa.core import pbcgrid as P

from modysa.io import trajectory as T
from modysa.io import index as NDX
import modysa.misc as M
import os.path, itertools, operator, time
import numpy as N

def buildMolecularSystem():

    lipidA = M.model.Molecule('lipid-A')
    water = M.model.Molecule('water')
    na_ion = M.model.Molecule('na+')
    mg_ion = M.model.Molecule('mg2+')
    ca_ion = M.model.Molecule('ca2+')
    cl_ion = M.model.Molecule('cl-')

    # setup from provided pdb files
    lipidA.fromPDB('data/lipidA.pdb')
    water.fromPDB('data/water.pdb')
    na_ion.fromPDB('data/na+.pdb')
    mg_ion.fromPDB('data/mg2+.pdb')
    ca_ion.fromPDB('data/ca2+.pdb')
    cl_ion.fromPDB('data/cl-.pdb')

    # system model description
    spec = ((lipidA, 72), (na_ion, 30), (mg_ion, 30), (ca_ion, 30), 
            (cl_ion, 6), (water, 6480))
    system = M.model.MolecularSystem('lipAH-2', *spec)

    return system


def getAtoms(system,molname,selected_atoms=None,verbose=1):
    # selected_atoms = m.hb_atoms.keys()
    molecules = [m for m in system.molecules if m.name == molname]
    if verbose: print len(molecules),' %s molecules found' % molname
    if selected_atoms is None:  # we take all atoms, numbering: 1..N
        natoms = len(molecules[0].pdb.atoms) ## FIXME, should be: self.atoms
        selected_atoms = range(1,natoms+1)
    atoms = []
    for m in molecules:
        atoms.extend([ai + m.atom_offset - 1 for ai in selected_atoms])
    return atoms

def make_ndx():
    system = buildMolecularSystem()
    lipidA = system.molecules[0]
    selected_atoms = system.molecules[0].hb_atoms.keys()
    selected_atoms.sort()
    atoms = getAtoms(system,'lipid-A',selected_atoms)
    ions = []
    for ion_name in ['na+','mg2+','ca2+']:
        ions.extend(getAtoms(system,ion_name))
    return atoms, ions

def read_ndx():
    ndx = NDX.IndexFile()
    filename = os.path.join(os.path.dirname(__file__), 'data/lipidA-ions.ndx')
    ndx.read(filename)
    ndx.flat('atoms','moieties')
    ndx.range('pairs')
    # test starts
    ndx.groups['pairs'].append((73,73))
    # test ends
    data = M.flatten(ndx.groups['pairs'])
    partners = {}
    for i,j in ndx.groups['pairs']:
        if partners.has_key(i): partners[i].append(j)
        else: partners[i] = [j]
        if i == j: continue # intramoiety pairs
        if partners.has_key(j): partners[j].append(i)
        else: partners[j] = [i]
    atom_count = {}
    for i,j in partners.items():
        atom_count[i] = sum([ndx.groups['moieties'][a-1] for a in j])
    atom_count = sorted(atom_count.items(), key=operator.itemgetter(1))
    npairs = len(ndx.groups['pairs'])
    result = {}
    i, j = atom_count.pop()
    result[i] = sorted(partners[i])
    done = len(partners[i])
    stored = [i]
    # collect unique pairs (they're cloned in partners due to symmetry
    # i.e. (a,b) is stored as (a,b) and (b,a)
    while done < npairs:
        i, j = atom_count.pop()
        # is sorted needed here?
        p = sorted(set(partners[i]).difference(stored))
        done += len(p)
        stored.append(i)
        result[i] = p
    # test starts
    x = result[73]
    del(x[12:15])
    del(x[45])
    # test ends
    natoms = [0]
    natoms.extend(N.add.accumulate(ndx.groups['moieties']))
    final = {}
    for k,v in result.items(): # (ge,lt)
        final[(natoms[k-1],natoms[k])] = \
            [(natoms[i[0]-1],natoms[i[-1]]) for i in M.slicer(v)]
    # TODO: omissions? split atom ranges
    print final
    # TODO: check for recurring values, and merge keys (i.e. ranges of atoms)
    # actually, this should be hardly be an improvement since we are
    # looping over atoms; for sloppy defs of groups it might, however,
    # do the trick

def test1():
    filename = os.path.join(os.path.dirname(__file__), 'data/lipidA-ions.ndx')
    rdf = RDF(filename)
    #traj_fn = 'lipAH-2a.part0010.xtc'
    traj_fn = '/datac1/murzyn/lipAH-2/xtc/lipAH-2.part0040.xtc'
    traj_fn = '/data/T/lipAH-2/lipAH-2.part0040.xtc'
    traj = T.Trajectory(traj_fn)
    nframes = 8000
    class Cell: pass
    start_time = time.time()
    for i in range(nframes):
        configuration = traj.getConfiguration()
        cell = traj.getCell()
        rdf(configuration, cell)
        traj.readFrame()
        t = traj.getTime()
        if not int(t) % 100: print t
    total_time = time.time() - start_time
    print ' Done in %.2f sec.' % total_time
    print ' %.4f sec. per configuration' % (total_time/nframes,)
    bins, data = rdf.get()
    f = open('test.dat','w')
    for x,y in zip(bins,data/nframes):
        f.write('%10.3f %10.5f\n' % (x,y))
    f.close()

def test2():
    filename = os.path.join(os.path.dirname(__file__), 'data/lipidA-ions.ndx')
    rdf = RDF(filename)
    traj_fn = '/datac1/murzyn/lipAH-2/xtc/lipAH-2.part0040.xtc'
    traj_fn = '/data/T/lipAH-2/lipAH-2.part0040.xtc'
    traj = T.Trajectory(traj_fn)
    nframes = 10
    grid = P.gridsearch()
    for i in range(nframes):
        configuration = traj.getConfiguration()
        cell = traj.getCell()
        print traj.getTime(),
        frame = configuration[rdf.atoms]
        grid.buildgrid(cell, frame, 1.0)
        for atom_range, intervals in rdf.pairs.items():
            for ai in range(*atom_range):
                neighbors = grid.neighbors(ai)
                if neighbors is None: continue
                n, s, d = neighbors
                if min(d) < 0.2: 
                    rai = rdf.atoms[ai]
                    raj = rdf.atoms[n[N.where(d==min(d))]]
                    print rai, raj, min(d)
                    print configuration[rai]
                    print configuration[raj]
                    break
            print
        traj.readFrame()
    return n, d, rdf
    

if __name__ == '__main__': test1()
