
import re, copy, itertools, os.path
import modysa.io.pdb as PDB
import modysa.io.index as NDX
import modysa.misc as M

# copy o demo.lipidA

class Molecule:

    def __init__(self, name, atom_offset=0, residue_offset=0):
        self.name = name
        self.number = 0
        self.residues = [] # FIXME: residues should be instances
        self.hb_atoms = {} # FIXME: this should be getHydrogenBondAtoms method
              # such props as hb_atoms should be shown as self.root.hb_atoms
              # self.root.name but self.system.number, 
              # self.system.atom_offset, self.system.residue_offset
        self.pdb = None    # FIXME: to be removed
        self.atom_offset = atom_offset       # FIXME: __atom_offset
        self.residue_offset = residue_offset # FIXME: __residue_offset

    def fromPDB(self, filename, max_distance=1.1):
        # atom indices start with 1 !!!
        self.pdb = PDB.PDB()
        self.pdb.read(filename) 
        self.residues = [(r.name, len(r)) for r in self.pdb.residues]
        ## find HBond donor (together with H atoms) and acceptor atoms
        # we assume that every atom whose name starts either with O or N
        # is a potential acceptor of HBond; we also include atoms whose
        # name starts with a single digit followed by O or N further 
        # followed by two arbitrary characters
        query_DA = {'atom_name': \
                        lambda n: re.match('[ON].*|\d[ON]..', n) is not None}
        donor_acceptors = [a.atom_number + self.atom_offset \
                               for a in self.pdb.filter(query_DA)]
        # no more than three H atoms expected
        query_H_after = lambda acc: \
            {'atom_number': lambda i: acc < i <= (acc+3),
             'atom_name': lambda n: re.match('[H].*|\dH..', n) is not None}
        query_H_before = lambda acc: \
            {'atom_number': lambda i: (acc-3) <= i < acc,
             'atom_name': lambda n: re.match('[H].*|\dH..', n) is not None}
        for at in donor_acceptors:
            atoms = [at]
            h = [a.atom_number + self.atom_offset \
                     for a in self.pdb.filter(query_H_after(at))]
            H_atoms = []
            acc_position = [getattr(self.pdb.atoms[at-1],which) \
                                for which in 'xyz']
            # no H atoms found *after* potential HB donor
            if not h or h[0]-at > 1: 
                h = [a.atom_number + self.atom_offset \
                         for a in self.pdb.filter(query_H_before(at))]
                atoms = h + [at]
                ### acc_name = self.pdb.atoms[at-1].atom_name
                for i in range(len(atoms)-1,0,-1):
                    ### aname = self.pdb.atoms[atoms[i-1]-1].atom_name
                    ### if atoms[i] - atoms[i-1] == 1 and acc_name[1:] in aname:
                    ###    H_atoms.append(atoms[i-1])
                    position = [getattr(self.pdb.atoms[atoms[i-1]-1],which) \
                                    for which in 'xyz']
                    distance = sum([(xi-xj)**2 for xi,xj in \
                                        zip(acc_position,position)])**.5
                    if distance < max_distance: H_atoms.append(atoms[i-1])
                    else: break
                H_atoms.sort()
            else: 
                atoms = [at] + h
                for i in range(1,len(atoms)):
                    position = [getattr(self.pdb.atoms[atoms[i]-1],which) \
                                    for which in 'xyz']
                    distance = sum([(xi-xj)**2 for xi,xj in \
                                        zip(acc_position,position)])**.5
                    if distance < max_distance: H_atoms.append(atoms[i])
                    else: break
                    ### if atoms[i] - atoms[i-1] == 1: H_atoms.append(atoms[i])
                    ### else: break
            # FIXME: take consequtive H atoms starting with the 1st element
            self.hb_atoms[at] = H_atoms
        # info about donors, acceptors and H_atoms/H_count in hb_atoms
        # {1: [2,3], 4: []}

    def fromITP(self, filename):
        # parse charge group information, identify positively/neg. charged
        # moieties
        pass
    

class Residue:

    pass

class Atom:

    __data = ('molecule_name', 'molecule_number',
              'residue_name', 'residue_number', 
              'name', 'number',
              'atom_offset', 'residue_offset')

    def __init__(self, data):
        for attr in self.__data:
            if data.has_key(attr): value = data[attr]
            else: value = None
            setattr(self, attr, value)

    def __str__(self):
        return '%s.%s:%s.%s:%s.%s' % tuple(getattr(self,attr) for attr in self.__data[:-2])

class MolecularSystem:


    def __init__(self, name, *spec):
        self.name = name
        self.molecules = []
        self.atoms = []
        atom_offset = 0 # 1
        residue_offset = 0 # 1
        molecules = {}
        for mol, count in spec: 
            nat = len(mol.pdb.atoms)
            nres = len(mol.residues)
            # same kinds of molecules might occur at various places
            # in the system molecule sequence
            # so we need to keep track of molecule_offsets, but 
            # atom_offset and residue_offset is handled consequtively
            if molecules.has_key(mol.name):
                molecule_offset = molecules[mol.name]
            else:
                molecules[mol.name] = count
                molecule_offset = 1
            for i in range(count):
                # we create a shalow copy here (i.e. copy of 
                # the toplevel object) changing the offset, i.e. 
                # atom_number for the first atom;
                # thus all data in m.pdb are the same, and in 
                # the toplevel object we might incrementally 
                # change m.atom_offset and m.residue_offset;
                # there's no need to use copy.deepcopy here
                # which would make also copies of all children object
                m = copy.copy(mol)
                m.atom_offset = atom_offset + nat*i
                m.residue_offset = residue_offset + nres*i
                # beware that each kind of molecules
                # is indexed independently, i.e.
                # the first lipid molecule, the first water mol.,
                # etc. starting with 1
                m.number = molecule_offset + i   
                self.addMolecule(m)
            if m: 
                atom_offset = m.atom_offset + nat
                residue_offset = m.residue_offset + nres
        # final info
        self.natoms = len(self.atoms)
            
        
    def addMolecule(self, molecule):
        self.molecules.append(molecule)
        data = dict.fromkeys(Atom._Atom__data)
        data.update({'molecule_name': molecule.name, 
                     'molecule_number': molecule.number,
                     'atom_offset': molecule.atom_offset,
                     'residue_offset': molecule.residue_offset})
        for atom in molecule.pdb.atoms:
            data.update(dict(name = atom.atom_name, 
                             number = atom.atom_number + molecule.atom_offset,
                             residue_name = atom.residue_name, 
                             residue_number = atom.residue_number + \
                                 molecule.residue_offset))
            a = Atom(data)
            self.atoms.append(a)

    def queryAtom(self, index, attr=None, fmt=None):
        assert 0 < index <= self.natoms, 'atoms are numbered 1..n'
        if attr is None: return self.atoms[index-1]
        elif type(attr) is type('') and attr in Atom._Atom__data: 
            return getattr(self.atoms[index-1],attr)
        elif type(attr) is type([]):
            fields = [getattr(self.atoms[index-1], a) for a in attr]
            if fmt is None: return ':'.join(map(str,fields))
            else: return fmt % tuple(fields)

    def queryResidue(self, index):
        pass
        
    def queryMolecule(self, index):
        pass

    def writePDB(self, configuration, cell=None):
        pass

    def summary(self):
        contents = {}
        for m in self.molecules:
            if contents.has_key(m.name): contents[m.name] += 1
            else: contents[m.name] = 1
        return contents

# --- auxillary functions

def hb_ndx(system, ndx_fn):
    """ given MolecularSystem object, write HB IndexFile to file"""
    assert not os.path.lexists(ndx_fn)
    donors = []
    acceptors = []
    H_atoms = []
    H_counts = []
    molecules = []
    md, ma = 0, 0
    for i in range(len(system.molecules)):
        m = system.molecules[i]
        hbat = m.hb_atoms.items()
        hbat.sort()
        offset = m.atom_offset
        da = [key+offset for key, value in hbat if len(value) > 0]
        # -NH3(+) group is not an acceptor..
        # FIMXE: apply some heuristics here
        # on the other hand, this flaw seems to be pretty harmless
        # since no HB will be detected (geom. criteria not met)
        aa = [key+offset for key, value in hbat]
        ha = M.flatten([[ih+offset for ih in value] \
                            for key, value in hbat if len(value) > 0])
        hc = [len(value) for key, value in hbat if len(value) > 0]
        # if no donors AND no acceptors: skip molecule
        if len(da) == 0 or len(aa) == 0: continue
        if len(da) > 0: 
            donors.append(da)
            H_atoms.append(ha)
            H_counts.append(hc)
            md = len(donors)
        else: md = -1
        if len(aa) > 0: 
            acceptors.append(aa)
            ma = len(acceptors)
        else: ma = -1
        molecules.append([md, ma])
    # pairs = list(itertools.combinations(range(1,73), 2))
    groups = dict(donor_atoms=donors, acceptor_atoms=acceptors, 
                  H_atoms=H_atoms, H_counts=H_counts, molecules=molecules)
    ndx = NDX.IndexFile(groups)
    ndx.write(ndx_fn)
    

