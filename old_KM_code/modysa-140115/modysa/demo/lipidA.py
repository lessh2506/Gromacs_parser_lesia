import re, copy, itertools
import modysa.io.pdb as PDB
import modysa.io.index as NDX
import modysa.misc as M

class Moiety:

    def __init__(self, name, atom_offset=0, residue_offset=0):
        self.name = name
        self.number = 0
        self.residues = []
        self.hb_atoms = {}
        self.pdb = None
        self.atom_offset = atom_offset
        self.residue_offset = residue_offset

    def fromPDB(self, filename):
        # atom indices start with 1 !!!
        self.pdb = PDB.PDB(filename) 
        self.residues = [(r.name, len(r)) for r in self.pdb.residues]
        ## find HBond donor (together with H atoms) and acceptor atoms
        # we assume that every atom whose name starts either with O or N
        # is a potential acceptor of HBond; we also include atoms whose
        # name starts with a single digit followed by O or N further 
        # followed by two arbitrary characters
        query_DA = {'atom_name': \
                        lambda n: re.match('[ON].*|\d[ON]..', n) is not None}
        acceptors = [a.atom_number + self.atom_offset \
                         for a in self.pdb.filter(query_DA)]
        # no more than three H atoms expected
        query_H = lambda acc: {'atom_number': lambda i: acc < i <= (acc+3),
                               'atom_name': lambda n: \
                                   re.match('[H].*|\dH..', n) is not None}
        for at in acceptors:
            atoms = [at] + \
                [a.atom_number + self.atom_offset \
                     for a in self.pdb.filter(query_H(at))]
            H_atoms = []
            for i in range(1,len(atoms)):
                if atoms[i]-atoms[i-1] == 1: H_atoms.append(atoms[i])
                else: break
            # FIXME: take consequtive H atoms starting with the 1st element
            self.hb_atoms[at] = H_atoms
        # info about donors, acceptors and H_atoms/H_count in hb_atoms
        # {1: [2,3], 4: []}


    

class System:


    def __init__(self, name, *spec):
        self.name = name
        self.molecules = []
        self.atoms = []
        for mol, count in spec: 
            nat = len(mol.pdb.atoms)
            nres = len(mol.residues)
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
                m.atom_offset = nat*i
                m.residue_offset = nres*i
                m.number = i
                self.addMolecule(m)
        
    def addMolecule(self, molecule):
        self.molecules.append(molecule)
        class Atom: pass
        data = {'name': None, 'number': None, 
                'residue_name': None, 'residue_number': None, 
                'molecule_name': molecule.name, 
                'molecule_number': molecule.number}
        for atom in molecule.pdb.atoms:
            data.update(dict(name=atom.atom_name, 
                             number=atom.atom_number + molecule.atom_offset,
                             residue_name=atom.residue_name, 
                             residue_number=atom.residue_number + \
                                 molecule.residue_offset))
            a = Atom()
            [setattr(a, attr, value) for attr, value in data.items()]
            self.atoms.append(a)

        
    def writePDB(self, configuration, cell=None):
        pass



# TODO: generate NDX file on fly

def do_hb():
    lipidA = Moiety('lipid-A')
    lipidA.fromPDB('data/lipidA.pdb')
    ion = Moiety('Na+')
    ion.fromPDB('data/na+.pdb')
    water = Moiety('water')
    water.fromPDB('data/water.pdb')
    spec = ((lipidA, 72), (ion, 72), (water, 5900))
    system = System('lipAH-1', *spec)
        
    ## create hb_wb NDX file from lipidA.hb_atoms
    # XX: in IndexFile, [ pairs ] group is not required, if missing
    # do all molecules-against-all
    # XX: in IndexFile, provide support for [ excluded_pairs ] group;
    # it makes sense only if [ pairs ] is missing; in that case, 
    # do all molecules-against-all but omit those pairs which are in
    # [ excluded_pairs ] group
    donors = []
    acceptors = []
    H_atoms = []
    H_counts = []
    molecules = []
    for i in range(72):
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
        ha = M.flat([[ih+offset for ih in value] \
                         for key, value in hbat if len(value) > 0])
        hc = [len(value) for key, value in hbat if len(value) > 0]
        # if no donors AND no acceptors: skip molecule
        md, ma = -1, -1
        if len(da) == 0 or len(aa) == 0: continue
        if len(da) > 0: 
            donors.append(da)
            H_atoms.append(ha)
            H_counts.append(hc)
            md = i+1
        if len(aa) > 0: 
            acceptors.append(aa)
            ma = i+1
        molecules.append([md, ma])
    # only intermolecule HBonds
    pairs = list(itertools.combinations(range(1,73), 2))
    # intramolecule HBonds
    # pairs.extend([(i,i) for i in range(1,73)])
    groups = dict(donor_atoms=donors, acceptor_atoms=acceptors, 
                  H_atoms=H_atoms, H_counts=H_counts, 
                  molecules=molecules, pairs=pairs)
    ndx = NDX.IndexFile(groups)
    #ndx.write('data/murutka.ndx')
    return system, ndx

if __name__ == '__main__': do_hb()
    
