
import string, os.path, re, math, string
import modysa.misc as M

# transfered from do_PmB6LipA66

class Atom:

    __fields = {'atom_number': int, 'atom_name': string.strip,
                'residue_name': string.strip, 'chain': string.strip,
                'residue_number': int, 
                'x': float, 'y': float, 'z': float}

    def __init__(self, **atom_data):
        for field in atom_data:
            if self.__fields.has_key(field): 
                setattr(self, field, self.__fields[field](atom_data[field]))

    def __str__(self):
        line = 'Atom %4s, %5d; Residue %4s, %5d; at position: %8.3f %8.3f %8.3f'
        return line % (self.atom_name, self.atom_number, self.residue_name,
                       self.residue_number, self.x, self.y, self.z)


class Residue:

    def __init__(self, name, number=1):
        self.name = name
        self.number = number
        self.atoms = []

    def addAtom(self, atom):
        self.atoms.append(atom)

    def __str__(self):
        line = 'Residue %4s, %5d; %d atoms'
        return line % (self.name, self.number, len(self))

    def __len__(self):
        return len(self.atoms)


class Molecule:

    def __init__(self, name):
        # represents chains as well
        self.name = name
        self.residues = []

    def addResidue(self, residue):
        self.residues.append(residue)

    def __len__(self):
        return len(self.residues)

class PDB:

    __fields = {'label': (0, 6), 'atom_number': (6, 11), 'atom_name': (12, 16),
              'residue_name': (17, 20), 'chain': (21, 22),
              'residue_number': (22, 26), 
              'x': (30, 38), 'y': (38, 46), 'z': (46, 54)}

    __records = {'ATOM': 'ATOM  %5s %-4s %-3s %c%4s    %8.3f%8.3f%8.3f',
                 'CRYST1': 'CRYST1%(a)9.3f%(b)9.3f%(c)9.3f' + \
                           '%(bc)7.2f%(ac)7.2f%(ab)7.2f P 1',
                 'HETATM': '',
                 'CONECT': ''}

    def __init__(self):
        self.atoms = []
        self.residues = []
        self.cell = []

    def read(self, filename, renumber=True):

        contents = open(filename,'r').read()
        cell = re.search('(?m)^CRYST1' \
                         ' *(?P<a>\d{1,5}\.\d{3}) *(?P<b>\d{1,5}\.\d{3})' \
                         ' *(?P<c>\d{1,5}\.\d{3}) *(?P<bc>\d{1,4}\.\d{2})' \
                         ' *(?P<ac>\d{1,4}\.\d{2}) *(?P<ab>\d{1,4}\.\d{2})',
                         ''.join(contents))
        if cell:
            self.cell = dict([(k, float(v)) \
                                  for k,v in cell.groupdict().items()])
        else: self.cell = None

        # FIXME: HETATM, CONECT
        data = re.findall('(?m)^ATOM  .+?$', contents)

        if renumber:
            atom_offset = 1 # renumber atoms starting with 1
            residue_offset = 1
        else: 
            atom_offset = None
            residue_offset = None

        def parsePdbLine(line):
            f = self.__fields
            fields = [(k, line[f[k][0]:f[k][1]]) for k,v in f.items()]
            return Atom(**dict(fields))

        self.atoms = [parsePdbLine(line) for line in data]
        self.residues = []
        residue = None
        # since we are not making deep copies of Atom instances
        # both self.atoms and self.residue.atoms contain the same Atom instances
        orig_residue_number, ires = None, 0
        for iat in range(len(self.atoms)): 
            at = self.atoms[iat]
            # renumbering atoms, starting with 1
            setattr(at, 'atom_number', 
                    atom_offset and atom_offset + iat or at.atom_number)
            # renumbering residues, starting with 1
            if residue is None:
                ires = residue_offset or at.residue_number
                orig_residue_number = at.residue_number
                residue = Residue(at.residue_name, ires)
            elif orig_residue_number != at.residue_number:
                self.residues.append(residue)
                ires += 1
                # FIXME: stop renumbering residues
                #ires = residue_offset and ires or at.residue_number
                orig_residue_number = at.residue_number
                residue = Residue(at.residue_name, ires)
            setattr(at, 'residue_number', ires)
            residue.addAtom(at)
        if residue: self.residues.append(residue)

    def filter(self, query):
        # query is a dictionary, e.g. {'x': lambda i: 12. <= i <= 15, 
        #       'atom_name': lambda i: re.match('C.{0,2}', i) is not None}
        result = []
        for atom in self.atoms:
            filtered = [v(getattr(atom,k)) for k,v in query.items() \
                            if hasattr(atom,k)]
            if sum(filtered) != len(filtered): # not all True
                continue
            result.append(atom)
        return result

    def renumberAtoms(self, offset):
        # do we really need this method? 
        # deleting residues/atoms should create new PDB instance
        for i in range(len(self.atoms)): 
            a = self.atoms[i]
            a.atom_number = i + offset

    def getMolecules(self):
        # return list of Molecule instances
        # join residues on the basis of known bonds between standard
        # residues or by chain_id 
        pass

    def getConfiguration(self, selection=None):
        """ return configuration [Ang] for given selection of atoms,
        if selection is None, configuration of all atoms is returned """
        if selection is None: selection = self.atoms
        return [(a.x, a.y, a.z) for a in selection]

    def setConfiguration(self, configuration, selection=None):
        if selection is None: selection = self.atoms
        assert configuration and M.isiterable(configuration) and \
            len(configuration) == len(selection) and \
            M.isiterable(configuration[0]) and len(configuration[0]) == 3, \
            'ERROR: a list of x,y,z coordinates expected'
        for atom, (x,y,z) in zip(selection, configuration):
            atom.x, atom.y, atom.z = map(float, (x, y, z))

    def getCell(self):
        # return tensor
        pass

    def write(self,filename,offsets={}):
        assert not os.path.lexists(filename)
        if offsets.has_key('atom_number'): an_off = offsets['atom_number']
        else: an_off = 0
        line = self.__records['ATOM'] # FIXME: check attr: hetatm = False
        check_rn = lambda chain, log_rn: \
           log_rn > 4 or string.ascii_uppercase[int(log_rn) - 4] and \
           chain or ' '
        contents = [line % (str(a.atom_number+an_off)[-5:], 
                            len(a.atom_name) == 4 and \
                                a.atom_name or ' '+a.atom_name,
                            a.residue_name, 
                            check_rn(a.chain, math.log10(a.residue_number)),
                            str(a.residue_number)[-4:], 
                            a.x, a.y, a.z) for a in self.atoms]
        fh = open(filename, 'w')
        if self.cell: # write CRYST1 too
            fh.write(self.__records['CRYST1'] % self.cell)
            fh.write('\n')
        fh.write('\n'.join(contents))
        fh.write('\n')
        fh.close()
                 
        
        
