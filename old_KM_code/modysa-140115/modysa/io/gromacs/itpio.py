# written by KM
# formerly misc_utils used together with pymol_utils

import re, operator, os.path as OP


class MoleculeType:

    __patre = {'comments': re.compile('^\s*;'),
               'sections': re.compile('^\[ +[a-z]+ +\].*$',re.M),
               'hashes': re.compile('^\s*#'),
               'newline': re.compile('[\n\r]+')}

    def __init__(self, ffname='ffoplsaa'):
        self.ff = ffname

    def read(self, filename):
        # ITP
        f = open(filename, 'r')
        # strip out comment lines
        content = ''.join([line for line in f.readlines() \
                  if re.match(self.__patre['comments'], line) is None])
        f.close()
        # get sections
        section_names = re.findall(self.__patre['sections'],content)
        section_data = [x for x in re.split(self.__patre['sections'],content) \
                        if x.strip()]
        # 
        self.data = {}
        for name, data in zip(section_names, section_data):
            kwd = name.split()[1]
            if self.data.has_key(kwd): self.data[kwd] += '\n'+data.strip()
            else: self.data[kwd] = data.strip()
        #
        def fetchEntries(data, fields):
            records = []
            ftypes = map(operator.itemgetter(1), fields)
            for line in re.split(self.__patre['newline'], data):
                if not line.strip(): continue
                # if line contains fewer entries, ftypes is truncated!
                # this applies in dihedrals section
                records.append([t(v) for t,v in zip(ftypes, line.split())])
            return records
            
        # parse atoms
        fields = (('atom_number', int), ('atom_type', str),
                 ('residue_number', int), ('residue_name', str),
                 ('atom_name', str), ('charge_group', int),
                 ('charge', float), ('mass', float))
        what = 'atoms'
        if self.data.has_key(what):
            self.data[what] = fetchEntries(self.data[what], fields)
        # parse bonds
        fields = (('iatom', int), ('jatom', int), ('ftype', int))
        what = 'bonds'
        if self.data.has_key(what):
            self.data[what] = fetchEntries(self.data[what], fields)
        # angles
        fields = (('iatom', int), ('jatom', int), ('katom', int), 
                  ('ftype', int))
        what = 'angles'
        if self.data.has_key(what):
            self.data[what] = fetchEntries(self.data[what], fields)
        # dihedrals
        fields = (('iatom', int), ('jatom', int), ('katom', int), 
                  ('latom', int), ('ftype', int), ('name', str))
        what = 'dihedrals'
        if self.data.has_key(what):
            self.data[what] = fetchEntries(self.data[what], fields)
        # parse pairs
        fields = (('iatom', int), ('jatom', int), ('ftype', int))
        what = 'pairs'
        if self.data.has_key(what):
            self.data[what] = fetchEntries(self.data[what], fields)

    def write(self, filename):
        assert not OP.lexists(filename)
        f = open(filename, 'w')
        #
        def dump(name, data):
            f.write('[ %s ]\n' % name)
            out = ['\t'.join(map(str, line)) for line in data]
            f.write('\n'.join(out))
            f.write('\n'*4)
        #
        dump('moleculetype', [[self.data['moleculetype']]])
        for item in ['atoms', 'bonds', 'angles', 'dihedrals', 'pairs']:
            dump(item, self.data[item])
        f.close()


def getGromacsForcefieldParameters(ff_filename,verbose=0):

    comment = re.compile('^\s*;')
    sections = re.compile('^\[ [a-z]+ \].*$',re.M)
    hashes = re.compile('^\s*#')
    
    # open file and strip out comment lines
    f = open(ff_filename,'r')
    ff_buf = ''.join([line for line in f.readlines() \
                      if re.match(comment,line) is None])
    # FIXME: hashes
    #   handling blocks: #ifdef..#else..#endif
    #   handling lines:  #define
    f.close()

    # sections
    section_names = re.findall(sections,ff_buf)
    section_data = [x for x in re.split(sections,ff_buf) if x.strip()]
    ff_database = {}
    for name, data in zip(section_names,section_data):
        if ff_database.has_key(name):
            ff_database[name] += '\n'+data.strip()
        else: ff_database[name] = data.strip()

    forcefield = {}

    what = 'defaults'
    if ff_database.has_key('[ %s ]' % what):
        fields = ff_database['[ %s ]' % what].split()
        forcefield[what] = {'nbfunc': int(fields[0]),
                            'comb-rule': int(fields[1]),
                            'gen-pairs': fields[2],
                            'fudgeLJ': float(fields[-2]),
                            'fudgeQQ': float(fields[-1])}
    
    def fetchEntries(what):
        entries = [line.strip() for line in \
                   ff_database['[ %s ]' % what].split('\n') \
                   if line.strip()]
        return entries
    
    data = {}
    atomtypes_line = None
    what = 'atomtypes'
    if ff_database.has_key('[ %s ]' % what):
        for entry in fetchEntries(what):
            if re.match(hashes,entry): continue
            if atomtypes_line is None:
                gj = entry.split()
                try:
                    fifth = int(gj[1])
                    # amber format found
                    atomtypes_line = (('name', str), ('atom.number', int), 
                                      ('mass', float), ('charge', float), 
                                      ('ptype', str), ('sigma', float), 
                                      ('epsilon', float))
                except ValueError:
                    # opls format found
                    atomtypes_line = (('nbtag', str), ('name', str), 
                                      ('atom.number', int), ('mass', float), 
                                      ('charge', float), ('ptype', str),
                                      ('sigma', float), ('epsilon', float))
                finally:
                    nfields = len(atomtypes_line)
            parsed = dict([(kind[0], kind[1](item)) for kind, item in \
                               zip(atomtypes_line, entry.split()[:nfields])])
            name = parsed['name']
            del(parsed['name'])
            if data.has_key(name) and verbose:
                print ' Warning: entry %s duplicated' % name
            data[name] = parsed
        forcefield[what] = data

    data = {}
    what = 'bondtypes'
    if ff_database.has_key('[ %s ]' % what):
        for entry in fetchEntries(what):
            if re.match(hashes,entry): continue
            iat, jat, ftype, b0, kb = entry.split()[:5]
            bid = [iat,jat]
            bid.sort()
            bid = tuple(bid)
            if data.has_key(bid) and verbose:
                print ' Warning: entry %s duplicated' % bid[0]
            data[bid] = {'ftype': int(ftype),
                         'b0': float(b0),
                         'kb': float(kb)}
        forcefield[what] = data

    data = {}
    what = 'angletypes'
    if ff_database.has_key('[ %s ]' % what):
        for entry in fetchEntries(what):
            if re.match(hashes,entry): continue
            items = entry.split()
            iat, jat, kat = items[:3]
            ftype = int(items[3])
            aid = [[iat,jat,kat],[kat,jat,iat]]
            aid.sort()
            aid = tuple(aid[0])
            if data.has_key(aid) and verbose:
                print ' Warning: entry %s duplicated' % aid[0]
            if ftype == 1:
                th0, cth = items[4:6]
                data[aid] = {'ftype': int(ftype),
                             'th0': float(th0),
                             'cth': float(cth)}
            elif ftype == 5: # charmm/stockholm
                th0, cth, ub0, cub = items[4:8]
                data[aid] = {'ftype': int(ftype),
                             'th0': float(th0),
                             'cth': float(cth),
                             'ub0': float(ub0),
                             'cub': float(cub)}
        forcefield[what] = data

    data = {}
    what = 'dihedraltypes'
    if ff_database.has_key('[ %s ]' % what):
        for entry in fetchEntries(what):
            if re.match(hashes,entry): continue
            iat, jat, kat, lat, ftype = entry.split()[:5]
            params = {}
            try: fid = int(ftype)
            except ValueError:
                print ' Warning: skiping "%s"' % entry
                continue # two central atoms (i.e. jat,kat) vs. four atoms 
            if fid == 3: # Ryckaert-Bellemans C0..C5
                params = dict([(name,value) for name, value \
                                   in zip(['C%d' % i for i in range(6)],
                                          map(float,entry.split()[5:11]))])
            elif fid in [1,4,9]: # amber/stockholm/charmm or impropers?
                params = dict([(name,fid is 9 and [value] or value) \
                                   for name, value \
                                   in zip(['phase','kd','pn'],
                                          map(float,entry.split()[5:8]))])
            elif fid == 2: # harmonic improper dihedral
                params = dict([(name,value) for name, value \
                                   in zip(['q0','cq'],
                                          map(float,entry.split()[5:7]))])
            assert params, ' empty params for dihedral ftype %d' % fid
            did = [[iat,jat,kat,lat],[lat,kat,jat,iat]]
            did.sort()
            did = tuple(did[0] + [fid]) 
            # BEWARE: dihedrals are given by atomtype quartet AND int(ftype)
            # this is necessery since certain quartets might be 
            # both proper (e.g. 1,3,5,8,9) and improper dihedrals (2,4)
            if data.has_key(did) and fid != 9 and verbose:
                # multiple entries possible *only* for ftype: 9
                print ' Warning: entry %s duplicated' % did
            if not data.has_key(did): data[did] = params
            elif fid == 9: # condition checked just in case
                for k in params.keys(): data[did][k].extend(params[k])
            else:
                # not storing in data dict !!!
                print ' Warning: entry %s duplicated' % did
                print '     >>>: %s', entry
        forcefield[what] = data

    data = {}
    what = 'pairtypes'
    if ff_database.has_key('[ %s ]' % what):
        for entry in fetchEntries(what):
            if re.match(hashes, entry): continue
            iat, jat, ftype, sigma14, epsilon14 = entry.split()[:5]
            pid = [[iat, jat], [jat, iat]]
            pid.sort()
            pid = tuple(pid[0])
            if data.has_key(pid) and verbose:
                print ' Warning: entry %s duplicated' % pid[0]
            data[pid] = {'ftype': int(ftype),
                         'sigma14': float(sigma14),
                         'epsilon14': float(epsilon14)}
        forcefield[what] = data

    return forcefield

def getParams(forcefield,typename,intername):
    if not forcefield.has_key(typename): return
    atoms = list(intername)
    intername.reverse()
    ratoms = tuple(atoms)
    atoms = tuple(intername)

    def printout(data):
        print '%6s'*len(atoms) % tuple(atoms),
        print '%5d' % data['func'],
        print data['params']

    if forcefield[typename].has_key(tuple(atoms)):
        printout(forcefield[typename][tuple(atoms)])
    elif forcefield[typename].has_key(tuple(ratoms)):
        printout(forcefield[typename][tuple(ratoms)])
