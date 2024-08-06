import os.path, re, string, sys, copy
from ruczaj import Misc

PAT_COMMENT = re.compile('^\s*;')
PAT_SECTION = re.compile('^\[ [a-z]+ \].*$',re.M)
PAT_IFDEF = re.compile('^\#ifdef.+?^\#endif',re.M+re.S)
PAT_IFDEF_ELSE_ENDIF = re.compile('\#ifdef\s+(?P<Flag>.+?)$'
                                  '(?P<True>.+)'
                                  '\#else(?P<False>.+)\#endif',
                                  re.M+re.S)
PAT_IFDEF_ENDIF = re.compile('\#ifdef\s+(?P<Flag>.+?)$(?P<True>.+)',re.M+re.S)


class Forcefield:

    def __init__(self,path,defs=None,verbose=0):
        # defs might be a list of flags recognised by #ifdef clauses
        # e.g. HEAVY_H, etc.        
        self.path = path
        self.verbose = verbose
        self.defs = defs is not None and defs or []

    def getEntries(self,data):
        return [line.strip() for line in data.split('\n') if line.strip()]

    def getSections(self,data,section=None):
        # if a file contains more than one section with a given name
        # all parts will be saved together
        if section is None: section = PAT_SECTION
        section_names = re.findall(section,data)
        section_data = [x for x in re.split(section,data) if x.strip()]
        ff_database = {}
        for name, data in zip(section_names,section_data):
            tag = re.match('\[\s+(?P<tag>.+)\s+\]',name).group('tag')
            # get rid of #ifdef .. #endif
            # - they will be handled later
            data = reduce(lambda x,y: x+y, re.split(PAT_IFDEF,data))
            if ff_database.has_key(tag):
                ff_database[tag] += '\n'+data.strip()
            else: ff_database[tag] = data.strip()
        return ff_database

    def printDihedralWildcards(self):
        keys = [x for x in self.dihedraltypes.keys() \
                if 'X' in x and type(x) == type(())]
        keys.sort()
        print '\n'.join(['%5s'*len(k) % k for k in keys])

    def read(self):
        # to be overwritten
        pass

class OPLSAA(Forcefield):

    PROPS = {'atom_mass': 2,   'bond_type': 0}
    ffsettings = {'bondfunc': 1, 'anglefunc': 1,
                  'dihedralfunc': 3, 'improperfunc': 1,
                  'pairfunc': 1}

    def getNonbonded(self, filename=None):
        
        # ----------------------
        # nonbonded interactions
        # ----------------------
        
        if filename is None: 
            filename = os.path.join(self.path,'ffnonbonded.itp')
        fh = open(filename,'r')
        data = ''.join([line for line in fh \
                        if re.match(PAT_COMMENT,line) is None])
        sections = self.getSections(data)
        if 'atomtypes' not in sections.keys():
            raise ValueError(' atomtypes not found in ffnonbonded.itp')

        # ---
        def parseEntry(entry):
            opls_name, atype, anumb, amass, acharge, ptype, sigma, epsilon = \
                       entry.split()[:8]
            is_opls = re.match('opls_(?P<tag>.+)',opls_name)
            if is_opls: tag = is_opls.group('tag')
            elif len(opls_name) <= 4: tag = opls_name
            else: tag = None
            return tag, (atype,int(anumb),float(amass),float(acharge),
                         ptype,float(sigma),float(epsilon))
        # ---
        
        if self.verbose: print 'ATOMTYPES'
        params = {}; itype = 'atomtypes'
        for entry in self.getEntries(sections['atomtypes']):
            tag, fields = parseEntry(entry)
            if self.verbose: print tag, fields
            if tag is None and self.verbose:
                print ' Warning! Entry: %s ... skipped' % entry[:40]
                continue
            params[tag] = fields
        if hasattr(self, itype): getattr(self, itype).update(params)
        else: setattr(self, itype, params)

        # check and add #ifdef ... #endif if present depending on
        # user provided condition
        for entry in re.findall(PAT_IFDEF,data):
            fields = re.match(PAT_IFDEF_ELSE_ENDIF,entry)
            if fields is not None:
                if fields.group('Flag') in self.defs: gj = 'True'
                else: gj = 'False'
            else:
                fields = re.match(PAT_IFDEF_ENDIF,entry)
                if fields is None: continue
                gj = 'True'
            for entry in self.getEntries(fields.group(gj)):
                tag, fields = parseEntry(entry)
                self.atomtypes[tag] = fields

    def getBonded(self, filename=None):

        # -------------------
        # bonded interactions
        # -------------------

        if filename is None: filename = os.path.join(self.path,'ffbonded.itp')
        fh = open(filename,'r')
        data = ''.join([line for line in fh \
                        if re.match(PAT_COMMENT,line) is None])
        sections = self.getSections(data)

        # ---
        def aUniqueKey(key, entry):
            if type(key) == type([]):
                if key[0] > key[-1]: key.reverse()
                key = tuple(key)
            if params.has_key(key) and self.verbose:
                print ' Warning! Entry for %s' % ' - '.join(key),
                print ' already seen'
                print '          The present definition: ',
                print ' '.join(map(str,params[key]))[:40]
                print '          is replaced with: %s ...' % \
                      " ".join(entry.split())[:50]
            return key
        # ---

        # bondtypes
        params = {}; itype = 'bondtypes'
        if sections.has_key(itype):
            for entry in self.getEntries(sections[itype]):
                atype_i, atype_j, func, b0, kb = entry.split()[:5]
                key = aUniqueKey([atype_i, atype_j], entry)
                params[key] = (int(func), float(b0), float(kb))
            if hasattr(self, itype): getattr(self,itype).update(params)
            else: setattr(self,itype,params)

        # angletypes
        params = {}; itype = 'angletypes'
        if sections.has_key(itype):
            for entry in self.getEntries(sections[itype]):
                atype_i, atype_j, atype_k, func, th0, cth0 = entry.split()[:6]
                key = aUniqueKey([atype_i, atype_j, atype_k], entry)
                params[key] = (int(func), float(th0), float(cth0))
            if hasattr(self, itype): getattr(self,itype).update(params)
            else: setattr(self,itype,params)

        # dihedraltypes
        params = {}; itype = 'dihedraltypes'
        if sections.has_key(itype):
            for entry in self.getEntries(sections[itype]):
                fields = entry.split(';')[0].split()[:11]
                if entry.startswith('#define'):
                    key = aUniqueKey(fields[1], entry)
                    params[key] = tuple(fields[2:])
                elif len(fields) == 11 and fields[4] == '3':
                    atype_i, atype_j, atype_k, atype_l = fields[:4]
                    key = aUniqueKey([atype_i, atype_j, atype_k, atype_l], entry)
                    params[key] = tuple([int(fields[4])] + map(float,fields[5:11]))
                elif len(fields) == 6 and fields[2] == '1':
                    # FIXME: check if several lines (fourier terms) possible
                    atype_i, atype_j = fields[:2]
                    key = aUniqueKey([atype_i,atype_j], entry)
                    params[key] = tuple([int(fields[2])] + map(float,fields[3:6]))
            if hasattr(self, itype): getattr(self,itype).update(params)
            else: setattr(self,itype,params)


    def read(self):
        self.getNonbonded()
        self.getBonded()



class PDBFile:

    PAT_HETATM = re.compile('^HETATM *?(?P<atnum>\d{1,5}) +?'
                            '(?P<atname>\w{1,4}).+?(?P<resname>\w{1,3}) +?'
                            '(?P<chain>[ A-Z]) *?(?P<resnum>\d{1,4}) +?'
                            '(?P<x>[\d\.-]{5,8}) *?(?P<y>[\d\.-]{5,8}) *?'
                            '(?P<z>[\d\.-]{5,8}) *?(?P<occ>[\d\.-]{4,6}) *?'
                            '(?P<fac>[\d\.-]{4,6}) +?(?P<segid>\w{1,4}) *?'
                            '(?P<el>\w{1,2})$',re.M)
    PAT_CONECT = re.compile('^CONECT *?(?P<atnum>\d+)(?P<data>[ \d]+)$',re.M)
    PAT_REMARK = re.compile("^REMARK (?P<tag>\w+) +(?P<data>.+)$",re.M)

    def __init__(self, filename):
        data = open(filename,'r').read()
        self.remark = re.findall(self.PAT_REMARK,data)
        self.hetatm = re.findall(self.PAT_HETATM,data)
        self.conect = re.findall(self.PAT_CONECT,data)
        self.dih = {}
        self.imp = {}
        self.parseRemarks()

    def parseRemarks(self):
        for tag, data in self.remark:
            if tag == 'NAME':
                self.molname = data.split()[0]
                self.nrexcl = len(self.hetatm) < 3 and 1 or 3
            elif tag == 'DESC':
                pass # FIXME
            elif tag in ['DIH','IMP']:
                db = getattr(self,tag.lower())
                res_id, ai, aj, ak, al, defname = data.split()[:6]
                # the allowed convention reads RESIDUE_NAME.RESIDUE_NUMBER
                if tag.lower() == 'dih': 
                    resname, resnumber = res_id.split('.')[:2]
                    atoms = []
                    for a in [ai, aj, ak, al]:
                        if '.' not in a:
                            # atoms from the named XXXX (DIH XXXX) residue
                            name = '%s.%s' % (a, resnumber)
                        else: name = a
                        atoms.append(name)
                    atoms = tuple(atoms)
                else: 
                    # improper definitions are generic
                    # but a hack is also needed here in cases improper definition
                    # span several residues 
                    atoms = tuple([ai,aj,ak,al])
                if db.has_key(res_id):
                    if db[res_id].has_key(atoms):
                        oldkey = db[res_id][atoms]
                        print ' Warning! Definition %s in %s for "%s" ' % \
                              (oldkey,resname,' '.join(atoms))
                        print '          replaced with %s' % defname
                    else: db[res_id][atoms] = defname
                else: db[res_id] = {atoms: defname}


class Pdb2Top:

    def __init__(self, pdbfilename, forcefield):
        """ forcefield is an instance """

        assert isinstance(forcefield, Forcefield)
        self.forcefield = forcefield
        self.pdb = PDBFile(pdbfilename)
        # if forcefield is 'oplsaa':
        #     self.forcefield = OPLS_aa('gromacs_3.3.1')
        #     self.forcefield.read()
        #     self.ffsettings = {'bondfunc': 1, 'anglefunc': 1,
        #                        'dihedralfunc': 3, 'improperfunc': 1,
        #                        'pairfunc': 1}
        # else:
        #     raise ValueError(' Forcefield "%s" unsupported' % forcefield)


    def addForceFieldParameters(self,itpfilename):
        # FIXME
        # add missing parameters to self.forcefield
        # redefince API: this should be a method of ForceField object
        # Pdb2Top should have method setForceField
        pass

    def __call__(self):
        self.output = []
        self.atoms = {}
        self.missing = {}
        self.molecule()
        self.atomList()
        self.bondList()
        self.angleList()
        self.dihedralList()
        self.pairList()
        self.improperList()
        return '\n'.join(self.output)

    def isDefinedDihedral(self,ai,aj,ak,al):
        """ checks if atoms given by numbers correspond to
        known definitions of dihedral. It's assumed that DIH definitions have
        no +/- signs. Signs +/i occur only in IMP definitions. This means
        that in DIH definition ALL FOUR ATOMS ARE IN THE SAME RESIDUE """
        
        A = self.atoms
        D = self.pdb.dih
        ##resnames = [A[x]['residue_name'] for x in [ai,aj,ak,al]]
        #
        # molecule might have several residues with the same residue_name
        # but they might be not exactly the same (i.e. KDO in lipid A)
        # due to links they make with other residues
        # 
        # the problem is that if a definition occurs with all residues
        # with a given name, DIH entries have to be in appropriate number
        # of copies
        # 
        # note that improper torsions are defined by name (number is neglected)
        resnames = ['%s.%d' % (A[x]['residue_name'],
                               A[x]['residue_number']) \
                        for x in [ai,aj,ak,al]]
        # if resnames.count(resnames[0]) != 4:
        #     # atoms from more than one residue
        #     # previously not allowed, now we provide support for atom_name.residue_number
        #     return
        #if resnames.count(resnames[0]) != 4:
        #    return # atoms from different residues

        # a sole condition to met is that the residue_name.residue_number need to be
        # read from PDB header (REMARK DIH field), matching four atom names with
        # PDB definitions is done just below
        #if resnames[0] not in D.keys():
        #    return # no dihedral definition known

        atnames = ['%s.%s' % (A[x]['atom_name'],A[x]['residue_number']) \
                       for x in [ai, aj, ak, al]]
        key = None
        for res_id in set(resnames):
            if not D.has_key(res_id): continue
            allowed_keys = D[res_id].keys()
            if tuple(atnames) in allowed_keys:
                key = tuple(atnames)
            reversed_key = tuple(reversed(atnames))
            if reversed_key in allowed_keys:
                key = reversed_key
            if key is not None: break
        if key is None: return  # no dihedral definition known
        else: return D[res_id][key]

    def molecule(self):
        output = ['[ moleculetype ]\n; name     nrexcl']
        output.append('%s  %d\n\n' % (self.pdb.molname, self.pdb.nrexcl))
        self.output.append('\n'.join(output))
    
    def atomList(self, charge_factor=0.1, chrgrp_factor=100):
        fields = ('nr','type','resnr','resn','atom','cgnr','charge','mass')
        output = ['[ atoms ]']
        output.append('; %5s  %-10s %4s %4s   %4s  %4s  %10s  %12s' %  fields)
        ATOMLINE = ' %6d  %-10s  %4d  %-4s  %-4s  %4d  %10.3f  %12.4f'
        FF = self.forcefield
        iamass = FF.PROPS['atom_mass']
        ibtype = FF.PROPS['bond_type']
        for atnum, atname, resname, chain, resnum, \
                x, y, z, occ, fac, segid, el in self.pdb.hetatm:
            if segid[0] in string.digits: # assuming opls_[0-9]{3}[A-Z]* 
                attype = 'opls_%s' % segid
            else: attype = segid # everything else, but 4 chars at most
            if not FF.atomtypes.has_key(segid):
                if self.missing.has_key(segid):
                    self.missing[(segid,)].append(atnum)
                else: self.missing[(segid,)] = [atnum]
                amass = 1.0
            else:
                amass = FF.atomtypes[segid][iamass]
            # int(float('0.29')*100) = 28
            # decimal.Decimal(occ)*chrgrp_factor is ok
            # but decimal is present in 2.4 so let's stick to built-in round
            chrgid = round(float(occ)*chrgrp_factor)
            iatom = int(atnum)
            if atname[0] in string.digits: atname = atname[1:] + atname[0]
            output.append(ATOMLINE % (iatom, attype, int(resnum), resname,
                                      atname, chrgid, charge_factor*float(fac),
                                      amass))
            if self.atoms.has_key(iatom):
                raise ValueError('Recurring atom number %d' % iatom)
            if not FF.atomtypes.has_key(segid):
                raise ValueError('Unrecognized atom type "%s"' % segid) 
            self.atoms[iatom] = {'opls_type': attype,
                                 'bond_type': FF.atomtypes[segid][ibtype],
                                 'residue_name': resname,
                                 'residue_number': int(resnum),
                                 'atom_name': atname }
        output.append('\n')
        if len(output) > 3: self.output.append('\n'.join(output))

    def bondList(self):

        bondfunc = self.forcefield.ffsettings['bondfunc']
        # find bonds defined in pdb.conect records
        self.bonds = {}
        for ai, di in self.pdb.conect:
            bi = int(ai)
            for bond in [x > bi and (bi,x) or (x,bi) \
                         for x in map(int,di.split())]: self.bonds[bond] = 0
        self.bonds = self.bonds.keys()
        self.bonds.sort()

        # save list of bonds in gromacs.top format
        fields = ('i','j','ftype')
        output = ['[ bonds ]']
        output.append('; %5s  %6s %4s' %  fields)
        BONDLINE = ' %6d  %6d  %4d'
        BT = self.forcefield.bondtypes
        for bi, bj in self.bonds:
            output.append(BONDLINE % (bi, bj, bondfunc))
            key = [self.atoms[x]['bond_type'] for x in [bi,bj]]
            key.sort()
            key = tuple(key)
            if not BT.has_key(key):
                if self.missing.has_key(key):
                    self.missing[key].append((bi,bj))
                else: self.missing[key] = [(bi,bj)]
        output.append('\n')
        if len(output) > 3: self.output.append('\n'.join(output))

    def angleList(self):

        anglefunc = self.forcefield.ffsettings['anglefunc']
        # find angles for given list of bonds
        centralatoms = {}
        for ai,aj in self.bonds:
            if centralatoms.has_key(ai): centralatoms[ai].append(aj)
            else: centralatoms[ai] = [aj]
            if centralatoms.has_key(aj): centralatoms[aj].append(ai)
            else: centralatoms[aj] = [ai]
        alists = [[x]+centralatoms[x] for x in centralatoms.keys() \
                  if len(centralatoms[x]) > 1]
        self.angles = {}
        for atoms in alists:
            for first,last in Misc.uniqueCombinations(atoms[1:],2):
                aid = first < last and (first,atoms[0],last) or \
                      (last,atoms[0],first)
                self.angles[aid] = 0
        self.angles = self.angles.keys()
        self.angles.sort()

        # save list of angles in gromacs.top format
        fields = ('i','j','k','ftype')
        output = ['[ angles ]']
        output.append('; %5s  %6s  %6s %5s' % fields)
        ANGLELINE = ' %6d  %6d  %6d  %4d'
        AT = self.forcefield.angletypes
        for ai,aj,ak in self.angles:
            output.append(ANGLELINE % (ai, aj, ak, anglefunc))
            key = [self.atoms[x]['bond_type'] for x in [ai, aj, ak]]
            key = key[0] < key[-1] and tuple(key) or tuple(key.__reversed__())
            if not AT.has_key(key):
                if self.missing.has_key(key):
                    self.missing[key].append((ai,aj,ak))
                else: self.missing[key] = [(ai,aj,ak)]
        output.append('\n')
        if len(output) > 3: self.output.append('\n'.join(output))

    def dihedralList(self):

        dihedralfunc = self.forcefield.ffsettings['dihedralfunc']
        # find dihedrals for given list of angles
        cbonds = {}
        # store bonds which make angles
        # each bond is saved with atom indices ordered from smaller to larger
        for ai, aj, ak in self.angles:
            b1 = ai < aj and (ai, aj) or (aj, ai)
            b2 = aj < ak and (aj, ak) or (ak, aj)
            if cbonds.has_key(b1): cbonds[b1].append(b2)
            else: cbonds[b1] = [b2]
            if cbonds.has_key(b2): cbonds[b2].append(b1)
            else: cbonds[b2] = [b1]
        self.dihedrals = {}
        for b in cbonds.keys():
            # let's consider every bond as a central bond in a torsion
            a1, a2 = b
            a1b = [x for x in cbonds[b] if a1 in x]
            a2b = [x for x in cbonds[b] if a2 in x]
            # if there're fewer then 2 other bonds, no torsion possible
            if len(cbonds[b]) < 2 or len(a1b) == 0 or len(a2b) == 0: continue
            # loop over bonds containing the first atom from the central bond
            for ai, aj in a1b:
                dihedral= [a1,a2]
                # pick the 3rd atom to place in a (still) potential
                # dihedral definition
                if ai == a1: dihedral.insert(0,aj)
                else: dihedral.insert(0,ai)
                for ak, al in a2b:
                    # several torsions possible for given three atoms
                    d = copy.copy(dihedral)
                    # complete the definition of a torsion
                    # by placing the 4th atom
                    if ak == a2: d.append(al)
                    else: d.append(ak)
                    # reorder atom indices if needed
                    did = d[0] < d[-1] and tuple(d) or tuple(d.__reversed__())
                    self.dihedrals[did] = 0
        self.dihedrals = self.dihedrals.keys()
        self.dihedrals.sort()

        # save list of dihedrals in gromacs.top format
        fields = ('i','j','k','l','ftype')
        output = ['[ dihedrals ]']
        output.append('; %5s  %6s  %6s  %6s %5s' % fields)
        DIHEDRALLINE = ' %6d  %6d  %6d  %6d  %4d'    
        DT = self.forcefield.dihedraltypes
        for ai,aj,ak,al in self.dihedrals:
            key = [self.atoms[x]['bond_type'] for x in [ai, aj, ak, al]]
            if key[0] == key[-1]:
                key = key[1] < key[2] and tuple(key) or \
                      tuple(key.__reversed__())
            else:
                key = key[0] < key[-1] and tuple(key) or \
                      tuple(key.__reversed__())
            keyX1 = ('X',key[1],key[2],'X')    # wildcard
            keyX2 = ('X',key[2],key[1],'X')    # wildcard
            keyX3 = (key[0],key[1],key[2],'X') # wildcard
            keyX4 = (key[3],key[2],key[1],'X') # wildcard
            # keyX5 = tuple(list(keyX3).__reversed__()) # just in case
            # keyX6 = tuple(list(keyX4).__reversed__()) # just in case
            if not (DT.has_key(key)   or \
                    DT.has_key(keyX1) or DT.has_key(keyX2) or \
                    DT.has_key(keyX3) or DT.has_key(keyX4)):
                if self.missing.has_key(key):
                    self.missing[key].append((ai,aj,ak,al))
                else: self.missing[key] = [(ai,aj,ak,al)]
            line = DIHEDRALLINE % (ai, aj, ak, al, dihedralfunc)
            # handle dihedral definitions
            defline = self.isDefinedDihedral(ai,aj,ak,al)
            if defline: line = '%s  %s' % (line,defline)
            output.append(line)
        output.append('\n')
        if len(output) > 3: self.output.append('\n'.join(output))

    def pairList(self):
        pairfunc = self.forcefield.ffsettings['pairfunc']
        fields = ('i','j','ftype')
        output = ['[ pairs ]']
        output.append('; %5s  %6s %5s' % fields)
        PAIRLINE = ' %6d  %6d  %4d'
        for ai, aj, ak, al in self.dihedrals:
            output.append(PAIRLINE % (ai, al, pairfunc))
        output.append('\n')
        if len(output) > 3: self.output.append('\n'.join(output))

    def improperList(self):
        """ four atoms involved in an improper definition might be in at most
        three consequtive residues. In that way they can be used to
        override some proper dihedrals. """

        improperfunc = self.forcefield.ffsettings['improperfunc']
        # collect information of consequtive residues
        residues = []
        rnumb = None
        atidx = self.atoms.keys()
        atidx.sort()
        for ia in atidx:
            rnu = self.atoms[ia]['residue_number']
            rna = self.atoms[ia]['residue_name']
            ana = self.atoms[ia]['atom_name']
            if rnu != rnumb:
                rnumb = rnu
                residues.append((rna,{}))
            residues[-1][1][ana] = ia
        # map improper definitions from atom names to atom numbers
        P = self.pdb.imp
        fields = ('i','j','k','l','ftype','name')
        output = ['[ dihedrals ]']
        output.append('; %5s  %6s  %6s  %6s %5s    %-6s' % fields)
        IMPROPERLINE = ' %6d  %6d  %6d  %6d   %3d    %-20s ; %s'
        for i in range(len(residues)):
            rname, atoms = residues[i]
            if not self.pdb.imp.has_key(rname): continue
            for imp in P[rname]:
                iname = P[rname][imp]
                atidx = []
                for ia in range(len(imp)):
                    aname = imp[ia]
                    if aname[0] == '+':
                        if i < (len(residues)-1):
                            rnu, aname = i + 1, aname[1:]
                        else: break
                    elif aname[0] == '-':
                        if i > 0:
                            rnu, aname = i - 1, aname[1:]
                        else: break
                    elif aname[0] == '*':
                        # Here we handle non-sequential crosslinks.
                        # imp[2] is a central atom in the imp def.
                        # To establish identity of imp[ia] atom
                        # we need to find a bond which connects it with
                        # imp[2] atom, and knowing its number we can
                        # find the residue it belongs to.
                        if ia == 2 or imp[2][0] in '+-*':
                            # FIXME: it's difficult to handle all possible
                            # improper definitions
                            break
                        atom_cross = residues[rnu][1][imp[2]]
                        bond = [b for b in self.bonds if atom_cross in b]
                        rnu = None
                        for ai,aj in bond:
                            if ai != atom_cross: a = self.atoms[ai]
                            else: a = self.atoms[aj]
                            # let's exclude atoms from the same residue
                            # there're no reason in marking'em with *
                            # in such a case
                            if (a['residue_number'] - 1) == i: continue
                            if a['atom_name'] == aname[1:]:
                                rnu = a['residue_number'] - 1
                                break
                        if rnu is None: break # ok, we have to give up here!
                        aname = aname[1:]
                    else: rnu = i
                    if residues[rnu][1].has_key(aname):
                        atidx.append(residues[rnu][1][aname])
                    else: 
                        print rnu, residues[rnu]
                        break
                if len(atidx) != 4:
                    print ' Warning! failed to match "%s" in %s' % \
                          (str(imp),rname)
                    print atidx
                    continue
                ai, aj, ak, al = atidx
                comment = '%s.%-4d:  %s' % (rname,i+1,'  '.join(imp))
                output.append(IMPROPERLINE % (ai,aj,ak,al,improperfunc,
                                              iname,comment))
        output.append('\n')
        if len(output) > 3: self.output.append('\n'.join(output))

    def saveTopologyFile(self, filename):
        open(filename,'w').write(self.get())
        self.validate()

    def validate(self,out=None,verbose=0):
        if out is None: out = sys.stdin
        labels = ['undefined','atom types','bond types',
                  'angle types', 'dihedral types']
        warnings = {'atom types': {}, 'bond types': {},
                    'angle types': {}, 'dihedral types': {}, 'undefined': {}}
        for k in self.missing.keys():
            warnings[labels[len(k) < 5 and len(k) or 0]][k] = self.missing[k]
        for label in labels:
            if not warnings[label]: continue
            print 'missing %s (%d):' % (label,len(warnings[label]))
            for k in warnings[label].keys():
                print ' -> ',k
                if verbose:
                    print '\t'+'\n\t'.join(map(str,warnings[label][k]))

    def queryTorsion(self,pattern=None):
        if pattern is None:
            pattern = re.compile('C\w+\s+C\w+\s+C\w+\s+C\w+')
        A = self.atoms
        for dihedral in self.dihedrals:
            btypes = [A[x]['bond_type'] for x in dihedral]
            names = [A[x]['atom_name'] for x in dihedral]
            resnumber = [A[x]['residue_number'] for x in dihedral]
            resname = [A[x]['residue_name'] for x in dihedral]
            hit = None
            if re.match(pattern, ' '.join(btypes)):
                hit = '%5s %5s %5s %5s' % tuple(names)
            elif re.match(pattern, ' '.join(btypes.__reversed__())): 
                hit = '%5s %5s %5s %5s' % tuple(names.__reversed__())
                btypes.reverse() ; resnumber.reverse() ; resname.reverse()
            if hit:
                hit += ' ::: '
                for btyp, rnam,rnum in zip(btypes,resname, resnumber):
                    hit += '%5s.%3s.%-4d' % (btyp,rnam,int(rnum))
                print hit

    def torsionREP(self, ai, aj, ak, al, pattern=None):
        if pattern is None:
            rep = re.compile('%s\s+%s\s+%s\s+%s' % (ai,aj,ak,al))
        else:
            rep = re.compile('%s%s\s+%s%s\s+%s%s\s+%s%s' % \
                             (ai,pattern,aj,pattern,ak,pattern,al,pattern))
        return rep

if __name__ == '__main__':

    p2t = Pdb2Top('molecules/pmb.pdb')
    print p2t()
    p2t = Pdb2Top('molecules/lipidA_03.pdb')
    x = p2t()
    pat_CT_CT_OH_HO = re.compile('C\w+\s+C\w+\s+O\w+\s+H\w+')
    p2t.queryTorsion(pat_CT_CT_OH_HO)

