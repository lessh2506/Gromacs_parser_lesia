from mdconv import *
from Scientific.IO import TextFile

# by MMTK.Biopolymers
# including terminal caps
VALID_AMINOACIDS = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS',
                    'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP',
                    'TYR','VAL','CYX','HSD','HSE','HSP','HID','HIE','HIP',
                    'NHE','NME','ACE']

TextFile.gzip = None # FIX, gzip.GzipFile.file return instance not file object (file.fileobj)

class AmberParmTopologyFile:

    def __init__(self,filename):
        self.__reset()
        file = TextFile.TextFile(filename)
        self.info = {'TITLE': file.readline()}
        header = getDataBlock(file,'I6',30)
        names = ['NATOM','NTYPES','NBONH','MBONA','NTHETH','MTHETA','NPHIH',
                 'NHPARM','NPARM', # these two are not used
                 'MPHIA','NEXT','NRES',
                 'NBONA','NTHETA','NPHIA','NUMBND','NUMANG','NPTRA', 'NATYP',
                 'NPHB','IFPERT','NBPER','NGPER','NDPER',
                 'MBPER','MGPER','MDPER','IFBOX','NMXRS','IFCAP']
        for i in range(len(names)):
            key = names[i]
            self.info[key] = header[i]

        self.igraph = getDataBlock(file,'A4',self.info['NATOM'])
                    # The user atoms names
        self.chrg = getDataBlock(file,'F16',self.info['NATOM'])
            # The atom charges
        self.amass = getDataBlock(file,'F16',self.info['NATOM'])
            # The atom masses
        self.iac = getDataBlock(file,'I6',self.info['NATOM'])
            # index for the atom types involved in Lennard Jones (6-12)
            # interactions. See ICO below.
        self.numex = getDataBlock(file,'I6',self.info['NATOM'])
            # total number of excluded atoms for atom 'i'.
            # See NATEX bellow.
        self.ico = getDataBlock(file,'I6',
                                self.info['NTYPES']*self.info['NTYPES'])
            # provides the index to the nonbon parameter arrays CN1,
            # CN2 and ASOL, BSOL. All possible 6-12 or 10-12 atoms type
            # interactions are represented.
            # NOTE: A particular atom type can have either a 10-12 or
            # a 6-12 interaction, but not both. The index is
            # calculated as follows:
            #
            #   index = ICO(NTYPES*(IAC(i)-1) + IAC(j))
            # 
            # If index is positive, this is an index into the 6-12
            # parameter arrays (CN1 and CN2) otherwise it is
            # an index into the 10-12 parameter arrays (ASOL & BSOL)
        self.labres = getDataBlock(file,'A4',self.info['NRES'])
            # the residue labels
        self.ipres = getDataBlock(file,'I6',self.info['NRES'])
            # atoms in each residue are listed for atom i in IPRES(i)
            # to IPRES(i+1)-1 ... (e.g. ipres[0]=ipres(1)=1,
            # ipres[1]=ipres(2)=9 ... 1st res has 8 atoms). In other words,
            # for each residue the index of first atom is given.
        self.rk = getDataBlock(file,'F16',self.info['NUMBND'])
            # force constant for the bonds of each type, kcal/mol
        self.req = getDataBlock(file,'F16',self.info['NUMBND'])
            # equilibrium bond length for the bonfs of each type,
            # angstroms
        self.tk = getDataBlock(file,'F16',self.info['NUMANG'])
            # force constant for angles of each type, kcal/mol A**2
        self.teq = getDataBlock(file,'F16',self.info['NUMANG'])
            # equilibrium angle for the angles of each type, degrees
        self.pk = getDataBlock(file,'F16',self.info['NPTRA'])
            # force constant for the dihedrals of each type, kcal/mol
        self.pn = getDataBlock(file,'F16',self.info['NPTRA'])
            # periodicity of the dihedral of a given type
        self.phase = getDataBlock(file,'F16',self.info['NPTRA'])
            # phase of the dihedral of a given type
        self.solty = getDataBlock(file,'F16',self.info['NATYP'])
            # currently unused Av5 (anyway must be read)
        self.cn1 = getDataBlock(file,'F16', \
                                self.info['NTYPES']*(self.info['NTYPES']+1)/2)
            # Lennard Jones r**12 terms for all possible atom type
            # interactions, indexed by ICO and IAC; for atom i and j where
            # i<j, the index into this array is as follows 
            # (assuming in index is positive):
            # cn1(ico(ntypes*iac(i)-1+iac(j)))
        self.cn2 = getDataBlock(file,'F16', \
                                self.info['NTYPES']*(self.info['NTYPES']+1)/2)
            # Lennard Jones r**6 terms for all possible atom type
            # interactions. Indexed like CN1 above.
            # -----------------------------------------------------------------
            # NOTE: the atom numbers in the arrays which follow that describe
            # bonds, angles, and dihedrals are obfuscated by the following 
            # formula (for runtime speed in indexing arrays and to make the 
            # user confused...). The true atom number equals the absolute value
            # of the number divided by three, plus one. In the case of
            # dihedrals, if the third atom is negative, this implies an
            # improper torsion and if the fourth atom is negative, this implies
            # that end group interactions are to be ignored.
            # End group interactions are ignored, for example,
            # in dihedrals of various ring systems (to prevent
            # double counting) and in multiterm dihedrals.
            # ----------------------------------------------------------------
        self.bonh = getDataBlock(file,'3I6',self.info['NBONH'])
            # indexes of atoms (one of them is hydrogen) forming a covalent 
            # bond; data are orgenised as follows: each single tuple
            # consists of indexes of atom i, atom j and an index into parameter
            # arrays RK and REQ
        self.bona = getDataBlock(file,'3I6',self.info['NBONA'])
            # as above but neither of atoms is hydrogen
        self.theth = getDataBlock(file,'4I6',self.info['NTHETH'])
            # tuple[i]: [0] atom involved in angle i, angle contains hydrogen;
            # [1] atom involved in angle i; [2] atom involved in angle i; 
            # [3] index into parameter arrays TK and TEQ for angle [0]-[1]-[2]
        self.theta = getDataBlock(file,'4I6',self.info['NTHETA'])
            # as above but neither of atoms is hydrogen
        self.phih = getDataBlock(file,'5I6',self.info['NPHIH'])
            # tuple[i]: [0]-[1]-[2]-[3] atoms involved in dihedral i,
            # angle contains hydrogen; [4] index into parameter arrays PK, PN,
            # and PHASE for this dihedral
        self.phia = getDataBlock(file,'5I6',self.info['NPHIA'])
            # as above but neither of atoms is hydrogen. NOTE: if the 
            # periodicity (PN) is negative, this implies the following entry
            # in the PK, PN, and PHASE arrays is another term in a multitermed
            # dihedral.
        self.natex = getDataBlock(file,'I6',self.info['NEXT'])
            # the excluded atom list. To get the excluded list for atom i
            # you need to traverse the NUMEX list, adding up all the previous
            # NUMEX values, since NUMEX(i) holds the number of excluded atoms
            # for atom i, not the index into the NATEX list.
            # Let IEXCL = SUM(NUMEX(j),j=1,i-1), then excluded atoms are
            # NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).
        self.asol = getDataBlock(file,'F16',self.info['NPHB'])
            # the value for r**12 term for hydrogen bonds of all possible
            # types. Index into these arrays is equivalent to the CN1 and
            # CN2 arrays, however the index is negative. For example,
            # for atoms i and j, with i<j, the index is:
            # -(NTYPES*IAC(i)-1+IAC(j))
        self.bsol = getDataBlock(file,'F16',self.info['NPHB'])
            # the value for r**10 term for hydrogen bonds of all possible
            # types. Indexed like ASOL.
        self.hbcut = getDataBlock(file,'F16',self.info['NPHB'])
            # no longer in use but read anyway
        self.isymbl = getDataBlock(file,'A4',self.info['NATOM'])
            # the AMBER atom types for each atom
        self.itree = getDataBlock(file,'A4',self.info['NATOM'])
            # the list of tree joining information, classified into five
            # types. M -- main chain, S -- side chain, B -- branch point,
            # 3 -- branch into three chains, E -- end of the chain.
            # M can have 1,2,3 or 4 atoms connected to it; E can have 1
            # atom connected; S - 2 atoms; B - 3 atoms; 3 - 4 atoms...
            # In A5 there are additionally types 4,5,6...
        self.join = getDataBlock(file,'I6',self.info['NATOM'])
            # tree joining information, potentially used in ancient
            # analysis programs (MMTK is quite young, though).
            # It can be used to extract information about, for example,
            # the number of molecules and the topology of each of them.
            # This is valid for nonperiodic systems. If IFBOX==1 this
            # this information is given explicitly below (NSPM,NSP).
        self.irotat = getDataBlock(file,'I6',self.info['NATOM'])
            # apparently the last atom that would move if atom i was
            # rotated, however the meaning has been lost over time.
        if self.info['IFBOX'] > 0:
            # 1 rectangular
            # 2 truncated octahedron
            self.info['IPTRES'], self.info['NSPM'], self.info['NSPSEL'] = \
                                 getDataBlock(file,'I6',3)
            self.nsp = getDataBlock(file,'I6',self.info['NSPM'])
            # the total number of atoms in each molecule,
            # necessary to correctly determine the pressure scaling
            box = getDataBlock(file,'F16',4)
            self.info['BETA'], self.info['BOX'] = box[0], tuple(box[1:])
            # BETA: periodic box, angle between the XY and YZ planes in degrees
            # BOX: the periodic box lengths in the X, Y, and Z directions
        if self.info['IFCAP'] > 0:
            self.info['NATCAP'] = getDataBlock(file,'I6',1)
            cap = getDataBlock(file,'F16',4)
            self.info['CUTCAP'], self.info['CAPXYZ'] = cap[0], tuple(cap[1:])
            # CUTCAP: the distance from the center of the cap to the outside
            # CAPXYZ: xyz coordinate for the center of the cap
        if self.info['IFPERT'] > 0:
            raise ' perturbation part of parmtop file is not implemented yet'
        file.close()

    def __reset(self):
        for name,value in [('_molinfo',None)]:
            setattr(self,name,value)

    def getMolInfo(self):

        if self._molinfo is None: 
            residues = self.getResidues()
            mol, order, bonds = self.getMolecules()
            self._molinfo = {'Residues': residues,
                             'Bonds': bonds,
                             'Molecules':mol,
                             'Order of molecules':order}
        return self._molinfo

    def getResidues(self):

        if self._molinfo is not None:
            return self._molinfo['Residues']
        #
        # Let's find the number of residues
        # with different labels
        #
        diffres = {}
        for i in self.labres:
            try: diffres[i] = diffres[i] + 1
            except KeyError:  diffres[i] = 1
        self.info['INUMRES'] = len(diffres)
        if diffres.has_key('WAT '):
            self.info['IWATMOL'] = diffres['WAT ']
        elif diffres.has_key('HOH '):
            self.info['IWATMOL'] = diffres['HOH ']
        else:
            self.info['IWATMOL'] = 0 
        #
        # Find the number of atoms in every unique residue
        # (terminal AA should be handled correctly);
        # in general it's possible that more then one residue
        # with the same name has diffrent number of atoms
        # (terminal AAs are a typical example). Basing on ATOMRES
        # we thus will update the number of REALLY different residues
        # in DIFFRES.
        #
        atomres = {}
        unique  = {}
        for i in range(len(self.ipres)):
            try: ires = tuple(range(self.ipres[i],self.ipres[i+1])) 
            except: 
                 ires = tuple(range(self.ipres[i],self.info['NATOM']+1))
            iname = self.labres[i]
            if not atomres.has_key(iname):
                atomres[iname] = [ires,]
                unique[iname]  = [iname,]
            else:
                hit = None
                for ia in range(len(unique[iname])):
                    if len(atomres[unique[iname][ia]][0]) == len(ires):
                        hit = unique[iname][ia]
                if not hit: hit = iname+str(len(unique[iname]))
                try: atomres[hit].append(ires)
                except KeyError:
                    atomres[hit] = [ires,]
                    unique[iname].append(hit)
        for ia in atomres.keys(): diffres[ia] = len(atomres[ia])

        return {'summary':diffres, 'atoms':atomres}

    def getResidueTopology(self,res_name,res_numb=0):
        """
        define topology for a residue with
        given name (res_name) and number (res_numb)
        (resnumb varies between 0 and the number of res_name
        residues in the system)
        """
        residues = self.getResidues()
        diffres = residues['summary']
        atomres = residues['atoms']
        try:
            bona = Numeric.array(self.bona)  # non-hydrogen atoms
            bona[:,0] = bona[:,0]/3 + 1      # the "true" indices
            bona[:,1] = bona[:,1]/3 + 1
        except IndexError: pass
        try:
            bonh = Numeric.array(self.bonh)  # hydrogen ones
            bonh[:,0] = bonh[:,0]/3 + 1
            bonh[:,1] = bonh[:,1]/3 + 1
        except IndexError: pass
        try: atomlst = atomres[res_name][res_numb]
        except IndexError: return 0
        intra_bonds = []
        inter_bonds = []
        for ib in atomlst:
            try:
                bona_hit = Numeric.where(Numeric.equal(bona,ib),bona,0)
                row_hita = Numeric.nonzero(bona_hit[:,0])
                for ic in range(len(row_hita)):
                    bond = bona[row_hita[ic]][:2]
                    if bond[1] <= atomlst[-1]: intra_bonds.append(bond)
                    else: inter_bonds.append(bond)
            except: pass
            try:
                bonh_hit = Numeric.where(Numeric.equal(bonh,ib),bonh,0)
                row_hith = Numeric.nonzero(bonh_hit[:,0])
                for ic in range(len(row_hith)):
                    bond = bonh[row_hith[ic]][:2]
                    if bond[1] <= atomlst[-1]: intra_bonds.append(bond)
                    else: inter_bonds.append(bond)
            except: pass
        intra_bonds = map(lambda v: v.tolist(),intra_bonds) # obvious...
        inter_bonds = map(lambda v: v.tolist(),inter_bonds) # with others
                                                            # residues
        return intra_bonds, inter_bonds

    def getAADefinition(self,res_name,bonds,res_numb=0):
        """ MMTK style: protein (for a backbone) and
        sidechain (for the rest)
        """
        pass

#     def getResidueDefinition(self,res_name,bonds,res_numb=0,
#                              diffres=None,atomres=None):
#         """
#         get the MMTK definition of a given residue and
#         return it as a list which can be printed or
#         easily transformed into evaluable Python code
#         """
#         definition = []
#         if diffres is None or atomres is None:
#             diffres, atomres = self.getResidues()
#         try: atoms = atomres[res_name][res_numb]
#         except IndexError: return 0
#         for ia in atoms:
#             line = string.strip(self.igraph[ia-1])+' = Atom(\''+\
#                    self.isymbl[ia-1][0]+'\')'
#             definition.append(line)
#         line = 'bonds = ['
#         for ia in bonds:
#             line = line + 'Bond('+string.strip(self.igraph[ia[0]-1])+\
#                    ','+string.strip(self.igraph[ia[1]-1])+'), '
#         line = line + ']'
#         definition.append(line)
#         line = 'name = \''+string.strip(res_name)+'\''
#         definition.append(line)
#         line = 'pdbmap = [(\''+res_name[:3]+'\', {'
#         for ia in atoms:
#             aname = string.strip(self.igraph[ia-1])
#             if len(aname)>3: aname = aname[-1] + aname[:3]
#             line = line + '\''+aname+'\': ' +\
#                    string.strip(self.igraph[ia-1])+', '
#         line = line + ' }, ), ]'
#         definition.append(line)
#         line = 'amber_atom_type = {'
#         for ia in atoms:
#             line = line + string.strip(self.igraph[ia-1])+': \''+\
#                    string.strip(self.isymbl[ia-1])+'\', '
#         line = line + ' }'
#         definition.append(line)
#         line = 'amber_charge = {'        # opls_charge ?
#         for ia in atoms:
#             line = line + string.strip(self.igraph[ia-1])+\
#                    ': '+str(round(self.chrg[ia-1]/18.2223,5))+', '
#         line = line + ' }'
#         definition.append(line)
#         return len(atoms),definition

    def getGroupDefinition(self,res_name,ff_name='opls',res_numb=0):
        """
        get the MMTK definition of a given residue and
        return it as a list which can be printed or
        easily transformed into evaluable Python code
        """
        if self._molinfo is None: gj = self.getMolInfo()
        residues = self.getResidues()
        atomres = residues['atoms']
        bonds = self._molinfo['Bonds']['intra'][res_name]
        
        definition = []
        try: atoms = atomres[res_name][res_numb]
        except IndexError: # return 0
            if not atomres.has_key(res_name):
                error = ' Residue \"%s\" not found ' % res_name
            else:
                error = ' There is only %d residues named \"%s\"' % \
                        (len(atomres),res_name)
        for ia in atoms:
            line = string.strip(self.igraph[ia-1])+' = Atom(\''+\
                   self.isymbl[ia-1][0]+'\')'
            definition.append(line)
        line = 'bonds = ['
        for ia in bonds:
            line = line + 'Bond('+string.strip(self.igraph[ia[0]-1])+\
                   ','+string.strip(self.igraph[ia[1]-1])+'), '
        line = line + ']'
        definition.append(line)
        line = 'name = \''+string.strip(res_name)+'\''
        definition.append(line)
        line = 'pdbmap = [(\''+res_name[:3]+'\', {'
        for ia in atoms:
            aname = string.strip(self.igraph[ia-1])
            if len(aname)>3: aname = aname[-1] + aname[:3]
            line = line + '\''+aname+'\': ' +\
                   string.strip(self.igraph[ia-1])+', '
        line = line + ' }, ), ]'
        definition.append(line)
        line = '%s_atom_type = {' % ff_name
        for ia in atoms:
            line = line + string.strip(self.igraph[ia-1])+': \''+\
                   string.strip(self.isymbl[ia-1])+'\', '
        line = line + ' }'
        definition.append(line)
        line = '%s_charge = {' % ff_name
        for ia in atoms:
            line = line + string.strip(self.igraph[ia-1])+\
                   ': '+str(round(self.chrg[ia-1]/18.2223,5))+', '
        line = line + ' }'
        definition.append(line)
        return definition

#     def getMoleculeDefinition(self,resseq,resdefnames,interbonds,\
#                               symbol,molname):
#         """
#         get the MMTK definition of a molecule which is defined
#         by a list of residues RESSEQ (filenames with Group
#         definitions which match these to be saved in a MMTK Database)
#         """
#         # the name scheme is that an underscore and the sequential number
#         # is added, this is however required only to avoid ambiguity,
#         # else stripped names of residues should be enough
#         gj = {}
#         for i in resseq:
#             try: gj[i] = gj[i]+1
#             except KeyError: gj[i] = 1
#         definition = []
#         bond = {}
#         for i in range(1,len(resseq)): # it's assumed ONE inter bond
#             name = resseq[i-1]
#             snameA = string.replace(string.strip(name[:4]),' ','_')
#             snameB = string.replace(string.strip(resseq[i][:4]),' ','_')
#             if len(gj) != len(resseq):
#                 # there are residues in the sequence
#                 # which have the same name
#                 snameA = snameA + "_" + str(i)
#                 snameB = snameB + "_" + str(i+1)
#             atomA = snameA + "." +\
#                     string.strip(self.igraph[interbonds[name][0][0]-1])
#             atomB = snameB + "." +\
#                     string.strip(self.igraph[interbonds[name][0][1]-1])
#              # (i,i+1) ? should be checked to which residue a given atom
#              # atomB belongs (on the other hand, we are processing
#              # SEQUENCE so the only information to be checked is
#              # which atom of the first residue is connected to which atom of
#              # the second residue
#             bond[name] = 'Bond(' + atomA + ',' + atomB +'), '
#             line = snameA + ' = Group(\'' + resdefnames[name] + '\')'
#             definition.append(line)
#         name = string.replace(string.strip(resseq[-1][:4])," ","_")
#         if len(gj) != len(resseq):
#             # there are residues in the sequence
#             # which has the same name
#             name = name + "_" + str(len(resseq))
#         definition.append(name+' = Group(\''+ resdefnames[resseq[-1]]+'\')')
#         line = 'bonds = ['
#         for i in bond.keys(): line = line + bond[i] # sequence unimportant?
#         line = line + ']'
#         definition.append(line)
#         definition.append('symbol = \''+symbol+'\'')
#         definition.append('name = \''+molname+'\'')
#         return len(resseq),definition

    def getMoleculeDefinition(self,molid,symbol='NN',molname='No-Name'):
        """
        get the MMTK definition of a molecule which is defined
        by a list of residues RESSEQ (filenames with Group
        definitions which match these to be saved in a MMTK Database)
        """
        # the name scheme is that an underscore and the sequential number
        # is added, this is however required only to avoid ambiguity,
        # else stripped names of residues should be enough
        if self._molinfo is None: gj = self.getMolInfo()
        moltypes = len(self._molinfo['Order of molecules'])
        if molid >= moltypes:
            raise NotImplementedError, \
                  ' Invalid molecule number (allowed 0..%d)' % moltypes-1
        resseq = self._molinfo['Order of molecules'][molid][0]
        interbonds = self._molinfo['Bonds']['inter']
        
        gj = {}
        for i in resseq:
            try: gj[i] = gj[i]+1
            except KeyError: gj[i] = 1
        definition = []
        bond = {}
        for i in range(1,len(resseq)): # it's assumed ONE inter bond
            name = resseq[i-1]
            print name
            snameA = string.replace(string.strip(name[:4]),' ','_')
            snameB = string.replace(string.strip(resseq[i][:4]),' ','_')
            if len(gj) != len(resseq):
                # there are residues in the sequence
                # which have the same name
                snameA = snameA + "_" + str(i)
                snameB = snameB + "_" + str(i+1)
            atomA = snameA + "." +\
                    string.strip(self.igraph[interbonds[name][0][0]-1])
            atomB = snameB + "." +\
                    string.strip(self.igraph[interbonds[name][0][1]-1])
             # (i,i+1) ? should be checked to which residue a given atom
             # atomB belongs (on the other hand, we are processing
             # SEQUENCE so the only information to be checked is
             # which atom of the first residue is connected to which atom of
             # the second residue
            bond[name] = 'Bond(' + atomA + ',' + atomB +'), '
            #line = snameA + ' = Group(\'' + resdefnames[name] + '\')'
            line = snameA + ' = Group(\'' + name.strip() + '\')'
            definition.append(line)
        name = string.replace(string.strip(resseq[-1][:4])," ","_")
        if len(gj) != len(resseq):
            # there are residues in the sequence
            # which has the same name
            name = name + "_" + str(len(resseq))
        #definition.append(name+' = Group(\''+ resdefnames[resseq[-1]]+'\')')
        definition.append(name+' = Group(\''+ resseq[-1].strip() +'\')')
        line = 'bonds = ['
        for i in bond.keys(): line = line + bond[i] # sequence unimportant?
        line = line + ']'
        definition.append(line)
        definition.append('symbol = \''+symbol+'\'')
        definition.append('name = \''+molname+'\'')
        return definition

    def getProteinDefinition(self,pid,model='polar',c_terminus=1,
                             n_terminus=1,circular=0):
        if self._molinfo is None: gj = self.getMolInfo()
        resseq = self._molinfo['Order of molecules'][pid][0]
        s3 = map(lambda x: string.strip(x[:min(len(x),4)]), resseq)
        otype = self.getObjectType(resseq)
        if otype != 'Protein':
            raise NotImplementedError, \
                  ' The following sequence contains residue '\
                  ' which has not been recognised as an aminoacid'\
                  ' (see MMTK.Biopolymers.defineAminoAcidResidue)\n'\
                  "-".join(s3)
        definition = []
        # FIX: inter residue connections
        definition.append('aaseq = %s' % str(s3))
        definition.append('chains = [ PeptideChain(aaseq,'\
                          'model=\'%s\',' % model + \
                          'c_terminus=%d,' % c_terminus + \
                          'n_terminus=%d,' % n_terminus + \
                          'circular=%d) ]' % circular)
        definition.append('del aaseq')
        return definition

    def getMolecules(self):
        """
        find the number of molecules and their topology,
        inter is a dict containing all inter_bonds
        (thus it requires prior evaluating getResidueTopology
        for all unique residues and (!!!) res_numb = 0)
        """
        residues = self.getResidues()
        diffres = residues['summary']
        atomres = residues['atoms']
        
        bter  = {}; btra = {}
        for i in atomres.keys():
            btra[i], bter[i] = self.getResidueTopology(i)
        molname = ""; crosslink = {}; molecules = {}; molorder = []; res1st = 0
        # --- mapping between atom numbers and residue numbers
        result = {}
        for res in range(self.info['NRES']-1):
            for at in range(self.ipres[res],self.ipres[res+1]):
                result[at] = [res,self.labres[res]]
        for at in range(self.ipres[self.info['NRES']-1],self.info['NATOM']+1):
                result[at] = [res,self.labres[res]]
        _atomresmap = result
        # ---
        for ia in range(self.info['NRES']):
            try: natom = self.ipres[ia+1] - self.ipres[ia]
            except: natom = self.info['NATOM'] - self.ipres[ia] + 1 # last
            name = self.labres[ia]
            ix = 0
            while len(atomres[name][0]) != natom:
               ix = ix+1
               name = '%s%d' % (name,ix)
               try: natom = len(atomres[name][0])
               except KeyError: return 0
            molname = molname + "&" + name
             # len(bter):
             #      0 - a terminal residue or an isolated molecule
             #      1 - a residue in a sequence
             #      2 - a cross-linked residue in a sequence
             #      what about 3 or more?
            if len(bter[name]) == 2:
                if crosslink.has_key(name): crosslink[name].append(ia)
                else: crosslink[name] = [ia]
            if len(bter[name]) == 0:
                # CYY -- PAL: PAL is after NHE (bter=0)
                # multichains proteins ?
                if len(crosslink) > 0:
                    for namres in crosslink.keys():
                        for ires in range(len(crosslink[namres])):
                            bondres = self.getResidueTopology(namres,ires,
                                                              diffres,
                                                              atomres)[1]
                            for at1, at2 in bondres:
                                gj = max(at1,at2)
                                print crosslink[namres][ires], namres,
                                if _atomresmap[gj][0] > ia:
                                    print 'cross (extra res)',
                                else: print 'cross (intra prot)',
                                print _atomresmap[gj][0],
                                print _atomresmap[gj][1]
                if 1:
                    try: molecules[molname] = molecules[molname] + 1
                    except KeyError: molecules[molname] = 1
                    mn_gj = tuple(string.split(molname[1:],"&"))
                    try:
                        if molorder[-1][0] != mn_gj: molorder.append([mn_gj,1])
                        else: molorder[-1][1] = molorder[-1][1] + 1
                    except IndexError: molorder.append([mn_gj,1])
                molname = ""; crosslink = {}; res1st = ia+1
        result = []
        for ia in molecules.keys():
            res = string.split(ia[1:],"&")
            result.append((res,molecules[ia]))
        return result, molorder, {'inter': bter, 'intra': btra} #, _atomresmap

    def getObjectType(self,reslist):
        """ based on the sequence of residues find the type of molecule;
        e.g. *Molecule* or *Protein* """

        s3 = map(lambda x: string.strip(x[:min(len(x),4)]),reslist)
        for r in s3:
            if VALID_AMINOACIDS.count(r) == 0: return 'Molecule'
        return 'Protein'

    def getMolMassAndCharge(self,reslist,resdict,parm):
        """ no coment """

        mass = charge = natom = 0
        for r in reslist:
            natom = natom + len(resdict[r][0])
            mass = mass + add.reduce(take(self.amass,array(resdict[r][0])-1))
            charge = charge + add.reduce(take(self.chrg,
                                         array(resdict[r][0])-1))/_factor
        if abs(charge) < _epsilon: charge = 0.0
        return mass, charge, natom


class AmberTrajectoryFile:

    def __init__(self,filename,natoms,pbc=''):

        self.__FIELDWIDTH = 8
        self.__FIELDTYPE = 'F'
        self.__NFIELDS = 10
        
        self.file = TextFile.TextFile(filename)
        self.info = {'TITLE': self.file.readline()}
        self.__nlines = natoms*3/self.__NFIELDS
        self.__pbc = pbc
        self.__natoms = natoms
        if natoms*3%self.__NFIELDS > 0:
            self.__nlines += 1
        if pbc: self.__nlines += 1

    def getChunk(self,skip=0):

        output = {'data': None, 'box': None}

        def cleanUp(self):
            try: self.file.close()
            except: pass
            self.file = None

        if skip != 0:
            for i in range(self.__nlines):
                line = self.file.readline()
            if not line: cleanUp(self)
            return output
        else:
            try: block = getDataBlock(self.file,'3F8',self.__natoms)
            except NotImplementedError:
                cleanUp(self)
                return output
            output['data'] = Numeric.reshape(block,(self.__natoms,3))

        if self.__pbc:
            try:
                # bez znaczenia czy brick/truncoct: zawsze 3F8
                output['box'] = getDataBlock(self.file,'3F8',1) 
            except NotImplementedError:
                cleanUp(self)
                return output
        else: output['box'] = None
        return output

def readAmberTrajectory(filename,
                        natoms,
                        frame_set=(0,None,1),
                        trajtype='xyz',
                        pbc = '',
                        verbose=0,
                        maxframes=100):

    FACTORS = {'xyz': Units.Ang,
               'vel': Numeric.sqrt(4.184) }
    
    output = {'data': [], 'box': [], 'eof': 0}
    fin1, finN, finS = frame_set
    if finN is None: finN = maxframes*finS + fin1

    opts = {'pbc': pbc}
    atraj = apply(AmberTrajectoryFile,[filename,natoms],opts)

    frame_it = 0 # first frame has index 0
    while 1:
        skip = (frame_it-fin1)%finS or frame_it<fin1
        if verbose:
            print frame_it, skip
        chunk = atraj.getChunk(skip=skip)
        if atraj.file is None:
            output['eof'] = 1
        if frame_it >= finN or atraj.file is None: break
        frame_it += 1
        if chunk['data'] is None: continue
        if trajtype == 'xyz':
            if chunk['box'] is not None:
                output['box'].append(\
                    Numeric.array(chunk['box'])*FACTORS['xyz'])
            else: pass # nonperiodic ?
            data = FACTORS['xyz']*Numeric.reshape(chunk['data'],(natoms,3))
        elif trajtype == 'vel':
            data = FACTORS['vel']*Numeric.reshape(chunk['data'],(natoms,3))
        else: data = []
        output['data'].append(data)
    output['frames'] = frame_it
    return output

        
def assignAmberIndices(universe,parm):

    def isdigit(char):
        try: return string.atoi(char)
        except ValueError: return None

    resseq = []
    logfile = []
    for i in universe:
        if i.__class__.__name__ == 'Protein':
            for ia in i.chains:
                for ib in ia[0].groups:
                    try: altmap = ib.pdb_alternative
                    except: altmap = None
                    for ic in ib.pdbmap:
                        atomlist = []
                        for id in ic[1].keys():
                            name = [id]
                            if altmap:
                                for ie in altmap.keys():
                                    if altmap[ie] == id: name.append(ie)
                            atomlist.append((name,ib.atomList()\
                                             [ic[1][id].number]))
                        resseq.append((ic[0],atomlist))
        elif i.__class__.__name__ == 'Molecule':
            for ia in i.pdbmap: # who cares for pdb_alternative...
                atomlist = []
                for ib in ia[1].keys():
                    atomlist.append(([ib],i.atomList()[ia[1][ib].number]))
                resseq.append((ia[0],atomlist))
        else:
            print i.__class__.__name__,' is unsupported ... yet'
            return # fix
    logfile.append(('C',
                    str(len(resseq))+' <-> in amber: '+str(len(parm.labres))))
    natom = parm.info['NATOM']
    for i in range(len(resseq)):
        logi = []
        first = parm.ipres[i]
        try: last = parm.ipres[i+1]
        except: last = natom + 1
        anum = last - first
        logi.append([i,resseq[i][0],first-1,last-2])
        names = map(string.strip,parm.igraph[first-1:last-1])
        names = map(string.upper,names)
        logi.append(names)
        for ia in resseq[i][1]:
            logi.append([ia[0]])
            index = None
            ibul = map(string.upper,ia[0])
            for ib in ibul: # ia[0]
                try: index = names.index(ib)
                except ValueError:
                    if isdigit(ib[-1]): gj = ib[-1] + ib[:-1]
                    else: gj = None
                    if names.count(gj) == 0:
                        if isdigit(ib[0]): gj = ib[1:] + ib[0]
                        else: gj = ib # added on 22.2.01
                    logi[-1].append([ib,gj])
                    try: # accept
                        index = names.index(gj)
                        logi[-1].append(gj)
                    except: pass # reject
                if index: break
            if index is None:
                print 'atom names do not match'
                logfile.append(('E',logi))
                return None,logfile
            index = index + first
            ia[1].setIndex(index-1)
            logi.append([ia[1],ia[1].index])
        if len(resseq[i][1]) != anum:
            logfile.append(('E',logi))
            universe = None
            break
    return universe,logfile

