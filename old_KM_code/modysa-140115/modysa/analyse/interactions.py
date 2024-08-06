# -*- coding: utf-8 -*-

from modysa.core import pbcgrid as P, geometry as G
import modysa.misc as M
import modysa.statistics as S
from modysa.io.index import IndexFile

import time, itertools, operator
import numpy as N


### zastosowanie: np. do_PmB6LipA66



class Interactions:

    def toFile(self):
        pass

    def fromFile(self):
        pass

    def number(self):
        return sum(map(len,self.data.values()))

    def reset(self):
        self.data = {}


class HydrogenBonds(Interactions):

    def __init__(self, donor_atoms, acceptor_atoms, h_counts, h_atoms, 
                 molecules, pairs, exclusions,
                 cutoff=0.325, angle=30., verbose=0):
        # angle is 35 in previous papers, 30 is more conservative/stronger HB
        self.cutoff = cutoff
        self.angle = angle
        self.verbose = verbose
        self.data = None
        
        # apply exclusions
        pairs = sorted(list(set(pairs) - set(exclusions)))
        
        ## collect all valid donors/acceptors
        # 1. extact unique mol id from [ pairs ]
        um_id = set(M.flatten(pairs))
        # 2. extract columns from [ molecules ]
        donors, acceptors = zip(*molecules) # vs. map(set, zip(*m))
        # 3. make lists of unique donors and acceptors
        # the new atom index corresponds to the list index
        # the value of the list element is the original trajectory atom index
        ###self.donor_atoms = M.flat([da[donors[i-1]-1] for i in um_id \
        ###                               if donors[i-1] > -1])
        # atom indices in NDX are 1..n, here we need array indices 0..(n-1)
        self.donor_atoms = \
            M.flatten([[idx-1 for idx in donor_atoms[donors[i-1]-1]] \
                           for i in um_id if donors[i-1] > -1])
        ###self.acceptor_atoms = M.flat([aa[acceptors[i-1]-1] for i in um_id \
        ###                                  if acceptors[i-1] > -1])
        # atom indices in NDX are 1..n, here we need array indices 0..(n-1)
        self.acceptor_atoms = \
            M.flatten([[idx-1 for idx in acceptor_atoms[acceptors[i-1]-1]] \
                           for i in um_id if acceptors[i-1] > -1])
        n_donors = len(self.donor_atoms)
        # 4. mapping for donors from old-to-new
        # since donor might be also acceptors, we need two maps not just one
        self.donor_map = dict(map(reversed, enumerate(self.donor_atoms)))
        self.acceptor_map = dict(map(reversed, 
                                     enumerate(self.acceptor_atoms, n_donors)))
        # let i atom be a donor and acceptor at the same time, then
        # self.donor_atoms[self.donor_map[i]] == \
        #    self.acceptor_atoms[self.acceptor_map[i]-len(self.donor_atoms)]

        ## map new donor atom indices to bound H atoms
        ###hatoms = M.flat([ha[donors[i-1]-1] for i in um_id])
        # atom indices in NDX are 1..n, here we need array indices 0..(n-1)
        hatoms = M.flatten([[idx-1 for idx in h_atoms[donors[i-1]-1]] \
                                   for i in um_id])
        hcounts = M.flatten([h_counts[donors[i-1]-1] for i in um_id])
        self.H_atoms = {}
        k = 0
        for i in range(len(self.donor_atoms)):
            self.H_atoms[i] = hatoms[k:k+hcounts[i]]
            k += hcounts[i]

        ## collecting decoding maps
        # 1. collect (donor atoms, acceptor atoms) for each molecule
        self.molecules = {}
        for i in um_id:
            self.molecules[i] = ([idx-1 for idx in donor_atoms[donors[i-1]-1]], 
                           [idx-1 for idx in acceptor_atoms[acceptors[i-1]-1]])
        # 2. map orig atom indices to molecule indices
        self.atoms = {}
        for m_id in self.molecules.keys():
            # unique donor and acceptor atoms of a molecule
            atoms = set(M.flatten(self.molecules[m_id])) 
            self.atoms.update(dict(zip(atoms,[m_id]*len(atoms))))
     
        # 3. collect allowed pairs as (lower, higher) entries
        self.pairs = {}
        for i in range(len(pairs)): self.pairs[self.__sorted_key(*pairs[i])] = i
     
    @staticmethod
    def __sorted_key(i, j):
        return i > j and (j, i) or (i, j)

    def __call__(self, conf): # , cell=None):
        # FIXME: __call__ more universal?
        # FIXME: cell is None: handle HB finding with no PBC
        cell = conf.cell.vectors
        configuration = conf.array
        if cell is None: raise ' PBC only'
        grid = P.gridsearch()
        ## extract valid atom trajectories
        atoms = self.donor_atoms + self.acceptor_atoms
        #frame = configuration[atoms] # in Numeric: N.take(configuration, atoms)
        frame = N.take(configuration, atoms, axis=0)
        grid.buildgrid(cell, frame, self.cutoff)
        ## find potential (distance criterion only) D-A pairs
        n_donors = len(self.donor_atoms)
        hb_donors = []
        hb_acceptors = []
        hb_acceptors_shift = []
        for i in range(n_donors):
            neighbors = grid.neighbors(i)
            if neighbors is None: continue
            n, s, d = neighbors
            # from list of atoms, get acceptors only;
            # since some donors might be on 'acceptor side', check
            # distances smaller than 1.0 \AA (0.1 nm) since they might indicate
            # the pair of the same donor atoms (actually, 
            # condition "> 0" would be enough)
            # N.nonzero: we only take the 1st dimension (no other present)
            acc_idx = N.nonzero(N.logical_and(N.greater_equal(n, n_donors),
                                              N.greater(d, 0.1)))[0]
            n_acceptors = len(acc_idx)
            if n_acceptors == 0: continue
            hb_donors.extend([i]* n_acceptors)
            hb_acceptors.extend(n[acc_idx])
            hb_acceptors_shift.extend(s[acc_idx])
        ## compute D-A vectors 
        xyzD = frame[hb_donors]
        xyzA = frame[hb_acceptors]
        xyzA_shifted = xyzA + N.inner(hb_acceptors_shift, cell)
        vDA = xyzD - xyzA_shifted
        # distances: N.sqrt(N.add.reduce(vDA*vDA,-1))
        ## compute D-H vectors
        hb_Hatoms = [self.H_atoms[i] for i in hb_donors]
        hb_Hcounts = M.flatten([[i]*len(hb_Hatoms[i]) \
                                    for i in range(len(hb_Hatoms))])
        #xyzH = configuration[M.flatten(hb_Hatoms)]
        xyzH = N.take(configuration, M.flatten(hb_Hatoms), axis=0)
        # duplicated donor entries (each donor might have diff. no of H atoms)
        #xyzD = xyzD[hb_Hcounts]
        xyzD = N.take(xyzD, hb_Hcounts, axis=0)
        vDH = xyzD - xyzH
        # length of vDH shouldn't be larger then, say, 1.2 \AA, i.e. D-H bond
        # since we assume that H atoms are together with donor atom
        # in charge group and as such don't require applying PBC
        ## compute angles between D-H and D-A
        #angles = G.angle(vDA[hb_Hcounts], vDH)/M.Units.deg
        angles = G.angle(N.take(vDA, hb_Hcounts, axis=0), vDH)/M.Units.deg
        hb = N.flatnonzero(N.less_equal(angles, self.angle))
        #di = N.array(hb_donors)[hb_Hcounts]
        di = N.take(hb_donors, hb_Hcounts)
        #ai = N.array(hb_acceptors)[hb_Hcounts]
        ai = N.take(hb_acceptors, hb_Hcounts)
        # return list of (donor_atom, acceptor_atom) indices
        #
        # recovering original atom numbering, i.e. 1..n
        # from array indices, i.e. 0..(n-1)
        # FIXME: HB(15, 26) found meaning intramolecular HB
        #        once HBs have been identified, filter them according
        #        to pairs group
        # the fact that intramolecular HBs are identified is a feature
        # not a bug because there is one list of donors and 
        # one list of acceptors
        ## wstępne rozwiązanie: 1) indeksy atomów na indeksy cząsteczek
        ## 2) wypełnij prostokątną macierz n x m (n cząsteczek oddziałujących
        ## z m cząsteczkami ???) - normalnie 0, tam gdzie są pairs wstaw 1 
        ## (masked array), 3) odzyskaj indeksy ???
        ## idx = [(2,4), (0,3), (1,1)]
        ## mask = N.zeros((3,5))
        ## x, y = zip(*idx)
        ## mask[x,y] = 1
        ## # masking.. N.ma
        ## which_ok = zip(*N.nonzero(x))
        ## http://www.scipy.org/Cookbook/Indexing
        ## x = N.ma.zeros((2000,2000),N.bool) # 4 MB
        self.data = []
        
        if len(di[hb]) > 1:
            da = zip(N.array(self.donor_atoms)[di[hb]], 
                     N.array(self.acceptor_atoms)[ai[hb]-n_donors])
            for d, a in da:
                # get molecule indices for donor/acceptor atoms
                allowed = self.__sorted_key(self.atoms[d], self.atoms[a])
                # recovering 'plus one' original numbering
                if self.pairs.has_key(allowed): self.data.append((d+1,a+1))
        elif len(di[hb]) == 1: 
            d, a = self.donor_atoms[di[hb]], \
                self.acceptor_atoms[ai[hb]-n_donors]
            
            allowed = self.__sorted_key(self.atoms[d], self.atoms[a])
            if self.pairs.has_key(allowed): self.data = [(d+1, a+1)]


# WaterBridges liczy WB dla przypadku setA - atomy jednego rodzaju
# czasteczek (np. same lipidy-A, POPE, PmB) oraz setB zawiera atomy
# wody (o tym ktory jest ktory decyduje zmienna `water`: 0 lub 1

# WaterBridges2 liczy WB między dwoma setami *różnych* cząsteczek,
# woda definiowana jest oddzielnie setW

# patrz: do_PmB6LipA66 (zastosowania WaterBridges2)

class WaterBridges(Interactions):

    def __init__(self, setA, setB, water=1, verbose=0):
        self.verbose = verbose
        self.nmolA, self.natA, self.subsetA, self.hA, self.offsetA = setA
        self.nmolB, self.natB, self.subsetB, self.hB, self.offsetB = setB
        assert water in [0, 1]
        self.water = water # FIXME: one of two sets contain water molecules
        self.reset()

    def reset(self):
        self.data = {}
        self.watbri = {}

    def find(self, hbonds):

        # fetch those from hbonds instance?
        w_offset = [self.offsetA, self.offsetB][self.water]
        w_nmol = [self.nmolA, self.nmolB][self.water]
        w_nat = [self.natA, self.natB][self.water]
        w_min, w_max = (w_offset, w_offset + w_nmol*w_nat) # FIXME: we assume here, water molecules are in one chunk

        other = [1, 0][self.water]
        o_offset = [self.offsetA, self.offsetB][other]
        o_nmol = [self.nmolA, self.nmolB][other]
        o_nat = [self.natA, self.natB][other]

        wb = {}
        for atA in hbonds.data.keys():
            for a, b in hbonds.data[atA]:
                if w_min <= a <= w_max:
                    if wb.has_key(a): wb[a].append(b)
                    else: wb[a] = [b]
                else:
                    if wb.has_key(b): wb[b].append(a)
                    else: wb[b] = [a]
        self.watbri = dict([(w, wb[w]) for w in wb.keys() if len(wb[w]) > 1])

        inter = []
        intra = []
        for w, at in self.watbri.items():
            # zip(*comb)
            wb = N.array(list(itertools.combinations(at, 2)))
            imol = ((wb-o_offset)/o_nat).astype(N.Int)
            diff = N.subtract.reduce(imol, -1)
            inter.extend([(a<b and a or b, w, a<b and b or a) for a,b \
                              in N.take(wb,N.nonzero(diff))])
            intra.extend([(a<b and a or b, w, a<b and b or a) for a,b in \
                              N.take(wb,N.nonzero(N.equal(diff,0)))])
        self.data = {'inter': inter, 'intra': intra, 
                     'water': self.watbri.keys()}
        
    def summary(self):
        if not self.watbri: return
        gj = map(len, self.watbri.values())
        max_order = max(gj)
        data = {}
        for i in range(2, max_order+1):
            data[i] = gj.count(i)
        return data

    def query(self, multiplicity):
        """ return water bridges with given multiplicity """
        return dict([(w, self.watbri[w]) for w in self.watbri.keys() \
                         if len(self.watbri[w]) == multiplicity])
        
    def multiple(self, which='inter'):
        if not self.watbri: return
        assert which in ['inter','intra']
        unique = {}
        for wbid in self.data[which]:
            a, w, b = wbid
            if unique.has_key((a,b)): unique[(a,b)].append(w)
            else: unique[(a,b)] = [w]
        # pairs of lipid atoms with
        # the number of bridging water molecules > 1
        return dict([(wbid, unique[wbid]) for wbid in unique.keys() \
                         if len(unique[wbid])>1])


class WaterBridges2(Interactions):

    def __init__(self, setA, setW, setB=None, verbose=0):
        self.atomsA, self.molA = setA # (first_atom, last_atom) we have to reorder atoms if necessery (pbcgrid)
        self.atomsW, self.molW = setW
        assert setA != setB
        if setB: self.atomsB, self.molB = setB
        self.verbose = verbose
        self.reset()

    def find(self, hbA, hbB = None):
        # hbA and hbB are dictionaries: HydrogenBonds.data

        assert hbA != hbB
        def process_HBonds(hbonds):
            wb = {}
            for atA in hbonds.keys():
                for a, b in hbonds[atA]:
                    aw = self.atomsW[0] <= a <= self.atomsW[1]
                    bw = self.atomsW[0] <= b <= self.atomsW[1]
                    assert not (aw and bw) # not both water molecules allowed
                    if aw:
                        if wb.has_key(a): wb[a].append(b)
                        else: wb[a] = [b]
                    elif bw:
                        if wb.has_key(b): wb[b].append(a)
                        else: wb[b] = [a]
            return wb

        wb = M.Results()
        wb.update(process_HBonds(hbA).items())
        if hbB: wb.update(process_HBonds(hbB).items())
        # filter out non-water bridge interactions
        # self.watbri contains all WB (including setB-water-setB, 
        #      setA-water-setA); we need setA-water-setB which are
        #      identified below and stored in self.data
        self.watbri = dict([(w, wb[w]) for w in wb.keys() if len(wb[w]) > 1])

        def validAtoms(a,b):
            aa = self.atomsA[0] <= a <= self.atomsA[1]
            ab = self.atomsB[0] <= a <= self.atomsB[1]
            ba = self.atomsA[0] <= b <= self.atomsA[1]
            bb = self.atomsB[0] <= b <= self.atomsB[1]
            if (aa and bb) or (ab and ba): return True

        self.data = []
        for w, at in self.watbri.items():
            wb = list(itertools.combinations(at, 2))
            # only allowed a from A and b from B or a from B and b from A
            if hbB: check_func = validAtoms
            else: check_func = lambda a,b: True
            wb = [(a<b and a or b,w, a<b and b or a) \
                      for a,b in wb if check_func(a,b)]
            self.data.extend(wb)


    def sort(self):
        sorted = {}
        for a,w,b in self.data:
            ia = a%self.molA.natoms
            ib = a%self.molB.natoms
            if sorted.has_key(ia): sorted[ia].append((a,w,b))
            else: sorted[ia] = [(a,w,b)]
        return sorted

    def reset(self):
        self.data = {}
        self.watbri = {}
        
    def number(self):
        return len(self.data)

    def summary(self):
        if not self.watbri: return
        gj = map(len, self.watbri.values())
        max_order = max(gj)
        data = {}
        for i in range(2, max_order+1):
            data[i] = gj.count(i)
        return data

    def query(self, multiplicity):
        """ return water bridges with given multiplicity """
        return dict([(w, self.watbri[w]) for w in self.watbri.keys() \
                         if len(self.watbri[w]) == multiplicity])
        
    def multiple(self):
        if not self.watbri: return
        unique = {}
        for wbid in self.data:
            a, w, b = wbid
            if unique.has_key((a,b)): unique[(a,b)].append(w)
            else: unique[(a,b)] = [w]
        # atom pairs with the number of bridging water molecules > 1
        return dict([(wbid, unique[wbid]) for wbid in unique.keys() \
                         if len(unique[wbid])>1])



class RadialDistributionFunction:

    def __init__(self, ndx_filename, r_max=2.0, width=.005):
        # FIXME: norm=(1,1,1) could be a vector which will allow
        #        to calculate 1-, 2-, or 3-dimensional RDFs, e.g.
        #        if one wants to calculate RDF in XY, norm reads (1,1,0)
        # REMARK: significance testing might be done with bootstrap test,
        #         e.g. integral as a difference between two RDF's
        #         calculated from generated bootstrap samples
        # lenghts/distances in nm!
        self.bins = N.arange(0, r_max+width, width)
        self.reset()
        # ** let's NDX **
        # atoms, moieties, pairs
        # FIXME: exclusions (moiety pairs), omissions (atom pairs)
        # omissions might be useful, if we're interested in intramolecular
        # RDF yet we'd like to skip certain atom pairs
        ndx = IndexFile()
        ndx.read(ndx_filename)
        # FIXME: *exclusions 1-72 1-72, intermolecular only
        assert ndx.valid('atoms', 'moieties', 'pairs')
        ndx.flat('atoms','moieties')
        ndx.range('pairs')
        # ** fetch atom indices **
        # atom indices in NDX are 1..n, here we need array indices 0..(n-1)
        # FIXME: do we really need this?
        self.atoms = [i-1 for i in ndx.groups['atoms']]
        # --- main index data processing starts ---
        # The heuristic strategy of speeding up calculations
        # will fail if moieties in the IndexFile are not given
        # sequentially since this leeds to presumably large number
        # of atom intervals (see comments below for explanation);
        # the remedy for this is to pre-sort moieties (i.e. we
        # should avoid situations lipids, ion, water, lipids, ion, water, etc.
        # instead we should have lipids, ion, water - i.e. continuous
        # sets of moieties of a given kind
        data = M.flatten(ndx.groups['pairs'])
        # 1. find which groups of atoms form the largest set of pairs
        # i.e. each atom vs paired atoms, the length of 'paierd atoms' list
        # is a criterion here; note that for sake of efficiency we are
        # dealing here with moieties *not* atoms
        partners = {}
        for i,j in ndx.groups['pairs']:
            if partners.has_key(i): partners[i].append(j)
            else: partners[i] = [j]
            if i == j: continue # intra-moiety pairs stored ONCE
            if partners.has_key(j): partners[j].append(i)
            else: partners[j] = [i]
        atom_count = {}
        for i,j in partners.items():
            atom_count[i] = sum([ndx.groups['moieties'][a-1] for a in j])
        atom_count = sorted(atom_count.items(), key=operator.itemgetter(1))
        # atom_count is now a list of (moiety_id, number_of_moieties)
        npairs = len(ndx.groups['pairs'])
        result = {}
        i, j = atom_count.pop()
        result[i] = sorted(partners[i])
        done = len(partners[i])
        stored = [i]
        # 2. collect unique pairs (they're cloned in partners due to symmetry
        # i.e. (a,b) is stored as (a,b) and (b,a); 
        # beware of the fact, that intra-moiety pairs (a,a) are stored once
        # see the comment above
        while done < npairs:
            i, j = atom_count.pop()
            # FIXME: is sorted needed here?
            p = sorted(set(partners[i]).difference(stored))
            done += len(p)
            stored.append(i)
            result[i] = p
        # 3. store final data in a form of dictionary
        # (atom_first, atom_last-but-one): list of intervals of atom indices
        natoms = [0]
        natoms.extend(N.add.accumulate(ndx.groups['moieties']))
        final = {}
        for k,v in result.items(): # (ge,lt)
            final[(natoms[k-1],natoms[k])] = \
                [(natoms[i[0]-1],natoms[i[-1]]) for i in M.slicer(v)]
        # TODO: omissions? split atom ranges
        # TODO: check for recurring values, and merge keys 
        # (i.e. ranges of atoms)
        # actually, this should be hardly be an improvement since we are
        # looping over atoms; for sloppy defs of groups it might, however,
        # do the trick
        self.pairs = final
        # --- main index data processing ends ---

    def reset(self):
        self.volume = []
        self.npairs = 0
        self.data = N.zeros(len(self.bins)-1)

    def __call__(self, configuration, cell=None):
        # FIXME: non-periodic systems/surface atoms in proteins
        assert cell is not None 
        grid = P.gridsearch()
        frame = configuration[self.atoms]
        r_max = self.bins[-1]
        grid.buildgrid(cell, frame, r_max)
        distances = []
        for atom_range, intervals in self.pairs.items():
            for ai in range(*atom_range):
                neighbors = grid.neighbors(ai)
                if neighbors is None: continue
                n, s, d = neighbors
                distances.extend(d[M.pick_items(n,intervals)])
                
        self.npairs += len(distances)
        self.volume.append( N.dot( N.cross(cell[0],cell[1]), cell[2]) )
        ##self.volume.append(cell.volume) # FIXME: cell is an instance not array!
        counts, bins = N.histogram(distances, bins=self.bins, new=True)
        self.data += counts

    def get(self):
        # returns bins, rdf; integral? (i.e. accum. sum of pairs vs distance)
        factor = 4./3 * N.pi * self.npairs * \
            (self.bins[1:]**3 - self.bins[:-1]**3) / S.mean(self.volume)
        rdf = self.data / factor
        return self.bins, rdf

    
class AtomContacts:
    # see demo.count_ions

    def __init__(self, atoms, moieties, pairs, 
                 r_cutoff=0.31, register=False, verbose=0):

        # r_cutoff: 0.31 for ion (Na+)/oxygen contacts
        #           0.75 for carbon/carbon contacts (nonpolar interactions)
        self.register = register
        self.r_cutoff = r_cutoff
        self.verbose = verbose

        self.atoms = [i-1 for i in atoms]

        ## see RadialDistributionFunction for other comments

        # transform list of pairs into lexicon mapping 
        # list of atoms (value) in contact with a given atom (key)
        partners = {}
        for i,j in pairs:
            if partners.has_key(i): partners[i].append(j)
            else: partners[i] = [j]
            if i == j: continue # intra-moiety pairs stored ONCE
            if partners.has_key(j): partners[j].append(i)
            else: partners[j] = [i]
            
        # transform partners lexicon into result lexicon
        # with keys: (moiety1_atom1, moiety1_atomN) and values
        # being list of atom ranges of other moieties;
        # NOTE performance tweak: keys are for smallest moieties
        # values are for potentially larger atom ranges
        items = sorted([(len(v),k) for k,v in partners.items()])
        # (smaller_index, larger_index)
        pairs = set([i > j and (j,i) or (i,j) for i, j in pairs])
        result = {}
        stored = []
        while pairs:
            count, k = items.pop()
            hits = [i > j and (j,i) or (i,j) \
                        for i in [k] for j in partners[k]]
            # c.f. RDF, counting pairs is probably faster
            pairs = pairs.difference(hits) 
            result[k] = sorted(set(partners[k]).difference(stored))
            stored.append(k)

        ## self.pairs: from [pairs] 1-72 73-75 (which corresponds to atom 
        ## ranges <0,1943> and <1944,2033>) to
        ## {(1944, 1974): [(0, 1944)], (2004, 2034): [(0, 1944)], 
        ##  (1974, 2004): [(0, 1944)]}

        natoms = [0]
        natoms.extend(N.add.accumulate(moieties))
        self.pairs = {}
        for k,v in result.items():
            self.pairs[(natoms[k-1],natoms[k])] = \
                [(natoms[i[0]-1],natoms[i[-1]]) for i in M.slicer(v)]        
        #
        self.data = []
        self.contacts = []

    def __call__(self, configuration):

        # TODO: implement no PBC case
        assert configuration.cell is not None, 'PBC only!'

        grid = P.gridsearch()
        frame = N.take(configuration.array,self.atoms,axis=0)
        grid.buildgrid(configuration.cell.vectors, frame, self.r_cutoff)
        # 
        self.data = []
        self.contacts = []
        for atom_range, intervals in self.pairs.items():
            for ai in range(*atom_range):
                neighbors = grid.neighbors(ai)
                if neighbors is None: continue
                n, s, d =  neighbors
                # filter out larger distances, no need? see neighbor.c: if(d < g->rcut)
                #if N.any(N.where(d>.28)[0]):
                #    print n, s, d
                #    raise
                #n = N.take(n,N.where(d<=self.r_cutoff)[0]) # N.where -> tuple
                # proceed..
                if self.verbose:
                    print ai, n, intervals
                    print self.atoms[ai], N.take(self.atoms,n)
                    print M.pick_items(n, intervals)
                    print N.take(n, M.pick_items(n, intervals), -1)
                    print '-'*20
                contacts = N.take(n, M.pick_items(n, intervals), -1)
                # the number of neighbors of given type
                num_n = contacts.shape[0] 
                if num_n == 0: continue
                else: self.data.append(num_n)
                if self.register:
                    # self.atoms contains atoms numbered 0..n-1
                    self.contacts.append(\
                        (self.atoms[ai]+1, 
                         (N.take(self.atoms, contacts)+1).tolist()))
                    # FIXME: how about counting unique atoms 
                    #        from setB (i.e. various atoms from setA 
                    #        might see the same atoms from setB)

    def count(self):
        """return the number of atom contacts found"""
        return sum(self.data)

    def unique(self, count=True):
        """return either the number or pairs (lexicon) of atoms in contact
        in the latter case, keys are setB atoms and values setA atoms"""
        assert self.register
        if count:
            result = len(set(M.flatten([x[1] for x in self.contacts])))
        else:
            result = {}
            for a, blist in self.contacts:
                for b in blist:
                    if result.has_key(b): result[b].append(a)
                    else: result[b] = [a]
        return result

    def occurrence(self):
        # histogram: (1,6), (2,4), (3,1) <- one ion bound to 6 setB atoms
        # two ions bound to 4 setB etoms etc.
        ##return [(i[0], len(list(i[1]))) \
        ##            for i in itertools.groupby(sorted(self.data))]
        return M.occurence(self.data)

