import string, os, re

class DSSP:

    def __init__(self,filename,**options):

        self.__cmd = 'dssp --'
        self.__result = ''
        self.__logfile = ''
        self.__summary = None
        self.__secst = None
        self.__chains = []
        f = open(filename,'r')
        self.__data = f.read()
        f.close()

    def run(self):

        stdin, stdout, stderr = os.popen3('%s' % self.__cmd)
        stdin.write(self.__data)
        stdin.close()
        self.__result = stdout.readlines()
        stdout.close()
        self.__logfile = stderr.read()
        stderr.close()

    def save(self,filename):

        if not self.__result: self.run()
        f = open(filename,'w')
        f.write(''.join(self.__result))
        f.close()

    __pat = [('num_res',"TOTAL NUMBER OF RESIDUES"),
             ('num_chains', "NUMBER OF CHAINS"),
             ('access', "ACCESSIBLE SURFACE OF PROTEIN"),
             ('hb_total',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)"),
             ('hb_par',
              "TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES"),
             ('hb_anti',
              "TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES"),
             ('hb_im5',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5)"),
             ('hb_im4',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4)"),
             ('hb_im3',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-3)"),
             ('hb_im2',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-2)"),
             ('hb_im1',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-1)"),
             ('hb_i0' ,
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+0)"),
             ('hb_ip1',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+1)"),
             ('hb_ip2',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2)"),
             ('hb_ip3',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3)"),
             ('hb_ip4',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4)"),
             ('hb_ip5',
              "TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5)"),
             ('res_aa_struct', "RESIDUE AA STRUCTURE")]

    __struc = {'bridge': "B",
               'extended': "E",
               'alpha-helix': "H",
               '3_10-helix': "G",
               'pi-helix': "I",
               'hturn': "T",
               'bend': "S"}


    def parseDsspOutput(self):

        if not self.__result: self.run()
        lines = self.__result
        result = {}
        lit = 0
        while(string.count(lines[lit],self.__pat[-1][1]) == 0):
            for id, name in self.__pat:
                if string.count(lines[lit],name) > 0:
                    result[id] = string.atof(string.split(lines[lit])[0])
                    break
            lit = lit + 1
        start = lit + 1

        secst = {}
        for i in range(start,len(lines)):
            chain = string.strip(lines[i][11:12])
            if not chain: continue
            resno = string.atoi(lines[i][:5])
            pdbresno = string.atoi(lines[i][5:10])
            amino = lines[i][13:14]
            struct = lines[i][16:17]
            access = string.atoi(lines[i][34:38])
            phi = string.atof(lines[i][103:109])
            psi = string.atof(lines[i][109:115])
            #
            if not secst.has_key(chain):
                secst[chain] = []
                self.__chains.append(chain)
            secst[chain].append((resno,pdbresno,amino,
                                 struct,access,phi,psi))
        self.__secst = secst
        self.__summary = result

    def secondaryStructure(self):
        if self.__secst is None: self.parseDsspOutput()
        return self.__secst

    def summary(self):
        if self.__summary is None: self.parseDsspOutput()
        return self.__summary

    def secondaryStructureSummary(self):

        if self.__secst is None: self.parseDsspOutput()
        self.__summary['num_chains'] = len(self.__secst)
        res = {}
        for chain in self.__chains:
            res[chain] = {}
            for resno,pdbresno,amino,struct,access,phi,psi \
                    in self.__secst[chain]:
                if res[chain].has_key(struct):
                    res[chain][struct] = res[chain][struct] + 1
                else: res[chain][struct] = 1
        return res

    def printSummary(self):

        if self.__summary is None: self.parseDsspOutput()
        for id, name in self.__pat:
            if not self.__summary.has_key(id): continue
            print '%-60s: %6d' % (name,self.__summary[id])
        
    def secondaryStructureRanges(self,skip_isolated=1):

        if self.__secst is None: self.parseDsspOutput()
        result = {}
        for chain in self.__chains:
            out = self.__secst[chain]
            sec = ''.join([x[3] for x in out])
            result[chain] = []
            for s in self.__struc.values():
                gj = re.findall('%s+' % s,sec)
                if not gj: continue
                end = 0
                for x in gj:
                    start = sec.index(x,end)
                    end = start+len(x)
                    if skip_isolated and start == end - 1:
                        continue
                    result[chain].append((start,end,s))
            # coils
            result[chain].sort()
        return result

def dssp_to_molscript(pdb_filename, mapping = {}):
    """ runs DSSP for a given protein PDB and produces
    valid molscript input describing protein secondary structure """

    default = {'H': 'helix', 'G': 'helix', 'I': 'helix',
               'E': 'strand', 'T': 'turn'}
    for k in default.keys():
        if mapping.has_key(k): continue
        mapping[k] = default[k]

    dssp = DSSP(pdb_filename)
    secranges = dssp.secondaryStructureRanges()
    full = dssp.secondaryStructure()
    chains = full.keys()
    result = []
    for chain in chains:
        out = [x for x in secranges[chain] \
               if mapping.has_key(x[-1])]
        previous = 0
        coils = []
        insert = len(result)
        for b, e, s in out:
            if b > previous:
                coils.append('%s from %c%d to %c%d;' % \
                             ('coil',
                              chain,full[chain][previous][1],
                              chain,full[chain][b][1]))
                previous = e
            key = mapping[s]
            result.append('%s from %c%d to %c%d;' % \
                          (key,
                           chain,full[chain][b][1],
                           chain,full[chain][e][1]))
            previous = e
        result.append('')
        result.insert(insert,'\n'.join(coils))
    return '\n'.join(result)


