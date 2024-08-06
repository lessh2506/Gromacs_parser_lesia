class ForceField:
    def __init__(self, data):
        parsed_data = self.text_to_ff(data)
        self.defaults = parsed_data['defaults']
        self.atomtypes = parsed_data['atomtypes']
        self.angletypes = parsed_data['angletypes']
        self.bondtypes = parsed_data['bondtypes']
        # self.constrainttypes = parsed_data['constrainttypes']
        self.dihedraltypes = parsed_data['dihedraltypes']
        self.pairtypes = parsed_data['pairtypes']

        
    def text_to_ff(self, data):
        defaults=[]     
        atomtypes=[]
        bondtypes=[]
        constrainttypes=[]
        angletypes=[]
        dihedraltypes=[]
        pairtypes=[]

        current_tag =""
        for line in data.split('\n'):
            #assign tag
            if "[" in line:
                current_tag=line.split('[')[1].split(']')[0].split()[0]
                continue
            elif len(line.split())>0:

                if current_tag=='defaults':
                    #check if line is correct
                    els=line.split()
                    if len(els)!=5:
                        assert print("incorect number of parameters")
                    defaults.append(Defaults( int(els[0]),int(els[1]), els[2],float(els[3]),float(els[4])))
                elif current_tag=='atomtypes':
                    els=line.split()
                    if len(els)!=8:
                        assert print("incorect number of parameters")
                    atomtypes.append(Atomtypes(els[0],els[1],int(els[2]),float(els[3]),float(els[4]),els[5],float(els[6]),float(els[7])))
                elif current_tag=='bondtypes':
                    els=line.split()
                    function_type=int(els[2])
                    if function_type == 1:
                        if len(els)!=5:
                            assert print("incorect number of parameters")
                        bondtypes.append(Bondtypes(els[0],els[1],int(els[2]),float(els[3]),float(els[4])))
                # elif current_tag=='constrainttypes':
                #     els=line.split()
                #     if len(els)!=4:
                #         assert print("incorect number of parameters")
                #     constrainttypes.append(Constrainttypes(els[0],els[1],int(els[2]),float(els[3])))
                
                elif current_tag=="angletypes":
                    els=line.split()
                    function_type=int(els[3])
                    if function_type == 5:
                        if len(els) != 8:
                            assert print("incorect number of parameters")
                        angletypes.append(Angletypes(els[0],els[1], els[2],int(els[3]),th0=float(els[4]),k0=float(els[5]),r13=float(els[6]),Kub=float(els[7])))
                    if function_type == 1:
                        if len(els) != 6:
                            assert print("incorect number of parameters")
                        angletypes.append(Angletypes(els[0],els[1], els[2],int(els[3]),th0=float(els[4]),k0=float(els[5])))
                elif current_tag=="dihedraltypes":
                    els=line.split()
                    function_type=int(els[4])
                    if function_type == 9:
                        if len(els) != 8:
                            assert print("incorect number of parameters")
                        dihedraltypes.append(Dihedraltypes(els[0],els[1],els[2],els[3],int(els[4]),float(els[5]),float(els[6]),int(els[7])))
                    if function_type == 2:
                        if len(els) != 7:
                            assert print("incorect number of parameters")
                        dihedraltypes.append(Dihedraltypes(els[0],els[1],els[2],els[3],int(els[4]),float(els[5]),float(els[6])))
                    if function_type == 3:
                        if len(els) != 11:
                            assert print("incorect number of parameters")
                        dihedraltypes.append(Dihedraltypes(els[0],els[1],els[2],els[3],int(els[4]),float(els[5]),float(els[6]),float(els[7]),float(els[8]),float(els[9]),float(els[10])))
                
                elif current_tag=="pairtypes":
                    els=line.split()
                    function_type=int(els[2])
                    if function_type == 1:
                        if len(els) != 5:
                            assert print("incorect number of parameters")
                        pairtypes.append(Pairtypes(els[0],els[1],int(els[2]),float(els[3]),float(els[4])))
                    if function_type == 2:
                        if len(els) != 7:
                            assert print("incorect number of parameters")
                        pairtypes.append(Pairtypes(els[0],els[1],int(els[2]),float(els[3]),float(els[4]),float(els[5]),float(els[6])))
                                
        return {'defaults':defaults}

#info to create classes
#https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#tab-topfile2

class Defaults:
    def __init__(self, nb_function_type, combination_rule, generate_pairs, fudge_LJ, fudge_QQ, ifdefs=[]):
        self.nb_function_type=nb_function_type
        self.combination_rule=combination_rule
        self.generate_pairs=generate_pairs
        self.fudge_LJ=fudge_LJ  
        self.fudge_QQ=fudge_QQ
        self.ifdefs=ifdefs

class Atomtypes:
    #bonded type and atomic number are optional
    def __init__(self, atom_type, bonded_type, atomic_number, m_u, q, particle_type, V, W, ifdefs=[]):
        self.atom_type=atom_type
        self.bonded_type=bonded_type
        self.atomic_number=atomic_number
        self.m_u=m_u
        self.q=q
        self.particle_type=particle_type
        self.V=V
        self.W=W
        self.ifdefs=ifdefs

class Bondtypes: 
    def __init__(self, bonded_type_i, bonded_type_j, function_type, b0="", kb="", D="", beta="", Ci="", bm="",table_number="",k="", low="", up="", kdr="", ifdefs=[]):
        self.bonded_type_i=bonded_type_i
        self.bonded_type_i=bonded_type_j
        self.ifdefs=ifdefs
        if function_type==1:
            #bond (harmonic)
            self.b0=b0
            self.kb=kb
        else:
            assert print(f"There is no code yet to handle this function type bondtypes:{function_type}")
        

            
class Constrainttypes:
    def __init__(self, bonded_type_i, bonded_type_j, function_type, b0=""):
        self.bonded_type_i=bonded_type_i
        self.bonded_type_j=bonded_type_j
        self.function_type=function_type
        self.b0=b0
        
        

class Angletypes:
    def __init__(self, bonded_type_i, bonded_type_j, bonded_type_k, function_type, th0="", k0="", r13="", Kub=""):
        self.bonded_type_i=bonded_type_i
        self.bonded_type_j=bonded_type_j
        self.bonded_type_k=bonded_type_k
        self.function_type=function_type
        if function_type==5:
            self.th0=th0
            self.k0=k0
            self.r13=r13
            self.Kub=Kub
        elif function_type==1:
            self.th0=th0
            self.k0=k0
        else:
            assert print(f"There is no code yet to handle this function type bondtypes:{function_type}")    
        

class Dihedraltypes:
    def __init__(self, bonded_type_i, bonded_type_j, bonded_type_k, bonded_type_l, function_type, phi0="", kp="", mult="", q0="", cq="", c0="", c1="", c2="",c3="", c4="",c5=""):
        self.bonded_type_i=bonded_type_i
        self.bonded_type_j=bonded_type_j
        self.bonded_type_k=bonded_type_k
        self.bonded_type_l=bonded_type_l
        self.function_type=function_type
        if function_type==9:
            self.phi0=phi0
            self.kp=kp
            self.mult=mult
        elif function_type==2:
            self.q0=q0
            self.cq=cq
        elif function_type==3:
            self.c0=c0
            self.c1=c1
            self.c2=c2
            self.c3=c3
            self.c4=c4
            self.c5=c5
        else:
            assert print(f"There is no code yet to handle this function type bondtypes:{function_type}")        
        

class Pairtypes:
    def __init__(self, bonded_type_i, bonded_type_j, function_type, fudgeQQ="", q1="", q2="", c6="", c12=""):
        self.bonded_type_i=bonded_type_i
        self.bonded_type_j=bonded_type_j
        self.function_type=function_type
        if function_type==1:
            self.c6=c6
            self.c12=c12
        if function_type==2:
            self.fudgeQQ=fudgeQQ
            self.q1=q1
            self.q2=q2
            self.c6=c6
            self.c12=c12
        else:
            assert print(f"There is no code yet to handle this function type bondtypes:{function_type}")                
        