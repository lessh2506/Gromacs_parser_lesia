import sys

class Atom:
    def __init__(self, residue_number, residue_name, atom_name, atom_number, x, y, z, v_x, v_y, v_z, structure, start_terminus=None, end_terminus=None, chain_id=None, PDB_record=None, occupancy=None, temp_factor=None):
        self.residue_number = residue_number
        # i don't understand what is self.residue_name_ff
        self.residue_name = residue_name
        self.residue_name_ff =residue_name
        self.atom_name = atom_name
        self.atom_number = atom_number
        self.chain_id = chain_id
        self.x = x  
        self.y = y
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        self.atoms_record=None
        self.start_terminus=start_terminus
        self.end_terminus=end_terminus
        self.structure=structure
        self.PDB_record=PDB_record
        self.occupancy=occupancy
        self.temp_factor=temp_factor

        def __repr__(self):
            return f"Chain ID: {self.chain_id}, Residue Name: {self.residue_name}, Residue Number: {self.residue_number}, Atom Name: {self.atom_name}"


class Atoms:
    def __init__(self, atoms_list, structure):
        self.atoms_list=atoms_list
        self.structure=structure



        @property
        def atom_length(self):
            return len(self.atoms_list)
        
        def __iter__(self):
            return iter(self.atoms_list)
        
        def pdb_text(self):
            pdb_lines = []
            for atom in self.atoms_list:
                pdb_lines.append(f"ATOM  {atom.atom_number:5d} {atom.atom_name:>4s} {atom.residue_name:>3s} {atom.chain_id} {atom.residue_number:4d}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy if atom.occupancy else 1.00:6.2f}{atom.temp_factor if atom.temp_factor else 0.00:6.2f}")
            return "\n".join(pdb_lines)
        
        def write(self, file_path):
            self.structure.write(file_path, self)

        def show(self, cartoon=True, sticks=False, highlight_atoms=None): 
            import nglview as nv
            # Create a PDB format string from the atoms data
            pdb_content=self.structure.pdb_text(self)

            # Create an NGLWidget and load the PDB structure
            view = nv.NGLWidget()
            view.add_component(pdb_content, ext='pdb')
            view.clear_representations()
            if cartoon:
                view.add_cartoon()
            if sticks:
                view.add_ball_and_stick()
            if highlight_atoms:
                view.add_ball_and_stick("sphere", selection=highlight_atoms, color="red")


            # Display the widget
            return view
        
        def __filter_by_residue(self, residue_name):
            return Atoms([atom for atom in self.atoms_list if atom.residue_name==residue_name], self.structure)
        

class StructureRead:
    def __init__(self, file_path="", atoms=None):
        super().__init__()
        if file_path != "":
            self.file_path = file_path
            self.atoms = self.read()
        else: ValueError("Please provide a file path")	


        def pdb_reader(self, file_path):
            with open(file_path, 'r') as file:
                lines = file.readlines()
                atoms = [] 
                for line in lines:
                    record_name = line[0:6].strip()

                #coordinates ATOM and HETATM
                    if record_name == 'ATOM' or record_name == "HETATM":
                        atom_number = int(line[6:11].strip())
                        atom_name = line[12:16]
                        alt_loc = line[16].strip()
                        residue_name = line[17:20].strip()
                        chain_id = line[21].strip()
                        residue_number = int(line[22:26].strip())
                        insertion_code = line[26].strip()
                        #devide by 10 to get nm
                        x = float(line[30:38].strip())/10
                        y = float(line[38:46].strip())/10
                        z = float(line[46:54].strip())/10
                        occupancy = float(line[54:60].strip()) if line[54:60].strip() != '' else None
                        temp_factor = float(line[60:66].strip()) if line[60:66].strip() != '' else None
                        segment_id = line[72:76].strip()
                        element_symbol = line[76:78].strip()
                        charge = line[78:80].strip()
                        #skip alt loc, gromacs does this too as gro file format can not handle alt loc
                        if not (alt_loc == '' or alt_loc == 'A'):
                            continue
                        atom = Atom(residue_number, residue_name, atom_name, atom_number, x,y,z,None,None,None, self,chain_id=chain_id, PDB_record=record_name, occupancy=occupancy, temp_factor=temp_factor)
                        atoms.append(atom)
