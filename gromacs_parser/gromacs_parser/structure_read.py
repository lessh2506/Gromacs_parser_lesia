import os
class Atom:
    def __init__(self, residue_number, residue_name, atom_name, atom_number, x, y, z, v_x, v_y, v_z, structure, start_terminus=None, end_terminus=None, chain_id=None, PDB_record=None, occupancy=None, temp_factor=None):
        self.residue_number = residue_number
        self.residue_name = residue_name
        self.residue_name_ff = residue_name
        self.atom_name = atom_name
        self.atom_number = atom_number
        self.chain_id = chain_id
        self.x = x  
        self.y = y
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
        self.atoms_record = None
        self.start_terminus = start_terminus
        self.end_terminus = end_terminus
        self.structure = structure
        self.PDB_record = PDB_record
        self.occupancy = occupancy
        self.temp_factor = temp_factor

    def __repr__(self):
        return f"Chain ID: {self.chain_id}, Residue Name: {self.residue_name}, Residue Number: {self.residue_number}, Atom Name: {self.atom_name}"

class Atoms:
    def __init__(self, atoms, structure):
        self.atoms = atoms
        self.structure = structure

    @property
    def atom_length(self):
        return len(self.atoms)
    
    def first_atom(self):
        return self.atoms[0]
        
    def __iter__(self):
        return iter(self.atoms)
        
    def pdb_text(self):
        pdb_lines = []
        for atom in self.atoms:
            pdb_lines.append(f"ATOM  {atom.atom_number:5d} {atom.atom_name:>4s} {atom.residue_name:>3s} {atom.chain_id} {atom.residue_number:4d}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy if atom.occupancy else 1.00:6.2f}{atom.temp_factor if atom.temp_factor else 0.00:6.2f}")
        return "\n".join(pdb_lines)
        
    def write(self, file_path):
        output_dir = os.path.dirname(file_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        

        with open(file_path, 'w') as file:
            file.write(self.pdb_text())

            
    def show(self, cartoon=True, sticks=False, highlight_atoms=None): 
        import nglview as nv
        pdb_content = self.pdb_text()



        view = nv.NGLWidget()
        view.add_component(pdb_content, ext='pdb')
        view.clear_representations()
        if cartoon:
            view.add_cartoon()
        if sticks:
            view.add_ball_and_stick()
        if highlight_atoms:
            view.add_ball_and_stick("sphere", selection=highlight_atoms, color="red")

        return view
        
    def __filter_by_residue(self, residue_name):
        return Atoms([atom for atom in self.atoms if atom.residue_name == residue_name], self.structure)

class StructureRead:
    def __init__(self, file_path="", atoms=None):
        super().__init__()
        if file_path != "":
            self.file_path = file_path
            self.atoms = self.read_pdb()
        else:
            raise ValueError("Please provide a file path")

    def read_pdb(self):
        atoms_list = []
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_number = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    residue_number = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    occupancy = float(line[54:60].strip()) if line[54:60].strip() else None
                    temp_factor = float(line[60:66].strip()) if line[60:66].strip() else None
                    atom = Atom(residue_number, residue_name, atom_name, atom_number, x, y, z, 0.0, 0.0, 0.0, self, chain_id=chain_id, occupancy=occupancy, temp_factor=temp_factor)
                    atoms_list.append(atom)
        return atoms_list
