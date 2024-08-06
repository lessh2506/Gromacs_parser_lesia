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
    def __init__(self, atoms_list, structure, reassign_structure_to_atoms=False):
        self.atoms_list=atoms_list
        self.structure=structure

        #why do we need this if?
        #if reassing structure to atoms
        if reassign_structure_to_atoms:
            for atom in self.atoms_list:
                atom.structure=structure


        @property
        def atom_length(self):
            return len(self.atoms_list)
        
        def __iter__(self):
            return iter(self.atoms_list)
        
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
        
