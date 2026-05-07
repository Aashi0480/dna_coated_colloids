import numpy as np

class Atom:
    def __init__(self, *values):
        field_names = [
            "type", "id", "x", "y", "z", "l", "na1", "na2", "nb1", "nb2",
            "2sbridge", "a1a2", "a2a1", "a1b1", "a2b2", "b1b2", "b2b1", "3sbridge",
            "a1a2b1", "a1a2b2", "b1b2a1", "b1b2a2", "sur_2sbrg", "sur_3sbrg",
            "tot_a1", "tot_a2", "tot_b1", "tot_b2", "nn"
        ]
        self.attributes = {field: self.convert(value) for field, value in zip(field_names, values)}
    
    def convert(self, value):
        try:
            return int(value) if '.' not in value else float(value)
        except ValueError:
            return value

    def __getitem__(self, key):
        return self.attributes[key]

    def __setitem__(self, key, value):
        self.attributes[key] = value


def read_trajectory(file_path, max_timesteps=None):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    all_timesteps = []
    current_timestep = None
    timestep_data = []
    
    lines_iter = iter(lines)
    for line in lines_iter:
        if line.startswith("ITEM: TIMESTEP"):
            if current_timestep is not None:
                current_timestep["atoms"] = read_atoms(timestep_data)
                all_timesteps.append(current_timestep)
            try:
                current_timestep = {"timestep": int(float(next(lines_iter).strip()))}
                timestep_data = []
            except (ValueError, StopIteration):
                current_timestep = None
        elif current_timestep is not None and line.startswith("ITEM: ATOMS"):
            atom_format = line.split()[2:]
        elif current_timestep is not None and not line.startswith("ITEM:"):
            timestep_data.append(line.strip())

    if current_timestep is not None:
        current_timestep["atoms"] = read_atoms(timestep_data)
        all_timesteps.append(current_timestep)
    
    if max_timesteps:
        return all_timesteps[-max_timesteps:]
    return all_timesteps


def read_atoms(lines):
    atoms = []
    for line in lines:
        values = line.split()
        if len(values) == 29:
            atom = Atom(*values)
            atoms.append(atom)
    return atoms


def group_atoms_by_z_range(atoms):
    ranges = {
        "z<7": [], "7<=z<15": [], "15<=z<23": [], "23<=z<31": [],
        "31<=z<39": [], "39<=z<47": [], "z>47": []
    }
    for atom in atoms:
        z = atom["z"]
        if z < 7:
            ranges["z<7"].append(atom)
        elif 7 <= z < 15:
            ranges["7<=z<15"].append(atom)
        elif 15 <= z < 23:
            ranges["15<=z<23"].append(atom)
        elif 23 <= z < 31:
            ranges["23<=z<31"].append(atom)
        elif 31 <= z < 39:
            ranges["31<=z<39"].append(atom)
        elif 39 <= z < 47:
            ranges["39<=z<47"].append(atom)
        else:
            ranges["z>47"].append(atom)
    return ranges


def calculate_averages(atoms, field_names):
    ranges = group_atoms_by_z_range(atoms)
    averages = []
    
    for z_range, atoms_in_range in ranges.items():
        if atoms_in_range:
            field_values = [sum([atom[field] for field in field_names]) for atom in atoms_in_range]
            z_values = [atom["z"] for atom in atoms_in_range]
            
            avg_field_value = np.mean(field_values)
            avg_z_value = np.mean(z_values)
            
            averages.append((avg_z_value, avg_field_value))
    
    return averages



def save_field_to_file(data, field_names, output_file):
    with open(output_file, 'w') as file:
        # Join field names as a string (e.g., 'l_na1_na2_nb1_nb2')
        field_name_str = "_".join(field_names)
        
        # Write the header
        file.write(f"{'timestep':<10}{'z':<15}{field_name_str:<15}\n")
        
        for timestep_data in data:
            timestep = timestep_data["timestep"]
            if not timestep_data["atoms"]:
                continue
            averages = calculate_averages(timestep_data["atoms"], field_names)
            for avg_z_value, avg_field_value in averages:
                if avg_field_value is not None:
                    file.write(f"{timestep:<10}{avg_z_value:<15.6f}{avg_field_value:<15.6f}\n")


def main():
    file_path = 'long_Gillespie_BD_negno_test18_AB_DGT0_DG_-20_rho_ideal_1.0e-05_1000.dump'
    #fields_to_save = ["a1a2", "a2a1", "a1b1", "a2b2", "b1b2", "b2b1",
                      #"a1a2b1", "a1a2b2", "b1b2a1", "b1b2a2"]
    
    fields_to_save = [["a1a2", "a2a1"]]
    
    max_timesteps = int(input("Enter the number of timesteps to process (or -1 for all): "))
    max_timesteps = None if max_timesteps == -1 else max_timesteps
    
    trajectory_data = read_trajectory(file_path, max_timesteps)
    
    for field in fields_to_save:
        #output_file = f"{field}_averages.txt"
        output_file = f"na1a2_na2a1_averages.txt"
        save_field_to_file(trajectory_data, field, output_file)
        print(f"Saved {field} averages to {output_file}")


if __name__ == "__main__":
    main()

