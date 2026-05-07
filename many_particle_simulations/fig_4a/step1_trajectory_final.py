import numpy as np

class Atom:
    def __init__(self, *values):
        # Use a dictionary to store attributes with original names
        field_names = [
            "type", "id", "x", "y", "z", "l", "na1", "na2", "nb1", "nb2",
            "2sbridge", "a1a2", "a2a1", "a1b1", "a2b2", "b1b2", "b2b1", "3sbridge",
            "a1a2b1", "a1a2b2", "b1b2a1", "b1b2a2", "sur_2sbrg", "sur_3sbrg",
            "tot_a1", "tot_a2", "tot_b1", "tot_b2", "nn"
        ]
        self.attributes = {field: self.convert(value) for field, value in zip(field_names, values)}

    def convert(self, value):
        """Convert values to appropriate types."""
        try:
            return int(value) if '.' not in value else float(value)
        except ValueError:
            return value

    def __getitem__(self, key):
        return self.attributes[key]

    def __setitem__(self, key, value):
        self.attributes[key] = value


def read_trajectory(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    timestep_data = {}
    current_timestep = None

    lines_iter = iter(lines)
    for line in lines_iter:
        if line.startswith("ITEM: TIMESTEP"):
            if current_timestep is not None:
                current_timestep["atoms"] = read_atoms(timestep_data)
                yield current_timestep
            try:
                # Convert scientific notation to float, then cast to int
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
        yield current_timestep


def read_atoms(lines):
    atoms = []
    for line in lines:
        values = line.split()
        if len(values) == 29:  # Expecting exactly 29 columns
            atom = Atom(*values)
            atoms.append(atom)
    return atoms


def calculate_averages(file_path, output_file):
    with open(output_file, 'w') as f:
        f.write("Timestep\tAverage_sur_2sbrg\tAverage_sur_3sbrg\n")

        for timestep_data in read_trajectory(file_path):
            timestep = timestep_data["timestep"]
            atoms = timestep_data["atoms"]

            if atoms:
                # Extract the sur_2sbrg and sur_3sbrg values
                sur_2sbrg_values = [atom["sur_2sbrg"] for atom in atoms]
                sur_3sbrg_values = [atom["sur_3sbrg"] for atom in atoms]

                # Calculate averages
                avg_sur_2sbrg = np.mean(sur_2sbrg_values)
                avg_sur_3sbrg = np.mean(sur_3sbrg_values)

                # Write the timestep and averages to the file
                f.write(f"{timestep}\t{avg_sur_2sbrg:.6f}\t{avg_sur_3sbrg:.6f}\n")


def main():
    file_path = 'long_Gillespie_BD_negno_test18_AB_DGT-2_DG_-20_rho_ideal_1.0e-05_1000_modified.dump'  # Replace with your file path
    output_file = 'averages_output.txt'  # Replace with your desired output file name
    calculate_averages(file_path, output_file)


if __name__ == "__main__":
    main()

