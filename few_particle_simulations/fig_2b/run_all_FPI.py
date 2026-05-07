import subprocess

# Generate values from 5.4 to 7.0 (step 0.1, 17 points)
values = [f"{5.4 + 0.1*i:.1f}" for i in range(17)]

results = []

for val in values:
    print(f"\nProcessing {val}...")

    clean_src = f"FPI_clean_{val}.c"
    clean_exe = f"FPI_clean_{val}"

    free_src = f"FPI_free_energy_all_system_{val}.c"
    free_exe = f"FPI_free_energy_all_system_{val}"

    # Compile
    subprocess.run(["g++", clean_src, "-o", clean_exe], check=True)
    subprocess.run(["g++", free_src, "-o", free_exe], check=True)

    # Run clean executable
    subprocess.run([f"./{clean_exe}"], check=True)

    # Run free energy executable and capture output
    result = subprocess.run(
        [f"./{free_exe}"],
        capture_output=True,
        text=True,
        check=True
    )

    # Extract last line
    last_line = result.stdout.strip().split("\n")[-1]

    # Parse distance and energy
    parts = last_line.split()
    distance = float(parts[0])
    energy = float(parts[1])

    print(f"Distance = {distance}, Energy = {energy}")

    results.append((distance, energy))

# ---- SORT in descending order of distance ----
results.sort(reverse=True, key=lambda x: x[0])

# ---- NORMALIZE (first entry energy = 0) ----
ref_energy = results[0][1]
normalized_results = [(d, e - ref_energy) for d, e in results]

# ---- SAVE ----
output_file = "FPIenergy_surB_-10-2.txt"

with open(output_file, "w") as out:
    for d, e in normalized_results:
        out.write(f"{d} {e}\n")

print(f"\nSaved normalized & sorted results to {output_file}")
