# normalize.py

x = []
y = []

# Read file
with open("FPI_distance_energy_sorted.txt", "r") as f:
    for line in f:
        if line.strip():
            xi, yi = map(float, line.split())
            x.append(xi)
            y.append(yi)

# Find maximum y
ymax = max(y)

# Normalize: y_new = y - max(y)
y_norm = [yi - ymax for yi in y]

# Write output
with open("normalized.dat", "w") as f:
    for xi, yi in zip(x, y_norm):
        f.write(f"{xi} {yi}\n")

print("Done. Output written to normalized.dat")
