

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import MDAnalysis as mda
from adjustText import adjust_text


### ADN ###

# FILE PATHS 
gro_file = "ADN_50ns_rep1.gro"
dx_file = "volmap_out_ADN.dx"

#  Load structure
u = mda.Universe(gro_file)
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname ADN")
ligand_coords = ligand.positions
lig_com = ligand.center_of_mass()

#  Load .dx file 
with open(dx_file, "r") as f:
    lines = f.readlines()

dimensions, origin, delta = None, None, []
data_start_index = 0

for i, line in enumerate(lines):
    if line.startswith("object 1 class gridpositions counts"):
        dimensions = list(map(int, line.strip().split()[-3:]))
    elif line.startswith("origin"):
        origin = list(map(float, line.strip().split()[1:]))
    elif line.startswith("delta"):
        delta.append(list(map(float, line.strip().split()[1:])))
    elif line.startswith("object 3 class array type double rank 0 items"):
        data_start_index = i + 1
        break

data = []
for line in lines[data_start_index:]:
    if line.strip().startswith("object"):
        break
    if line.strip():
        data.extend(map(float, line.strip().split()))

volmap = np.array(data).reshape(dimensions)

# Select Z-slice 
ligand_z_coords = ligand_coords[:, 2]
grid_z = np.array([origin[2] + i * delta[2][2] for i in range(dimensions[2])])
ligand_z_counts = [
    np.sum(np.abs(ligand_z_coords - z_val) < 3.0)
    for z_val in grid_z
]
z_index = int(np.argmax(ligand_z_counts))
print(f"Auto-selected Z-slice with most ligand atoms: {z_index} (count: {ligand_z_counts[z_index]})")

z_slice = volmap[:, :, z_index]
z_threshold = origin[2] + delta[2][2] * z_index

# Protein atoms near Z-slice 
protein_coords = protein.positions
protein_slice_coords = protein_coords[np.abs(protein_coords[:, 2] - z_threshold) < 1.5]
x_idx = ((protein_slice_coords[:, 0] - origin[0]) / delta[0][0]).astype(int)
y_idx = ((protein_slice_coords[:, 1] - origin[1]) / delta[1][1]).astype(int)

# Ligand atoms near Z-slice 
ligand_slice_coords = ligand_coords[np.abs(ligand_coords[:, 2] - z_threshold) < 3.0]
lig_x_idx = ((ligand_slice_coords[:, 0] - origin[0]) / delta[0][0]).astype(int)
lig_y_idx = ((ligand_slice_coords[:, 1] - origin[1]) / delta[1][1]).astype(int)

lig_com_x = int((lig_com[0] - origin[0]) / delta[0][0])
lig_com_y = int((lig_com[1] - origin[1]) / delta[1][1])

print(f"Ligand atoms in slice: {len(lig_x_idx)}")

# Manually select and label specific residues 
target_residues = [
    ("SER", 277),
    ("HSD", 278),
    ("LEU", 249),
    ("TRP", 246),
    ("THR", 88),
    ("GLU", 169),
]

protein_CAs = u.select_atoms("protein and name CA")
closest_residues = [res for res in protein_CAs if (res.resname, res.resid) in target_residues]

print(" Manually labeling residues:", [f"{r.resname}{r.resid}" for r in closest_residues])

# PLOT 
plt.figure(figsize=(10, 8))
plt.imshow(z_slice.T, origin='lower', cmap='viridis')
plt.colorbar(label='Water Occupancy')

# Protein points
plt.scatter(x_idx, y_idx, color='red', s=15, label='Protein')

# Ligand
if len(lig_x_idx) > 0:
    plt.scatter(lig_x_idx, lig_y_idx, color='deepskyblue', s=70, marker='*',
                edgecolor='black', linewidth=0.5, label='Ligand Atoms')
    plt.scatter(lig_com_x, lig_com_y, color='white', s=80, marker='X',
                edgecolor='black', linewidth=1.2, label='Ligand COM')

# Residue labels
texts = []

for residue in closest_residues:
    res_x = int((residue.position[0] - origin[0]) / delta[0][0])
    res_y = int((residue.position[1] - origin[1]) / delta[1][1])
    label = f"{residue.resname}{residue.resid}"

    texts.append(
        plt.text(
            res_x, res_y, label,
            fontsize=7, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )
    )

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

# Safe zoom only if ligand atoms exist
if len(lig_x_idx) > 0 and len(lig_y_idx) > 0:
    xmin = min(np.min(x_idx), np.min(lig_x_idx), lig_com_x) - 10
    xmax = max(np.max(x_idx), np.max(lig_x_idx), lig_com_x) + 10
    ymin = min(np.min(y_idx), np.min(lig_y_idx), lig_com_y) - 10
    ymax = max(np.max(y_idx), np.max(lig_y_idx), lig_com_y) + 10
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
else:
    print("Warning: No ligand atoms in this Z-slice — skipping zoom.")

# Labels and grid
plt.title(f"Water Occupancy Map ADN (Z = {z_index})")
plt.xlabel("X Grid")
plt.ylabel("Y Grid")
plt.grid(visible=True, color='gray', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("hydration_map_residuesADN.png", dpi=300, bbox_inches='tight')








### MOD1 ###

# FILE PATHS 
gro_file = "mod1_50ns.gro"
dx_file = "volmap_out_mod1.dx"

#  Load structure 
u = mda.Universe(gro_file)
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname UNK")
ligand_coords = ligand.positions
lig_com = ligand.center_of_mass()

# Load .dx file 
with open(dx_file, "r") as f:
    lines = f.readlines()

dimensions, origin, delta = None, None, []
data_start_index = 0

for i, line in enumerate(lines):
    if line.startswith("object 1 class gridpositions counts"):
        dimensions = list(map(int, line.strip().split()[-3:]))
    elif line.startswith("origin"):
        origin = list(map(float, line.strip().split()[1:]))
    elif line.startswith("delta"):
        delta.append(list(map(float, line.strip().split()[1:])))
    elif line.startswith("object 3 class array type double rank 0 items"):
        data_start_index = i + 1
        break

data = []
for line in lines[data_start_index:]:
    if line.strip().startswith("object"):
        break
    if line.strip():
        data.extend(map(float, line.strip().split()))

volmap = np.array(data).reshape(dimensions)

# Select Z-slice 
ligand_z_coords = ligand_coords[:, 2]
grid_z = np.array([origin[2] + i * delta[2][2] for i in range(dimensions[2])])

ligand_z_counts = [
    np.sum(np.abs(ligand_z_coords - z_val) < 3.0)
    for z_val in grid_z
]

z_index = int(np.argmax(ligand_z_counts))
print(f"Auto-selected Z-slice with most ligand atoms: {z_index} (count: {ligand_z_counts[z_index]})")

z_slice = volmap[:, :, z_index]
z_threshold = origin[2] + delta[2][2] * z_index

# Protein atoms near Z-slice 
protein_coords = protein.positions
protein_slice_coords = protein_coords[np.abs(protein_coords[:, 2] - z_threshold) < 1.5]
x_idx = ((protein_slice_coords[:, 0] - origin[0]) / delta[0][0]).astype(int)
y_idx = ((protein_slice_coords[:, 1] - origin[1]) / delta[1][1]).astype(int)

#  Ligand atoms near Z-slice 
ligand_slice_coords = ligand_coords[np.abs(ligand_coords[:, 2] - z_threshold) < 3.0]
lig_x_idx = ((ligand_slice_coords[:, 0] - origin[0]) / delta[0][0]).astype(int)
lig_y_idx = ((ligand_slice_coords[:, 1] - origin[1]) / delta[1][1]).astype(int)

lig_com_x = int((lig_com[0] - origin[0]) / delta[0][0])
lig_com_y = int((lig_com[1] - origin[1]) / delta[1][1])

print(f"Ligand atoms in slice: {len(lig_x_idx)}")

# Get yellow (high occupancy) voxels in this Z slice
yellow_thresh = 0.4  
yellow_voxels = np.argwhere(z_slice > yellow_thresh)

# Convert them to real 3D coords
yellow_coords_real = np.array([
    [
        origin[0] + x * delta[0][0],
        origin[1] + y * delta[1][1],
        z_threshold
    ]
    for x, y in yellow_voxels
])

# Manually select residues to label
protein_CAs = u.select_atoms("protein and name CA")

target_residues = [
    ("MET", 270),
    ("SER", 67),
    ("TYR", 271),
    ("ILE", 64),
    ("ILE", 274),
    ("LEU", 267),
    ("ALA", 63),
]

closest_residues = [res for res in protein_CAs if (res.resname, res.resid) in target_residues]

print(" Manually labeling residues:", [f"{r.resname}{r.resid}" for r in closest_residues])

# PLOT 
plt.figure(figsize=(10, 8))
plt.imshow(z_slice.T, origin='lower', cmap='viridis')
plt.colorbar(label='Water Occupancy')

# Protein points
plt.scatter(x_idx, y_idx, color='red', s=15, label='Protein')

# Ligand
if len(lig_x_idx) > 0:
    plt.scatter(lig_x_idx, lig_y_idx, color='deepskyblue', s=70, marker='*',
                edgecolor='black', linewidth=0.5, label='Ligand Atoms')
    plt.scatter(lig_com_x, lig_com_y, color='white', s=80, marker='X',
                edgecolor='black', linewidth=1.2, label='Ligand COM')

# Residue labels
texts = []

for residue in closest_residues:
    res_x = int((residue.position[0] - origin[0]) / delta[0][0])
    res_y = int((residue.position[1] - origin[1]) / delta[1][1])
    label = f"{residue.resname}{residue.resid}"

    texts.append(
        plt.text(
            res_x, res_y, label,
            fontsize=7, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )
    )

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

# Safe zoom only if ligand atoms exist
if len(lig_x_idx) > 0 and len(lig_y_idx) > 0:
    xmin = min(np.min(x_idx), np.min(lig_x_idx), lig_com_x) - 10
    xmax = max(np.max(x_idx), np.max(lig_x_idx), lig_com_x) + 10
    ymin = min(np.min(y_idx), np.min(lig_y_idx), lig_com_y) - 10
    ymax = max(np.max(y_idx), np.max(lig_y_idx), lig_com_y) + 10
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
else:
    print("Warning: No ligand atoms in this Z-slice — skipping zoom.")

# === Labels and grid
plt.title(f"Water Occupancy Map MOD1 (Z = {z_index})")
plt.xlabel("X Grid")
plt.ylabel("Y Grid")
plt.grid(visible=True, color='gray', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("hydration_map_residuesMOD1.png", dpi=300, bbox_inches='tight')









### MOD2 ###

# FILE PATHS 
gro_file = "mod2_50ns.gro"
dx_file = "volmap_out_mod2.dx"

# Load structure 
u = mda.Universe(gro_file)
protein = u.select_atoms("protein")
ligand = u.select_atoms("resname UNK")
ligand_coords = ligand.positions
lig_com = ligand.center_of_mass()

#  Load .dx file 
with open(dx_file, "r") as f:
    lines = f.readlines()

dimensions, origin, delta = None, None, []
data_start_index = 0

for i, line in enumerate(lines):
    if line.startswith("object 1 class gridpositions counts"):
        dimensions = list(map(int, line.strip().split()[-3:]))
    elif line.startswith("origin"):
        origin = list(map(float, line.strip().split()[1:]))
    elif line.startswith("delta"):
        delta.append(list(map(float, line.strip().split()[1:])))
    elif line.startswith("object 3 class array type double rank 0 items"):
        data_start_index = i + 1
        break

data = []
for line in lines[data_start_index:]:
    if line.strip().startswith("object"):
        break
    if line.strip():
        data.extend(map(float, line.strip().split()))

volmap = np.array(data).reshape(dimensions)

# Select Z-slice 
ligand_z_coords = ligand_coords[:, 2]
grid_z = np.array([origin[2] + i * delta[2][2] for i in range(dimensions[2])])

ligand_z_counts = [
    np.sum(np.abs(ligand_z_coords - z_val) < 3.0)
    for z_val in grid_z
]

z_index = int(np.argmax(ligand_z_counts))
print(f"Auto-selected Z-slice with most ligand atoms: {z_index} (count: {ligand_z_counts[z_index]})")

z_slice = volmap[:, :, z_index]
z_threshold = origin[2] + delta[2][2] * z_index

# Protein atoms near Z-slice 
protein_coords = protein.positions
protein_slice_coords = protein_coords[np.abs(protein_coords[:, 2] - z_threshold) < 1.5]
x_idx = ((protein_slice_coords[:, 0] - origin[0]) / delta[0][0]).astype(int)
y_idx = ((protein_slice_coords[:, 1] - origin[1]) / delta[1][1]).astype(int)

#  Ligand atoms near Z-slice 
ligand_slice_coords = ligand_coords[np.abs(ligand_coords[:, 2] - z_threshold) < 3.0]
lig_x_idx = ((ligand_slice_coords[:, 0] - origin[0]) / delta[0][0]).astype(int)
lig_y_idx = ((ligand_slice_coords[:, 1] - origin[1]) / delta[1][1]).astype(int)

lig_com_x = int((lig_com[0] - origin[0]) / delta[0][0])
lig_com_y = int((lig_com[1] - origin[1]) / delta[1][1])

print(f"Ligand atoms in slice: {len(lig_x_idx)}")

# Get yellow (high occupancy) voxels in this Z slice
yellow_thresh = 0.4  
yellow_voxels = np.argwhere(z_slice > yellow_thresh)

# Convert them to real 3D coords
yellow_coords_real = np.array([
    [
        origin[0] + x * delta[0][0],
        origin[1] + y * delta[1][1],
        z_threshold
    ]
    for x, y in yellow_voxels
])

# =Manually select residues to label
protein_CAs = u.select_atoms("protein and name CA")

target_residues = [
    ("TRP", 246),
    ("SER", 277),
    ("PHE", 168),
    ("LEU", 85),
    ("LEU", 249),
    ("ILE", 274),
    ("THR", 88)
]

closest_residues = [res for res in protein_CAs if (res.resname, res.resid) in target_residues]

print("Manually labeling residues:", [f"{r.resname}{r.resid}" for r in closest_residues])

# PLOT 
plt.figure(figsize=(10, 8))
plt.imshow(z_slice.T, origin='lower', cmap='viridis')
plt.colorbar(label='Water Occupancy')

# Protein points
plt.scatter(x_idx, y_idx, color='red', s=15, label='Protein')

# Ligand
if len(lig_x_idx) > 0:
    plt.scatter(lig_x_idx, lig_y_idx, color='deepskyblue', s=70, marker='*',
                edgecolor='black', linewidth=0.5, label='Ligand Atoms')
    plt.scatter(lig_com_x, lig_com_y, color='white', s=80, marker='X',
                edgecolor='black', linewidth=1.2, label='Ligand COM')

# Residue labels
texts = []

for residue in closest_residues:
    res_x = int((residue.position[0] - origin[0]) / delta[0][0])
    res_y = int((residue.position[1] - origin[1]) / delta[1][1])
    label = f"{residue.resname}{residue.resid}"

    texts.append(
        plt.text(
            res_x, res_y, label,
            fontsize=7, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )
    )

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

# Safe zoom only if ligand atoms exist
if len(lig_x_idx) > 0 and len(lig_y_idx) > 0:
    xmin = min(np.min(x_idx), np.min(lig_x_idx), lig_com_x) - 10
    xmax = max(np.max(x_idx), np.max(lig_x_idx), lig_com_x) + 10
    ymin = min(np.min(y_idx), np.min(lig_y_idx), lig_com_y) - 10
    ymax = max(np.max(y_idx), np.max(lig_y_idx), lig_com_y) + 10
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
else:
    print("Warning: No ligand atoms in this Z-slice — skipping zoom.")

# Labels and grid
plt.title(f"Water Occupancy Map MOD1 (Z = {z_index})")
plt.xlabel("X Grid")
plt.ylabel("Y Grid")
plt.grid(visible=True, color='gray', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("hydration_map_residuesMOD2.png", dpi=300, bbox_inches='tight')
