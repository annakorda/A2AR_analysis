# A2AR_analysis

## Ligand Binding Analysis and Steered MD
This repository contains scripts, plots, and data files from our Master's project in Molecular Dynamics Simulations class. The main goal of the project is to analyze ligand binding modes, protein-ligand interactions, solvent density around an adenosine receptor (A2AR) in the presence of different ligands and also ligand unbinding. Due to file size constraints, only a subset of the project is hosted here; the complete dataset including trajectory files and heavy outputs is available via a shared Drive folder, available upon request.

## Project Overview
We investigated the interaction of the A2A adenosine receptor with three ligands:

ADN (Adenosine)
MOD1 and MOD2 (Ligand modifications)

For each ligand system, we performed:

1. Molecular dynamics simulations
2. Ligand RMSD and protein RMSD analysis
3. Contact frequency calculations
4. Water occupancy mapping 
5. Ligand binding pose comparison and superimposition
6. Steered Molecular Dynamics

## Repository Contents
### A2AR-ADN/
Contains data and plots specifically for the adenosine-bound receptor simulation.

contacts.txt, final_contacts.tsv: Raw and processed contact files.

lig_rmsd.dat, protein_rmsd.dat: Ligand and protein RMSD values.

Plots_rmsd_contacts.R: R script for generating analysis plots.

plots/: Includes various visualization outputs like RMSD curves and contact maps.

watermap.dx: Water occupancy map for solvent density.

### Ligand-Modifications/
Main folder for comparative analysis between ADN, MOD1, and MOD2.

Charmm-gui-input/: Finalized ligand structures and CHARMM parameters used in simulation.

RMSD_Contacts/: RMSD and contact analysis for all three ligands.

Superimposition/: Images comparing receptor conformations across ligand states.

Water_map/: Includes .dx volumetric maps and gro files used for water occupancy calculations, along with a Python script for plotting.

Final_plots_all/: Summary plots combining data from all systems for quick comparison.

### Steered_Molecular_Dynamics/
Simulations and analysis of ligand unbinding via Steered MD (SMD).

A2AR-ADN/namd/ and A2AR-MOD2/namd/: NAMD SMD input files and scripts

results/: Output logs, force data, and contact results

Results_analysis/:

plots/: Heatmaps, unbinding curves, PMF and VMD snapshots

scripts/: Analysis notebooks for contacts, unbinding patter, and PMF profile.

Unbinding_graph.xlsx: Processed ligand_smd.out files of both ligands ready for plotting.

## LICENSE
Open-source license for the code and data shared here.

## Accessing Full Results
Due to GitHub's file size limitations, large files such as trajectory outputs are provided via a Google Drive folder, available upon request.
Contents:
Full GROMACS trajectories, CHARMM-GUI output directories, Steered MD trajectories

## How to Reproduce the Analysis
Contact Analysis: Use the provided .sh scripts to extract and reformat contacts from simulation output.

Plotting: Run the R script (Plots_rmsd_contacts.R) to regenerate RMSD and contact plots.

Water Maps: Use water_plots_code.py to visualize .dx files. Requires MDAnalysis, matplotlib, and numpy.

SMD Visualization: Use Unbinding_graph.ipynb anf PMF.ipynb to explore force-time relationships and PMF.

