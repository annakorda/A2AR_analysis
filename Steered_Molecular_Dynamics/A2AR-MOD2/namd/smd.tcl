# List of ligand atoms (IDs obtained from VMD)
set ligand_atoms {4975 4976 4977 4978 4979 4980 4981 4982 4983 4984 4985 4986 4987 4988 4989 4990 4991 4992 4993 4994 4995 4996 4997 4998 4999 5000 5001 5002 5003 5004 5005 5006 5007 5008 5009}

# Create a group with all ligand atoms
set a2 [addgroup $ligand_atoms]

# Simulation parameters
set Tclfreq 50
set t 0
set k 7.2        ;# Force constant (kcal/mol/Å²)
set v 0.00002    ;# Pulling velocity (Å/timestep)
set target_z 90.0

# Output file
set outfilename ligand_smd.out
open $outfilename w

# Function to calculate COM Z using the "coor" array (valid in NAMD Tcl)
proc calcCOMZ_simple { group } {
    loadcoords coor
    set sum_z 0.0
    set n 0

    foreach atom $group {
        set z [lindex $coor($atom) 2]
        set sum_z [expr {$sum_z + $z}]
        incr n
    }

    return [expr {$sum_z / double($n)}]
}

# Main function to apply the force
proc calcforces {} {

  global Tclfreq t k v a2 target_z outfilename

  # Update coordinates and dynamically calculate COM Z
  set comz [calcCOMZ_simple $a2]

  # Target Z coordinate for this step
  set target_current_z [expr {$comz + ($v * $t)}]

  # Calculate force
  set fz [expr {$k * ($target_current_z - $comz)}]

  # Avoid NaN values
  if { $fz != $fz } {
      puts "ERROR: fz is NaN, ending simulation"
      exit
  }

  # Apply force
  set force [list 0.0 0.0 $fz]
  addforce $a2 $force

  # Save results every Tclfreq steps
  if { [expr {$t % $Tclfreq}] == 0 } {
      set outfile [open $outfilename a]
      set time [expr {$t*2/1000.0}]  ;# Convert timesteps to ps
      puts $outfile "$time $comz $fz"
      close $outfile
  }

  # Stop if target Z is reached
  if { $comz >= $target_z } {
      puts "SMD finished: Ligand reached Z = $target_z Å"
      exit
  }

  incr t
  return
}

