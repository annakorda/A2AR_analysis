set ref [atomselect top "resname UNK" frame 0]
set all [atomselect top "resname UNK"]
set nframes [molinfo top get numframes]
set outfile [open "MOD2_ligand_rmsd.dat" w]
for {set i 0} {$i < $nframes} {incr i} {
	$all frame $i
	set rmsd [measure rmsd $all $ref]
	puts $outfile "$i $rmsd"
}
close $outfile

