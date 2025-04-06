set ref [atomselect top "resname ADN" frame 0]
set all [atomselect top "resname ADN"]
set nframes [molinfo top get numframes]
set outfile [open "ADN_ligand_rmsd.dat" w]
for {set i 0} {$i < $nframes} {incr i} {
	$all frame $i
	set rmsd [measure rmsd $all $ref]
	puts $outfile "$i $rmsd"
}
close $outfile

