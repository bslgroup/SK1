# Load structure and trajectory
mol new step5_input.psf
mol addfile step5_input.pdb
mol addfile dcd/sk1.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

# Define selection text strings
set seltext1  "protein and resid 146 to 157 and alpha"
set seltext2  "protein and resid 350 to 362 and alpha"

# Output labels
set names {TM1 TM2}

# Get number of frames
set num_steps [molinfo top get numframes]

# Loop through all 14 selections
for {set i 1} {$i <= 2} {incr i} {
    set seltext  [set seltext$i]
    set selname  [lindex $names [expr $i - 1]]

    set sel      [atomselect top $seltext]
    set sel_ref  [atomselect top $seltext frame 0]

    set outfile [open "${selname}_rmsd.txt" "w"]
    puts $outfile "Frame RMSD"

    for {set frame 1} {$frame < $num_steps} {incr frame} {
        $sel frame $frame
        set trans [measure fit $sel $sel_ref]
        $sel move $trans
        set rmsd [measure rmsd $sel $sel_ref]
        puts $outfile "$frame $rmsd"
    }

    close $outfile
}
quit
