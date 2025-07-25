# Load topology + trajectory
mol new step5_input.psf
mol addfile step5_input.pdb
mol addfile dcd/sk1.dcd waitfor all

# Analysis parameters
set cutoff       4.0
set residues     {27 29 186 194 197 198}
set lipid_types  {POPC POPE POPS SAPI CHL1}
set nf           [molinfo top get numframes]

# Initialize counters: contact_counts(resid, lipid)
array set contact_counts {}
foreach resid $residues {
    foreach lip $lipid_types {
        set contact_counts($resid,$lip) 0
    }
}

# Loop over frames
for {set i 0} {$i < $nf} {incr i} {
    foreach resid $residues {
        # 1) select only the Cα atom of this residue
        set sel_ca [atomselect top "resid $resid and name CA" frame $i]
        if {[$sel_ca num] != 1} {
            # if it's zero (or >1) something's off—skip
            $sel_ca delete
            continue
        }

        # 2) for each lipid type, check Cα → lipid within cutoff
        foreach lip $lipid_types {
            set sel_near [atomselect top \
              "resname $lip and within $cutoff of (resid $resid and name CA)" \
              frame $i]
            if {[$sel_near num] > 0} {
                incr contact_counts($resid,$lip)
            }
            $sel_near delete
        }

        $sel_ca delete
    }
}

# Write out results
set out [open "lipid_contact_CA_summary.txt" w]
puts $out "Residue\tPOPC\tPOPE\tPOPS\tSAPI\tCHL1"
foreach resid $residues {
    puts -nonewline $out "$resid"
    foreach lip $lipid_types {
        puts -nonewline $out "\t$contact_counts($resid,$lip)"
    }
    puts $out ""
}
close $out

puts "✅ Written lipid_contact_CA_summary.txt"
quit
