# Load structure and trajectory
mol new step5_input.psf
mol addfile step5_input.pdb
mol addfile dcd/sk1.dcd waitfor all

# Parameters
set lipid_sel_text "resname POPC POPE POPS SAPI CHL1"
set cutoff 4.0
set nf [molinfo top get numframes]

# Residues of interest (name + resid)
set residues_of_interest {
    "LYS 27"
    "LYS 29"
    "ARG 186"
    "LEU 194"
    "PHE 197"
    "LEU 198"
}

# Generate residue keys like LYS27
set res_keys {}
foreach res $residues_of_interest {
    set parts [split $res " "]
    lappend res_keys "[lindex $parts 0][lindex $parts 1]"
}

# Initialize contact counter
array set res_contacts {}
foreach key $res_keys {
    set res_contacts($key) 0
}

# Loop through frames
for {set i 0} {$i < $nf} {incr i} {
    set lipid_sel [atomselect top "$lipid_sel_text" frame $i]
    set contact_sel [atomselect top "name CA and protein and within $cutoff of ($lipid_sel_text)" frame $i]

    set resnames [$contact_sel get resname]
    set resids   [$contact_sel get resid]
    set seen {}

    for {set j 0} {$j < [llength $resids]} {incr j} {
        set key "[lindex $resnames $j][lindex $resids $j]"
        if {[lsearch -exact $seen $key] == -1} {
            lappend seen $key
            if {[info exists res_contacts($key)]} {
                incr res_contacts($key)
            }
        }
    }

    $lipid_sel delete
    $contact_sel delete
}

# Output results
set out [open "lipid_contacts_percentage.txt" w]
puts $out "Residue\tContact_Frames\tContact_Percentage"

foreach key $res_keys {
    set contact_frames $res_contacts($key)
    set percent [format "%.2f" [expr {100.0 * $contact_frames / $nf}]]
    puts $out "$key\t$contact_frames\t$percent"
}
close $out

puts "Analysis complete. Output written to lipid_contacts_percentage.txt"

quit
