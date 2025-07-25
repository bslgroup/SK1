#!/bin/bash

cat > RMSF.tcl << EOF
#Writes a RMSF scrpits for a single chain/segment you can combine this scrpit if you have one segment/chain in the system

mol new step5_input.psf
mol addfile step5_input.pdb
mol addfile dcd/sk1.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

source macros_sk1.tcl

set all [atomselect top "allprot"]
set sel0 [\$all num]
set sel [atomselect top "proa1"]
set ref [atomselect top "proa1" frame 0]

set num_steps [molinfo 0 get numframes]
for {set frame 0} {\$frame < \$num_steps} {incr frame} {
    \$all frame \$frame
    \$sel frame \$frame

set trans [measure fit \$sel \$ref]
\$all move \$trans
}

for {set i 0} {\$i < [\$sel num]} {incr i} {
#set rmsf [measure rmsf \$sel first 0 last -1 step 5]
     set rmsf [measure rmsf \$sel first 0 last -1]
     puts stderr "[expr {\$i+5}] \t [lindex \$rmsf \$i]"
#prints overall RMSF for each residue in the chain/segment
}

quit
EOF

vmd -dispdev text -e RMSF.tcl 2> rmsf.txt
rm RMSF.tcl
