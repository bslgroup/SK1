#!/bin/bash

cat > RMSD.tcl << EOF #writes a file for RMSD analysis

mol new step5_input.psf
mol addfile step5_input.pdb
mol addfile dcd/sk1.dcd type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

source macros_sk1.tcl

set allprot [atomselect top "allprot"]
#atom sections from macros file, should list all the selections required for respective analysis.
set allprota [atomselect top "allprota"]
set PROA1 [atomselect top "proa1"]
#set PROB1 [atomselect top "prob1"]
#set PROA2 [atomselect top "proa2"]
#set PROB2 [atomselect top "prob2"]


set allR [atomselect top "allprot" frame 0 ]
#include the respective reference (frame 0 which is pdb crystal structure) for each selection of above for rmsd calculation.
set allaR [atomselect top "allprota" frame 0 ]

set PROA1R [atomselect top "proa1" frame 0 ]
#include the respective reference (frame 0 which is pdb crystal structure) for each selection of above for rmsd calculation.
#set PROB1R [atomselect top "prob1" frame 0 ]
#set PROA2R [atomselect top "proa2" frame 0 ]
#set PROB2R [atomselect top "prob2" frame 0 ]


set num_steps [molinfo 0 get numframes]
#calculates the number of frames in the trajectory
for {set frame 1} {\$frame < \$num_steps} {incr frame} {
#loops through each frame from 0 to n for measuring the analysis
    \$allprot frame \$frame
    \$allprota frame \$frame
    \$PROA1 frame \$frame
#    \$PROB1 frame \$frame
#    \$PROA2 frame \$frame
#    \$PROB2 frame \$frame

    set trans [measure fit \$allprot  \$allR]
 #alignment of protein with intial structure
    \$allprot move \$trans
    set rmsdh0 [measure rmsd \$allprot \$allR]
    set transa [measure fit \$allprota  \$allaR]
    \$allprota move \$transa
    set rmsdh1 [measure rmsd \$allprota \$allaR]
 #calculation of rmsd with respective to intial structure
    set trans [measure fit \$PROA1 \$PROA1R]
    \$allprot move \$trans
    set rmsdh2 [measure rmsd \$PROA1 \$PROA1R]
#    set trans [measure fit \$PROB1 \$PROB1R]
#    \$allprot move \$trans
#    set rmsdh3 [measure rmsd \$PROB1 \$PROB1R]
#    set trans [measure fit \$PROA2 \$PROA2R]
#    \$allprot move \$trans
#    set rmsdh5 [measure rmsd \$PROA2 \$PROA2R]
#    set trans [measure fit \$PROB2 \$PROB2R]
#    \$allprot move \$trans
#    set rmsdh6 [measure rmsd \$PROB2 \$PROB2R]

     puts stderr "\$frame \$rmsdh0 \$rmsdh1 \$rmsdh2"
 # prints the value of each caluction from above part and prints in order (colum 0 = frame number number, colum 1 = values of rmsh1 ......)
}

quit

EOF

vmd -dispdev text -e  RMSD.tcl 2> rmsd.txt
rm RMSD.tcl
