##!init
set pdb structure.filtered.pdb
set xtc traj.filtered.xtc

package require pbc

##!load
set t [mol new $pdb]
animate delete all
mol addfile $xtc waitfor all

##!frame
# Number of frames
molinfo top get numframes

set t [atomselect top all]
$t num;         # Number of atoms

$t frame 0
$t get {x y z}; # Coordinates

pbc get;        # Unit cell

##!select
set y72_oeta [atomselect top "resid 72 and name OH and chain 0"]
set w288_chi1 [atomselect top "resid 288 and name N CA CB CG and chain 0"]

# Access the "occupancy" property
# of a single atom
$y72_oeta get occupancy

##!filter
set bb [atomselect top backbone]

# Write backbone frame 0
animate write pdb bb_frame0.pdb               beg 0 end 0 sel $bb

##!geometry
# Center of mass
$bb frame 0
measure center $bb weight mass

# W288, chi 1, first frame
measure dihed [$w288_chi1 get index]

# All frames
measure dihed [$w288_chi1 get index] first 0 last 40


##!align
# Convenience functions from
# github.com/tonigi/vmd_extensions
source ~/VMDextensions.tcl

# t: trajectory; r: reference;
# alg: alignment; meas: measurement
set meas_t [atomselect top protein]
set meas_r [atomselect top protein
	    frame 0]
set alg_t [atomselect top backbone]
set alg_r [atomselect top backbone
	   frame 0]

set rmsd_traj [rmsdOf $meas_t $meas_r
	              $alg_t  $alg_r]

##!finalize
set f [open "rmsd_vmd.out" w]
foreach i $rmsd_traj {
    puts $f [format "%.3f" $i]
}
close $f
quit
