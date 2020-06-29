##!init
pdb = "structure.filtered.pdb"
xtc = "traj.filtered.xtc"

##!load
from htmd.ui import *
t=Molecule(pdb)
t.read(xtc)

##!frame
t.numFrames
t.numAtoms

# Atoms by 3 by frames
t.coords

# Unit cell
t.box[:,0]


##!select
y72_oeta = t.atomselect("resid 72 and name OH and chain 0")
w288_chi1 = t.atomselect("resid 288 and name N CA CB CG and chain 0")

t.occupancy[y72_oeta]


##!filter
bb = t.copy()
bb.filter("backbone")
bb.dropFrames(keep=0)
bb.write("bb_frame0.pdb")


##!geometry
# Center of geometry: add "weights=bb.masses" if available in topology
np.average(bb.coords, axis=0).T

# First frame
htmd.molecule.util.dihedralAngle(
             t.coords[w288_chi1,:,0])

# All frames
htmd.molecule.util.dihedralAngle(
             t.coords[w288_chi1,:,:])


##!align
meas_r = t.copy()
meas_r.dropFrames(keep=0)
meas_t = t.copy()

meas_t.align("backbone",meas_r)

meas_set = meas_t.atomselect("protein")
rmsd_traj = htmd.molecule.util.molRMSD(
    meas_t,meas_r,meas_set,meas_set)

##!finalize
import numpy as np
np.savetxt("rmsd_htmd.out",rmsd_traj,fmt="%.3f")
