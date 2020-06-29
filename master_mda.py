##!init
pdb = "structure.filtered.pdb"
xtc = "traj.filtered.xtc"

##!load
import MDAnalysis as mda
t = mda.Universe(pdb, xtc)

##!frame
# Self-explanatory
t.atoms.n_atoms
t.trajectory.n_frames

# Atoms by 3
t.atoms.positions

# Unit cell
t.atoms.dimensions
 

##!select
y72_oeta = t.select_atoms("resid 72 and name OH and segid 0")
w288_chi1 = t.select_atoms("resid 288 and name N CA CB CG and segid 0")

y72_oeta.occupancies
y72_oeta[0].occupancy   # also

##!filter
bb = t.select_atoms("backbone")

with mda.Writer("bb_frame0.pdb") as w:
    t.trajectory[0]
    w.write(bb)

# Also bb.write() for single frames

##!geometry
# Self-explanatory
bb.center_of_mass()

# Current frame
w288_chi1.dihedral.value()

# All frames (iterator)
[w288_chi1.dihedral.value()
   for f in t.trajectory]


##!align
from MDAnalysis.analysis.rms import RMSD
R = RMSD(atomgroup=t, 
         reference=t, 
         select="backbone",  # align set
         groupselections=["protein"]) 
R.run()

# Measures found in column 4 and on
rmsd_traj = R.rmsd[:,3]      

##!finalize
import numpy as np
np.savetxt("rmsd_mda.out",rmsd_traj,fmt="%.3f")
