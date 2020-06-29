##!init
import numpy as np
pdb = "structure.filtered.pdb"
xtc = "traj.filtered.xtc"


##!load
import mdtraj as mdt
t = mdt.load(xtc, top=pdb)

##!frame
# Number of frames
len(t)                         

# Frames by Atoms by 3
t.xyz.shape
# Coordinates in frame 0
t.xyz[0]                        

# Unit cell
t.unitcell_lengths[0,:]


##!select
y72_oeta = t.topology.select("residue 72 and name OH and chainid 0")
w288_chi1 = t.topology.select("residue 288 and name N CA CB CG and chainid 0")

t.atom_slice(y72_oeta).topology.\
    atom(0).element


##!filter
bb = t.topology.select("backbone")
t_bb = t.atom_slice(bb)   # Subset

# Select frame 0 (also [0])
bb_ref = t_bb.slice(0)

# Write to file
bb_ref.save("bb_frame0.pdb")


##!geometry
# Self-explanatory
mdt.compute_center_of_mass(bb_ref)

# First frame (to degrees)
mdt.compute_dihedrals(t[0],
    [w288_chi1])*180.0/np.pi

# All frames
mdt.compute_dihedrals(t,
    [w288_chi1])*180.0/np.pi




##!align
fitted = t[:]            # Copy

# Align
alg_set = t.topology.select("backbone")
fitted.superpose(reference=t,
                 frame=0, 
                 atom_indices=alg_set)

meas_set = t.topology.select("protein")
meas_r = t.atom_slice(meas_set)
meas_t = fitted.atom_slice(meas_set)

d = meas_t.xyz-meas_r.xyz[0]
rmsd_traj = 10 * np.sqrt( np.mean(
        np.sum(d**2,axis=2),axis=1))

##!finalize
import numpy as np
np.savetxt("rmsd_mdt.out",rmsd_traj,fmt="%.3f")
