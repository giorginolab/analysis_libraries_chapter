##!init
pdb <- "structure.filtered.pdb"
xtc <- "traj.filtered.xtc"
dcd <- "traj.filtered.dcd"


##!load
library(bio3d)
tp <- read.pdb(pdb) 
tp$xyz <- read.dcd(dcd)

##!frame
nrow(tp$xyz)        # 40 frames
nrow(tp$atom)       # 28799 atoms

## Accessing coordinates in frame 0
## reshaped for convenience
xyz <- tp$xyz[1,]
xyz <- matrix(xyz, ncol=3, byrow=T)

## Or: array(xyz,c(40,3,28799))

##!select
pdb <- tp$atom
y72_oeta <- pdb[pdb$resno == 72   &
                pdb$elety == "OH" &
                pdb$chain == "0" , ]

w288_chi1 <- atom.select(tp,
         elety=c("N","CA","CB","CG"),
         resno=288, chain="0")

y72_oeta$o


##!filter
bb <- atom.select(tp, "backbone")
tp_bb <- trim(tp, bb)   

## Select frame 1 (i.e. 0) only
bb_ref <- tp_bb
bb_ref$xyz <- trim(bb_ref$xyz, 1)

write.pdb(bb_ref, "bb_frame0.pdb")

##!geometry
# Center of mass
com(tp_bb)

# Torsion, first frame
tmp <- tp$xyz[1, w288_chi1$xyz]
torsion.xyz(c(t(tmp)))

# All frames (reshape as 1D vector)
tmp <- tp$xyz[ , w288_chi1$xyz]
torsion.xyz(c(t(tmp)))


##!align
meas.r <- trim(tp$xyz,1)
meas.t <- tp$xyz
alg.set <- bb$xyz

fitted <- fit.xyz(fixed = meas.r,
                  mobile = meas.t,
                  fixed.inds = alg.set,
                  mobile.inds = alg.set)

meas.set <- atom.select(tp,
                        "protein")$xyz
rmsd.traj <- rmsd(a = meas.r,
                  b = fitted,
                  a.inds = meas.set,
                  b.inds = meas.set)

##!finalize
write.table(rmsd.traj,"rmsd_bio3d.out")

