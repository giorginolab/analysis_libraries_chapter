Errata
======

This file collects known errata arising since the publication, mostly due to changes in underlying software.


Listing 12: Geometry operations
-------------------

Several functions from `htmd` were moved to the related `moleculekit` package. This affects `dihedralAngle`:

    htmd.molecule.util.dihedralAngle(...

is now (with the appropriate import)

    import moleculekit.dihedral
    moleculekit.dihedral.dihedralAngle(...
    
    
Listing 13: RMSD-based alignments
-------------------

    htmd.molecule.util.molRMSD(...
    
is now

    import moleculekit.util
    moleculekit.util.molRMSD(...
    

Listing 14: Initialization code
-------------------

The OPM files are no longer available at the URL. You may download similar files e.g. from
the [GPCRmd](https://submission.gpcrmd.org/home/) database.

