.. index:: fix adres/region

fix adres/region command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID adres/region region-atomistic region-transition region-CG

* ID, group-ID are documented in :doc:`fix <fix>` command
* adres/region = style name of this fix command
* region-atomistic = ID of region where atoms are treated atomistically
* region-transition = ID of region where atoms transition between resolutions
* region-CG = ID of region where atoms are treated as coarse-grained

Examples
""""""""

.. code-block:: LAMMPS

   region atomistic block -10.0 10.0 -10.0 10.0 -10.0 10.0
   region transition block -12.0 -10.0 -10.0 10.0 -10.0 10.0
   region cg block -20.0 -12.0 -10.0 10.0 -10.0 10.0
   fix 1 all adres/region atomistic transition cg

Description
"""""""""""

This fix manages the resolution regions for Adaptive Resolution (AdRes)
simulations. It tracks which region each atom belongs to and calculates
the lambda switching parameter for atoms in the transition region.

The fix defines three regions:
* **atomistic region**: Atoms are treated with full atomistic detail
* **transition region**: Atoms transition between atomistic and CG representations
* **CG region**: Atoms are treated as coarse-grained particles

The fix updates the resolution state of each atom based on its current
position. For atoms in the transition region, it calculates a lambda
value (0 to 1) that determines the mixing between atomistic and CG
interactions.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes per-atom data (resolution state and lambda) to binary
restart files, so that simulations can be restarted with the same
resolution assignments.

Restrictions
""""""""""""

This fix is part of the USER-ADRES package. It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`fix adres/thermo <fix_adres_thermo>`, :doc:`pair adres <pair_adres>`,
:doc:`compute adres/stats <compute_adres_stats>`

Default
"""""""

none

