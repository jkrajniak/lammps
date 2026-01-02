.. index:: fix adress/region

fix adress/region command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID adress/region region-atomistic region-transition region-CG axis

* ID, group-ID are documented in :doc:`fix <fix>` command
* adress/region = style name of this fix command
* region-atomistic = ID of region where atoms are treated atomistically
* region-transition = ID of region where atoms transition between resolutions
* region-CG = ID of region where atoms are treated as coarse-grained
* axis = *x* or *y* or *z* (optional, default: *x*)
  
  Direction along which the transition occurs. Lambda is calculated as the squared normalized distance from the CG region boundary along this axis.

Examples
""""""""

.. code-block:: LAMMPS

   # Transition along x-axis (default)
   region atomistic block 15.0 25.0 0.0 20.0 0.0 20.0
   region transition block 10.0 15.0 0.0 20.0 0.0 20.0
   region cg block 0.0 10.0 0.0 20.0 0.0 20.0
   fix 1 all adress/region atomistic transition cg x

   # Transition along y-axis
   region atomistic block 0.0 20.0 15.0 25.0 0.0 20.0
   region transition block 0.0 20.0 10.0 15.0 0.0 20.0
   region cg block 0.0 20.0 0.0 10.0 0.0 20.0
   fix 2 all adress/region atomistic transition cg y

Description
"""""""""""

This fix manages the resolution regions for Adaptive Resolution (AdResS)
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

Lambda Calculation
""""""""""""""""""

Lambda is calculated using a squared distance formula along the specified
transition axis:

.. math::

   d = \frac{pos - CG_{boundary}}{transition_{width}}

   \lambda = d^2 \quad \text{(clamped to [0, 1])}

Where:
* :math:`pos` = atom position along transition axis
* :math:`CG_{boundary}` = coordinate where CG region ends
* :math:`transition_{width}` = width of transition region along axis
* :math:`d` = normalized distance (0 to 1)
* :math:`\lambda` = lambda value (0 to 1)

The regions should be arranged along the transition axis:
CG Region → Transition Region → Atomistic Region

At the CG boundary, :math:`\lambda = 0`. At the atomistic boundary,
:math:`\lambda = 1`. The squared formula provides a smooth transition
between the two regions.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes per-atom data (resolution state and lambda) to binary
restart files, so that simulations can be restarted with the same
resolution assignments.

Restrictions
""""""""""""

This fix is part of the USER-ADRESS package. It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`fix adress/thermo <fix_adres_thermo>`, :doc:`pair adress <pair_adres>`,
:doc:`compute adress/stats <compute_adres_stats>`

Default
"""""""

none

