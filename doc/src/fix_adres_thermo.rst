.. index:: fix adres/thermo

fix adres/thermo command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID adres/thermo region-ID transition-width keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* adres/thermo = style name of this fix command
* region-ID = ID of region where thermodynamic corrections are applied
* transition-width = width of transition region for force calculations
* zero or more keyword/value pairs may be appended
* keyword = *fix* or *region*

  .. parsed-literal::

     *fix* value = fix-ID
       fix-ID = ID of fix adres/region command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all adres/thermo transition 2.0 fix adres_region
   fix 1 all adres/thermo transition 2.0

Description
"""""""""""

This fix applies thermodynamic corrections (thermodynamic force) to
maintain proper energy and pressure in transition regions of Adaptive
Resolution (AdRes) simulations.

The thermodynamic force is calculated from density gradients and applied
as a correction to maintain proper thermodynamics when atoms transition
between different resolution levels. The force is applied in the
transition region to compensate for the change in resolution.

The *transition-width* parameter determines the width of the region
over which the thermodynamic force is calculated and applied.

If the *fix* keyword is used, the fix will use the specified
:doc:`fix adres/region <fix_adres_region>` to determine which atoms
are in transition regions and their lambda values.

Restrictions
""""""""""""

This fix is part of the USER-ADRES package. It is only enabled if LAMMPS
was built with that package. See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`fix adres/region <fix_adres_region>`, :doc:`pair adres <pair_adres>`,
:doc:`compute adres/stats <compute_adres_stats>`

Default
"""""""

none

