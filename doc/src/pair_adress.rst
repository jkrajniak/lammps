.. index:: pair adres

pair_style adress command
========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style adress cutoff

* adress = style name of this pair style
* cutoff = global cutoff for interactions (distance units)

.. code-block:: LAMMPS

   pair_coeff N M cutoff fix-ID

* N, M = atom types (see :doc:`pair_coeff <pair_coeff>`)
* cutoff = cutoff for this type pair (distance units)
* fix-ID = ID of fix adress/region command (optional)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style adress 10.0
   pair_coeff * * 10.0 adres_region

Description
"""""""""""

The *adres* pair style implements hybrid interactions for Adaptive
Resolution (AdResS) simulations. It applies switching functions based
on the lambda parameter from :doc:`fix adress/region <fix_adres_region>`
to smoothly transition between atomistic and coarse-grained interactions.

The pair style uses the lambda values from the fix adress/region command
to determine how to weight interactions. For pairs of atoms, the effective
lambda is the average of the two atoms' lambda values.

The switching function linearly interpolates between full atomistic
interactions (lambda = 1) and full CG interactions (lambda = 0) based
on the lambda value.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing, shift, table, or tail
corrections. It does support restart and rRESPA.

Restrictions
""""""""""""

This pair style is part of the USER-ADRESS package. It is only enabled
if LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix adress/region <fix_adres_region>`, :doc:`fix adress/thermo <fix_adres_thermo>`,
:doc:`compute adress/stats <compute_adres_stats>`

Default
"""""""

none

