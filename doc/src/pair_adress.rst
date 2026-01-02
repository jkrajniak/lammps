.. index:: pair adres

pair_style adress command
========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style adress cutoff atomistic_style atomistic_args ... cg_style cg_args ...

* adress = style name of this pair style
* cutoff = global cutoff for interactions (distance units)
* atomistic_style = pair style for atomistic interactions (e.g., lj/cut)
* atomistic_args = arguments for atomistic pair style
* cg_style = pair style for coarse-grained interactions (e.g., lj/cut)
* cg_args = arguments for CG pair style

.. code-block:: LAMMPS

   pair_coeff * * args ... fix fix-ID

* * * = wildcard for all atom types (see :doc:`pair_coeff <pair_coeff>`)
* args = arguments for both atomistic and CG pair styles
* fix = keyword to specify fix ID
* fix-ID = ID of fix adress/region command (optional)

Examples
""""""""

.. code-block:: LAMMPS

   # Atomistic: lj/cut with cutoff 2.5, CG: lj/cut with cutoff 5.0
   pair_style adress 10.0 lj/cut 2.5 lj/cut 5.0
   pair_coeff * * 1.0 1.0 2.5 fix adres_region

   # Different pair styles for atomistic and CG
   pair_style adress 10.0 lj/cut 2.5 lj/cut/coul/cut 5.0
   pair_coeff * * 1.0 1.0 2.5 fix adres_region

Description
"""""""""""

The *adress* pair style implements hybrid interactions for Adaptive
Resolution (AdResS) simulations. It interpolates between atomistic and
coarse-grained pair interactions based on the lambda parameter from
:doc:`fix adress/region <fix_adres_region>`.

The pair style requires two sub-styles:
* **atomistic_style**: The pair style used for full atomistic interactions
* **cg_style**: The pair style used for coarse-grained interactions

Forces are interpolated using the formula:
:math:`F_i = \lambda_i \cdot F_{at,i} + (1 - \lambda_i) \cdot F_{cg,i}`

Where:
* :math:`\lambda_i` = lambda value for atom i (from fix adress/region)
* :math:`F_{at,i}` = force on atom i from atomistic pair style
* :math:`F_{cg,i}` = force on atom i from CG pair style

The same coefficients are applied to both sub-styles via the
:doc:`pair_coeff <pair_coeff>` command. The global cutoff should be
set to the maximum of the atomistic and CG cutoffs.

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

