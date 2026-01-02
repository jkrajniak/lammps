.. index:: compute adress/stats

compute adress/stats command
============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID adress/stats fix-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* adress/stats = style name of this compute command
* fix-ID = ID of fix adress/region command

Examples
""""""""

.. code-block:: LAMMPS

   compute 1 all adress/stats adres_region
   thermo_style custom step c_1[1] c_1[2] c_1[3]

Description
"""""""""""

Define a computation that calculates statistics about the Adaptive
Resolution (AdResS) simulation regions. The compute returns a vector
with three values:

* c_ID[1] = number of atoms in atomistic region
* c_ID[2] = number of atoms in transition region
* c_ID[3] = number of atoms in CG region

The compute requires a :doc:`fix adress/region <fix_adres_region>`
command to determine which region each atom belongs to.

Output info
"""""""""""

This compute calculates a global vector of length 3, which can be
accessed by indices 1-3. These values can be used by any command that
uses global vector values from a compute.

The vector values are "intensive".

Restrictions
""""""""""""

This compute is part of the USER-ADRESS package. It is only enabled if
LAMMPS was built with that package. See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix adress/region <fix_adres_region>`, :doc:`fix adress/thermo <fix_adres_thermo>`,
:doc:`pair adress <pair_adres>`

Default
"""""""

none

