==============
Bitmask Module
==============

The central part of this module is the ``bitmasks_module.f90`` file. It contains
the constants that will be used to define on which kind of integer the bitmasks
will be defined.

In the program, when an integer ``X`` is used to represent a bit string (like a determinant
for example), it should be defined as, for example:

.. code-block:: fortran

  use bitmasks
  integer(bit_kind)  :: X


The ``bitmasks_routines.irp.f`` contains helper routines to manipulate bitmassk, like
transforming a bit string to a list of integers for example.

Assumptions
===========

.. Do not edit this section. It was auto-generated from the
.. NEEDED_MODULES_CHILDREN file by the `update_README.py` script.

``bit_kind_shift``, ``bit_kind_size`` and ``bit_kind`` are coherent:

.. code_block:: fortran

  2**bit_kind_shift = bit_kind_size
  bit_kind = bit_kind_size / 8




Needed Modules
==============

.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.

.. image:: tree_dependency.png

* `MO_Basis <http://github.com/LCPQ/quantum_package/tree/master/src/MO_Basis>`_

Needed Modules
==============
.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.


.. image:: tree_dependency.png

* `MO_Basis <http://github.com/LCPQ/quantum_package/tree/master/src/MO_Basis>`_

Documentation
=============
.. Do not edit this section It was auto-generated
.. by the `update_README.py` script.


`bitstring_to_hexa <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L98>`_
  Transform a bit string to a string in hexadecimal format for printing


`bitstring_to_list <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L1>`_
  Gives the inidices(+1) of the bits set to 1 in the bit string


`bitstring_to_str <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L65>`_
  Transform a bit string to a string for printing


`cas_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L220>`_
  Bitmasks for CAS reference determinants. (N_int, alpha/beta, CAS reference)


`cis_ijkl_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L32>`_
  Bitmask to include all possible single excitations from Hartree-Fock


`closed_shell_ref_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L389>`_
  Undocumented


`core_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L374>`_
  Reunion of the inactive, active and virtual bitmasks


`debug_det <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L131>`_
  Subroutine to print the content of a determinant in '+-' notation and
  hexadecimal representation.


`debug_spindet <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L166>`_
  Subroutine to print the content of a determinant in '+-' notation and
  hexadecimal representation.


`full_ijkl_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L12>`_
  Bitmask to include all possible MOs


`generators_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L147>`_
  Bitmasks for generator determinants.
  (N_int, alpha/beta, hole/particle, generator).
  .br
  3rd index is :
  .br
  * 1 : hole     for single exc
  .br
  * 2 : particle for single exc
  .br
  * 3 : hole     for 1st exc of double
  .br
  * 4 : particle for 1st exc of double
  .br
  * 5 : hole     for 2nd exc of double
  .br
  * 6 : particle for 2nd exc of double
  .br


`generators_bitmask_restart <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L103>`_
  Bitmasks for generator determinants.
  (N_int, alpha/beta, hole/particle, generator).
  .br
  3rd index is :
  .br
  * 1 : hole     for single exc
  .br
  * 2 : particle for single exc
  .br
  * 3 : hole     for 1st exc of double
  .br
  * 4 : particle for 1st exc of double
  .br
  * 5 : hole     for 2nd exc of double
  .br
  * 6 : particle for 2nd exc of double
  .br


`hf_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L44>`_
  Hartree Fock bit mask


`i_bitmask_gen <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L400>`_
  Current bitmask for the generators


`inact_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L254>`_
  Bitmasks for the inactive orbitals that are excited in post CAS method


`inact_virt_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L362>`_
  Reunion of the inactive and virtual bitmasks


`index_holes_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L237>`_
  Index of the holes in the generators_bitmasks


`index_particl_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L248>`_
  Index of the holes in the generators_bitmasks


`initialize_bitmask_to_restart_ones <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L3>`_
  Initialization of the generators_bitmask to the restart bitmask


`is_a_1h1p <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L407>`_
  Undocumented


`is_a_2p <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L418>`_
  Undocumented


`is_a_two_holes_two_particles <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L213>`_
  Undocumented


`is_the_hole_in_det <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/find_hole.irp.f#L1>`_
  Undocumented


`is_the_particl_in_det <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/find_hole.irp.f#L29>`_
  Undocumented


`list_act <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L433>`_
  list of active orbitals


`list_inact <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L304>`_
  Undocumented


`list_to_bitstring <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L29>`_
  Returns the physical string "string(N_int,2)" from the array of
  occupations "list(N_int*bit_kind_size,2)


`list_virt <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L305>`_
  Undocumented


`modify_bitmasks_for_hole <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L25>`_
  modify the generators_bitmask in order that one can only excite
  the electrons occupying i_hole


`modify_bitmasks_for_particl <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L60>`_
  modify the generators_bitmask in order that one can only excite
  the electrons to the orbital i_part


`n_act_orb <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L421>`_
  number of active orbitals


`n_cas_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L190>`_
  Number of bitmasks for CAS


`n_core_orb <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L375>`_
  Reunion of the inactive, active and virtual bitmasks


`n_generators_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L70>`_
  Number of bitmasks for generators


`n_inact_orb <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L256>`_
  Bitmasks for the inactive orbitals that are excited in post CAS method


`n_int <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L3>`_
  Number of 64-bit integers needed to represent determinants as binary strings


`n_virt_orb <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L257>`_
  Bitmasks for the inactive orbitals that are excited in post CAS method


`number_of_holes <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L4>`_
  Undocumented


`number_of_holes_verbose <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L432>`_
  Undocumented


`number_of_particles <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L110>`_
  Undocumented


`number_of_particles_verbose <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmask_cas_routines.irp.f#L460>`_
  Undocumented


`print_det <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L149>`_
  Subroutine to print the content of a determinant using the '+-' notation


`print_generators_bitmasks_holes <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L146>`_
  Undocumented


`print_generators_bitmasks_holes_for_one_generator <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L190>`_
  Undocumented


`print_generators_bitmasks_particles <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L168>`_
  Undocumented


`print_generators_bitmasks_particles_for_one_generator <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L213>`_
  Undocumented


`print_spindet <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L182>`_
  Subroutine to print the content of a determinant using the '+-' notation


`print_string <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks_routines.irp.f#L120>`_
  Undocumented


`ref_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L62>`_
  Reference bit mask, used in Slater rules, chosen as Hartree-Fock bitmask


`reunion_of_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L325>`_
  Reunion of the inactive, active and virtual bitmasks


`reunion_of_cas_inact_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L349>`_
  Reunion of the inactive, active and virtual bitmasks


`reunion_of_core_inact_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L337>`_
  Reunion of the inactive, active and virtual bitmasks


`set_bitmask_hole_as_input <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L121>`_
  set the generators_bitmask for the holes
  as the input_bimask


`set_bitmask_particl_as_input <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/modify_bitmasks.irp.f#L96>`_
  set the generators_bitmask for the particles
  as the input_bimask


`unpaired_alpha_electrons <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L409>`_
  Bitmask reprenting the unpaired alpha electrons in the HF_bitmask


`virt_bitmask <http://github.com/LCPQ/quantum_package/tree/master/src/Bitmask/bitmasks.irp.f#L255>`_
  Bitmasks for the inactive orbitals that are excited in post CAS method

