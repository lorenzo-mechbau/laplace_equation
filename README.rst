
================
Laplace Equation
================

Test cases for new implementation of domain mappings of elements, nodes and dofs that avoid the GLOBAL_TO_LOCAL_MAP. See Opencmiss iron issue #111.

dark green: internal elements/nodes
light green: ghost elements/nodes
red: boundary nodes

2 Processes, 1 field variable with 2 components, 2x node based dofs

|p2|

2 Processes, 1 field variable with 3 components, 2x node based dofs, 1x element based dofs

|p23|

3 Processes, 1 field variable with 2 components, 2x node based dofs

|p3|

3 Processes, 1 field variable with 3 components, 2x node based dofs, 1x element based dofs

|p33|

4 Processes

|p4|

.. |p2| image:: ./docs/images/PartitionedMesh_4x4b.svg
   :align: middle
   
.. |p23| image:: ./docs/images/PartitionedMesh_4x4c.svg
   :align: middle
   
.. |p3| image:: ./docs/images/PartitionedMesh_4x4-3.svg
   :align: middle
   
.. |p33| image:: ./docs/images/PartitionedMesh_4x4-3a.svg
   :align: middle
   
.. |p4| image:: ./docs/images/PartitionedMesh_4x4-2c.svg
   :align: middle
   
Laplace Equation
================

This example solves the weak form of the standard Laplace equation, 

|standard_laplace_equation| 

using the Galerkin Finite Element method. |phi| is a spatially varying scalar field (e.g. potential). |nabla| and |laplacian| are 'nabla' (divergence operator) and 'laplacian' (laplace operator) respectively and given by,

|nabla_equation|

|laplacian_equation|.

The generalised Laplace equation,

|generalised_laplace_equation|

where |conductivity_tensor| is the rank two, positive definite and symmetric conductivity tensor, can be reduced to the standard Laplace equation by setting |equation| (identity tensor).


Only Dirichlet boundary conditions are prescribed in this example.

.. |standard_laplace_equation| image:: ./docs/images/standard_laplace_equation.svg
   :align: middle
   
.. |phi| image:: ./docs/images/phi.svg
   :align: middle   
   
.. |nabla| image:: ./docs/images/nabla.svg
   :align: middle
   
.. |laplacian| image:: ./docs/images/laplacian.svg
   :align: middle
   
.. |nabla_equation| image:: ./docs/images/nabla_equation.svg
   :align: middle
   
.. |laplacian_equation| image:: ./docs/images/laplacian_equation.svg
   :align: middle
   
.. |generalised_laplace_equation| image:: ./docs/images/generalised_laplace_equation.svg
   :align: middle
   
.. |conductivity_tensor| image:: ./docs/images/conductivity_tensor.svg
   :align: middle
   
.. |equation| image:: ./docs/images/equation.svg 
   :align: middle  
   
.. |dphi_dn| image:: ./docs/images/dphi_dn.svg 
   :align: middle  
  
Building the example
====================

The fortran version of the example can be configured and built with CMake::

  git clone https://github.com/OpenCMISS-Examples/laplace_equation
  mkdir laplace_equation-build
  cd laplace_equation-build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../laplace_equation
  make

This will create the example executable "laplace_equation" in ./src/fortran/ directory.

Running the example
===================

Fortran version::

  cd ./src/fortran/
  ./laplace_equation

Python version::

  cd ./laplace_equation/src/python
  source  /path/to/opencmisslibs/install/.../.../virtualenvironments/oclibs_pyXY_release/bin/activate
  python laplace_equation.py
  
  Note: If the above fails, try the following.
  cd ./laplace_equation/src/python
  PYTHONPATH="/path/to/opencmisslibs/install/.../.../Release/opencmiss.iron" python laplace_equation.py  

Verifying the example
=====================

Results can be visualised by running `visualise.cmgui <./src/fortran/visualise.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.

The following figure shows the finite element mesh (computational domain) and solutions (|phi| and |dphi_dn|).  
 
.. |figure1a| image:: ./docs/images/mesh.svg
   :align: middle
   :width: 250
   :scale: 100
   
.. |figure1b| image:: ./docs/images/solution_phi.svg
   :align: middle  
   :width: 250
   :scale: 100
   
.. |figure1c| image:: ./docs/images/solution_dphi_dn.svg
   :align: middle  
   :width: 250
   :scale: 100   
   
    
|figure1a|  |figure1b|  |figure1c|

Figure 1. (a) Finite element mesh (b) Primary variable solution (c) Secondary variable solution
  
The expected results from this example are available in `expected_results <./src/fortran/expected_results>`_ folder.  
  
Prerequisites
=============

There are no additional input files required for this example as it is self-contained.

License
=======

License applicable to this example is described in `LICENSE <./LICENSE>`_.
