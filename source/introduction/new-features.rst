Version 4: new features
=======================


|
|


Hash Matrix
-----------

A new internal managment of matrix inside the FreeFEM core have been introduced in FreeFEM 4.0, for better performance.


===============

Surface Finite Element
----------------------
.. role:: freefem(code)
   :language: freefem
   
The surface finite element method is available since FreeFEM V4.2.1. Some examples in examples/3dSurf. 
Abstract about Surface FEM in FreeFEM.

* new **meshS** type, refer to the section :ref:`The type meshS in 3 dimension <meshStype>`  
  
  - new type of surface mesh: :freefem:`meshS`
  - the functionalities on the :freefem:`meshS` type, it is necessary to load the adapted plugin (for example: load ”msh3”, load ”tetgen”, load ”medit”) 
  - generator of meshS :freefem:`square3`, :freefem:`sphere`, :freefem:`ellipsoide` or composition of FreeFEM function like existing generator of 2D mesh and :freefem:`movemesh23`.
  - the old surface :freefem:`mesh3` object is removed and remplaced by :freefem:`meshS` type 
  - two possibility to define a meshS:    
    
	+ the considered domain is a 3D surface, so naturally the FreeFEM type is a :freefem:`meshS`
	+ let be Th3 a volume mesh (:freefem:`mesh3`) and its border a surface :math:`\Gamma`. FreeFEM allows to define the volume part with a :freefem:`mesh3` type, and can extract and generate the border surface domain :freefem:`meshS ThS=Th3.Gamma` (the surface domain at sens of the full border domain)
 
  
It is possible to build a meshS from a mesh3 with the function :freefem:`buildSurface`. Let suppose a declared :freefem:`mesh3`, use :freefem:`Th3=buildSurface(Th3)` to build the surface domain associated to Th3. 
  - link with the :freefem:`mesh3` type
  
    + by the command meshS ThS=Th3. :freefem:`Gamma`,
    + operator on :freefem:`meshS` type such as :freefem:`movemesh`, :freefem:`trunc`, :freefem:`change`...
    + operator in relation mesh3 / meshS such as :freefem:`extract`, :freefem:`buildSurface`, gluing of meshS with the operator :freefem:`+` 
	+ :freefem:`tetg` allows to tetrahedralize the interior of the surface mesh with tetgen

* new FESpace with surface finite element type,  refer to the section :ref:`surface Lagrangian Finite Elements <surfacePkLagrange>`
 
 - :freefem:`FESpace` :freefem:`P0` :freefem:`P1`, :freefem:`P2`, :freefem:`P1b` Lagrange finite elements


* define the variational problem associated with surface PDE

  - :freefem:`problem` 
  - :freefem:`varf` to access to matrix and RHS vector
  - available operators are :freefem:`int1d`, :freefem:`int2d`, :freefem:`on`


* visualisation tools 

  - plot with :freefem:`plot` of ffglut, :freefem:`medit` meshes meshS and surface solutions
  - loading, saving of meshes and solution at freefem's format ;
    + .mesh, .meshb: mesh format file of Medit (P. Frey LJLL) 
    + .msh for mesh and .sol data solution at freefem format
    + .msh data file of Gmsh (Mesh generator) (load  "gmsh")
    + vtk format for meshes and solutions (load "iovtk")

.. note::
   If you want use a surface mesh, the freefem type is :freefem:`meshS`. For a FreeFEM V3 script working with surface meshes, try to change :freefem:`mesh3` by :freefem:`meshS`.  
   

===============

CMake
-----

A compilation process using CMake is under development      


.. is available in FreeFEM 4.1  , see the :ref:`compilation process <cmake>`.
