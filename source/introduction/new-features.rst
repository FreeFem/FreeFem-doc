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
  - the functionalities on the :freefem:`meshS` type, it is necessary to load the plugin ”msh3”. 
  - generator of meshS :freefem:`square3`, :freefem:`sphere`, :freefem:`ellipsoide` or from a :freefem:`mesh` using the command :freefem:`movemesh23`.
  - the old surface :freefem:`mesh3` object is removed and replaced by :freefem:`meshS` type 
  - two possibilities to define a :freefem:`meshS` :    
    
	+ the considered domain is a 3D surface, so naturally the FreeFEM type is a :freefem:`meshS`
	+ let Th3 be a volume mesh (:freefem:`mesh3`) and its border :math:`\Gamma`. FreeFEM allows to define the volume part with a :freefem:`mesh3` type, and can extract and generate the entire border surface domain with :freefem:`meshS ThS=Th3.Gamma` or specific borders with :freefem:`meshS ThS=extract(Th3,label=llabs)`.  
 
  
    It is possible to build the underlying :freefem:`meshS` from a :freefem:`mesh3` with the function :freefem:`buildSurface`: :freefem:`Th3=buildSurface(Th3)` builds the surface domain associated to the :freefem:`mesh3` Th3. 
  - link with the :freefem:`mesh3` type
  
    + by the command :freefem:`meshS` ThS=Th3. :freefem:`Gamma`
    + operator on :freefem:`meshS` type such as :freefem:`movemeshS`, :freefem:`trunc`, :freefem:`change`...
    + operator in relation :freefem:`mesh3` / :freefem:`meshS` such as :freefem:`extract`, :freefem:`buildSurface`, gluing of meshS with the operator :freefem:`+` 
    + :freefem:`tetg` allows to tetrahedralize the interior of the surface mesh with tetgen

* new FESpace with surface finite element type, see the section :ref:`surface Lagrangian Finite Elements <surfacePkLagrange>`
 
 - :freefem:`FESpace` :freefem:`P0` :freefem:`P1`, :freefem:`P2`, :freefem:`P1b` Lagrange finite elements


* as in the standard 2d or 3d case, the variational problem associated to surface PDE can be defined by using the keywords

  - :freefem:`problem` 
  - :freefem:`varf` to access to matrix and RHS vector
  - available operators are :freefem:`int1d`, :freefem:`int2d`, :freefem:`on`


* visualisation tools 

  - plot with :freefem:`plot` of ffglut, :freefem:`medit` meshes meshS and surface solutions
  - loading, saving of meshes and solution at freefem's format
    
    + ".mesh"  mesh format file of Medit (P. Frey LJLL) 
    + ".msh" for mesh and ".sol" data solution at freefem format
    + ".msh" data file of Gmsh (Mesh generator) (load  "gmsh")
    + vtk format for meshes and solutions (load "iovtk")


.. warning::
   Since the release 4.2.1, the surface :freefem:`mesh3` object (list of vertices and border elements, without tetahedra elements) is remplaced by :freefem:`meshS` type.  For a FreeFEM V3 script working with surface meshes, try to change :freefem:`mesh3` by :freefem:`meshS`. 
   
===============

CMake
-----

A compilation process using CMake is under development      


.. is available in FreeFEM 4.1  , see the :ref:`compilation process <cmake>`.
