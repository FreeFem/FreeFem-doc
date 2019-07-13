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
   
The release version of the surface finite element is available in FreeFEM 4.2.1. with some examples in examples/3dSurf. 
Abstract about Surface FEM in FreeFEM

* new **meshS** type,  refer to the section :ref:`The type meshS in 3 dimension <meshStype>`  
  
  - the functionalities on the meshS type, it is necessary to load the adapted plugin (for example: load ”msh3”, load ”tetgen”, load ”medit”) 
  - new type of surface mesh: **meshS**
  - the old surface mesh3 object is removed and remplaced by this new type 
  - two possibility to use a meshS: 
    + the domain is described by a 3D surface and defined a meshS then 
    + let a volume mesh (mesh3) and its border domain is a surface that FreeFEM can extract and generate the border surface domain (the surface domain at sens of the full border domain)
	
  - generator of meshS :freefem:`square3`, :freefem:`sphere`, :freefem:`ellipsoide` or composition of FreeFEM function like existing generator of 2D mesh and :freefem:`movemesh23`.
  
It is possible to build a meshS from a mesh3 with the function :freefem:`buildSurface`. Let suppose a declared mesh3, use Th3=buildSurface(Th3) to build the surface domain associated to Th3. 
  - link with the mesh3 type
  
    + by the command meshS ThS=Th3. :freefem:`Gamma`, you can build and manipulate the border mesh independly of the volume mesh Th3 ; such as the surface is describeb by triangle elements and edges border elements in 3d.
    + operator on **meshS** type such as :freefem:`movemesh`, :freefem:`trunc`, :freefem:`change`
    + operator in relation mesh3 / meshS such as :freefem:`extract`, :freefem:`buildSurface` gluing of meshS with the operator :freefem:`+` 
	+ tetg allows to tetrahedralize the interior of the surface mesh with tetgen

* new FESpace with surface finite element type,  refer to the section :ref:`surface Lagrangian Finite Elements <surfacePkLagrange>`
 
 - :freefem:`FESpace` :freefem:`P0` :freefem:`P1` and :freefem:`P2` Lagrange finite elements for the moment


* new psoosibility of definition for a variational form defined by a 

  - :freefem:`problem` 
  - :freefem:`varf` to access to matrix and RHS vector


* visualisation tools 

  - :freefem:`plot` of ffglut, :freefem:`medit`, export at gmsh and vtk format for meshS and surface solutions
  - loading, saving of meshes and solution at freefem's format ;
    + .mesh, .meshb: Mesh Format File of Medit (P. Frey LJLL) 
    + .msh for mesh and .sol data solution at freefem format
    + .msh data file of Gmsh (Mesh generator) (load  "gmsh")
    + vtk format for meshes and solutions (load "iovtk")

.. note::
   **Warning now the surface mesh type is meshS and not mesh3 as before version 4.2.1 **.
   If you want use a surface mesh, the freefem type is **meshS**. For a FreeFEM V3 script working with surface meshes, change mesh3 by meshS.  


===============

CMake
-----

A compilation process using CMake is under development      


.. is available in FreeFEM 4.1  , see the :ref:`compilation process <cmake>`.
