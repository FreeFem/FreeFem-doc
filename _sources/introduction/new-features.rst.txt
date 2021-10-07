Version 4.5: new features
=========================

|
|

Release, binaries packages 
--------------------------

* Since the version 4.5, the FreeFEM binary packages provides with a compiled PETSc library.
* FreeFEM is now interfaced with ParMmg.

New meshes and FEM border 
-------------------------
.. role:: freefem(code)
   :language: freefem

After Surface FEM, Line FEM is possible with a new mesh type :freefem:`meshL`, :freefem:`P0` :freefem:`P1` :freefem:`P2` :freefem:`P1dc` FE, basic FEM, mesh generation.
This new development allows to treat a 1d problem, such as a problem described on a 3d curve.

=======

Abstract about Line FEM in FreeFEM.

* new **meshL** type, refer to the section :ref:`The type meshL in 3 dimension <meshStype>`  
  
  - new type of surface mesh: :freefem:`meshL`
  - the functionalities on the :freefem:`meshL` type, it is necessary to load the plugin ”msh3”. 
  - generator of meshL :freefem:`segment`, define multi :freefem:`border` and :freefem:`buildmesh` function.
  - basic transformation are avalaible: :freefem:`movemesh`, :freefem:`trunc`, :freefem:`extract`, :freefem:`checkmesh`, :freefem:`change`, :freefem:`AddLayers`, glue of :freefem:`meshL`.
  
    It is possible to build the underlying :freefem:`meshL` from a :freefem:`meshS` with the function :freefem:`buildBdMesh`: :freefem:`ThS=buildBdMesh(ThS)` builds the boundary domain associated to the :freefem:`meshS` ThS and extract it by the command :freefem:`meshL` ThL=ThS. :freefem:`Gamma`. 
  

* new finite element space with curve finite element type
 
 - :freefem:`FESpace` :freefem:`P0` :freefem:`P1`, :freefem:`P2`, :freefem:`P1dc` Lagrange finite elements and possible to add a custumed finite element with the classical method (like a plugin).

* as in the standard 2d, 3d, surface 3d case, the variational problem associated to surface PDE can be defined by using the keywords

  - :freefem:`problem` 
  - :freefem:`varf` to access to matrix and RHS vector
  - available operators are :freefem:`int1d`, :freefem:`on` and the operator :freefem:`int0d` to define a Neumann boundary condition 


* visualisation tools 

  - plot with :freefem:`plot` of ffglut, :freefem:`medit` meshes meshL and solutions
  - 2d or 3d view, with in 3d the option to visualize the elememt Normals at element (touch 'T') and the deformed domain according to it (touch '2').
  - loading, saving of meshes and solution at FreeFEM's format
    
    + ".mesh"  mesh format file of Medit (P. Frey LJLL) 
    + ".msh" for mesh and ".sol" data solution at freefem format
    + ".msh" data file of Gmsh (Mesh generator) (load  "gmsh")
    + vtk format for meshes and solutions (load "iovtk" and use the ".vtu" extension)

 
===============


Boundary Element Method
-----------------------

Allows to define and solve a 2d/3d BEM formulation and rebuild the associated potential.
The document is in construction.
