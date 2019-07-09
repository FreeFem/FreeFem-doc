New features
============

Hash Matrix
-----------

A new internal managment of matrix inside the FreeFEM core have been introduced in FreeFEM 4.0, for better performance.

Surface Finite Element
----------------------

The release version of the surface finite element is available in FreeFEM 4.2.1. with some examples in examples/3dSurf. 
Abstract about Surface FEM in FreeFEM

**meshS** type:
- new type of surface mesh: **meshS**, the old surface mesh3 object is removed and remplaced by this new type
- generator of meshS   square3 sphere ellipsoide   
- a volume mesh is a mesh3 and can contain a meshS describing is surface domain (the surface domain at sens of the full border domain)
It is possible to build a meshS from a mesh3 with the function buildSurface(). Let suppose a declared mesh3, use Th3=buildSurface(Th3) to build the surface domain associated to Th3. 
By the command meshS ThS=Th3.gamma() you can manipulate the created surface mesh independly of the volume mesh




.. note::
   If you want use a surface mesh, the freefem type is **mesh**. For a FreeFEM V3 script working with surface mesh, change mesh3 by meshS.  

CMake
-----

A compilation process using CMake is under development      


.. is available in FreeFEM 4.1  , see the :ref:`compilation process <cmake>`.
