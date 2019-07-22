History 
=======


  The project has evolved from MacFem, PCfem, written in Pascal.
  The first C version lead to ``freefem 3.4``; it offered mesh adaptivity on a single mesh only.

  A thorough rewriting in ``C++`` led to ``freefem+`` (``freefem+`` 1.2.10 was its last release), which included interpolation over multiple meshes (functions defined on one mesh can be used on any other mesh); this software is no longer maintained but is still in use because it handles a problem description using the strong form of the PDEs.
  Implementing the interpolation from one unstructured mesh to another was not easy because it had to be fast and non-diffusive; for each point, one had to find the containing triangle.
  This is one of the basic problems of computational geometry (see [PREPARATA1985]_ for example).
  Doing it in a minimum number of operations was the challenge.
  Our implementation is :math:`\mathcal{O}(n log n)` and based on a quadtree.
  This version also grew out of hand because of the evolution of the template syntax in ``C++``.

  We have been working for a few years now on **FreeFEM** , entirely re-written again in ``C++`` with a thorough usage of template and generic programming for coupled systems of unknown size at compile time.
  Like all versions of ``freefem``, it has a high level user friendly input language which is not too far from the mathematical writing of the problems.

  The ``freefem`` language allows for a quick specification of any partial differential system of equations.
  The language syntax of **FreeFEM** is the result of a new design which makes use of the STL [STROUSTRUP2000]_, templates, and bison for its implementation; more details can be found in [HECHT2002]_.
  The outcome is a versatile software in which any new finite elements can be included in a few hours; but a recompilation is then necessary.
  Therefore the library of finite elements available in **FreeFEM** will grow with the version number and with the number of users who program more new elements.
  So far we have discontinuous :math:`P_0` elements,linear :math:`P_1` and quadratic :math:`P_2` Lagrangian elements, discontinuous :math:`P_1` and Raviart-Thomas elements and a few others like bubble elements.

|
|
|

.. rst-class:: fake-title

   The development of FreeFEM through more than 30 years

|

**1987**

  MacFem/PCFem the old ones (O. Pironneau in Pascal) no free.
  

**1992** 

  `FreeFem <http://www3.freefem.org/ff++/freefem/fraold.htm>`__ rewrite in C++ (P1,P0 one mesh ) O. Pironneau, D. Bernardi, F.Hecht (mesh adaptation , bamg) , C. Prudhomme .
  

**1996** 

  `FreeFem+ <http://www3.freefem.org/ff++/freefem/index.html>`__ rewrite in C++ (P1,P0 more mesh) O. Pironneau, D. Bernardi, F.Hecht (algebra of function).
  

**1998**

  FreeFem++ rewrite with an other finite element kernel and an new language F. Hecht, O. Pironneau, K.Ohtsuka.  

**1999**

  FreeFem 3d (S. Del Pino),  a fist 3d version base on fictitious domaine method.
  

**2008**

  FreeFem++ v3 use a new finite element kernel multidimensionnels: 1d,2d,3d...
  
**2014** 

  FreeFem++ v3.34 parallel version
  
**2017** 

  FreeFem++ v3.57 parallel version  

**2018**

  FreeFem++ v4: New matrix type, Surface element, New Parallel tools ...
  

  
  
  
  
  