Documentation
=============

The fruit of a long maturing process, **freefem**, in its last avatar, **FreeFEM** , is a high level integrated development environment (IDE) for numerically solving partial differential equations (PDE) in dimension 2 and 3.
It is the ideal tool for teaching the finite element method but it is also perfect for research to quickly test new ideas or multi-physics and complex applications.

**FreeFEM** has an advanced automatic mesh generator, capable of a posteriori mesh adaptation; it has a general purpose elliptic solver interfaced with fast algorithms, such as the multi-frontal method UMFPACK, SuperLU, MUMPS.
Hyperbolic and parabolic problems are solved by iterative algorithms prescribed by the user with the high level language of **FreeFEM**.
It has several triangular finite elements, including discontinuous elements.
Everything is there in **FreeFEM** to prepare research quality reports with online color display, zooming and other features as well as postscript printouts.

This manual is meant for students at a Masters level, for researchers at any level, and for engineers (including financial engineering) with some understanding of variational methods for partial differential equations.

.. rst-class:: fake-title

   Introduction

A partial differential equation is a relation between a function of several variables and its (partial) derivatives.
Many problems in physics, engineering, mathematics and even banking are modeled by one or several partial differential equations.

**FreeFEM** is a software to solve these equations numerically.
As its name implies, it is a free software (see the copyrights for full detail) based on the Finite Element Method; it is not a package, it is an integrated product with its own high level programming language.
This software runs on all UNIX OS (with g++ 3.3 or later, and OpenGL), on Window XP, Vista and 7, 8, 10 and on MacOS 10 intel.

Moreover **FreeFEM** is highly adaptive.
Many phenomena involve several coupled systems.
Fluid-structure interactions, Lorentz forces for aluminum casting and ocean-atmosphere problems are three such systems.
These require different finite element approximations and polynomial degrees, possibly on different meshes.
Some algorithms like the Schwarz’ domain decomposition method also requires data interpolation on multiple meshes within one program.
**FreeFEM** can handle these difficulties, i.e. arbitrary finite element spaces on arbitrary unstructured and adapted bi-dimensional meshes.

The characteristics of **FreeFEM** are:

-  Problem description (real or complex valued) by their variational formulations, with access to the internal vectors and matrices if needed.
-  Multi-variables, multi-equations, bi-dimensional and three-dimensional static or time dependent, linear or nonlinear coupled systems; however the user is required to describe the iterative procedures which reduce the problem to a set of linear problems.
-  Easy geometric input by analytic description of boundaries by pieces; however this part is not a CAD system; for instance when two boundaries intersect, the user must specify the intersection points.
-  Automatic mesh generator, based on the Delaunay-Voronoi algorithm; the inner point density is proportional to the density of points on the boundaries [GEORGE1996]_.
-  Metric-based anisotropic mesh adaptation.
   The metric can be computed automatically from the Hessian of any **FreeFEM** function [HECHT1998]_.
-  High level user friendly typed input language with an algebra of analytic and finite element functions.
-  Multiple finite element meshes within one application with automatic interpolation of data on different meshes and possible storage of the interpolation matrices.
-  A large variety of triangular finite elements: linear, quadratic Lagrangian elements and more, discontinuous P1 and Raviart-Thomas elements, elements of a non-scalar type, the mini-element,... (but no quadrangles).
-  Tools to define discontinuous Galerkin finite element formulations P0, P1dc, P2dc and keywords: jump, mean, intalledges.
-  A large variety of linear direct and iterative solvers (LU, Cholesky, Crout, CG, GMRES, UMFPACK, MUMPS, SuperLU, …) and eigenvalue and eigenvector solvers (ARPARK) .
-  Near optimal execution speed (compared with compiled ``C++`` implementations programmed directly).
-  Online graphics, generation of ,.txt,.eps,.gnu, mesh files for further manipulations of input and output data.
-  Many examples and tutorials: elliptic, parabolic and hyperbolic problems, Navier-Stokes flows, elasticity, fluid structure interactions, Schwarz’s domain decomposition method, eigenvalue problem, residual error indicator, …
-  A parallel version using MPI

.. rst-class:: fake-title

   History

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

.. toctree::

   notations
   meshGeneration
   finiteElement
   visualization
   algorithmsOptimization
   parallelization
   plugins
   developers
   ffddm/index
