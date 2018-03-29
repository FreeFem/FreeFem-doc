Fruit of a long maturing process, freefem, in its last avatar, FreeFem++ , is a high level integrated development environment (IDE) for numerically solving partial differential equations (PDE) in dimension 2 and 3. It is the ideal tool for teaching the finite element method
but it is also perfect for research to quickly test new ideas or multi-physics and complex applications.

FreeFem++ has an advanced automatic mesh generator, capable of a posteriori mesh adaptation; it has a general purpose elliptic solver interfaced with fast algorithms such as the multi-frontal method UMFPACK, SuperLU, MUMPS . Hyperbolic and parabolic problems are solved by iterative algorithms prescribed by the user with the high level language of FreeFem++ . It has several triangular finite elements, including discontinuous elements. Finally everything is there in FreeFem++ to prepare research quality reports: color display online with zooming and other features and postscript printouts.

This manual is meant for students at Master level, for researchers at any level, and for engineers (including financial engineering) with some understanding of variational methods for partial differential equations.

## Introduction

A partial differential equation is a relation between a function of several variables and its (partial) derivatives. Many problems in physics, engineering, mathematics and even banking are modeled by one or several partial differential equations.

FreeFem++ is a software to solve these equations numerically. As its name implies, it is a free software (see the copyrights for full detail) based on the Finite Element Method; it is not a package, it is an integrated product with its own high level programming language. This software runs on all UNIX OS (with g++ 3.3 or later, and OpenGL) , on Window XP,
Vista and 7,8,10 and on MacOS 10 intel.

Moreover FreeFem++ is highly adaptive. Many phenomena involve several coupled systems, for example: fluid-structure interactions, Lorentz forces for aluminum casting and ocean-atmosphere problems are three such systems. These require different finite element approximations and polynomial degrees, possibly on different meshes. Some algorithms like Schwarz’ domain decomposition method also require data interpolation on multiple meshes within one program. FreeFem++ can handle these difficulties, i.e. arbitrary finite element spaces on arbitrary unstructured and adapted bi-dimensional meshes.

The characteristics of FreeFem++ are:

 * Problem description (real or complex valued) by their variational formulations, with access to the internal vectors and matrices if needed.
 * Multi-variables, multi-equations, bi-dimensional and three-dimensional static or time dependent, linear or nonlinear coupled systems; however the user is required to describe the iterative procedures which reduce the problem to a set of linear problems.
 * Easy geometric input by analytic description of boundaries by pieces; however this part is not a CAD system; for instance when two boundaries intersect, the user must specify the intersection points.
 * Automatic mesh generator, based on the Delaunay-Voronoi algorithm; the inner point density is proportional to the density of points on the boundaries [GEORGE1996](#GEORGE1996).
 * Metric-based anisotropic mesh adaptation. The metric can be computed automatically from the Hessian of any FreeFem++ function [HECHT1998](#HECHT1998).
 * High level user friendly typed input language with an algebra of analytic and finite element functions.
 * Multiple finite element meshes within one application with automatic interpolation of data on different meshes and possible storage of the interpolation matrices.
 * A large variety of triangular finite elements : linear, quadratic Lagrangian elements and more, discontinuous P1 and Raviart-Thomas elements, elements of a non-scalar type, the mini-element,. . . (but no quadrangles).
 * Tools to define discontinuous Galerkin finite element formulations P0, P1dc, P2dc and keywords: jump, mean, intalledges.
 * A large variety of linear direct and iterative solvers (LU, Cholesky, Crout, CG, GMRES, UMFPACK, MUMPS, SuperLU, ...) and eigenvalue and eigenvector solvers (ARPARK) .
 * Near optimal execution speed (compared with compiled `C++` implementations programmed directly).
 * Online graphics, generation of ,.txt,.eps,.gnu, mesh files for further manipulations of input and output data.
 * Many examples and tutorials: elliptic, parabolic and hyperbolic problems, Navier-Stokes flows, elasticity, Fluid structure interactions, Schwarz’s domain decomposition method, eigenvalue problem, residual error indicator, ...
 * A parallel version using MPI

## History

The project has evolved from MacFem, PCfem, written in Pascal. The first C version lead to freefem 3.4; it offered mesh adaptivity on a single mesh only.

A thorough rewriting in `C++` led to freefem+ (freefem+ 1.2.10 was its last release), which included interpolation over multiple meshes (functions defined on one mesh can be used on any other mesh); this software is no longer maintained but still in use because it handles a problem description using the strong form of the PDEs. Implementing the interpolation from one unstructured mesh to another was not easy because it had to be fast and non-diffusive; for each point, one had to find the containing triangle. This is one of the basic problems of computational geometry (see [PREPARATA1985](#PREPARATA1985) for example). Doing it in a minimum number of operations was the challenge. Our implementation is $\mathcal{O}(n log n)$ and based on a quadtree. This version also grew out of hand because of the evolution of the template syntax in `C++`.

We have been working for a few years now on FreeFem++ , entirely re-written again in `C++` with a thorough usage of template and generic programming for coupled systems of unknown size at compile time. Like all versions of freefem it has a high level user friendly input language which is not too far from the mathematical writing of the problems.

The freefem language allows for a quick specification of any partial differential system of equations. The language syntax of FreeFem++ is the result of a new design which makes use of the STL [STROUSTRUP2000](#STROUSTRUP2000), templates and bison for its implementation; more detail can be found in [HECHT2002](#HECHT2002). The outcome is a versatile software in which any new finite element can be included in a few hours; but a recompilation is then necessary. Therefore the library of finite elements available in FreeFem++ will grow with the version number and with the number of users who program more new elements. So far we have discontinuous $P_0$ elements,linear $P_1$ and quadratic $P_2$ Lagrangian elements, discontinuous $P_1$ and Raviart-Thomas elements and a few others like bubble elements.

## References

<a name="GEORGE1996">[GEORGE1996]</a> GEORGE, P. L. et BOROUCHAKI, H. Automatic triangulation. 1996.

<a name="HECHT1998">[HECHT1998]</a> HECHT, F. The mesh adapting software: bamg. INRIA report, 1998, vol. 250, p. 252.

<a name="PREPARATA1985">[PREPARATA1985]</a> PREPARATA, F. P. et SHAMOS, M. I. Computational Geometry Springer-Verlag. New York, 1985.

<a name="STROUSTRUP2000">[STROUSTRUP2000]</a> STROUSTRUP, Bjarne. The C++ programming language. Pearson Education India, 2000.

<a name="HECHT2002">[HECHT2002]</a> HECHT, Frédéric. C++ Tools to construct our user-level language. ESAIM: Mathematical Modelling and Numerical Analysis, 2002, vol. 36, no 5, p. 809-836.
