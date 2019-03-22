.. role:: freefem(code)
  :language: freefem

.. _globalVariables:

Global variables
================

area
----

Area of the current triangle.

.. code-block:: freefem
   :linenos:

   fespace Vh0(Th, P0);
   Vh0 A = area;

ARGV
----

Array that contains all the command line arguments.

.. code-block:: freefem
   :linenos:

   for (int i = 0; i < ARGV.n; i++)
       cout << ARGV[i] << endl;

See :ref:`Command line arguments example <exampleCommandLineArguments>` for a complete example.

BoundaryEdge
------------

Return 1 if the current edge is on a boundary, 0 otherwise.

.. code-block:: freefem
   :linenos:

   real B = int2d(Th)(BoundaryEdge);

CG
--

Conjugate gradient solver.

Usable in :ref:`problem <typeProblem>` and :ref:`solve <typeSolve>` definition

.. code-block:: freefem
   :linenos:

   problem Laplacian (U, V, solver=CG) = ...

Or in :ref:`matrix <typeMatrix>` construction

.. code-block:: freefem
   :linenos:

   matrix A = vLaplacian(Uh, Uh, solver=CG);

Or in :ref:`set function <functionSet>`

.. code-block:: freefem
   :linenos:

   set(A, solver=CG);

Cholesky
--------

Cholesky solver.

Crout
-----

Crout solver.

edgeOrientation
---------------

Sign of :math:`i-j` if the current edge is :math:`[q_i, q_j]`.

.. code-block:: freefem
   :linenos:

   real S = int1d(Th, 1)(edgeOrientation);

false
-----

False boolean value.

.. code-block:: freefem
   :linenos:

   bool b = false;

.. _globalVariablesGMRES:

GMRES
-----

GMRES solver (Generalized minimal residual method).

hTriangle
---------

Size of the current triangle.

.. code-block:: freefem
   :linenos:

   fespace Vh(Th, P0);
   Vh h = hTriangle;

include
-------

Include an :ref:`external library <externalLibraries>`.

.. code-block:: freefem
   :linenos:

   include "iovtk"

InternalEdge
------------

Return 0 if the current edge is on a boundary, 1 otherwise.

.. code-block:: freefem
   :linenos:

   real I = int2d(Th)(InternalEdge);

label
-----

Label number of a boundary if the current point is on a boundary, 0 otherwise.

.. code-block:: freefem
   :linenos:

   int L = Th(xB, yB).label;

lenEdge
-------

Length of the current edge.

For an edge :math:`[q_i, g_j]`, return :math:`|q_i-q_j|`.

.. code-block:: freefem
   :linenos:

   real L = int1d(Th, 1)(lenEdge);

load
----

Load a script.

.. code-block:: freefem
   :linenos:

   load "Element_P3"

LU
--

LU solver.

N
-

Outward unit normal at the current point if it is on a curve defined by a border.
:freefem:`N.x, N.y, N.z` are respectively the :math:`x`, :math:`y` and :math:`z` components of the normal.

.. code-block:: freefem
   :linenos:

   func Nx = N.x;
   func Ny = N.y;
   func Nz = N.z;

nTonEdge
--------

Number of adjacent triangles of the current edge.

.. code-block:: freefem
   :linenos:

   real nTE = int2d(Th)(nTonEdge);

nuEdge
------

Index of the current edge in the triangle.

.. code-block:: freefem
   :linenos:

   real nE = int2d(Th)(nuEdge);

nuTriangle
----------

Index of the current triangle.

.. code-block:: freefem
   :linenos:

   fespace Vh(Th, P0);
   Vh n = nuTriangle;

P
-

Current point.

.. code-block:: freefem
   :linenos:

   real cx = P.x;
   real cy = P.y;
   real cz = P.z;

pi
--

Pi = 3.14159.

.. code-block:: freefem
   :linenos:

   real Pi = pi;

This is a real value.

region
------

Region number of the current point. If the point is outside, then :freefem:`region == notaregion` where :freefem:`notaregion` is a **FreeFEM** integer constant.

.. code-block:: freefem
   :linenos:

   int R = Th(xR, yR).region;

sparsesolver
------------

Sparse matrix solver.

true
----

True boolean value.

.. code-block:: freefem
   :linenos:

   bool b = true;

verbosity
---------

Verbosity level.

.. code-block:: freefem
   :linenos:

   int Verbosity = verbosity;
   verbosity = 0;

0 = nothing, 1 = little information, 10 = a lot of information, …

This is an integer value.

version
-------

**FreeFEM** version.

.. code-block:: freefem
   :linenos:

   cout << version << endl;

volume
------

Volume of the current tetrahedra.

.. code-block:: freefem
   :linenos:

   fespace Vh0(Th, P0);
   Vh0 V = volume;

x
-

The :math:`x` coordinate at the current point.

.. code-block:: freefem
   :linenos:

   real CurrentX = x;

This is a real value.

y
-

The :math:`y` coordinate at the current point.

.. code-block:: freefem
   :linenos:

   real CurrentY = y;

This is a real value.

z
-

The :math:`z` coordinate at the current point.

.. code-block:: freefem
   :linenos:

   real CurrentZ = z;

This is a real value.
