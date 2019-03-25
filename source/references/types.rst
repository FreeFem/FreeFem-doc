.. role:: freefem(code)
  :language: freefem

.. role:: cpp(code)
 :language: cpp

.. _referenceTypes:

Types
=====

Standard types
--------------

int
~~~

Integer value (equivalent to :cpp:`long` in ``C++``).

.. code-block:: freefem
   :linenos:

   int i = 0;

bool
~~~~

Boolean value.

.. code-block:: freefem
   :linenos:

   bool b = true;

.. tip:: The result of a comparison is a boolean

   :freefem:`bool b = (1 < 2);`

real
~~~~

Real value (equivalent to :cpp:`double` in ``C++``).

.. code-block:: freefem
   :linenos:

   real r = 0.;

complex
~~~~~~~

Complex value (equivalent to two :cpp:`double` or :cpp:`complex<double>` in ``C++``).

.. code-block:: freefem
   :linenos:

   complex c = 0. + 1i;

The imaginary number :math:`i` is defined as ``1i``

.. tip:: Example

   .. code-block:: freefem
      :linenos:

      complex a = 1i, b = 2 + 3i;
      cout << "a + b = " << a + b << endl;
      cout << "a - b = " << a - b << endl;
      cout << "a*b = " << a*b << endl;
      cout << "a/b = " << a/b << endl;

   The output of this script is:

   .. code-block:: bash

      a + b = (2,4)
      a - b = (-2,-2)
      a*b = (-3,2)
      a/b = (0.230769,0.153846)

.. note:: See :ref:`Complex example <exampleComplex>` for a detailed example.

string
~~~~~~

String value.

.. code-block:: freefem
   :linenos:

   string s = "this is a string";

.. note:: :freefem:`string` value is enclosed within double quotes.

Other types can be concatenate to a string, like:

.. code-block:: freefem
   :linenos:

   int i = 1;
   real r = 1.;
   string s = "the int i = " + i +", the real r = " + r + ", the complex z = " + (1. + 1i);

To append a string in a string at position 4:

.. code-block:: freefem
   :linenos:

   s(4:3) = "+++";

To copy a substring in an other string:

.. code-block:: freefem
   :linenos:

   string s2 = s1(5:10);

See :ref:`String Example <exampleString>` for a complete example.

Mesh design
-----------

border
~~~~~~

Border type.

.. code-block:: freefem
   :linenos:

   border b(t=0., 1.){x=cos(2.*pi*t); y=sin(2.*pi*t);}

Define the 2D geometrical border in parametric coordinates.

.. note:: Label

   A label can be defined with the border:

   .. code-block:: freefem
      :linenos:

      border b(t=0., 1.){x=cos(2.*pi*t); y=sin(2.*pi*t); label=1;}

.. note:: Inner variable

   An inner variable can be defined inside a border:

   .. code-block:: freefem
      :linenos:

      border b(t=0., 1.){real tt=2.*pi*t; x=cos(tt); y=sin(tt);}

.. note:: From vector

   A border can be defined from two vectors using :freefem`P.x` and :freefem:`P.y`:

   .. code-block:: freefem
      :linenos:

      border b(t=0, vectorX.n-1){x=vectorX[t]; P.x=vectorY[t];}

mesh
~~~~

2D Mesh type (see :ref:`Mesh Generation <meshGeneration>`).

.. code-block:: freefem
   :linenos:

   mesh Th;

mesh3
~~~~~

3D mesh type (see :ref:`Mesh Generation <meshGeneration>`).

.. code-block:: freefem
   :linenos:

   mesh3 Th;

Finite element space design
---------------------------

fespace
~~~~~~~

Finite element space type (see :ref:`Finite Element <finiteElement>`).

.. code-block:: freefem
   :linenos:

   fespace Uh(Th, P1);
   fespace UPh(Th, [P2, P2, P1]);

A finite element space is based on a mesh (``Th``) with an element definition, scalar (:freefem:`P1`) or vector (:freefem:`[P2, P2, P1]`).

**Available finite element space:**

Generic:

-  :freefem:`P0 / P03d`
-  :freefem:`P0Edge`
-  :freefem:`P1 / P13d`
-  :freefem:`P1dc`
-  :freefem:`P1b / P1b3d`
-  :freefem:`P1bl / P1bl3d`
-  :freefem:`P1nc`
-  :freefem:`P2 / P23d`
-  :freefem:`P2b`
-  :freefem:`P2dc`
-  :freefem:`P2h`
-  :freefem:`RT0 / RT03d`
-  :freefem:`RT0Ortho`
-  :freefem:`Edge03d`

Using *Element_P3*:

-  :freefem:`P3`

Using *Element_P3dc*:

-  :freefem:`P3dc`

Using *Element_P4*:

-  :freefem:`P4`

Using *Element_P4dc*:

-  :freefem:`P4dc`

Using *Element_PkEdge*:

-  :freefem:`P1Edge`
-  :freefem:`P2Edge`
-  :freefem:`P3Edge`
-  :freefem:`P4Edge`
-  :freefem:`P5Edge`

Using *Morlay*:

-  :freefem:`P2Morley`

Using *HCT*:

-  :freefem:`HCT`

Using *BernardiRaugel*:

-  :freefem:`P2BR`

Using *Element_Mixte*:

-  :freefem:`RT1`
-  :freefem:`RT1Ortho`
-  :freefem:`RT2`
-  :freefem:`RT2Ortho`
-  :freefem:`BDM1`
-  :freefem:`BDM1Ortho`

Using *Element_Mixte3d*:

-  :freefem:`Edge13d`
-  :freefem:`Edge23d`

Using *Element_QF*:

-  :freefem:`FEQF`

A finite element function is defined as follow:

.. code-block:: freefem
   :linenos:

   fespace Uh(Th, P1);
   Uh u;

   fespace UPh(Th, [P2, P2, P1]);
   UPh [Ux, Uy, p];

.. _typeMacroDesign:

Macro design
------------

macro
~~~~~

Macro type.

.. code-block:: freefem
   :linenos:

   macro vU() [Ux, Uy] //
   macro grad(u) [dx(u), dy(u)] //

Macro ends with :freefem:`//`.

.. note:: Macro concatenation

   You can use the C concatenation operator :cpp:`##` inside a macro using :freefem:`#`.

   If `Ux` and `Uy` are defined as finite element function, you can define:

   .. code-block:: freefem
      :linenos:

      macro Grad(U) [grad(U#x), grad(U#y)] // End of macro

See :ref:`Macro example <exampleMacro>`

NewMacro / EndMacro
~~~~~~~~~~~~~~~~~~~

.. warning:: In developement - Not tested

Set and end a macro

.. code-block:: freefem
   :linenos:

   NewMacro grad(u) [dx(u), dy(u)] EndMacro

IFMACRO
~~~~~~~

Check if a macro exists and check its value.

.. code-block:: freefem
   :linenos:

   IFMACRO(AA) //check if macro AA exists
   ...
   ENDIFMACRO

   IFMACRO(AA, tt) //check if amcro exists and is equall to tt
   ...
   ENDIFMACRO

ENDIFMACRO
~~~~~~~~~~

Functions design
----------------

func
~~~~

Function type.

Function without parameters (:math:`x`, :math:`y` and :math:`z` are implicitly considered):

.. code-block:: freefem
   :linenos:

   func f = x^2 + y^2;

.. note:: Function’s type is defined by the expression’s type.

Function with parameters:

.. code-block:: freefem
   :linenos:

   func real f (real var){
       return x^2 + y^2 + var^2;
   }

.. _typeElementaryFunctions:

Elementary functions
~~~~~~~~~~~~~~~~~~~~

Class of basic functions (polynomials, exponential, logarithmic, trigonometric, circular) and the functions obtained from those by the four arithmetic operations

.. math::
   f(x) + g(x),\, f(x) - g(x),\, f(x)g(x),\, f(x)/g(x)

and by composition :math:`f(g(x))`, each applied a finite number of times.

In **FreeFEM**, all elementary functions can thus be created.
The derivative of an elementary function is also an elementary function; however, the indefinite integral of an elementary function cannot always be expressed in terms of elementary functions.

See :ref:`Elementary function example <exampleElementaryFunction>` for a complete example.

Random functions
~~~~~~~~~~~~~~~~

**FreeFEM** includes the `Mersenne Twister <http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html>`__ random number generator.
It is a very fast and accurate random number generator of period :math:`2^{219937}-1`.

See :ref:`randint32() <referenceRandInt32>`, :ref:`randint31() <referenceRandInt31>`, :ref:`randreal1() <referenceRandReal1>`, :ref:`randreal2() <referenceRandReal2>`, :ref:`randreal3() <referenceRandReal3>`, :ref:`randres53() <referenceRandRes53>`, :ref:`randinit(seed) <referenceRandInit>`.

In addition, the ``ffrandom`` plugin interface :freefem:`random`, :freefem:`srandom` and :freefem:`srandomdev` functions of the Unix ``libc`` library.
The range is :math:`0 -- 2^{31}-1`.

.. note:: If :freefem:`srandomdev` is not defined, a seed based on the current time is used.

:freefem:`gsl` plugin equally allows usage of all random functions of the ``gsllib``, see :ref:`gsl external library <referenceFFGSLAWK>`.

FE-functions
~~~~~~~~~~~~

Finite element functions are also constructed like elementary functions by an arithmetic formula involving elementary functions.

The difference is that they are evaluated at declaration time and **FreeFEM** stores the array of its values at the places associated with he degree of freedom of the finite element type.
By opposition, elementary functions are evaluated only when needed.
Hence FE-functions are not defined only by their formula but also by the mesh and the finite element which enter in their definitions.

If the value of a FE-function is requested at a point which is not a degree of freedom, an interpolation is used, leading to an interpolation error, while by contrast, an elementary function can be evaluated at any point exactly.

.. code-block:: freefem
   :linenos:

   func f = x^2*(1+y)^3 + y^2;
   mesh Th = square(20, 20, [-2+4*x, -2+4*y]); // ]-2, 2[^2
   fespace Vh(Th, P1);
   Vh fh=f; //fh is the projection of f to Vh (real value)
   func zf = (x^2*(1+y)^3 + y^2)*exp(x + 1i*y);
   Vh<complex> zh = zf; //zh is the projection of zf to complex value Vh space

The construction of :freefem:`fh = f` is explained in :ref:`Finite Element <finiteElement>`.

.. warning:: The :freefem:`plot` command only works for real or complex FE-functions, not for elementary functions.

.. _typeProblemDesign:

Problem design
--------------

.. _typeProblem:

problem
~~~~~~~

Problem type.

.. code-block:: freefem
   :linenos:

   problem Laplacian (u, uh) = ...

**FreeFEM** needs the variational form in the problem definition.

In order to solve the problem, just call:

.. code-block:: freefem
   :linenos:

   Laplacian;

.. note:: Solver

   A solver can be specified in the problem definition:

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, solver=CG) = ...

   The default solver is :freefem:`sparsesolver` or :freefem:`LU` if any direct sparse solver is available.

   Solvers are listed in the :ref:`Global variables <globalVariables>` section.

.. note:: Stop test

   A criterion :math:`\varepsilon` can be defined for iterative methods, like CG for example:

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, solver=CG, eps=1.e-6) = ...

   If :math:`\varepsilon>0`, the stop test is:

   .. math::
      ||Ax-b|| < \varepsilon

   Else, the stop test is:

   .. math::
      ||Ax-b|| < \frac{|\varepsilon|}{||Ax_0-b||}

.. note:: Reconstruction

   The keyword :freefem:`init` controls the reconstruction of the internal problem matrix.

   If :freefem:`init` is set to :freefem:`false` or :freefem:`0`, the matrix is reconstructed et each problem calls (or after a mesh modification), else the previously constructed matrix is used.

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, init=1) = ...

.. note:: Preconditioning

   A preconditioner can be specified in the problem definition:

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, precon=P) = ...

   The preconditioning function must have a prototype like:

   .. code-block:: freefem
      :linenos:

      func real[int] P(real[int] &xx);

.. note:: "*Très grande valeur*"

   The “*Très grand valeur*" :freefem:`tgv` (or *Terrible giant value*) used to implement the Dirichlet conditions can be modified in the problem definition:

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, tgv=1e30) = ...

   Refere to :ref:`Problem definition <problemDefinition>` for a description of the Dirichlet condition implementation.

.. note:: Pivot tolerance

   The tolerance of the pivot in :freefem:`UMFPACK`, :freefem:`LU`, :freefem:`Crout`, :freefem:`Cholesky` factorization can be modified in the problem definition:

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, solver=LU, tolpivot=1e-20) = ...

.. note:: :freefem:`UMFPACK`

   Two specific parameters for the :freefem:`UMFPACK` can be modifed:

   * Tolerance of the pivot sym
   * strategy

   .. code-block:: freefem
      :linenos:

      problem Laplacian(u, uh, solver=LU, tolpivotsym=1e-1, strategy=0) = ...

   Refer to the `UMFPACK website <http://faculty.cse.tamu.edu/davis/research.html>`__ for more informations.

.. note:: :freefem:`dimKrylov`

   Dimension of the Krylov space

Usage of :freefem:`problem` is detailled in the :ref:`tutorials <tutorial>`.

.. _typeSolve:

solve
~~~~~

Solve type.

Identical to :ref:`problem <typeProblem>` but automatically solved.

Usage of :freefem:`solve` is detailled in the :ref:`tutorials <tutorial>`.

.. _typeVarf:

varf
~~~~

Variational form type.

.. code-block:: freefem
   :linenos:

   varf vLaplacian (u, uh) = ...

Directly define a variational form.

This is the other way to define a problem in order to directly manage matrix and right hang side.

Usage of :freefem:`varf` is detailed in the :ref:`tutorial <tutorial>`.

Array
-----

An array stores multiple objects, and there are 2 kinds of arrays:

-  the first is similar to vector, i.e. array with integer indices
-  the second is array with string indices

In the first case, the size of the array must be known at execution time, and implementation is done with the :cpp:`KN<>` class and all the vector operator of :cpp:`KN<>` are implemented.

Arrays can be set like in Matlab or Scilab with the operator ``::``, the array generator of ``a:c`` is equivalent to ``a:1:c``, and the array set by ``a:b:c`` is set to size :math:`\lfloor |(b-a)/c|+1 \rfloor` and the value :math:`i` is set by :math:`a + i (b-a)/c`.

There are :freefem:`int,real, complex` array with, in the third case, two operators (:freefem:`.im`, :freefem:`.re`) to generate the real and imaginary real array from the complex array (without copy).

.. note:: Quantiles are points taken at regular intervals from the cumulative distribution function of a random variable.
   Here the array values are random.

   This statistical function :freefem:`a.quantile(q)` computes :math:`v` from an array :math:`a` of size :math:`n` for a given number :math:`q\in ]0,1[` such that:

   .. math::
      \#\{ i / a[i] < v \} \sim q*n

   it is equivalent to :math:`v = a[q*n]` when the array :math:`a` is sorted.

For example, to declare, fill and display an array of :freefem:`real` of size ``n``:

.. code-block:: freefem
   :linenos:

   int n = 5;
   real[int] Ai(n);
   for (int i = 0; i < n; i++)
       Ai[i] = i;
   cout << Ai << endl;

The output of this script is:

.. code-block:: bash

   5
         0   1   2   3   4

See the :ref:`Array example <exampleArray>` for a complete example.

Array index
~~~~~~~~~~~

Array index can be int or string:

.. code-block:: freefem
   :linenos:

   real[int] Ai = [1, 1, 0, 0];
   real[string] As = [1, 1, 0, 0];

Array size
~~~~~~~~~~

The size of an array is obtained using the keyword ``n``:

.. code-block:: freefem
   :linenos:

   int ArraySize = Ai.n;

Array sort
~~~~~~~~~~

To sort an array:

.. code-block:: freefem
   :linenos:

   Ai.sort;

Double array
~~~~~~~~~~~~

A double array (matrix) can be defined using two indexes:

.. code-block:: freefem
   :linenos:

   real[int, int] Aii = [[1, 1], [0, 0]];

The two sizes are obtained using the keywords ``n`` and ``m``:

.. code-block:: freefem
   :linenos:

   int ArraySize1 = Aii.n;
   int ArraySize2 = Aii.m;

The minimum and maximum values of an array (simple or double) can be obtained using:

.. code-block:: freefem
   :linenos:

   real ArrayMin = Aii.min;
   real ArrayMax = Aii.max;

Th minimum and maximum position of an array can be obtained using:

.. code-block:: freefem
   :linenos:

   int mini = Aii.imin;
   int minj = Aii.jmin;

   int maxi = Aii.imax;
   int maxj = Aii.jmax;

.. tip:: An array can be obtained from a finite element function using:

   .. code-block:: freefem
      :linenos:

      real[int] aU = U[];

   where :freefem:`U` is a finite element function.

Array of FE functions
~~~~~~~~~~~~~~~~~~~~~

It is also possible to make an array of FE functions, with the same syntax, and we can treat them as vector valued function if we need them.

The syntax for space or vector finite function is

.. code-block:: freefem
   :linenos:

   int n = 100; //size of the array.
   Vh[int] wh(n); //real scalar case
   Wh[int] [uh,vh](n); //real vectorial case
   Vh<complex>[int] cwh(n); //complex scalar case
   Wh<complex>[int] [cuh, cvh](n); //complex vectorial case
   [cuh[2], cvh[2]] = [x, y]; //set interpolation of index 2

   // Array of Array
   real [int][int] V(10);
   matrix[int] B(10);
   real [int, int][int] A(10);

.. tip:: Example

   In the following example, Poisson’s equation is solved for 3 different given functions :math:`f=1,\, \sin(\pi x)\cos(\pi y),\, |x-1||y-1|`, whose solutions are stored in an array of FE function.

   .. code-block:: freefem
      :linenos:

      // Mesh
      mesh Th = square(20, 20, [2*x, 2*y]);

      // Fespace
      fespace Vh(Th, P1);
      Vh u, v, f;

      // Problem
      problem Poisson (u, v)
         = int2d(Th)(
              dx(u)*dx(v)
            + dy(u)*dy(v)
         )
         + int2d(Th)(
            - f*v
         )
         + on(1, 2, 3, 4, u=0)
         ;

      Vh[int] uu(3); //an array of FE function
      // Solve problem 1
      f = 1;
      Poisson;
      uu[0] = u;
      // Solve problem 2
      f = sin(pi*x)*cos(pi*y);
      Poisson;
      uu[1] = u;
      // Solve problem 3
      f = abs(x-1)*abs(y-1);
      Poisson;
      uu[2] = u;

      // Plot
      for (int i = 0; i < 3; i++)
         plot(uu[i], wait=true);

   See :ref:`FE array example <exampleFEArray>`.

Map arrays
~~~~~~~~~~

.. code-block:: freefem
   :linenos:

   real[string] map; //a dynamic array

   map["1"] = 2.0;
   map[2] = 3.0; //2 is automatically cast to the string "2"

   cout << "map[\"1\"] = " << map["1"] << endl;
   cout << "map[2] = " << map[2] << endl;

It is just a map of the standard template library so no operations on vector are allowed, except the selection of an item.

.. _typeMatrix:

matrix
------

Defines a sparse matrix.

Matrices can be defined like vectors:

.. code-block:: freefem
   :linenos:

   matrix A = [[1, 2, 3],
               [4, 5, 6],
               [7, 8, 9]];

or using a variational form type (see :ref:`Finite Element <variationalFormSparseMatrixPDE>`):

.. code-block:: freefem
   :linenos:

   matrix Laplacian = vLaplacian(Uh, Uh);

or from block of matrices:

.. code-block:: freefem
   :linenos:

   matrix A1, ..., An;
   matrix A = [[A1, ...], ..., [..., An]];

or using sparse matrix set:

.. code-block:: freefem
   :linenos:

   A = [I, J, C];

.. note:: ``I`` and ``J`` are :freefem:`int[int]` and ``C`` is :freefem:`real[int]`. The matrix is defined as:

   .. math::
      A = \sum_k{C[k]M_{I[k], J[k]}}


   where :math:`M_{a, b} = \left(\delta_{ia}\delta_{jb}\right)_{ij}`

   `I`, `J` and `C` can be retrieved using `[I, J, C] = A` (arrays are automatically resized).

   The size of the matrix is :freefem:`n = I.max;`, :freefem:`m = J.max;`.

Matrices are designed using templates, so they can be real or complex:

.. code-block:: freefem
   :linenos:

   matrix<real> A = ...
   matrix<complex> Ai = ...

.. note:: Solver

   See :ref:`problem <typeProblem>`.

   The default solver is :ref:`GMRES <globalVariablesGMRES>`.

   .. code-block:: freefem
      :linenos:

      matrix A = vLaplacian(Uh, Uh, solver=sparsesolver);

   or

   .. code-block:: freefem
      :linenos:

      set(A , solver=sparsesolver);

.. note:: Factorize

   If :freefem:`true`, the factorization is done for :freefem:`LU`, :freefem:`Cholesky` or :freefem:`Crout`.

   .. code-block:: freefem
      :linenos:

      matrix A = vLaplacian(Uh, Uh, solver=LU, factorize=1);

   or

   .. code-block:: freefem
      :linenos:

      set(A , solver=LU, factorize=1);

.. note:: Stop test

   See :ref:`problem <typeProblem>`.

.. note:: *Très grande valeur*

   See :ref:`problem <typeProblem>`.

.. note:: Preconditioning

   See :ref:`problem <typeProblem>`.

.. note:: Pivot tolerance

   See :ref:`problem <typeProblem>`.

.. note::  :freefem:`UMFPACK`

   See :ref:`problem <typeProblem>`.

.. note:: :freefem:`dimKrylov`

   See :ref:`problem <typeProblem>`.

.. note:: datafilename

   Name of the file containing solver parameters, see :ref:`Parallel sparse solvers <parallelSparseSolvers>`

.. note:: :freefem:`lparams`

   Vector of integer parameters for the solver, see :ref:`Parallel sparse solvers <parallelSparseSolvers>`

.. note:: :freefem:`dparams`

   Vector of real parameters for the solver, see :ref:`Parallel sparse solvers <parallelSparseSolvers>`

.. note::  :freefem:`sparams`

   String parameters for the solver, see :ref:`Parallel sparse solvers <parallelSparseSolvers>`

.. tip:: To modify the :freefem:`solver`, the stop test,… after the matrix construction, use the :ref:`set keyword <functionSet>`.

Matrix size
~~~~~~~~~~~

The size of a matrix is obtain using:

.. code-block:: freefem
   :linenos:

   int NRows = A.n;
   int NColumns = A.m;

Matrix resize
~~~~~~~~~~~~~

To resize a matrix, use:

.. code-block:: freefem
   :linenos:

   A.resize(n, m);

.. warning:: When resizing, all new terms are set to zero.

Matrix diagonal
~~~~~~~~~~~~~~~

The diagonal of the matrix is obtained using:

.. code-block:: freefem
   :linenos:

   real[int] Aii = A.diag;

Matrix renumbering
~~~~~~~~~~~~~~~~~~

.. code-block:: freefem
   :linenos:

   int[int] I(15, J(15);
   matrix B = A;
   B = A(I, J);
   B = A(I^-1, J^-1);

Complex matrix
~~~~~~~~~~~~~~

Use ``.im`` and ``.re`` to get the imaginary and real part of a complex matrix, respectively:

.. code-block:: freefem
   :linenos:

   matrix<complex> C = ...
   matrix R = C.re;
   matrix I = C.im;

Dot product / Outer product
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dot product of two matrices is realized using:

.. code-block:: freefem
   :linenos:

   real d = A' * B;

The outer product of two matrices is realized using:

.. code-block:: freefem
   :linenos:

   matrix C = A * B'

See :ref:`Matrix operations example <exampleMatrixOperations>` for a complete example.

Matrix inversion
~~~~~~~~~~~~~~~~

See :ref:`Matrix inversion example <exampleMatrixInversion>`.
