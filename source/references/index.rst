.. role:: freefem(code)
  :language: freefem

Language references
===================

In essence **FreeFEM** is a compiler: its language is typed, polymorphic, with exception and reentrant.
Every variable must be declared of a certain type, in a declarative statement; each statement are separated from the next by a semicolon ``;``.

The language allows the manipulation of basic types integers (:freefem:`int`), reals (:freefem:`real`), strings (:freefem:`string`), arrays (example: :freefem:`real[int]`), bi-dimensional (2D) finite element meshes (:freefem:`mesh`), 2D finite element spaces (:freefem:`fespace`), analytical functions (:freefem:`func`), arrays of finite element functions (:freefem:`func[basic_type]`), linear and bilinear operators, sparse matrices, vectors , etc.
For example:

.. code-block:: freefem
    :linenos:

    int i, n = 20; //i, n are integer
    real[int] xx(n), yy(n); //two array of size n
    for (i = 0; i < n; i++){ //which can be used in statements such as
        xx[i] = cos(i*pi/10);
        yy[i] = sin(i*pi/10);
    }

The life of a variable is the current block ``{...}``, except the :freefem:`fespace` variable, and the variables local to a block are destroyed at the end of the block as follows.

.. tip:: Example

    .. code-block:: freefem
        :linenos:

        real r = 0.01;
        mesh Th = square(10, 10); //unit square mesh
        fespace Vh(Th, P1); //P1 Lagrange finite element space
        Vh u = x + exp(y);
        func f = z*x + r*log(y);
        plot(u, wait=true);
        { // new block
            real r = 2; //not the same r
            fespace Vh(Th, P1); //error because Vh is a global name
        }// end of block
        //here r back to 0.01

The type declarations are mandatory in **FreeFEM**; in the end this feature is an asset because it is easy to make bugs in a language with many implicit types.

The variable name is just an alphanumeric string, the underscore character ``_`` is not allowed, because it will be used as an operator in the future.

.. toctree::

   types
   global-variables
   quadrature-formulae
   operators
   loops
   IO
   functions
   external-libraries
