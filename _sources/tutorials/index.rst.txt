.. role:: freefem(code)
   :language: freefem

.. _tutorial:

Learning by Examples
====================

|
|

The **FreeFEM** language is :ref:`typed <referenceTypes>`, polymorphic and reentrant with :ref:`macro generation <typeMacroDesign>`.

Every variable must be typed and declared in a statement, that is separated from the next by a semicolon ``;``.

The **FreeFEM** language is a C++ idiom with something that is more akin to LaTeX.


For the specialist, one key guideline is that **FreeFEM** rarely generates an internal finite element array, this was adopted for speed and consequently **FreeFEM** could be hard to beat in terms of execution speed, except for the time lost in the interpretation of the language (which can be reduced by a systematic usage of :freefem:`varf` and :freefem:`matrix` instead of :freefem:`problem`).

.. rst-class:: fake-title

   The Development Cycle: Edit–Run/Visualize–Revise

Many examples and tutorials are given there after and in the :ref:`examples section <examples>`.
It is better to study them and learn by example.

If you are a beginner in the finite element method, you may also have to read a book on variational formulations.

The development cycle includes the following steps:

**Modeling:** From strong forms of PDE to weak forms, one must know the variational formulation to use **FreeFEM**; one should also have an eye on the reusability of the variational formulation so as to keep the same internal matrices; a typical example is the time dependent heat equation with an implicit time scheme: the internal matrix can be factorized only once and **FreeFEM** can be taught to do so.

**Programming:** Write the code in **FreeFEM** language using a text editor such as the one provided in your integrated environment.

**Run:** Run the code (here written in file ``mycode.edp``).
That can also be done in terminal mode by :

.. code-block:: bash
   :linenos:

   FreeFem++ mycode.edp

**Visualization:** Use the keyword :freefem:`plot` directly in ``mycode.edp`` to display functions while **FreeFEM** is running.
Use the plot-parameter :freefem:`wait=1` to stop the program at each plot.

**Debugging:** A global variable ``debug`` (for example) can help as in :freefem:`wait=true` to :freefem:`wait=false`.

.. code-block:: freefem
   :linenos:

   bool debug = true;

   border a(t=0, 2.*pi){x=cos(t); y=sin(t); label=1;};
   border b(t=0, 2.*pi){x=0.8+0.3*cos(t); y=0.3*sin(t); label=2;};

   plot(a(50) + b(-30), wait=debug); //plot the borders to see the intersection
   //so change (0.8 in 0.3 in b)
   //if debug == true, press Enter to continue

   mesh Th = buildmesh(a(50) + b(-30));
   plot(Th, wait=debug); //plot Th then press Enter

   fespace Vh(Th,P2);
   Vh f = sin(pi*x)*cos(pi*y);
   Vh g = sin(pi*x + cos(pi*y));

   plot(f, wait=debug); //plot the function f
   plot(g, wait=debug); //plot the function g

Changing debug to false will make the plots flow continuously.
Watching the flow of graphs on the screen (while drinking coffee) can then become a pleasant experience.

.. rst-class:: fake-title

   Error management

Error messages are displayed in the console window.
They are not always very explicit because of the template structure of the C++ code (we did our best!).
Nevertheless they are displayed at the right place.
For example, if you forget parenthesis as in:

.. code-block:: freefem
   :linenos:

   bool debug = true;
   mesh Th = square(10,10;
   plot(Th);

then you will get the following message from **FreeFEM**:

.. code-block:: bash
   :linenos:

   2 : mesh Th = square(10,10;
   Error line number 2, in file bb.edp, before  token ;
   parse error
   current line = 2
   syntax error
   current line = 2
   Compile error : syntax error
   line number :2, ;
   error Compile error : syntax error
   line number :2, ;
   code = 1 mpirank: 0

If you use the same symbol twice as in:

.. code-block:: freefem
   :linenos:

   real aaa = 1;
   real aaa;

then you will get the message:

.. code-block:: bash
       :linenos:

       2 : real aaa; The identifier aaa exists
             the existing type is <Pd>
             the new  type is <Pd>

If you find that the program isn’t doing what you want you may also use :freefem:`cout` to display in text format on the console window the value of variables, just as you would do in C++.

The following example works:

.. code-block:: freefem
   :linenos:

   ...
   fespace Vh(Th, P1);
   Vh u;
   cout << u;
   matrix A = a(Vh, Vh);
   cout << A;

Another trick is to *comment in and out* by using :freefem:`//` as in C++.
For example:

.. code-block:: freefem
   :linenos:

   real aaa =1;
   // real aaa;

.. toctree::

   poisson
   equationsClassification
   membrane
   heatExchanger
   acoustics
   thermalConduction
   fanBlade
   rotatingHill
   elasticity
   stokes
   navierStokesProjection
   navierStokesNewton
   aLargeFluidProblem
   complexNumbers
   optimalControl
   flowWithShocks
   heatEquationOptimization
   timeDependentStokes
   wifiPropagation
   matlabOctavePlot
