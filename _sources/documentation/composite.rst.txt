.. role:: freefem(code)
  :language: freefem

.. role:: console(code)
  :language: console

.. highlight:: freefem
  :linenothreshold: 1

.. raw:: html

  <style> .red {color:red; font-style: italic; font-weight: bold;} </style>

.. role:: red

.. _composite:

Composite finite element spaces :red:`NEW!`
===========================================

As of `FreeFEM v4.13 <https://github.com/FreeFem/FreeFem-sources/releases/tag/v4.13>`__, we introduce the notion of *composite* finite element spaces in the language in order to generalize the definition of a variational form for coupled problems. This means that coupled variables are allowed to be defined on different meshes or even mesh types. This can be useful for example for domain coupling, surface-volume coupling, FEM-BEM coupling, etc. 

.. _baseFEspaces:

Limitations of the current syntax
---------------------------------

In Freefem, you can define scalar and vector finite element spaces:

.. code-block:: freefem

  fespace Vh(Th, P1); // scalar space
  fespace Uh(Th, [P2,P2,P1], periodic=[[2,y],[4,y]]); // vector space

The current definition of FE spaces and variational problems with multiple variables or components is subject to two limitations:

- When defining periodic vector spaces, all components are considered periodic (see :freefem:`Uh` above)

- All variables/components have to be defined on the same mesh

The introduction of composite spaces aims at lifting these limitations and facilitates the definition of coupled problems, allowing to easily define and solve more general problems mixing different meshes or mesh types.

Definition of composite spaces
------------------------------

A composite FE space is defined as a cartesian product of two or more "standard" FE spaces:

.. math ::
  :label: eq_def_FEcomposite

  X_{h} = U_{h}^{1}(T_{h}^{1}, {FE}^{1})  \times U_{h}^{2} (T_{h}^{2}, {FE}^{2}) \times \cdots \times U_{h}^{n} (T_{h}^{n}, {FE}^{n}),

where :math:`T_{h}^{i}` is the mesh of the :math:`i\text{th}` FE space and :math:`{FE}^{i}` is the type of finite element:

- scalar FE element :freefem:`P1`, :freefem:`P2` :math:`\dots`

- vector FE element :freefem:`[P1,P1]`, :freefem:`RT0` :math:`\dots`

- FE element with periodic boundary :freefem:`P1,periodic=[[2,y],[4,y]]`

Each FE space :math:`U_{h}^{i}` can be defined on a different mesh :math:`T_{h}^{i}` ; meshes :math:`\left( T_{h}^{i} \right)_{i=1,n}` can even be of different types (:freefem:`mesh`, :freefem:`mesh3`, :freefem:`meshL`, :freefem:`meshS`).


First examples
~~~~~~~~~~~~~~

Defining a composite space as a product of FE spaces can be done by writing the product directly with :freefem:`*` or by using the angular bracket syntax :freefem:`<` and :freefem:`>`. For example:

- vector FE space with a periodic component in one direction: 

.. code-block:: freefem

  fespace Uh1(Th,P2);
  fespace Uh2(Th,P2,periodic=[[1,x],[3,x]]);  
  fespace Ph(Th,P1);

  // definition of the composite space Xh = Uh1 X Uh2 X Ph
  // any of the three following lines works:
  fespace Xh(<Uh1,Uh2,Ph>);
  fespace Xh(Uh1*Uh2*Ph);
  fespace Xh=Uh1*Uh2*Ph;

- composite space with the first component defined on a triangular 2D mesh :freefem:`Th` and the second component defined on a curve mesh :freefem:`ThL` (here :freefem:`ThL` is the boundary of :freefem:`Th`). This can be useful for example for volume-surface coupling or FEM-BEM coupling:

.. code-block:: freefem

  load "msh3"
  mesh Th = square(50,50); // Th is a 2d volume mesh
  meshL ThL = extract(Th); // ThL is a 1d curve mesh
  fespace Uh(Th,P1);
  fespace UhL(ThL,P1);

  // definition of the composite space Xh = Uh X UhL
  // any of the three following lines works:
  fespace Xh(<Uh,UhL>);
  fespace Xh(Uh*UhL);
  fespace Xh=Uh*UhL;


Stokes with P2-iso-P1 elements
------------------------------

Let us illustrate the use of composite spaces with the following 2D Stokes problem:

.. math::
    :label: eq_stokes_2d_BEM

    \left \{
    \begin{array}{rcl}
        -\Delta \mathbf{u} + \nabla p &=& \mathbf{f} \quad \text{in} \quad \Omega, \\
        \nabla \cdot \mathbf{u} &=& 0 \quad \text{in} \quad \Omega, \\
        \mathbf{u} &=& g \quad \text{on} \quad \Gamma, 
    \end{array}
    \right .

where :math:`\mathbf{u}=(u_1,u_2)` is the fluid velocity and :math:`p` the pressure.

In order to define the variational form, we multiply the first equation (resp. the second equation) of :eq:`eq_stokes_2d_BEM` by a test function :math:`\mathbf{v}` (resp. :math:`q`) and integrate on :math:`\Omega`:

.. math::
  
  \forall (\mathbf{v}, \: q ), \quad \left \{
  \begin{aligned}
  \int_{\Omega} \nabla \mathbf{u} . \nabla \mathbf{v} \: d\boldsymbol{x} + \int_{\Omega} \nabla p \cdot \mathbf{v} d\boldsymbol{x} &= \int_{\Omega} \mathbf{f} \cdot \mathbf{v} \: d\boldsymbol{x}, \\
  - \int_{\Omega} div(\mathbf{u}) q d\boldsymbol{x} &= 0.
  \end{aligned}
  \right . 

Using Green's formula in the second integral, we obtain:

.. math::
  :label: eq_varf_stokes_2d_BEM_dTh

  \forall (\mathbf{v}, \: q ), \left \{
  \begin{aligned}
  \int_{\Omega} \nabla \mathbf{u} . \nabla \mathbf{v} \:  d\boldsymbol{x} - \int_{\Omega} p div(\mathbf{v})  d\boldsymbol{x}  &= \int_{\Omega} \mathbf{f} \cdot \mathbf{v} \: d\boldsymbol{x}, \\
  - \int_{\Omega} div(\mathbf{u}) q  d\boldsymbol{x} \color{blue}{-\int_{\Omega} \epsilon p q d\boldsymbol{x}} &= 0.
  \end{aligned}
  \right . 
  
The stabilization term in blue is added to the variational form to fix the constant for the pressure (we take for example :math:`\epsilon=10^{-10}`).

We choose the **P2-iso-P1** finite element for the velocity and pressure. This element requires two different meshes: a coarser mesh for pressure and a finer mesh for velocity obtained by splitting each triangle of the pressure mesh into four triangles.

using *solve* or *problem*
~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to solve this problem, we can make use of the new composite spaces as in the script below:

.. code-block:: freefem

  int nn = 30; // number of edges in each direction
  mesh ThP = square(nn,nn,[2*pi*x,2*pi*y],flags=3); // Pressure mesh
  mesh ThU = trunc(ThP,1,split=2);  // Velocity mesh

  fespace Uh(ThU,[P1,P1]); // Velocity space
  fespace Ph(ThP,P1);      // Pressure space

  macro grad(u) [dx(u),dy(u)] //
  macro Grad(u1,u2) [grad(u1), grad(u2)] //
  macro div(u1,u2) (dx(u1)+dy(u2)) //

  // definition of the boundary condition
  func g1 = sin(x)*cos(y);
  func g2 = -cos(x)*sin(y);

  // definition of the right-hand side
  func f1 = 0;
  func f2 = -4*cos(x)*sin(y);

  Uh [u1,u2],[v1,v2];
  Ph p,q;

  solve Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>) = int2d(ThU)((Grad(u1,u2):Grad(v1,v2)))
  + int2d(ThU)(-div(u1,u2)*q -div(v1,v2)*p)
  + int2d(ThP)(-1e-10*p*q)
  - int2d(ThU)([f1,f2]'*[v1,v2])
  + on(1,2,3,4, u1=g1, u2=g2);

  plot([u1,u2], cmm="u");
  plot(p, cmm="p");

You can also find this example in the **FreeFEM** distribution `here <https://github.com/FreeFem/FreeFem-sources/blob/master/examples/examples/stokes_composite.edp>`__.  

Note that with the :freefem:`problem` or :freefem:`solve` syntax, the composite nature of the FE space has to be indicated directly in the :freefem:`problem/solve` instruction using angular brackets :freefem:`<` and :freefem:`>`: 

.. code-block:: freefem

  solve Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>) = ...

The explicit definition of the composite :freefem:`fespace Xh=Uh*Ph` is optional.  

Remark that if you omit to indicate the composite nature of the problem with :freefem:`<`, :freefem:`>` and write

.. code-block:: freefem

  solve Stokes ([u1,u2,p], [v1,v2,q]) = ...

**FreeFEM** falls back to the "standard" evaluation of the variational form and outputs the following error message: :console:`Exec error : all the finite element spaces must be defined on the same mesh in solve`.

using *varf* and *matrix*
~~~~~~~~~~~~~~~~~~~~~~~~~

Composite FE spaces can also be used with the :freefem:`varf`/:freefem:`matrix` syntax. We can replace the :freefem:`solve` instruction

.. code-block:: freefem

  solve Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>) = int2d(ThU)((Grad(u1,u2):Grad(v1,v2)))
  + int2d(ThU)(-div(u1,u2)*q -div(v1,v2)*p)
  + int2d(ThP)(-1e-10*p*q)
  - int2d(ThU)([f1,f2]'*[v1,v2])
  + on(1,2,3,4, u1=g1, u2=g2);

by

.. code-block:: freefem

  fespace Xh=Uh*Ph;

  varf Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>) = int2d(ThU)((Grad(u1,u2):Grad(v1,v2)))
  + int2d(ThU)(-div(u1,u2)*q -div(v1,v2)*p)
  + int2d(ThP)(-1e-10*p*q)
  + int2d(ThU)([f1,f2]'*[v1,v2])
  + on(1,2,3,4, u1=g1, u2=g2);

  matrix M = Stokes(Xh,Xh);
  real[int] b = Stokes(0,Xh);
  real[int] sol = M^-1*b;

  [u1[],p[]] = sol; // dispatch the solution

.. note::
  The sign of the linear form is flipped because with the :freefem:`varf`/:freefem:`matrix` syntax we are solving :math:`M x = b` whereas with :freefem:`problem` or :freefem:`solve` we are solving :math:`M x - b = 0`.

Note that in this case the explicit definition of the composite :freefem:`fespace Xh=Uh*Ph` is mandatory, as it is used in the instantiation of the variational form in :freefem:`Stokes(Xh,Xh)`. However, this time the :freefem:`<`, :freefem:`>` notation in the :freefem:`varf` is optional: we can also write

.. code-block:: freefem

  varf Stokes ([u1,u2,p], [v1,v2,q]) = ...

In the last instruction, the solution vector obtained when solving the linear system is dispatched to the velocity and pressure FE functions:

.. code-block:: freefem

  [u1[],p[]] = sol; // dispatch the solution

Under the hood, the matrix and right-hand side of a composite problem are assembled in a contiguous way with respect to the different components. This means that the linear system we are solving has the expected structure

.. math::
    \left(
    \begin{array}{cc}
    \mathbf{A}&-B^T\\
    -B&-\epsilon I
    \end{array}
    \right)
    \left(
    \begin{array}{cc}
    \mathbf{u}\\
    p
    \end{array}
    \right)
    =
    \left(
    \begin{array}{cc}
    \mathbf{f}\\
    0
    \end{array}
    \right)

Parallel assembly and solution
------------------------------

With the definition of composite problems, we also introduce a way to parallelize the assembly and solution of the linear system in a transparent way. There are two ways to make use of this functionality: with the parallel direct solver *MUMPS* and with *PETSc*.

with *MUMPS*
~~~~~~~~~~~~

With minimal changes in the script, we can parallelize the assembly of the linear system on multiple MPI processes and solve it in parallel using the *MUMPS* solver.  

First, we need to load the :freefem:`MUMPS` and :freefem:`bem` plugins (this functionality is implemented in the :freefem:`bem` plugin for now but this will be changed in future releases):

.. code-block:: freefem

  load "MUMPS"
  load "bem"

Now we just need to specify that we want to solve the problem with a sparse direct solver (*MUMPS*), with the linear system being distributed over all MPI processes (with :freefem:`master=-1`):

.. code-block:: freefem

  solve Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>, solver=sparsesolver, master=-1) = ...

It works the same way when using the :freefem:`varf`/:freefem:`matrix` syntax:

.. code-block:: freefem
  :emphasize-lines: 5

  fespace Xh=Uh*Ph;

  varf Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>) = ...

  matrix M = Stokes(Xh,Xh, solver=sparsesolver, master=-1);
  real[int] b = Stokes(0,Xh);
  real[int] sol = M^-1*b;

  [u1[],p[]] = sol; // dispatch the solution

At the end of the solution phase, all MPI processes hold the global solution of the problem in (:freefem:`u1,u2,p`).

In order to use multiple MPI processes, we run the script in parallel on e.g. 4 computing cores with

.. code-block:: console

  ff-mpirun -np 4 script.edp -wg

with *PETSc*
~~~~~~~~~~~~

Alternatively, we can use the *PETSc* library (see :ref:`PETSc and SLEPc <petscslepc>`) to solve the linear system. In contrast to the usual fully distributed *PETSc* framework (see :ref:`PETSc examples <petscExamples>`), with the composite syntax the parallel data distribution is hidden to the user and the solution output is given in the global space. This is part of the continuous efforts to hide the difficulties associated with the parallelization of a user script, making it more transparent to the user.  

In order to use *PETSc* as a solver, we simply load the :freefem:`PETSc` plugin and define a *PETSc* :freefem:`Mat` instead of a :freefem:`matrix` for the linear system matrix:

.. code-block:: freefem
  :emphasize-lines: 1, 12

  load "PETSc"
  load "bem"

  fespace Xh=Uh*Ph;

  varf Stokes (<[u1,u2],[p]>, <[v1,v2],[q]>) = int2d(ThU)((Grad(u1,u2):Grad(v1,v2)))
  + int2d(ThU)(-div(u1,u2)*q -div(v1,v2)*p)
  + int2d(ThP)(-1e-10*p*q)
  + int2d(ThU)([f1,f2]'*[v1,v2])
  + on(1,2,3,4, u1=g1, u2=g2);

  Mat M = Stokes(Xh,Xh);
  real[int] b = Stokes(0,Xh);
  real[int] sol = M^-1*b;

  [u1[],p[]] = sol; // dispatch the solution

We can then indicate what type of *PETSc* solver we want to use with *PETSc*-specific command-line arguments. We can start by using a simple distributed LU solver:

.. code-block:: console

  ff-mpirun -np 4 script.edp -wg -pc_type lu

We can also specify the solver by putting the corresponding *PETSc* flags directly in the script, with the :freefem:`set()` instruction:

.. code-block:: freefem
  :emphasize-lines: 2

  Mat M = Stokes(Xh,Xh);
  set(M, sparams = "-pc_type lu");

An interesting feature of using *PETSc* to solve a composite problem is that we can easily define a *fieldsplit* preconditioner combining separate preconditioners defined for each variable of the composite space, as the underlying composite *PETSc* operator already carries the block structure corresponding to the different variables.  

This means that we can reference each variable -- or *split* -- right out of the box in order to define our fieldsplit preconditioner. They come in the same order as in the definition of the composite space ; in our Stokes example, the velocity space is split 0 and the pressure space is split 1. As an illustration, we can use the Schur complement method to solve our Stokes problem, using a multigrid *HYPRE* preconditioner for the velocity block: 

.. code-block:: freefem
  :emphasize-lines: 2

  Mat M = Stokes(Xh,Xh);
  set(M, sparams = "-pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_pc_type hypre");
