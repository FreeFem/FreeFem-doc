.. role:: freefem(code)
  :language: freefem

Heat Exchanger
==============

**Summary:**
*Here we shall learn more about geometry input and triangulation files, as well as read and write operations.*

**The problem** Let :math:`\{C_{i}\}_{1,2}`, be 2 thermal conductors within an enclosure :math:`C_0` (see :numref:`figHeatGeo`).

.. figure:: images/heat_exchangerGeo.png
  :alt: HeatExchangerGeo
  :width: 50%
  :name: figHeatGeo

  Heat exchanger geometry

The first one is held at a constant temperature :math:`{u} _{1}` the other one has a given thermal conductivity :math:`\kappa_2` 3 times larger than the one of :math:`C_0`.

We assume that the border of enclosure :math:`C_0` is held at temperature :math:`20^\circ C` and that we have waited long enough for thermal equilibrium.

In order to know :math:`{u} (x)` at any point :math:`x` of the domain :math:`\Omega`, we must solve:

.. math::
   \nabla\cdot(\kappa\nabla{u}) = 0 \hbox{ in } \Omega,
   \quad {u}_{|\Gamma} = g

where :math:`\Omega` is the interior of :math:`C_0` minus the conductor :math:`C_1` and :math:`\Gamma` is the boundary of :math:`\Omega`, that is :math:`C_0\cup C_1`.

Here :math:`g` is any function of :math:`x` equal to :math:`{u}_i` on :math:`C_i`.

The second equation is a reduced form for:

.. math::
   {u} ={u} _{i} \hbox{ on } C_{i}, \quad i=0,1.

The variational formulation for this problem is in the subspace :math:`H^1_0(\Omega) \subset H^1(\Omega)` of functions which have zero traces on :math:`\Gamma`.

.. math::
   u-g\in H^1_0(\Omega): \int_{\Omega}{\nabla u \nabla v} = 0\forall v\in H^1_0(\Omega)

Let us assume that :math:`C_0` is a circle of radius 5 centered at the origin, :math:`C_i` are rectangles, :math:`C_1` being at the constant temperature :math:`u_1=60^\circ C` (so we can only consider its boundary).

.. code-block:: freefem
   :linenos:

   // Parameters
   int C1=99;
   int C2=98; //could be anything such that !=0 and C1!=C2

   // Mesh
   border C0(t=0., 2.*pi){x=5.*cos(t); y=5.*sin(t);}

   border C11(t=0., 1.){x=1.+t; y=3.; label=C1;}
   border C12(t=0., 1.){x=2.; y=3.-6.*t; label=C1;}
   border C13(t=0., 1.){x=2.-t; y=-3.; label=C1;}
   border C14(t=0., 1.){x=1.; y=-3.+6.*t; label=C1;}

   border C21(t=0., 1.){x=-2.+t; y=3.; label=C2;}
   border C22(t=0., 1.){x=-1.; y=3.-6.*t; label=C2;}
   border C23(t=0., 1.){x=-1.-t; y=-3.; label=C2;}
   border C24(t=0., 1.){x=-2.; y=-3.+6.*t; label=C2;}

   plot(   C0(50) //to see the border of the domain
       + C11(5)+C12(20)+C13(5)+C14(20)
       + C21(-5)+C22(-20)+C23(-5)+C24(-20),
       wait=true, ps="heatexb.eps");

   mesh Th=buildmesh(C0(50)
       + C11(5)+C12(20)+C13(5)+C14(20)
       + C21(-5)+C22(-20)+C23(-5)+C24(-20));

   plot(Th,wait=1);

   // Fespace
   fespace Vh(Th, P1);
   Vh u, v;
   Vh kappa=1 + 2*(x<-1)*(x>-2)*(y<3)*(y>-3);

   // Solve
   solve a(u, v)
       = int2d(Th)(
             kappa*(
                 dx(u)*dx(v)
               + dy(u)*dy(v)
           )
       )
       +on(C0, u=20)
       +on(C1, u=60)
       ;

   // Plot
   plot(u, wait=true, value=true, fill=true, ps="HeatExchanger.eps");

Note the following:

-  ``C0`` is oriented counterclockwise by :math:`t`, while ``C1`` is oriented clockwise and ``C2`` is oriented counterclockwise.
   This is why ``C1`` is viewed as a hole by :freefem:`buildmesh`.
-  ``C1`` and ``C2`` are built by joining pieces of straight lines.
   To group them in the same logical unit to input the boundary conditions in a readable way we assigned a label on the boundaries.
   As said earlier, borders have an internal number corresponding to their order in the program (check it by adding a :freefem:`cout << C22;` above).
   This is essential to understand how a mesh can be output to a file and re-read (see below).
-  As usual the mesh density is controlled by the number of vertices assigned to each boundary.
   It is not possible to change the (uniform) distribution of vertices but a piece of boundary can always be cut in two or more parts, for instance ``C12`` could be replaced by ``C121+C122``:

.. code-block:: freefem
   :linenos:

   // border C12(t=0.,1.){x=2.; y=3.-6.*t; label=C1;}
   border C121(t=0.,0.7){x=2.; y=3.-6.*t; label=C1;}
   border C122(t=0.7,1.){x=2.; y=3.-6.*t; label=C1;}
   ...
   buildmesh(.../*+ C12(20) */ + C121(12) + C122(8) + ...);

.. subfigstart::

.. _figHeatMesh:

.. figure:: images/heat_exchangerTh.png
   :alt: HeatExchangerTh
   :width: 90%

   Heat exchanger mesh

.. _figHeatSolution:

.. figure:: images/heat_exchanger.png
   :alt: HeatExchanger
   :width: 90%

   Heat exchanger solution

.. subfigend::
   :width: 0.49
   :alt: HeatExchanger
   :label: HeatExchanger

   Heat exchanger

.. tip:: **Exercise :**

   Use the symmetry of the problem with respect to the axes.

   Triangulate only one half of the domain, and set Dirichlet conditions on the vertical axis, and Neumann conditions on the horizontal axis.

**Writing and reading triangulation files** Suppose that at the end of the previous program we added the line

.. code-block:: freefem
   :linenos:

   savemesh(Th, "condensor.msh");

and then later on we write a similar program but we wish to read the mesh from that file.
Then this is how the condenser should be computed:

.. code-block:: freefem
   :linenos:

   // Mesh
   mesh Sh = readmesh("condensor.msh");

   // Fespace
   fespace Wh(Sh, P1);
   Wh us, vs;

   // Solve
   solve b(us, vs)
       = int2d(Sh)(
             dx(us)*dx(vs)
           + dy(us)*dy(vs)
       )
       +on(1, us=0)
       +on(99, us=1)
       +on(98, us=-1)
       ;

   // Plot
   plot(us);

Note that the names of the boundaries are lost but either their internal number (in the case of ``C0``) or their label number (for ``C1`` and ``C2``) are kept.
