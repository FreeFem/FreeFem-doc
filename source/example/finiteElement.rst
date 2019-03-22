.. _exampleFinteElement:

Finite Element
==============

.. _examplePeriodic3D:

Periodic 3D
-----------

.. code-block:: freefem
   :linenos:

   load "msh3"
   load "medit"

   // Parameters
   searchMethod=1; // More safe seach algo
   real a = 1, d = 0.5, h = 0.5;
   int nnb = 7, nni = 10;
   int nz = 3;
   func zmin = 0;
   func zmax = h;

   // Mesh 2D
   border b1(t=0.5, -0.5){x=a*t; y=-a/2; label=1;}
   border b2(t=0.5, -0.5){x=a/2; y=a*t; label=2;}
   border b3(t=0.5, -0.5){x=a*t; y=a/2; label=3;}
   border b4(t=0.5, -0.5){x=-a/2; y=a*t; label=4;}
   border i1(t=0, 2.*pi){x=d/2*cos(t); y=-d/2*sin(t); label=7;}
   mesh Th = buildmesh(b1(-nnb) + b3(nnb) + b2(-nnb) + b4(nnb) + i1(nni));

   { // Cleaning the memory correctly
       int[int] old2new(0:Th.nv-1);
       fespace Vh2(Th, P1);
       Vh2 sorder = x + y;
       sort(sorder[], old2new);
       int[int] new2old = old2new^-1; // Inverse permutation
       Th = change(Th, renumv=new2old);
       sorder[] = 0:Th.nv-1;
   }
   {
       fespace Vh2(Th, P1);
       Vh2 nu;
       nu[] = 0:Th.nv-1;
       plot(nu, cmm="nu=", wait=true);
   }

   // Mesh 3D
   int[int] rup = [0, 5], rlow = [0, 6], rmid = [1, 1, 2, 2, 3, 3, 4, 4, 7, 7], rtet = [0, 41];
   mesh3 Th3 = buildlayers(Th, nz, zbound=[zmin, zmax],
       reftet=rtet, reffacemid=rmid, reffaceup=rup, reffacelow=rlow);
   for(int i = 1; i <= 6; ++i)
       cout << " int " << i << " : " << int2d(Th3,i)(1.) << " " << int2d(Th3,i)(1./area) << endl;

   plot(Th3, wait=true);
   medit("Th3", Th3);

   fespace Vh(Th3, P2, periodic=[[1, x, z], [3, x, z], [2, y, z], [4, y, z], [5, x, y], [6, x, y]]);

.. figure:: images/Periodic.jpg
   :width: 50%

   Periodic mesh

.. _exampleLagrangeMultipliers:

Lagrange multipliers
--------------------

.. code-block:: freefem
   :linenos:

   // Parameters
   func f = 1 + x - y;

   // Mesh
   mesh Th = square(10, 10);

   // Fespace
   fespace Vh(Th, P1);
   int n = Vh.ndof;
   int n1 = n+1;
   Vh uh, vh;

   // Problem
   varf va (uh, vh)
       = int2d(Th)(
             dx(uh)*dx(vh)
           + dy(uh)*dy(vh)
       )
       ;

   varf vL (uh, vh) = int2d(Th)(f*vh);
   varf vb (uh, vh) = int2d(Th)(1.*vh);

   matrix A = va(Vh, Vh);
   real[int] b = vL(0, Vh);
   real[int] B = vb(0, Vh);

   // Block matrix
   matrix AA = [ [ A, B ], [ B', 0 ] ];
   set(AA, solver=sparsesolver);

   real[int] bb(n+1), xx(n+1), b1(1), l(1);
   b1 = 0;
   // Builds the right hand side block
   bb = [b, b1];

   // Solve
   xx = AA^-1 * bb;

   // Set values
   [uh[],l] = xx;

   // Display
   cout << " l = " << l(0) << " , b(u, 1) =" << B'*uh[] << endl;

   // Plot
   plot(uh);

.. figure:: images/LagrangeMultipliers.jpg
   :width: 50%

   Result
