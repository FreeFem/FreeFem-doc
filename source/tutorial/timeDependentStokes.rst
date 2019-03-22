.. role:: freefem(code)
  :language: freefem

Tutorial to write a transient Stokes solver in matrix form
==========================================================

Consider the following script to solve a time dependent Stokes problem in a cavity

.. code-block:: freefem
   :linenos:

   // Parameters
   real nu = 0.1;
   real T=1.;
   real dt = 0.1;

   // Mesh
   mesh Th = square(10, 10);

   // Fespace
   fespace Vh(Th, P2)
   Vh u, v;
   Vh uu, vv;
   Vh uold=0, vold=0;

   fespace Qh(Th, P1);
   Qh p;
   Qh pp;

   // Problem
   problem stokes (u, v, p, uu, vv, pp)
       = int2d(Th)(
             (u*uu+v*vv)/dt
           + nu*(dx(u)*dx(uu) + dy(u)*dy(uu) + dx(v)*dx(vv) + dy(v)*dy(vv))
           - p*pp*1.e-6
           - p*(dx(uu) + dy(vv))
           - pp*(dx(u) + dy(v))
       )
       - int2d(Th)(
             (uold*uu+vold*vv)/dt
       )
       + on(1, 2, 4, u=0, v=0)
       + on(3, u=1, v=0)
       ;

   // Time loop
   int m, M = T/dt;
   for(m = 0; m < M; m++){
       stokes;
       uold = u;
       vold = v;
   }

   // Plot
   plot(p, [u, v], value=true, wait=true, cmm="t="+m*dt);

Every iteration is in fact of the form :math:`A[u,v,p] = B[uold,vold,pold] + b` where :math:`A,B` are matrices and :math:`b` is a vector containing the boundary conditions.
:math:`A,B,b` are constructed by:

.. code-block:: freefem
   :linenos:

   fespace Xh(Th, [P2, P2, P1]);
   varf aa ([u, v, p], [uu, vv, pp])
       = int2d(Th)(
             (u*uu+v*vv)/dt
           + nu*(dx(u)*dx(uu) + dy(u)*dy(uu) + dx(v)*dx(vv) + dy(v)*dy(vv))
           - p*pp*1.e-6
           - p*(dx(uu) + dy(vv))
           - pp*(dx(u) + dy(v))
       )
       + on(1, 2, 4, u=0, v=0)
       + on(3, u=1, v=0)
       ;

   varf bb ([uold, vold, pold], [uu, vv, pp])
       = int2d(Th)(
             (uold*uu+vold*vv)/dt
       )
       //+ on(1, 2, 4, uold=0, vold=0)
       //+ on(3, uold=1, vold=0)
       ;

   varf bcl ([uold, vold, pold], [uu, vv, pp])
       = on(1, 2, 4, uold=0, vold=0)
       + on(3, uold=1, vold=0)
       ;

   matrix A = aa(Xh, Xh, solver=UMFPACK);
   matrix B = bb(Xh, Xh);
   real[int] b = bcl(0, Xh);

Note that the boundary conditions are not specified in :math:`bb`.
Removing the comment :freefem:`//` would cause the compiler to multiply the diagonal terms corresponding to a Dirichlet degree of freedom by a very large term (:freefem:`tgv`); if so :math:`b` would not be needed, on the condition that :math:`uold=1` on boundary 3 initially.
Note also that b has a tgv on the Dirichlet nodes, by construction, and so does A.

The loop will then be:

.. code-block:: freefem
   :linenos:

   real[int] sol(Xh.ndof), aux(Xh.ndof);
   for (m = 0; m < M; m++){
       aux = B*sol; aux += b;
       sol = A^-1 * aux;
   }

There is yet a difficulty with the initialization of :freefem:`sol` and with the solution from :freefem:`sol`.
For this we need a temporary vector in :math:`X_h` and here is a solution:

.. code-block:: freefem
   :linenos:

   Xh [w1, w2, wp] = [uold, vold, pp];
   sol = w1[]; //cause also the copy of w2 and wp
   for (m = 0; m < M; m++){
       aux = B*sol; aux += b;
       sol = A^-1 * aux;
   }
   w1[]=sol; u=w1; v= w2; p=wp;
   plot(p, [u, v], value=true, wait=true, cmm="t="+m*dt);

The freefem team agrees that the line :freefem:`sol=w1[];` is mysterious as it copies also w2 and wp into sol.
Structured data such as vectors of :math:`X_h` here cannot be written component by component.
Hence ``w1=u`` is not allowed.
