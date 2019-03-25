Time dependent schema optimization for heat equations
=====================================================

First, it is possible to define variational forms, and use this forms to build matrix and vector to make very fast script (4 times faster here).

For example solve the :ref:`ThermalConduction <thermalConduction>` problem, we must solve the temperature equation in :math:`\Omega` in a time interval (0,T).

.. math::
   \begin{array}{rcll}
       \partial_t u -\nabla\cdot(\kappa\nabla u) &=& 0 &\hbox{ in } \Omega\times(0,T)\\
       u(x,y,0) &=& u_0 + x u_1\\
       u &=& 30 &\hbox{ on } \Gamma_{24}\times(0,T)\\
       \kappa\frac{\partial u}{\partial n} + \alpha(u-u_e) &=& 0 &\hbox{ on } \Gamma\times(0,T)
   \end{array}

The variational formulation is in :math:`L^2(0,T;H^1(\Omega))`; we shall seek :math:`u^n` satisfying:

.. math::
   \forall w \in V_{0};\ \int_\Omega \frac{u^n-u^{n-1}}{\delta t} w + \kappa\nabla u^n\nabla w) +\int_\Gamma\alpha(u^n-u_{ue})w=0

where :math:`V_0 = \{w\in H^1(\Omega)/ w_{|\Gamma_{24}}=0\}`.

So, to code the method with the matrices :math:`A=(A_{ij})`, :math:`M=(M_{ij})`, and the vectors :math:`u^n, b^n, b',b", b_{cl}` (notation if :math:`w` is a vector then :math:`w_i` is a component of the vector).

.. math::
   u^n = A^{-1} b^n,
       \quad b' = b_0 + M u^{n-1},
       \quad b"= \frac{1}{\varepsilon} \; b_{cl} ,
       \quad b^n_i = \left\{
           \begin{array}{cl}
               b''_i & \mbox{if }\ i \in \Gamma_{24} \\
               b'_i & \mbox{else }
           \end{array}\right.
       \label{eq tgv}

Where with :math:`\frac{1}{\varepsilon} = \mathtt{tgv} = 10^{30}`:

.. math::
    \begin{array}{rcl}
        A_{ij} &=&
          \left\{\begin{array}{cl}
          \frac{1}{\varepsilon} & \mbox{if } i \in \Gamma_{24}, \mbox{and}\quad j=i\\
          \displaystyle{\int_{\Omega} w_j w_i / dt + k (\nabla w_j. \nabla w_i ) + \int_{\Gamma_{13}} \alpha w_j w_i} & \mbox{else}
          \end{array}\right.\\
        M_{ij} &=&
          \left\{\begin{array}{cl}
          \frac{1}{\varepsilon} & \mbox{if } i \in \Gamma_{24}, \mbox{and}\quad j=i \\
          \displaystyle n{\int_{\Omega} w_j w_i / dt} & \mbox{else}
          \end{array}\right. \\
        b_{0,i} &=& n{\int_{\Gamma_{13}} \alpha u_{ue} w_i } \\
        b_{cl} &=& u^{0} \quad \mbox{the initial data}
    \end{array}

The Fast version script:

.. code-block:: freefem
    :linenos:

    ...
    Vh u0=fu0, u=u0;

Create three variational formulation, and build the matrices :math:`A`,\ :math:`M`.

.. code-block:: freefem
   :linenos:

   varf vthermic (u, v)
       = int2d(Th)(
             u*v/dt
           + k*(dx(u)*dx(v) + dy(u)*dy(v))
       )
       + int1d(Th, 1, 3)(
             alpha*u*v
       )
       + on(2,4,u=1)
       ;

   varf vthermic0 (u, v)
       = int1d(Th, 1, 3)(
             alpha*ue*v
       )
       ;
   varf vMass (u,v)
       = int2d(Th)(
             u*v/dt
       )
       + on(2, 4, u=1)
       ;

   real tgv = 1e30;
   matrix A = vthermic(Vh, Vh, tgv=tgv, solver=CG);
   matrix M = vMass(Vh, Vh);

Now, to build the right hand size; we need 4 vectors.

.. code-block:: freefem
   :linenos:

   real[int] b0 = vthermic0(0,Vh); //constant part of RHS
   real[int] bcn = vthermic(0,Vh); //tgv on Dirichlet part
   real[int] bcl = tgv*u0[];   //the Dirichlet B.C. part

   // The fast loop
   for(real t = 0; t < T; t += dt){
       real[int] b = b0;   //the RHS
       b += M*u[]; //add the the time dependent part
       b = bcn ? bcl : b; //do $\forall i$: b[i] = bcn[i] ? bcl[i] : b[i];
       u[] = A^-1*b; //solve linear problem
       plot(u);
   }
