Acoustics
=========

**Summary :**
*Here we go to grip with ill posed problems and eigenvalue problems*

Pressure variations in air at rest are governed by the wave equation:

.. math::
   {\partial^2 u \over \partial t^2} - c^2 \Delta u =0

When the solution wave is monochromatic (and that depends on the boundary and initial conditions), :math:`u` is of the form :math:`u(x,t)=Re(v(x) e^{ik t})` where :math:`v` is a solution of Helmholtz’s equation:

.. math::
    \begin{array}{rcl}
        k^{2}v + c^{2}\Delta v &= 0 &\hbox{ in } \Omega\\
        \frac{\partial v}{\partial\boldsymbol{n}}|_\Gamma &= g &
    \end{array}

where :math:`g` is the source.

Note the “+" sign in front of the Laplace operator and that :math:`k > 0` is real.
This sign may make the problem ill posed for some values of :math:`\frac c k`, a phenomenon called "resonance".

At resonance there are non-zero solutions even when :math:`g=0`.
So the following program may or may not work:

.. code-block:: freefem
   :linenos:

   // Parameters
   real kc2 = 1.;
   func g = y*(1.-y);

   // Mesh
   border a0(t=0., 1.){x=5.; y=1.+2.*t;}
   border a1(t=0., 1.){x=5.-2.*t; y=3.;}
   border a2(t=0., 1.){x=3.-2.*t; y=3.-2.*t;}
   border a3(t=0., 1.){x=1.-t; y=1.;}
   border a4(t=0., 1.){x=0.; y=1.-t;}
   border a5(t=0., 1.){x=t; y=0.;}
   border a6(t=0., 1.){x=1.+4.*t; y=t;}

   mesh Th = buildmesh(a0(20) + a1(20) + a2(20)
       + a3(20) + a4(20) + a5(20) + a6(20));

   // Fespace
   fespace Vh(Th, P1);
   Vh u, v;

   // Solve
   solve sound(u, v)
      = int2d(Th)(
           u*v * kc2
         - dx(u)*dx(v)
         - dy(u)*dy(v)
      )
      - int1d(Th, a4)(
           g * v
      )
      ;

   // Plot
   plot(u, wait=1, ps="Sound.eps");

Results are on :numref:`figAcoustics1`. But when :math:`kc2` is an eigenvalue of the problem, then the solution is not unique:

-  if :math:`u_e \neq 0` is an eigen state, then for any given solution :math:`u+u_e` is **another** solution.

To find all the :math:`u_e` one can do the following :

.. code-block:: freefem
   :linenos:

   // Parameters
   real sigma = 20; //value of the shift

   // Problem
   // OP = A - sigma B ; // The shifted matrix
   varf op(u1, u2)
      = int2d(Th)(
           dx(u1)*dx(u2)
         + dy(u1)*dy(u2)
         - sigma* u1*u2
      )
      ;

   varf b([u1], [u2])
      = int2d(Th)(
           u1*u2
      )
      ; // No Boundary condition see note \ref{note BC EV}

   matrix OP = op(Vh, Vh, solver=Crout, factorize=1);
   matrix B = b(Vh, Vh, solver=CG, eps=1e-20);

   // Eigen values
   int nev=2; // Number of requested eigenvalues near sigma

   real[int] ev(nev);  // To store the nev eigenvalue
   Vh[int] eV(nev);    // To store the nev eigenvector

   int k=EigenValue(OP, B, sym=true, sigma=sigma, value=ev, vector=eV,
      tol=1e-10, maxit=0, ncv=0);

   cout << ev(0) << " 2 eigen values " << ev(1) << endl;
   v = eV[0];
   plot(v, wait=true, ps="eigen.eps");

.. subfigstart::

.. _figAcoustics1:

.. figure:: images/acoustics_0.png
   :alt: Amplitude of an acoustic signal coming from the left vertical wall
   :width: 90%
   :align: center

   Amplitude of an acoustic signal coming from the left vertical wall.

.. _figAcoustics2:

.. figure:: images/acoustics.png
   :alt: First eigen state
   :width: 90%
   :align: center

   First eigen state (:math:`\lambda=(k/c)^2=19.4256`) close to :math:`20` of eigenvalue problem: :math:`-\Delta \varphi = \lambda\varphi` and :math:`\frac{\partial \varphi}{\partial\boldsymbol{n}} = 0` on :math:`\Gamma`}

.. subfigend::
   :width: 0.49
   :alt: Acoustics
   :label: figAcoustics0

   Acoustics
