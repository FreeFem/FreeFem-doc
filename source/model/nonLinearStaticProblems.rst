Non-linear static problems
==========================

Here we propose to solve the following non-linear academic problem of minimization of a functional:

.. math::
    J(u) = \int_\Omega \frac{1}{2} f(|\nabla u|^2) - u*b

where :math:`u` is function of :math:`H^1_0(\Omega)` and :math:`f` defined by:

.. math::
    f(x) = a*x + x-ln(1+x),\ f'(x) = a+\frac{x}{1+x},\ f''(x) = \frac{1}{(1+x)^2}

Newton-Raphson algorithm
------------------------

Now, we solve the Euler problem :math:`\nabla J (u) = 0` with Newton-Raphson algorithm, that is:

.. math::
    u^{n+1} = u^n - ( \nabla^2 J (u^{n}))^{-1}*\nabla J(u^n)

.. code-block:: freefem
    :linenos:

    // Parameters
    real a = 0.001;
    func b = 1.;

    // Mesh
    mesh Th = square(10, 10);
    Th = adaptmesh(Th, 0.05, IsMetric=1, splitpbedge=1);
    plot(Th, wait=true);

    // Fespace
    fespace Vh(Th, P1);
    Vh u=0;
    Vh v, w;

    fespace Ph(Th, P1dc);
    Ph alpha; //to store |nabla u|^2
    Ph dalpha ; //to store 2f''(|nabla u|^2)

    // Function
    func real f (real u){
        return u*a + u - log(1.+u);
    }
    func real df (real u){
        return a +u/(1.+u);
    }
    func real ddf (real u){
        return 1. / ((1.+u)*(1.+u));
    }

    // Problem
    //the variational form of evaluate dJ = nabla J
    //dJ = f'()*(dx(u)*dx(vh) + dy(u)*dy(vh))
    varf vdJ (uh, vh)
        = int2d(Th)(
              alpha*(dx(u)*dx(vh) + dy(u)*dy(vh))
            - b*vh
        )
        + on(1, 2, 3, 4, uh=0)
        ;

    //the variational form of evaluate ddJ = nabla^2 J
    //hJ(uh,vh) = f'()*(dx(uh)*dx(vh) + dy(uh)*dy(vh))
    //  + 2*f''()(dx(u)*dx(uh) + dy(u)*dy(uh)) * (dx(u)*dx(vh) + dy(u)*dy(vh))
    varf vhJ (uh, vh)
        = int2d(Th)(
              alpha*(dx(uh)*dx(vh) + dy(uh)*dy(vh))
            + dalpha*(dx(u)*dx(vh) + dy(u)*dy(vh))*(dx(u)*dx(uh) + dy(u)*dy(uh))
        )
        + on(1, 2, 3, 4, uh=0)
        ;

    // Newton algorithm
    for (int i = 0; i < 100; i++){
        // Compute f' and f''
        alpha = df(dx(u)*dx(u) + dy(u)*dy(u));
        dalpha = 2*ddf(dx(u)*dx(u) + dy(u)*dy(u));

        // nabla J
        v[]= vdJ(0, Vh);

        // Residual
        real res = v[]'*v[];
        cout << i << " residu^2 = " << res << endl;
        if( res < 1e-12) break;

        // HJ
        matrix H = vhJ(Vh, Vh, factorize=1, solver=LU);

        // Newton
        w[] = H^-1*v[];
        u[] -= w[];
    }

    // Plot
    plot (u, wait=true, cmm="Solution with Newton-Raphson");
