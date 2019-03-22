.. role:: freefem(code)
  :language: freefem

Evolution problems
==================

**FreeFEM** also solves evolution problems such as the heat equation:

.. math::
    \begin{array}{rcll}
        \frac{\p u}{\p t}-\mu\Delta u &=& f & \textrm{ in }\Omega\times ]0,T[\\
        u(\mathbf{x},0) &=& u_0(\mathbf{x}) & \textrm{ in }\Omega\\
        \left(\p u/\p n\right)(\mathbf{x},t) &=& 0 & \textrm{ on }\p\Omega\times ]0,T[
    \end{array}
    :label: eqn::heatequation

with a positive viscosity coefficient :math:`\mu` and homogeneous Neumann boundary conditions.

We solve :eq:`eqn::heatequation` by FEM in space and finite differences in time.

We use the definition of the partial derivative of the solution in the time derivative:

.. math::
    \frac{\p u}{\p t}(x,y,t) = \lim_{\tau \to 0}\frac{u(x,y,t)-u(x,y,t-\tau )}{\tau }

which indicates that :math:`u^m(x,y)=u(x,y,m\tau )` will satisfy approximatively:

.. math::
    \frac{\p u}{\p t}(x,y,m\tau )\simeq \frac{u^m(x,y)-u^{m-1}(x,y)}{\tau }

The time discretization of heat equation :eq:`eqn::heatequation` is as follows, :math:`\forall m=0,\cdots,[T/\tau ]`:

.. math::
    \begin{array}{rcll}
        \frac{u^{m+1}-u^{m}}{\tau }-\mu\Delta u^{m+1} &=& f^{m+1} & \textrm{ in }\Omega\\
        u^0(\mathbf{x}) &=& u_0(\mathbf{x}) & \textrm{ in }\Omega\\
        \p u^{m+1}/\p n(\mathbf{x}) &=& 0 & \textrm{ on }\p\Omega
    \end{array}

which is so-called *backward Euler method* for :eq:`eqn::heatequation`.

To obtain the variational formulation, multiply with the test function :math:`v` both sides of the equation:

.. math::
    \int_{\Omega}\{u^{m+1}v-\tau \Delta u^{m+1}v\}=\int_{\Omega}\{u^m+\tau f^{m+1}\}v

By the divergence theorem, we have:

.. math::
    \int_{\Omega}\{u^{m+1}v+\tau\nabla u^{m+1}\cdot \nabla v\}
    -\int_{\p\Omega} \tau \left( \p u^{m+1}/\p n\right) v
    =\int_{\Omega }\{u^mv+\tau f^{m+1}v\}

By the boundary condition :math:`\p u^{m+1}/\p n=0`, it follows that:

.. math::
    \int_{\Omega} \{u^{m+1}v+\tau \nabla u^{m+1}\cdot \nabla v\}
    -\int_{\Omega }\{u^mv+\tau f^{m+1}v\}
    =0
    :label: eqn::heatequationBWE

Using the identity just above, we can calculate the finite element approximation :math:`u_h^m` of :math:`u^m` in a step-by-step manner with respect to :math:`t`.

.. tip:: Example

    We now solve the following example with the exact solution :math:`u(x,y,t)=tx^4`, :math:`\Omega = ]0,1[^2`.

    .. math::
        \begin{array}{rcll}
            \frac{{\p u}}{{\p t}} - \mu \Delta u &=& x^4 - \mu 12tx^2 & \textrm{ in }\Omega\times ]0,3[\\
            u(x,y,0) &=& 0 & \textrm{ on }\Omega\\
            \left. u \right|_{\p\Omega} &=& t*x^4
        \end{array}

    .. code-block:: freefem
        :linenos:

        // Parameters
        real dt = 0.1;
        real mu = 0.01;

        // Mesh
        mesh Th = square(16, 16);

        // Fespace
        fespace Vh(Th, P1);
        Vh u, v, uu, f, g;

        // Problem
        problem dHeat (u, v)
            = int2d(Th)(
                  u*v
                + dt*mu*(dx(u)*dx(v) + dy(u)*dy(v))
            )
            + int2d(Th)(
                - uu*v
                - dt*f*v
            )
            + on(1, 2, 3, 4, u=g)
            ;

        // Time loop
        real t = 0;
        uu = 0;
        for (int m = 0; m <= 3/dt; m++){
            // Update
            t = t+dt;
            f = x^4 - mu*t*12*x^2;
            g = t*x^4;
            uu = u;

            // Solve
            dHeat;

            // Plot
            plot(u, wait=true);
            cout << "t=" << t << " - L^2-Error=" << sqrt(int2d(Th)((u-t*x^4)^2)) << endl;
        }

    In the last statement, the :math:`L^2`-error :math:`\left(\int_{\Omega}\left| u-tx^4\right|^2\right)^{1/2}` is calculated at :math:`t=m\tau, \tau =0.1`. At :math:`t=0.1`, the error is 0.000213269. The errors increase with :math:`m` and 0.00628589 at :math:`t=3`.

    The iteration of the backward Euler :eq:`eqn::heatequationBWE` is made by :ref:`for loop<loopFor>`.

    .. note:: The stiffness matrix in the loop is used over and over again.
        **FreeFEM** support reuses of stiffness matrix.

Mathematical Theory on Time Difference Approximations.
------------------------------------------------------

In this section, we show the advantage of implicit schemes.
Let :math:`V, H` be separable Hilbert space and :math:`V` is dense in :math:`H`.
Let :math:`a` be a continuous bilinear form over :math:`V \times V` with coercivity and symmetry.

Then :math:`\sqrt{a(v,v)}` become equivalent to the norm :math:`\| v\|` of :math:`V`.

**Problem Ev(f,\Omega)**: For a given :math:`f\in L^2(0,T;V'),\, u^0\in H`

.. math::
    \begin{array}{rcl}
        \frac{d}{dt}(u(t),v)+a(u(t),v)&=&( f(t),v)\qquad \forall v\in V,\quad a.e. \, t\in [0,T]\\
        u(0)&=&u^0\nonumber
    \end{array}

where :math:`V'` is the dual space of :math:`V`.

Then, there is an unique solution :math:`u\in L^{\infty}(0,T;H)\cap L^2(0,T;V)`.

Let us denote the time step by :math:`\tau>0`, :math:`N_T=[T/\tau]`.
For the discretization, we put :math:`u^n = u(n\tau)` and consider the time difference for each :math:`\theta\in [0,1]`

.. math::
    \begin{array}{rcl}
        \frac{1}{\tau}\left( u_h^{n+1}-u_h^n,\phi_i\right) +a\left( u_h^{n+\theta},\phi_i\right)&=&\langle f^{n+\theta},\phi_i\rangle\\
        i=1,&\cdots&, m,\quad n=0,\cdots, N_T\nonumber\\
        u_h^{n+\theta}&=&\theta u_h^{n+1}+(1-\theta)u_h^n,\\
        f^{n+\theta}&=&\theta f^{n+1}+(1-\theta)f^n\nonumber
    \end{array}
    :label: eqn::t-method

Formula :eq:`eqn::t-method` is the *forward Euler scheme* if :math:`\theta=0`, *Crank-Nicolson scheme* if :math:`\theta=1/2`, the *backward Euler scheme* if :math:`\theta=1`.

Unknown vectors :math:`u^n=(u_h^1,\cdots,u_h^M)^T` in

.. math::
    u_h^n(x)=u^n_1\phi_1(x)+\cdots+u^n_m\phi_m(x),\quad u^n_1,\cdots,u^n_m\in \R

are obtained from solving the matrix

.. math::
    (M+\theta\tau A)u^{n+1}=\{M-(1-\theta)\tau A\}u^n
    +\tau\left\{\theta f^{n+1}+(1-\theta)f^n\right\}\\
    M=(m_{ij}),\quad m_{ij}=(\phi_j,\phi_i),\qquad
    A=(a_{ij}),\quad a_{ij}=a(\phi_j,\phi_i)\nonumber
    :label: eqn::Evolution-1

Refer [TABATA1994]_, pp.70–75 for solvability of :eq:`eqn::Evolution-1`. The stability of :eq:`eqn::Evolution-1` is in [TABATA1994]_, Theorem 2.13:

Let :math:`\{\mathcal{T}_h\}_{h\downarrow 0}` be regular triangulations (see :ref:`Regular Triangulation <meshRegularTriangulation>`).
Then there is a number :math:`c_0>0` independent of :math:`h` such that,

.. math::
    |u_h^n|^2\le
    \left\{
    \begin{array}{lr}
    \frac{1}{\delta}\left\{
    |u^0_h|^2+\tau \sum_{k=0}^{n-1}\|f^{k+\theta}\|^2_{V_h'}
    \right\}&\theta\in [0,1/2)\\
    |u^0_h|^2+\tau \sum_{k=0}^{n-1}\|f^{k+\theta}\|^2_{V_h'}&\theta\in [1/2,1]
    \end{array}
    \right.

if the following are satisfied:

1. When :math:`\theta\in [0,1/2)`, then we can take a time step :math:`\tau` in such a way that

    .. math::
        \tau <\frac{2(1-\delta)}{(1-2\theta)c_0^2}h^2

    for arbitrary :math:`\delta\in (0,1)`.

2. When :math:`1/2\leq \theta\leq 1`, we can take :math:`\tau` arbitrary.

.. tip:: Example

    .. code-block:: freefem
        :linenos:

        // Parameters
        real tau = 0.1; real
        theta = 0.;

        // Mesh
        mesh Th = square(12, 12);

        // Fespace
        fespace Vh(Th, P1);
        Vh u, v, oldU;
        Vh f1, f0;

        fespace Ph(Th, P0);
        Ph h = hTriangle; // mesh sizes for each triangle

        // Function
        func real f (real t){
            return x^2*(x-1)^2 + t*(-2 + 12*x - 11*x^2 - 2*x^3 + x^4);
        }

        // File
        ofstream out("err02.csv"); //file to store calculations
        out << "mesh size = " << h[].max << ", time step = " << tau << endl;
        for (int n = 0; n < 5/tau; n++)
            out << n*tau << ",";
        out << endl;

        // Problem
        problem aTau (u, v)
            = int2d(Th)(
                  u*v
                + theta*tau*(dx(u)*dx(v) + dy(u)*dy(v) + u*v)
            )
            - int2d(Th)(
                  oldU*v
                - (1-theta)*tau*(dx(oldU)*dx(v) + dy(oldU)*dy(v) + oldU*v)
            )
            - int2d(Th)(
                  tau*(theta*f1 + (1-theta)*f0)*v
            )
            ;

        // Theta loop
        while (theta <= 1.0){
            real t = 0;
            real T = 3;
            oldU = 0;
            out << theta << ",";
            for (int n = 0; n < T/tau; n++){
                // Update
                t = t + tau;
                f0 = f(n*tau);
                f1 = f((n+1)*tau);

                // Solve
                aTau;
                oldU = u;

                // Plot
                plot(u);

                // Error
                Vh uex = t*x^2*(1-x)^2; //exact solution = tx^2(1-x)^2
                Vh err = u - uex; // err = FE-sol - exact
                out << abs(err[].max)/abs(uex[].max) << ",";
            }
            out << endl;
            theta = theta + 0.1;
        }

    .. figure:: images/EvolutionProblems_TimeDifference.png
        :name: figEvolutionTimeDiff

        :math:`\max_{x\in\Omega}\vert u_h^n(\theta)-u_{ex}(n\tau)\vert\max_{x\in\Omega}\vert u_{ex}(n\tau)\vert at n=0,1,\cdots,29`

    We can see in :numref:`figEvolutionTimeDiff` that :math:`u_h^n(\theta)` become unstable at :math:`\theta=0.4`, and figures are omitted in the case :math:`\theta<0.4`.

Convection
----------

The hyperbolic equation

.. math::
    \p_t u +\mathbf{\alpha} \cdot \nabla u=f;\ \textrm{ for a vector-valued function }\mathbf{\alpha}
    :label: eqn::conv

appears frequently in scientific problems, for example in the Navier-Stokes equations, in the Convection-Diffusion equation, etc.

In the case of 1-dimensional space, we can easily find the general solution :math:`(x,t)\mapsto u(x,t)=u^0(x-\alpha t)` of the following equation, if :math:`\alpha` is constant,

.. math::
    \begin{array}{rcl}
        \p_t u +\alpha\p_x u &=& 0\\
        u(x,0) &=& u^0(x),
    \end{array}
    :label: eqn::conv0

because :math:`\p_t u +\alpha\p_x u=-\alpha\dot{u}^0+a\dot{u}^0=0`, where :math:`\dot{u}^0=du^0(x)/dx`.

Even if :math:`\alpha` is not constant, the construction works on similar principles.
One begins with the ordinary differential equation (with the convention that :math:`\alpha` is prolonged by zero apart from :math:`(0,L)\times (0,T)`):

.. math::
    \dot{X}(\tau )=+\alpha(X(\tau ),\tau ),\ \tau \in (0,t)\quad X(t)=x

In this equation :math:`\tau` is the variable and :math:`x,t` are parameters, and we denote the solution by :math:`X_{x,t}(\tau )`.
Then it is noticed that :math:`(x,t)\rightarrow v(X(\tau),\tau)` in :math:`\tau=t` satisfies the equation

.. math::
    \p _{t}v+\alpha\p _{x}v=\p _{t}X\dot{v}+a\p _{x}X\dot{v}%
    =0

and by the definition :math:`\p _{t}X=\dot{X}=+\alpha` and :math:`\p_{x}X=\p _{x}x` in :math:`\tau=t`, because if :math:`\tau =t` we have :math:`X(\tau )=x`.

The general solution of :eq:`eqn::conv0` is thus the value of the boundary condition in :math:`X_{x, t}(0)`, that is to say :math:`u(x,t)=u^{0}(X_{x,t}(0))` where :math:`X_{x,t}(0)` is on the :math:`x` axis, :math:`u(x,t)=u^{0}(X_{x,t}(0))` if :math:`X_{x,t}(0)` is on the axis of :math:`t`.

In higher dimension :math:`\Omega \subset R^{d},~d=2,3`, the equation for the convection is written

.. math::
    \p _{t}u+\mathbf{\alpha}\cdot \nabla u=0\hbox{ in }\Omega \times (0,T)

where :math:`\mathbf{a}(x,t)\in \R^{d}`.

**FreeFEM** implements the Characteristic-Galerkin method for convection operators.
Recall that the equation :eq:`eqn::conv` can be discretized as

.. math::
    \frac{Du}{Dt} = f\;\;\textrm{i.e. }\frac{du}{dt}\left( {X(t),t} \right) = f\left(X( t ),t \right)\textrm{ where }\frac{dX}{dt}( t ) = \mathbf{\alpha}( {X(t),t})

where :math:`D` is the total derivative operator.
So a good scheme is one step of backward convection by the method of Characteristics-Galerkin

.. math::
    \frac{1}{{\tau }}\left(u^{m + 1}(x) - u^m(X^m(x))\right) = f^m (x)
    :label: eqn::Charac

where :math:`X^m (x)` is an approximation of the solution at :math:`t = m\tau` of the ordinary differential equation

.. math::
    \frac{d\mathbf{X}}{dt}(t) = \mathbf{\alpha}^m(\mathbf{X}(t)), \mathbf{X}((m + 1)\tau) = x.

where :math:`\mathbf{\alpha}^m(x)=(\alpha_1(x,m\tau ),\alpha_2(x,m\tau))`.
Because, by Taylor’s expansion, we have

.. math::
    \begin{array}{rcl}
        u^m(\mathbf{X}(m\tau ))&=&
        u^m(\mathbf{X}((m+1)\tau )) -
        \tau \sum_{i=1}^d \frac{\p u^m}{\p x_i}(\mathbf{X}((m+1)\tau ))
        \frac{\p X_i}{\p t}((m+1)\tau )
        +o(\tau )\nonumber\\
        &=&u^m(x)-\tau \mathbf{\alpha}^m(x)\cdot \nabla u^m(x)+o(\tau )
    \end{array}
    :label: eqn::conv1

where :math:`X_i(t)` are the i-th component of :math:`\mathbf{X}(t)`, :math:`u^m(x)=u(x,m\tau )` and we used the chain rule and :math:`x=\mathbf{X}((m+1)\tau )`.
From :eq:`eqn::conv1`, it follows that

.. math::
    u^m(X^m(x))=u^m(x)-\tau \mathbf{\alpha}^m(x)\cdot \nabla u^m(x)+o(\tau )

Also we apply Taylor’s expansion for :math:`t \rightarrow u^m(x-\mathbf{\alpha}^m(x)t),0\le t\le \tau`, then

.. math::
    u^m(x-\mathbf{\alpha}\tau )=u^m(x)-\tau \mathbf{\alpha}^m(x)\cdot \nabla u^m(x)+o(\tau ).

Putting

:freefem:`convect`:math:`\left( {\mathbf{\alpha},-\tau ,u^m } \right)\approx u^m \left(x - \mathbf{\alpha}^m\tau \right)`

we can get the approximation

:math:`u^m \left( {X^m( x )} \right) \approx` :freefem:`convect` :math:`\left( {[a_1^m ,a_2^m],-\tau ,u^m } \right)` by :math:`X^m \approx x \mapsto x- \tau [a_1^m(x) ,a_2^m(x)]`

A classical convection problem is that of the “rotating bell" (quoted from [LUCQUIN1998]_, p.16).

Let :math:`\Omega` be the unit disk centered at 0, with its center rotating with speed :math:`\alpha_1 = y,\, \alpha_2 = -x`.
We consider the problem :eq:`eqn::conv` with :math:`f=0` and the initial condition :math:`u(x,0)=u^0(x)`, that is, from :eq:`eqn::Charac`

:math:`u^{m + 1}(x) = u^m(X^m(x))\approx` :freefem:`convect`\ :math:`(\mathbf{\alpha},-\tau ,u^m)`

The exact solution is :math:`u(x, t) = u(\mathbf{X}(t))` where :math:`\mathbf{X}` equals :math:`x` rotated around the origin by an angle :math:`\theta = -t` (rotate in clockwise).
So, if :math:`u^0` in a 3D perspective looks like a bell, then :math:`u` will have exactly the same shape, but rotated by the same amount.
The program consists in solving the equation until :math:`T = 2\pi`, that is for a full revolution and to compare the final solution with the initial one; they should be equal.

.. tip:: Convect

    .. code-block:: freefem
        :linenos:

        // Parameters
        real dt = 0.17;

        // Mesh
        border C(t=0, 2*pi){x=cos(t); y=sin(t);}
        mesh Th = buildmesh(C(70));

        // Fespace
        fespace Vh(Th, P1);
        Vh u0;
        Vh a1 = -y, a2 = x; //rotation velocity
        Vh u;

        // Initialization
        u = exp(-10*((x-0.3)^2 +(y-0.3)^2));

        // Time loop
        real t = 0.;
        for (int m = 0; m < 2*pi/dt; m++){
            // Update
            t += dt;
            u0 = u;

            // Convect
            u = convect([a1, a2], -dt, u0); //u^{m+1}=u^m(X^m(x))

            // Plot
            plot(u, cmm=" t="+t+", min="+u[].min+", max="+u[].max);
        }

    .. note:: The scheme :freefem:`convect` is unconditionally stable, then the bell become lower and lower (the maximum of :math:`u^{37}` is :math:`0.406` as shown in :numref:`figEvolutionConvect`.

    .. subfigstart::

    .. _figEvolutionConvect:

    .. figure:: images/EvolutionProblem_Convect.png
        :width: 90%
        :alt: EvolutionProblem_Convect

        :math:`u^0=e^{-10((x-0.3)^2 +(y-0.3)^2)}`

    .. _figEvolutionConvect2:

    .. figure:: images/EvolutionProblem_Convect2.png
        :width: 90%
        :alt: EvolutionProblem_Convect2

        The bell at :math:`t=6.29`

    .. subfigend::
       :width: 0.49
       :alt: EvolutionProblem_Convect
       :label: EvolutionProblem_Convect

2D Black-Scholes equation for an European Put option
----------------------------------------------------

In mathematical finance, an option on two assets is modeled by a Black-Scholes equations in two space variables, (see for example [WILMOTT1995]_ or [ACHDOU2005]_).

.. math::
    \begin{array}{rcl}
        \p _t u &+& \frac{{\left( {\sigma _1 x } \right)^2 }}{2}\frac{{\p ^2 u}}{{\p x^2 }} + \frac{{\left( {\sigma _2 y } \right)^2 }}{2}\frac{{\p ^2 u}}{{\p y^2 }} \\
        &&{\rm{ }} + \rho x y \frac{{\p ^2 u}}{{\p x \p y }} + rS_1 \frac{{\p u}}{{\p x }} + rS_2 \frac{{\p u}}{{\p y }} - rP = 0 \nonumber
    \end{array}

which is to be integrated in :math:`\left( {0,T} \right) \times \R^ + \times \R^ +` subject to, in the case of a put

.. math::
    u\left( {x , y ,T} \right) = \left( {K - \max \left( {x ,y } \right)} \right)^+

Boundary conditions for this problem may not be so easy to device.
As in the one dimensional case the PDE contains boundary conditions on the axis :math:`x_1 = 0` and on the axis :math:`x_2 = 0`, namely two one dimensional Black-Scholes equations driven respectively by the data :math:`u\left( {0, + \infty ,T} \right)` and :math:`u\left( { + \infty ,0,T} \right)`.
These will be automatically accounted for because they are embedded in the PDE.
So if we do nothing in the variational form (i.e. if we take a Neumann boundary condition at these two axis in the strong form) there will be no disturbance to these.
At infinity in one of the variable, as in 1D, it makes sense to impose :math:`u=0`.
We take

.. math::
    \sigma _1  = 0.3,\;\;\sigma _2  = 0.3,\;\;\rho  = 0.3,\;\;r = 0.05,\;\;K = 40,\;\;T = 0.5

An implicit Euler scheme is used and a mesh adaptation is done every 10 time steps.
To have an unconditionally stable scheme, the first order terms are treated by the Characteristic Galerkin method, which, roughly, approximates

.. math::
    \frac{{\p u}}{{\p t}} + a_1 \frac{{\p u}}{{\p x}} + a_2 \frac{{\p u}}{{\p y}} \approx \frac{1}{{\tau }}\left( {u^{n + 1} \left( x \right) - u^n \left( {x - \mathbf{\alpha}\tau } \right)} \right)

.. tip:: Black-Scholes

    .. code-block:: freefem
        :linenos:

        // Parameters
        int m = 30; int L = 80; int LL = 80; int j = 100; real sigx = 0.3; real sigy = 0.3; real rho = 0.3; real r = 0.05; real K = 40; real dt = 0.01;

        // Mesh
        mesh th = square(m, m, [L*x, LL*y]);

        // Fespace
        fespace Vh(th, P1);
        Vh u = max(K-max(x,y),0.);
        Vh xveloc, yveloc, v, uold;

        // Time loop
        for (int n = 0; n*dt <= 1.0; n++){
            // Mesh adaptation
            if (j > 20){
                th = adaptmesh(th, u, verbosity=1, abserror=1, nbjacoby=2,
                err=0.001, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
                splitpbedge=1, maxsubdiv=5, rescaling=1);
                j = 0;
                xveloc = -x*r + x*sigx^2 + x*rho*sigx*sigy/2;
                yveloc = -y*r + y*sigy^2 + y*rho*sigx*sigy/2;
                u = u;
            }

            // Update
            uold = u;

            // Solve
            solve eq1(u, v, init=j, solver=LU)
                = int2d(th)(
                      u*v*(r+1/dt)
                    + dx(u)*dx(v)*(x*sigx)^2/2
                    + dy(u)*dy(v)*(y*sigy)^2/2
                    + (dy(u)*dx(v) + dx(u)*dy(v))*rho*sigx*sigy*x*y/2
                )
                - int2d(th)(
                      v*convect([xveloc, yveloc], dt, uold)/dt
                )
                + on(2, 3, u=0)
                ;

            // Update
            j = j+1;
        };

        // Plot
        plot(u, wait=true, value=true);

    Results are shown on :numref:`figEvolutionBlackSholes1` and :numref:`figEvolutionBlackSholes2`.

    .. subfigstart::

    .. _figEvolutionBlackSholes1:

    .. figure:: images/EvolutionProblems_BlackSholes.png
        :width: 90%
        :alt: EvolutionProblems_BlackSholes

        The adapted triangulation

    .. _figEvolutionBlackSholes2:

    .. figure:: images/EvolutionProblems_BlackSholes2.png
        :width: 90%
        :alt: EvolutionProblems_BlackSholes2

        The level line of the European basquet put option

    .. subfigend::
        :width: 0.49
        :alt: EvolutionProblems_BlackSholes
        :label: EvolutionProblems_BlackSholes
