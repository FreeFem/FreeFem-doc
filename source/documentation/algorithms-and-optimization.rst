.. role:: freefem(code)
   :language: freefem

Algorithms & Optimization
=========================

Conjugate Gradient/GMRES
------------------------

Suppose we want to solve the Euler problem (here :math:`x` has nothing to do with the reserved variable for the first coordinate in **FreeFEM**):

find :math:`x\in \mathbb{R}^n` such that

.. math::
   \nabla J(x) = \left(\frac{\partial J}{\partial x_i} (\mathbf{x})\right) = 0
   :label: eqndJ

where :math:`J` is a function (to minimize for example) from :math:`\mathbb{R}^n` to :math:`\mathbb{R}`.

If the function is convex we can use the conjugate gradient algorithm to solve the problem, and we just need the function (named :freefem:`dJ` for example) which computes :math:`\nabla J`, so the parameters are the name of that function with prototype :freefem:`func real[int] dJ(real[int] &xx);` which computes :math:`\nabla J`, and a vector :freefem:`x` of type (of course the number 20 can be changed) :freefem:`real[int] x(20);` to initialize the process and get the result.

Given an initial value :math:`\mathbf{x}^{(0)}`, a maximum number :math:`i_{\max}` of iterations, and an error tolerance :math:`0<\epsilon<1`:

Put :math:`\mathbf{x}=\mathbf{x}^{(0)}` and write

.. code-block:: freefem
   :linenos:

   NLCG(dJ, x, precon=M, nbiter=imax, eps=epsilon, stop=stopfunc);

will give the solution of :math:`\mathbf{x}` of :math:`\nabla J(\mathbf{x})=0`.
We can omit parameters :freefem:`precon, nbiter, eps, stop`.
Here :math:`M` is the preconditioner whose default is the identity matrix.

The stopping test is

.. math::
   \|\nabla J(\mathbf{x})\|_P\le \epsilon\| \nabla J(\mathbf{x}^{(0)})\|_P

Writing the minus value in :freefem:`eps=`, i.e.,

.. code-block:: freefem
   :linenos:

   NLCG(dJ, x, precon=M, nbiter=imax, eps=-epsilon);

We can use the stopping test:

.. math::

   \| \nabla J(\mathbf{x})\|_P^2\le \epsilon

The parameters of these three functions are:

-  :freefem:`nbiter=` set the number of iteration (by default 100)
-  :freefem:`precon=` set the preconditioner function (:freefem:`P` for example) by default it is the identity, note the prototype is :freefem:`func real[int] P(real[int] &x)`.
-  :freefem:`eps=` set the value of the stop test :math:`\varepsilon` (:math:`=10^{-6}` by default) if positive then relative test :math:`||\nabla J(x)||_P\leq \varepsilon||\nabla J(x_0)||_P`, otherwise the absolute test is :math:`||\nabla J(x)||_P^2\leq |\varepsilon|`.
-  :freefem:`veps=` set and return the value of the stop test, if positive, then relative test is :math:`||\nabla J(x)||_P\leq \varepsilon||\nabla J(x_0)||_P`, otherwise the absolute test is :math:`||\nabla J(x)||_P^2\leq |\varepsilon|`.
   The return value is minus the real stop test (remark: it is useful in loop).
-  :freefem:`stop=` :freefem:`stopfunc` add your test function to stop before the :freefem:`eps` criterion. The prototype for the function :freefem:`stopfunc` is

   .. code-block:: freefem
      :linenos:

      func bool stopfunc(int iter, real[int] u, real[int] g)

   where :freefem:`u` is the current solution, and :freefem:`g`, the current gradient, is not preconditioned.

.. tip:: :ref:`Algorithms.edp <exampleAlgorithms>`

   For a given function :math:`b`, let us find the minimizer :math:`u` of the function

   .. math::
       \begin{array}{rcl}
          J(u) &=& \frac{1}{2}\int_{\Omega} f(|\nabla u|^2) - \int_{\Omega} u b \\
          f(x) &=& ax + x-\ln(1+x), \quad f'(x) = a+\frac{x}{1+x}, \quad f''(x) = \frac{1}{(1+x)^2}
       \end{array}

   under the boundary condition :math:`u=0` on :math:`\partial\Omega`.

   .. code-block:: freefem
      :linenos:

      fespace Ph(Th, P0);
      Ph alpha; //store df(|nabla u|^2)

      // The functionn J
      //J(u) = 1/2 int_Omega f(|nabla u|^2) - int_Omega u b
      func real J (real[int] & u){
         Vh w;
         w[] = u;
         real r = int2d(Th)(0.5*f(dx(w)*dx(w) + dy(w)*dy(w)) - b*w);
         cout << "J(u) = " << r << " " << u.min << " " << u.max << endl;
         return r;
      }

      // The gradiant of J
      func real[int] dJ (real[int] & u){
         Vh w;
         w[] = u;
         alpha = df(dx(w)*dx(w) + dy(w)*dy(w));
         varf au (uh, vh)
            = int2d(Th)(
                 alpha*(dx(w)*dx(vh) + dy(w)*dy(vh))
               - b*vh
            )
            + on(1, 2, 3, 4, uh=0)
            ;

         u = au(0, Vh);
         return u; //warning: no return of local array
      }

   We also want to construct a preconditioner :math:`C` with solving the problem:

   find :math:`u_h \in V_{0h}` such that:

   .. math::
      \forall v_h \in V_{0h}, \quad \int_\Omega \alpha \nabla u_h . \nabla v_h = \int_\Omega b v_h

   where :math:`\alpha=f'(|\nabla u|^2)`.

   .. code-block:: freefem
      :linenos:

      alpha = df(dx(u)*dx(u) + dy(u)*dy(u));
      varf alap (uh, vh)
         = int2d(Th)(
              alpha*(dx(uh)*dx(vh) + dy(uh)*dy(vh))
         )
         + on(1, 2, 3, 4, uh=0)
         ;

      varf amass(uh, vh)
         = int2d(Th)(
              uh*vh
         )
         + on(1, 2, 3, 4, uh=0)
         ;

      matrix Amass = amass(Vh, Vh, solver=CG);
      matrix Alap= alap(Vh, Vh, solver=Cholesky, factorize=1);

      // Preconditionner
      func real[int] C(real[int] & u){
         real[int] w = u;
         u = Alap^-1*w;
         return u; //warning: no return of local array variable
      }

   To solve the problem, we make 10 iterations of the conjugate gradient, recompute the preconditioner and restart the conjugate gradient:

   .. code-block:: freefem
      :linenos:

      int conv=0;
      for(int i = 0; i < 20; i++){
         conv = NLCG(dJ, u[], nbiter=10, precon=C, veps=eps, verbosity=5);
         if (conv) break;

         alpha = df(dx(u)*dx(u) + dy(u)*dy(u));
         Alap = alap(Vh, Vh, solver=Cholesky, factorize=1);
         cout << "Restart with new preconditionner " << conv << ", eps =" << eps << endl;
      }

      // Plot
      plot (u, wait=true, cmm="solution with NLCG");

For a given symmetric positive matrix :math:`A`, consider the quadratic form

.. math::

   J(\mathbf{x})=\frac{1}{2}\mathbf{x}^TA\mathbf{x}-\mathbf{b}^T\mathbf{x}

then :math:`J(\mathbf{x})` is minimized by the solution :math:`\mathbf{x}` of :math:`A\mathbf{x}=\mathbf{b}`.
In this case, we can use the function :freefem:`AffineCG`

.. code-block:: freefem
   :linenos:

   AffineCG(A, x, precon=M, nbiter=imax, eps=epsilon, stop=stp);

If :math:`A` is not symmetric, we can use GMRES(Generalized Minimum Residual) algorithm by

.. code-block:: freefem
   :linenos:

   AffineGMRES(A, x, precon=M, nbiter=imax, eps=epsilon);

Also, we can use the non-linear version of GMRES algorithm (the function :math:`J` is just convex)

.. code-block:: freefem
   :linenos:

   AffineGMRES(dJ, x, precon=M, nbiter=imax, eps=epsilon);

For the details of these algorithms, refer to [PIRONNEAU1998]_, Chapter IV, 1.3.

Algorithms for Unconstrained Optimization
-----------------------------------------

Two algorithms of COOOL package are interfaced with the Newton Raphson method (called :freefem:`Newton`) and the :freefem:`BFGS` method.
These two are directly available in **FreeFEM** (no dynamical link to load).
Be careful with these algorithms, because their implementation uses full matrices.
We also provide several optimization algorithms from the `NLopt library <https://nlopt.readthedocs.io/en/latest/>`__ as well as an interface for Hansen’s implementation of CMAES (a MPI version of this one is also available).

Example of usage for BFGS or CMAES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. tip:: BFGS

   .. code-block:: freefem
      :linenos:

      real[int] b(10), u(10);

      //J
      func real J (real[int] & u){
         real s = 0;
         for (int i = 0; i < u.n; i++)
            s += (i+1)*u[i]*u[i]*0.5 - b[i]*u[i];
         if (debugJ)
            cout << "J = " << s << ", u = " << u[0] << " " << u[1] << endl;
         return s;
      }

      //the gradiant of J (this is a affine version (the RHS is in)
      func real[int] DJ (real[int] &u){
         for (int i = 0; i < u.n; i++)
            u[i] = (i+1)*u[i];
         if (debugdJ)
            cout << "dJ: u = " << u[0] << " " << u[1] << " " << u[2] << endl;
         u -= b;
         if (debugdJ)
            cout << "dJ-b: u = " << u[0] << " " << u[1] << " " << u[2] << endl;
         return u; //return of global variable ok
      }

      b=1;
      u=2;
      BFGS(J, DJ, u, eps=1.e-6, nbiter=20, nbiterline=20);
      cout << "BFGS: J(u) = " << J(u) << ", err = " << error(u, b) << endl;

It is almost the same a using the CMA evolution strategy except, that since it is a derivative free optimizer, the :freefem:`dJ` argument is omitted and there are some other named parameters to control the behavior of the algorithm.
With the same objective function as above, an example of utilization would be (see :ref:`CMAES Variational inequality <exampleCMAESVariationalInequality>` for a complete example):

.. code-block:: freefem
   :linenos:

   load "ff-cmaes"
   //define J, u, ...
   real min = cmaes(J, u, stopTolFun=1e-6, stopMaxIter=3000);
   cout << "minimum value is " << min << " for u = " << u << endl;

This algorithm works with a normal multivariate distribution in the parameters space and tries to adapt its covariance matrix using the information provided by the successive function evaluations (see `NLopt documentation <https://nlopt.readthedocs.io/en/latest/>`__ for more details).
Therefore, some specific parameters can be passed to control the starting distribution, size of the sample generations, etc…
Named parameters for this are the following:

-  :freefem:`seed=` Seed for random number generator (:freefem:`val` is an integer).
   No specified value will lead to a clock based seed initialization.
-  :freefem:`initialStdDev=` Value for the standard deviations of the initial covariance matrix ( :freefem:`val` is a real).
   If the value :math:`\sigma` is passed, the initial covariance matrix will be set to :math:`\sigma I`.
   The expected initial distance between initial :math:`X` and the :math:`argmin` should be roughly initialStdDev. Default is 0.3.
-  :freefem:`initialStdDevs=` Same as above except that the argument is an array allowing to set a value of the initial standard deviation for each parameter.
   Entries differing by several orders of magnitude should be avoided (if it can’t be, try rescaling the problem).
-  :freefem:`stopTolFun=` Stops the algorithm if function value differences are smaller than the passed one, default is :math:`10^{-12}`.
-  :freefem:`stopTolFunHist=` Stops the algorithm if function value differences from the best values are smaller than the passed one, default is 0 (unused).
-  :freefem:`stopTolX=` Stopping criteria is triggered if step sizes in the parameters space are smaller than this real value, default is 0.
-  :freefem:`stopTolXFactor=` Stopping criteria is triggered when the standard deviation increases more than this value. The default value is :math:`10^{3}`.
-  :freefem:`stopMaxFunEval=` Stops the algorithm when :freefem:`stopMaxFunEval` function evaluations have been done.
   Set to :math:`900(n+3)^{2}` by default, where :math:`n` is the parameters space dimension.
-  :freefem:`stopMaxIter=` Integer stopping the search when :freefem:`stopMaxIter` generations have been sampled.
   Unused by default.
-  :freefem:`popsize=` Integer value used to change the sample size.
   The default value is :math:`4+ \lfloor 3\ln (n) \rfloor`.
   Increasing the population size usually improves the global search capabilities at the cost of, at most, a linear reduction of the convergence speed with respect to :freefem:`popsize`.
-  :freefem:`paramFile=` This :freefem:`string` type parameter allows the user to pass all the parameters using an extern file, as in Hansen’s original code.
   More parameters related to the CMA-ES algorithm can be changed with this file.
   Note that the parameters passed to the CMAES function in the **FreeFEM** script will be ignored if an input parameters file is given.

IPOPT
-----

The :freefem:`ff-Ipopt` package is an interface for the `IPOPT <https://projects.coin-or.org/Ipopt>`__ [WÄCHTER2006]_ optimizer.
IPOPT is a software library for large scale, non-linear, constrained optimization.
It implements a primal-dual interior point method along with filter method based line searches.

IPOPT needs a direct sparse symmetric linear solver.
If your version of **FreeFEM** has been compiled with the :freefem:`--enable-downlad` tag, it will automatically be linked with a sequential version of MUMPS.
An alternative to MUMPS would be to download the HSL subroutines (see `Compiling and Installing the Java Interface JIPOPT <https://projects.coin-or.org/Ipopt/wiki/JavaInterface>`__) and place them in the :freefem:`/ipopt/Ipopt-3.10.2/ThirdParty/HSL` directory of the **FreeFEM** downloads folder before compiling.

Short description of the algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section, we give a very brief glimpse at the underlying mathematics of IPOPT.
For a deeper introduction on interior methods for nonlinear smooth optimization, one may consult [FORSGREN2002]_, or [WÄCHTER2006]_ for more IPOPT specific elements.
IPOPT is designed to perform optimization for both equality and inequality constrained problems.
However, nonlinear inequalities are rearranged before the beginning of the optimization process in order to restrict the panel of nonlinear constraints to those of the equality kind.
Each nonlinear inequality is transformed into a pair of simple bound inequalities and nonlinear equality constraints by the introduction of as many slack variables as is needed : :math:`c_{i}(x)\leq 0` becomes :math:`c_{i}(x) + s_{i} = 0` and :math:`s_{i}\leq 0`, where :math:`s_{i}` is added to the initial variables of the problems :math:`x_{i}`.
Thus, for convenience, we will assume that the minimization problem does not contain any nonlinear inequality constraint.
It means that, given a function :math:`f:\mathbb{R}^{n}\mapsto\mathbb{R}`, we want to find:

.. math::
   x_{0} = \underset{x\in V}{\operatorname{argmin}} f(x) \\
   \mathrm{with}\ V = \left\lbrace x\in\mathbb{R}^{n}\ \vert\ c(x)= 0 \ \text{and}\ x_{l}\leq x\leq x_{u}\right\rbrace
   :label: minimproblem

Where :math:`c:\mathbb{R}^{n}\rightarrow\mathbb{R}^{m}` and :math:`x_{l},x_{u}\in\mathbb{R}^{n}` and inequalities hold componentwise.
The :math:`f` function as well as the constraints :math:`c` should be twice-continuously differentiable.

As a barrier method, interior points algorithms try to find a Karush-Kuhn-Tucker point for :eq:`minimproblem` by solving a sequence of problems, unconstrained with respect to the inequality constraints, of the form:

.. math::
   \mathrm{for\ a\ given\ }\mu > 0,\ \mathrm{find}\ x_{\mu} = \underset{x\in\mathbb{R}^{n}\ \vert\ c(x)=0}{\operatorname{argmin}}\ B(x,\mu)
   :label: barrier

Where :math:`\mu` is a positive real number and

.. math::
   B(x,\mu) = f(x) - \displaystyle{\mu\sum_{i=1}^{n} \ln (x_{u,i}-x_{i})} - \displaystyle{\mu\sum_{i=1}^{m} \ln(x_{i}-x_{l,i})}

The remaining equality constraints are handled with the usual Lagrange multipliers method.
If the sequence of barrier parameters :math:`\mu` converge to 0, intuition suggests that the sequence of minimizers of :eq:`barrier` converge to a local constrained minimizer of :eq:`minimproblem`.
For a given :math:`\mu`, :eq:`barrier` is solved by finding :math:`(x_{\mu},\lambda_{\mu})\in\mathbb{R}^{n}\times\mathbb{R}^{m}` such that:

.. math::
    \nabla B(x_{\mu},\mu) + \displaystyle{\sum_{i=1}^{m}\lambda_{\mu,i}\nabla c_{i}(x_{\mu})}= \nabla B(x_{\mu},\mu) + J_{c}(x_{\mu})^{T}\lambda_{\mu}&= 0\\
    c(x_{\mu}) &= 0
    :label: muproblem

The derivations for :math:`\nabla B` only holds for the :math:`x` variables, so that:

.. math::
   \nabla B(x,\mu) = \nabla f(x) + \left(\begin{matrix}\mu/(x_{u,1}-x_{1}) \\ \vdots \\ \mu/(x_{u,n}-x_{n})\end{matrix}\right) - \left(\begin{matrix}\mu/(x_{1}-x_{l,1}) \\ \vdots \\ \mu/(x_{n}-x_{l,n})\end{matrix}\right)

If we respectively call :math:`z_{u}(x,\mu) = \left(\mu/(x_{u,1}-x_{1}),\dots, \mu/(x_{u,n}-x_{n})\right)` and :math:`z_{l}(x,\mu)` the other vector appearing in the above equation, then the optimum :math:`(x_{\mu},\lambda_{\mu})` satisfies:

.. math::
   \nabla f(x_{\mu}) + J_{c}(x_{\mu})^{T}\lambda_{\mu}+ z_{u}(x_{\mu},\mu) - z_{l}(x_{\mu},\mu) = 0 \quad \text{and} \quad c(x_{\mu}) = 0
   :label: muproblemlambda

In this equation, the :math:`z_l` and :math:`z_u` vectors seem to play the role of Lagrange multipliers for the simple bound inequalities, and indeed, when :math:`\mu\rightarrow 0`, they converge toward some suitable Lagrange multipliers for the KKT conditions, provided some technical assumptions are fulfilled (see [FORSGREN2002]_).

Equation :eq:`muproblemlambda` is solved by performing a Newton method in order to find a solution of :eq:`muproblem` for each of the decreasing values of :math:`\mu`.
Some order 2 conditions are also taken into account to avoid convergence to local maximizers, see [FORSGREN2002]_ for details about them.
In the most classic IP algorithms, the Newton method is directly applied to :eq:`muproblem`.
This is in most case inefficient due to frequent computation of infeasible points.
These difficulties are avoided in Primal-Dual interior point methods where :eq:`muproblem` is transformed into an extended system where :math:`z_u` and :math:`z_l` are treated as unknowns and the barrier problems are finding :math:`(x,\lambda,z_u,z_l)\in\mathbb{R}^n\times\mathbb{R}^m\times\mathbb{R}^n\times\mathbb{R}^n` such that:

.. math::
   \left\lbrace
   \begin{array}{rcl}
      \nabla f(x) + J_{c}(x)^{T}\lambda+ z_{u} - z_{l} & = & 0 \\
      c(x) & = & 0 \\
      (X_u - X) z_u - \mu e & = & 0 \\
      (X - X_l) z_l - \mu e & = & 0
   \end{array}
   \right.
   :label: PrimalDualIPBarrierProblem

Where if :math:`a` is a vector of :math:`\mathbb{R}^n`, :math:`A` denotes the diagonal matrix :math:`A=(a_i \delta_{ij})_{1\leq i,j\leq n}` and :math:`e\in\mathbb{R}^{n} = (1,1,\dots,1)`.
Solving this nonlinear system by the Newton method is known as being the *primal-dual* interior point method.
Here again, more details are available in [FORSGREN2002]_.
Most actual implementations introduce features in order to globalize the convergence capability of the method, essentially by adding some line-search steps to the Newton algorithm, or by using trust regions.
For the purpose of IPOPT, this is achieved by a *filter line search* methods, the details of which can be found in [WÄCHTER2006]_.

More IPOPT specific features or implementation details can be found in [WÄCHTER2006]_.
We will just retain that IPOPT is a smart Newton method for solving constrained optimization problems, with global convergence capabilities due to a robust line search method (in the sense that the algorithm will converge no matter the initializer).
Due to the underlying Newton method, the optimization process requires expressions of all derivatives up to the order 2 of the fitness function as well as those of the constraints.
For problems whose Hessian matrices are difficult to compute or lead to high dimensional dense matrices, it is possible to use a BFGS approximation of these objects at the cost of a much slower convergence rate.

IPOPT in **FreeFEM**
~~~~~~~~~~~~~~~~~~~~~~~~~~

Calling the IPOPT optimizer in a **FreeFEM** script is done with the :freefem:`IPOPT` function included in the :freefem:`ff-Ipopt` dynamic library.
IPOPT is designed to solve constrained minimization problems in the form:

.. math::
   \mathrm{find} & x_{0} = \underset{x\in\mathbb{R}^{n}}{\operatorname{argmin}} f(x) \\
   \mathrm{s.t.} & \left\lbrace
   \begin{array}{l r}
      \forall i\leq n,\ x_{i}^{\mathrm{lb}}\leq x_{i}\leq x_{i}^{\mathrm{ub}} & \mathrm{\ (simple\ bounds)} \\
      \forall i\leq m,\ c_{i}^{\mathrm{lb}}\leq c_{i}(x)\leq c_{i}^{\mathrm{ub}} & \mathrm{(constraints\ functions)}
   \end{array}
   \right.

Where :math:`\mathrm{ub}` and :math:`\mathrm{lb}` stand for "upper bound" and "lower bound".
If for some :math:`i, 1\leq i\leq m` we have :math:`c_{i}^{\mathrm{lb}} = c_{i}^{\mathrm{ub}}`, it means that :math:`c_{i}` is an equality constraint, and an inequality one if :math:`c_{i}^{\mathrm{lb}} < c_{i}^{\mathrm{ub}}`.

There are different ways to pass the fitness function and constraints.
The more general one is to define the functions using the keyword :freefem:`func`.
Any returned matrix must be a sparse one (type :freefem:`matrix`, not a :freefem:`real[int,int]`):

.. code-block:: freefem
   :linenos:

   func real J (real[int] &X) {...} //Fitness Function, returns a scalar
   func real[int] gradJ (real[int] &X) {...} //Gradient is a vector

   func real[int] C (real[int] &X) {...} //Constraints
   func matrix jacC (real[int] &X) {...} //Constraints Jacobian

.. warning:: In the current version of **FreeFEM**, returning a :freefem:`matrix` object that is local to a function block leads to undefined results.
   For each sparse matrix returning function you define, an extern matrix object has to be declared, whose associated function will overwrite and return on each call.
   Here is an example for :freefem:`jacC`:

   .. code-block:: freefem
      :linenos:

      matrix jacCBuffer; //just declare, no need to define yet
      func matrix jacC (real[int] &X){
         ...//fill jacCBuffer
         return jacCBuffer;
      }

.. warning:: IPOPT requires the structure of each matrix at the initialization of the algorithm.
   Some errors may occur if the matrices are not constant and are built with the :freefem:`matrix A = [I,J,C]` syntax, or with an intermediary full matrix (:freefem:`real[int,int]`), because any null coefficient is discarded during the construction of the sparse matrix.
   It is also the case when making matrices linear combinations, for which any zero coefficient will result in the suppression of the matrix from the combination.
   Some controls are available to avoid such problems.
   Check the named parameter descriptions (:freefem:`checkindex`, :freefem:`structhess` and :freefem:`structjac` can help).
   We strongly advice to use :freefem:`varf` as much as possible for the matrix forging.

The Hessian returning function is somewhat different because it has to
be the Hessian of the Lagrangian function:

.. math::
   (x,\sigma_{f},\lambda)\mapsto\sigma_{f}\nabla^{2}f(x)+\displaystyle{\sum_{i=1}^{m}\lambda_{i}\nabla^{2}c_{i}(x)}\ \mathrm{ where }\ \lambda\in\mathbb{R}^{m}\ \mathrm{ and }\ \sigma\in\mathbb{R}

Your Hessian function should then have the following prototype:

.. code-block:: freefem
   :linenos:

   matrix hessianLBuffer; //Just to keep it in mind
   func matrix hessianL (real[int] &X, real sigma, real[int] &lambda){...}

If the constraints functions are all affine, or if there are only simple bound constraints, or no constraint at all, the Lagrangian Hessian is equal to the fitness function Hessian, one can then omit the :freefem:`sigma` and :freefem:`lambda` parameters:

.. code-block:: freefem
   :linenos:

   matrix hessianJBuffer;
   func matrix hessianJ (real[int] &X){...} //Hessian prototype when constraints are affine

When these functions are defined, IPOPT is called this way:

.. code-block:: freefem
   :linenos:

   real[int] Xi = ... ; //starting point
   IPOPT(J, gradJ, hessianL, C, jacC, Xi, /*some named parameters*/);

If the Hessian is omitted, the interface will tell IPOPT to use the (L)BFGS approximation (it can also be enabled with a named parameter, see further).
Simple bound or unconstrained problems do not require the constraints part, so the following expressions are valid:

.. code-block:: freefem
   :linenos:

   IPOPT(J, gradJ, C, jacC, Xi, ... ); //IPOPT with BFGS
   IPOPT(J, gradJ, hessianJ, Xi, ... ); //Newton IPOPT without constraints
   IPOPT(J, gradJ, Xi, ... ); //BFGS, no constraints

Simple bounds are passed using the :freefem:`lb` and :freefem:`ub` named parameters, while constraint bounds are passed with the :freefem:`clb` and :freefem:`cub` ones.
Unboundedness in some directions can be achieved by using the :math:`1e^{19}` and :math:`-1e^{19}` values that IPOPT recognizes as :math:`+\infty` and :math:`-\infty`:

.. code-block:: freefem
   :linenos:

   real[int] xlb(n), xub(n), clb(m), cub(m);
   //fill the arrays...
   IPOPT(J, gradJ, hessianL, C, jacC, Xi, lb=xlb, ub=xub, clb=clb, cub=cub, /*some other named parameters*/);

**P2 fitness function and affine constraints function :** In the case where the fitness function or constraints function can be expressed respectively in the following forms:

.. math::
   \begin{array}{c c}
       \forall x\in\mathbb{R}^{n},\ f(x) = \frac{1}{2}\left\langle Ax,x \right\rangle + \left\langle b,x\right\rangle & (A,b)\in\mathcal{M}_{n,n}(\mathbb{R})\times\mathbb{R}^{n} \\
       \mathrm{or} ,\ C(x) = Ax + b & (A,b)\in\mathcal{M}_{m,n}(\mathbb{R})\times\mathbb{R}^{m}
   \end{array}

where :math:`A` and :math:`b` are constant, it is possible to directly pass the :math:`(A,b)` pair instead of defining 3 (or 2) functions.
It also indicates to IPOPT that some objects are constant and that they have to be evaluated only once, thus avoiding multiple copies of the same matrix.
The syntax is:

.. code-block:: freefem
   :linenos:

   // Affine constraints with "standard" fitness function
   matrix A = ... ; //linear part of the constraints
   real[int] b = ... ; //constant part of constraints
   IPOPT(J, gradJ, hessianJ, [A, b], Xi, /*bounds and named parameters*/);
   //[b, A] would work as well.

Note that if you define the constraints in this way, they don’t contribute to the Hessian, so the Hessian should only take one :freefem:`real[int]` as an argument.

.. code-block:: freefem
   :linenos:

   // Affine constraints and P2 fitness func
   matrix A = ... ; //bilinear form matrix
   real[int] b = ... ; //linear contribution to f
   matrix Ac = ... ; //linear part of the constraints
   real[int] bc = ... ; //constant part of constraints
   IPOPT([A, b], [Ac, bc], Xi, /*bounds and named parameters*/);

If both objective and constraint functions are given this way, it automatically activates the IPOPT ``mehrotra_algorithm`` option (better for linear and quadratic programming according to the documentation).
Otherwise, this option can only be set through the option file (see the named parameters section).

A false case is the one of defining :math:`f` in this manner while using standard functions for the constraints:

.. code-block:: freefem
   :linenos:

   matrix A = ... ; //bilinear form matrix
   real[int] b = ... ; //linear contribution to f
   func real[int] C(real[int] &X){...} //constraints
   func matrix jacC(real[int] &X){...} //constraints Jacobian
   IPOPT([A, b], C, jacC, Xi, /*bounds and named parameters*/);

Indeed, when passing :freefem:`[A, b]` in order to define :math:`f`, the Lagrangian Hessian is automatically built and has the constant :math:`x \mapsto A` function, with no way to add possible constraint contributions, leading to incorrect second order derivatives.
So, a problem should be defined like that in only two cases:

1. constraints are nonlinear but you want to use the BFGS mode (then add :freefem:`bfgs=1` to the named parameter),
2. constraints are affine, but in this case, compatible to pass in the same way

Here are some other valid definitions of the problem (cases when :math:`f` is a pure quadratic or linear form, or :math:`C` a pure linear function, etc…):

.. code-block:: freefem
   :linenos:

   // Pure quadratic f - A is a matrix
   IPOPT(A, /*constraints arguments*/, Xi, /*bound and named parameters*/);
   // Pure linear f - b is a real[int]
   IPOPT(b, /*constraints arguments*/, Xi, /*bound and named parameters*/);
   // Linear constraints - Ac is a matrix
   IPOPT(/*fitness function arguments*/, Ac, Xi, /*bound and named parameters*/);

**Returned Value :** The :freefem:`IPOPT` function returns an error code of type :freefem:`int`.
A zero value is obtained when the algorithm succeeds and positive values reflect the fact that IPOPT encounters minor troubles.
Negative values reveal more problematic cases.
The associated IPOPT return tags are listed in the table below.
The `IPOPT pdf documentation <https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw>`__ provides a more accurate description of these return statuses:

+------------------------------------------+------------------------------------+
| Success                                  | Failures                           |
+==========================================+====================================+
| 0 ``Solve_Succeeded``                    |                                    |
+------------------------------------------+------------------------------------+
| 1 ``Solved_To_Acceptable_Level``         | -1 ``Maximum_Iterations_Exceeded`` |
+------------------------------------------+------------------------------------+
| 2 ``Infeasible_Problem_Detected``        | -2 ``Restoration_Failed``          |
+------------------------------------------+------------------------------------+
| 3 ``Search_Direction_Becomes_Too_Small`` | -3 ``Error_In_Step_Computation``   |
+------------------------------------------+------------------------------------+
| 4 ``Diverging_Iterates``                 | -4 ``Maximum_CpuTime_Exceeded``    |
+------------------------------------------+------------------------------------+
| 5 ``User_Requested_Stop``                |                                    |
+------------------------------------------+------------------------------------+
| 6 ``Feasible_Point_Found``               |                                    |
+------------------------------------------+------------------------------------+

+---------------------------------------+------------------------------------+
| Problem definition issues             | Critical errors                    |
+=======================================+====================================+
| -10 ``Not_Enough_Degrees_Of_Freedom`` | -100 ``Unrecoverable_Exception``   |
+---------------------------------------+------------------------------------+
| -11 ``Invalid_Problem_Definition``    | -101 ``NonIpopt_Exception_Thrown`` |
+---------------------------------------+------------------------------------+
| -12 ``Invalid_Option``                | -102 ``Insufficient_Memory``       |
+---------------------------------------+------------------------------------+
| -13 ``Invalid_Number_Detected``       | -199 ``Internal_Error``            |
+---------------------------------------+------------------------------------+

**Named Parameters :** The available named parameters in this interface are those we thought to be the most subject to variations from one optimization to another, plus a few that are interface specific.
Though, as one could see at `IPOPT Linear solver <https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_MA97_Linear_Solver>`__, there are many parameters that can be changed within IPOPT, affecting the algorithm behavior.
These parameters can still be controlled by placing an option file in the execution directory.
Note that `IPOPT’s pdf documentation <https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw>`__ may provides more information than the previously mentioned online version for certain parameters.
The in-script available parameters are:

-  :freefem:`lb`, :freefem:`ub` : :freefem:`real[int]` for lower and upper simple bounds upon the search variables must be of size :math:`n` (search space dimension).
   If two components of the same index in these arrays are equal then the corresponding search variable is fixed.
   By default IPOPT will remove any fixed variable from the optimization process and always use the fixed value when calling functions.
   It can be changed using the :freefem:`fixedvar` parameter.
-  :freefem:`clb`, :freefem:`cub` : :freefem:`real[int]` of size :math:`m` (number of constraints) for lower and upper constraints bounds.
   Equality between two components of the same index :math:`i` in :freefem:`clb` and :freefem:`cub` reflect an equality constraint.
-  :freefem:`structjacc` : To pass the greatest possible structure (indexes of non null coefficients) of the constraint Jacobians under the form :freefem:`[I,J]` where :freefem:`I` and :freefem:`J` are two integer arrays.
   If not defined, the structure of the constraint Jacobians, evaluated in :freefem:`Xi`, is used (no issue if the Jacobian is constant or always defined with the same :freefem:`varf`, hazardous if it is with a triplet array or if a full matrix is involved).
-  :freefem:`structhess` : Same as above but for the Hessian function (unused if :math:`f` is P2 or less and constraints are affine).
   Here again, keep in mind that it is the Hessian of the Lagrangian function (which is equal to the Hessian of :math:`f` only if constraints are affine).
   If no structure is given with this parameter, the Lagrangian Hessian is evaluated on the starting point, with :math:`\sigma=1` and :math:`\lambda = (1,1,\dots,1)` (it is safe if all the constraints and fitness function Hessians are constant or build with :freefem:`varf`, and here again it is less reliable if built with a triplet array or a full matrix).
-  :freefem:`checkindex` : A :freefem:`bool` that triggers a dichotomic index search when matrices are copied from **FreeFEM** functions to IPOPT arrays.
   It is used to avoid wrong index matching when some null coefficients are removed from the matrices by **FreeFEM**.
   It will not solve the problems arising when a too small structure has been given at the initialization of the algorithm.
   Enabled by default (except in cases where all matrices are obviously constant).
-  :freefem:`warmstart` : If set to :freefem:`true`, the constraints dual variables :math:`\lambda`, and simple bound dual variables are initialized with the values of the arrays passed to :freefem:`lm`, :freefem:`lz` and :freefem:`uz` named parameters (see below).
-  :freefem:`lm` : :freefem:`real[int]` of size :math:`m`, which is used to get the final values of the constraints dual variables :math:`\lambda` and/or initialize them in case of a warm start (the passed array is also updated to the last dual variables values at the end of the algorithm).
-  :freefem:`lz`, :freefem:`uz` : :freefem:`real[int]` of size :math:`n` to get the final values and/or initialize (in case of a warm start) the dual variables associated to simple bounds.
-  :freefem:`tol` : :freefem:`real`, convergence tolerance for the algorithm, the default value is :math:`10^{-8}`.
-  :freefem:`maxiter` : :freefem:`int`, maximum number of iterations with 3000 as default value.
-  :freefem:`maxcputime` : :freefem:`real` value, maximum runtime duration. Default is :math:`10^{6}` (almost 11 and a halfdays).
-  :freefem:`bfgs` : :freefem:`bool` enabling or not the (low-storage) BFGS approximation of the Lagrangian Hessian.
   It is set to false by default, unless there is no way to compute the Hessian with the functions that have been passed to IPOPT.
-  :freefem:`derivativetest` : Used to perform a comparison of the derivatives given to IPOPT with finite differences computation.
   The possible :freefem:`string` values are : :freefem:`"none"` (default), :freefem:`"first-order"`, :freefem:`"second-order"` and :freefem:`"only-second-order"`.
   The associated derivative error tolerance can be changed via the option file.
   One should not care about any error given by it before having tried, and failed, to perform a first optimization.
-  :freefem:`dth` : Perturbation parameter for the derivative test computations with finite differences.
   Set by default to :math:`10^{-8}`.
-  :freefem:`dttol` : Tolerance value for the derivative test error detection (default value unknown yet, maybe :math:`10^{-5}`).
-  :freefem:`optfile` : :freefem:`string` parameter to specify the IPOPT option file name.
   IPOPT will look for a :freefem:`ipopt.opt` file by default.
   Options set in the file will overwrite those defined in the **FreeFEM** script.
-  :freefem:`printlevel` : An :freefem:`int` to control IPOPT output print level, set to 5 by default, the possible values are from 0 to 12.
   A description of the output information is available in the `PDF documentation <https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw>`__ of IPOPT.
-  :freefem:`fixedvar` : :freefem:`string` for the definition of simple bound equality constraints treatment : use :freefem:`"make_parameter"` (default value) to simply remove them from the optimization process (the functions will always be evaluated with the fixed value for those variables), :freefem:`"make_constraint"` to treat them as any other constraint or :freefem:`"relax_bounds"` to relax fixing bound constraints.
-  :freefem:`mustrategy` : a :freefem:`string` to choose the update strategy for the barrier parameter :math:`\mu`.
   The two possible tags are :freefem:`"monotone"`, to use the monotone (Fiacco-McCormick) strategy, or :freefem:`"adaptive"` (default setting).
-  :freefem:`muinit` : :freefem:`real` positive value for the barrier parameter initialization.
   It is only relevant when :freefem:`mustrategy` has been set to :freefem:`monotone`.
-  :freefem:`pivtol` : :freefem:`real` value to set the pivot tolerance for the linear solver. A smaller number pivots for sparsity, a larger number pivots for stability.
   The value has to be in the :math:`[0,1]` interval and is set to :math:`10^{-6}` by default.
-  :freefem:`brf` : Bound relax factor: before starting the optimization, the bounds given by the user are relaxed.
   This option sets the factor for this relaxation.
   If it is set to zero, then the bound relaxation is disabled.
   This :freefem:`real` has to be positive and its default value is :math:`10^{-8}`.
-  :freefem:`objvalue` : An identifier to a :freefem:`real` type variable to get the last value of the objective function (best value in case of success).
-  :freefem:`mumin` : minimum value for the barrier parameter :math:`\mu`, a :freefem:`real` with :math:`10^{-11}` as default value.
-  :freefem:`linesearch` : A boolean which disables the line search when set to :freefem:`false`.
   The line search is activated by default.
   When disabled, the method becomes a standard Newton algorithm instead of a primal-dual system.
   The global convergence is then no longer assured, meaning that many initializers could lead to diverging iterates.
   But on the other hand, it can be useful when trying to catch a precise local minimum without having some out of control process making the iterate caught by some other near optimum.

Some short examples using IPOPT
-------------------------------

.. tip:: Ipopt variational inequality
   A very simple example consisting of, given two functions :math:`f` and :math:`g` (defined on :math:`\Omega\subset\mathbb{R}^{2}`), minimizing :math:`J(u) = \displaystyle{\frac{1}{2}\int_{\Omega} \vert\nabla u\vert^{2} - \int_{\Omega}fu}\ `, with :math:`u\leq g` almost everywhere:

   .. code-block:: freefem
      :linenos:

      // Solve
      //- Delta u = f
      //u < g
      //u = 0 on Gamma
      load "ff-Ipopt";

      // Parameters
      int nn = 20;
      func f = 1.; //rhs function
      real r = 0.03, s = 0.1;
      func g = r - r/2*exp(-0.5*(square(x-0.5) + square(y-0.5))/square(s));

      // Mesh
      mesh Th = square(nn, nn);

      // Fespace
      fespace Vh(Th, P2);
      Vh u = 0;
      Vh lb = -1.e19;
      Vh ub = g;

      // Macro
      macro Grad(u) [dx(u), dy(u)] //

      // Problem
      varf vP (u, v)
         = int2d(Th)(
              Grad(u)'*Grad(v)
         )
         - int2d(Th)(
              f*v
         )
         ;

   Here we build the matrix and second member associated to the function to fully and finally minimize it.
   The :freefem:`[A,b]` syntax for the fitness function is then used to pass it to IPOPT.

   .. code-block:: freefem
      :linenos:

      matrix A = vP(Vh, Vh, solver=CG);
      real[int] b = vP(0, Vh);

   We use simple bounds to impose the boundary condition :math:`u=0` on :math:`\partial\Omega`, as well as the :math:`u\leq g` condition.

   .. code-block:: freefem
      :linenos:

      varf vGamma (u, v) = on(1, 2, 3, 4, u=1);
      real[int] onGamma = vGamma(0, Vh);

      //warning: the boundary conditions are given with lb and ub on border
      ub[] = onGamma ? 0. : ub[];
      lb[] = onGamma ? 0. : lb[];

      // Solve
      IPOPT([A, b], u[], lb=lb[], ub=ub[]);

      // Plot
      plot(u);

.. tip:: Ipopt variational inequality 2

   Let :math:`\Omega` be a domain of :math:`\mathbb{R}^{2}`.
   :math:`f_{1}, f_{2}\in L^{2}(\Omega)` and :math:`g_{1}, g_{2} \in L^{2}(\partial\Omega)` four given functions with :math:`g_{1}\leq g_{2}` almost everywhere.
   We define the space:

   .. math::
      V = \left\lbrace (v_{1},v_{2})\in H^{1}(\Omega)^{2} ; v_{1}\vert_{\partial\Omega}=g_{1}, v_{2}\vert_{\partial\Omega}=g_{2}, v_{1}\leq v_{2}\ \mathrm{a.e.}\ \right\rbrace

   as well as the function :math:`J:H^{1}(\Omega)^{2}\longrightarrow \mathbb{R}`:

   .. math::
      J(v_{1},v_{2}) = \displaystyle{\frac{1}{2}\int_{\Omega}\vert\nabla v_{1}\vert^{2} - \int_{\Omega} f_{1}v_{1} + \frac{1}{2}\int_{\Omega}\vert\nabla v_{2}\vert^{2} - \int_{\Omega} f_{2}v_{2}}

   The problem entails finding (numerically) two functions :math:`(u_{1},u_{2}) = \underset{(v_{1},v_{2})\in V}{\operatorname{argmin}} J(v_{1},v_{2})`.

   .. code-block:: freefem
      :linenos:

      load "ff-Ipopt";

      // Parameters
      int nn = 10;
      func f1 = 10;//right hand side
      func f2 = -15;
      func g1 = -0.1;//Boundary condition functions
      func g2 = 0.1;

      // Mesh
      mesh Th = square(nn, nn);

      // Fespace
      fespace Vh(Th, [P1, P1]);
      Vh [uz, uz2] = [1, 1];
      Vh [lz, lz2] = [1, 1];
      Vh [u1, u2] = [0, 0]; //starting point

      fespace Wh(Th, [P1]);
      Wh lm=1.;

      // Macro
      macro Grad(u) [dx(u), dy(u)] //

      // Loop
      int iter=0;
      while (++iter){
         // Problem
         varf vP ([u1, u2], [v1, v2])
            = int2d(Th)(
                 Grad(u1)'*Grad(v1)
               + Grad(u2)'*Grad(v2)
            )
            - int2d(Th)(
                 f1*v1
               + f2*v2
            )
            ;

         matrix A = vP(Vh, Vh); //fitness function matrix
         real[int] b = vP(0, Vh); //and linear form

         int[int] II1 = [0], II2 = [1];//Constraints matrix
         matrix C1 = interpolate (Wh, Vh, U2Vc=II1);
         matrix C2 = interpolate (Wh, Vh, U2Vc=II2);
         matrix CC = -1*C1 + C2; // u2 - u1 > 0
         Wh cl = 0; //constraints lower bounds (no upper bounds)

         //Boundary conditions
         varf vGamma ([u1, u2], [v1, v2]) = on(1, 2, 3, 4, u1=1, u2=1);
         real[int] onGamma = vGamma(0, Vh);
         Vh [ub1, ub2] = [g1, g2];
         Vh [lb1, lb2] = [g1, g2];
         ub1[] = onGamma ? ub1[] : 1e19; //Unbounded in interior
         lb1[] = onGamma ? lb1[] : -1e19;

         Vh [uzi, uzi2] = [uz, uz2], [lzi, lzi2] = [lz, lz2];
         Wh lmi = lm;
         Vh [ui1, ui2] = [u1, u2];

         // Solve
         IPOPT([b, A], CC, ui1[], lb=lb1[], clb=cl[], ub=ub1[], warmstart=iter>1, uz=uzi[], lz=lzi[], lm=lmi[]);

         // Plot
         plot(ui1, ui2, wait=true, nbiso=60, dim=3);

         if(iter > 1) break;

         // Mesh adpatation
         Th = adaptmesh(Th, [ui1, ui2], err=0.004, nbvx=100000);
         [uz, uz2] = [uzi, uzi2];
         [lz, lz2] = [lzi, lzi2];
         [u1, u2] = [ui1, ui2];
         lm = lmi;
      }

   .. subfigstart::

   .. _figAlgoVarineqFill:

   .. figure:: images/VarIneqFill.jpg
      :alt: VarIneqFill
      :width: 90%

      Numerical Approximation of the Variational Inequality

   .. _figAlgoVarineqIso:

   .. figure:: images/VarIneqIso.jpg
      :width: 90%

      Numerical Approximation of the Variational Inequality

   .. subfigend::
      :width: 0.49
      :alt: VariationalInequality
      :label: VariationalInequality

      Variational inequality

3D constrained minimum surface with IPOPT
-----------------------------------------

Area and volume expressions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example is aimed at numerically solving some constrained minimum surface problems with the IPOPT algorithm.
We restrain to :math:`C^{k}` (:math:`k\geq 1`), closed, spherically parametrizable surfaces, i.e. surfaces :math:`S` such that:

.. math::
   \exists \rho \in C^{k}([0,2\pi ]\times[0,\pi] ) \vert
   S = \left\lbrace
   X = \left(
   \begin{array} {c}
    \rho(\theta,\phi) \\
    0 \\
    0
   \end{array}
   \right)
   , (\theta,\phi) \in [0,2\pi ]\times[0,\pi]
    \right\rbrace

Where the components are expressed in the spherical coordinate system.
Let’s call :math:`\Omega` the :math:`[0,2\pi ]\times[0,\pi]` angular parameters set.
In order to exclude self crossing and opened shapes, the following assumptions upon :math:`\rho` are made:

.. math::
   \rho \geq 0\ \ \mathrm{and}\ \ \forall \phi, \rho(0,\phi) = \rho(2\pi,\phi)

For a given function :math:`\rho` the first fundamental form (the metric) of the defined surface has the following matrix representation:

.. math::
   G =
   \left(
   \begin{array}{c c}
       \rho^{2}\sin^{2}(\phi) + (\partial_{\theta}\rho)^{2} &\partial_{\theta}\rho\partial_{\phi}\rho \\
       \partial_{\theta}\rho\partial_{\phi}\rho & \rho^{2} + (\partial_{\phi}\rho)^{2} \\
   \end{array}
   \right)
   :label: msfff

This metric is used to express the area of the surface.
Let :math:`g=\det(G)`, then we have:

.. math::
    \begin{array}{ll}
        \mathcal{A}(\rho) &= \int{\Omega}{\left\| \partial_{\theta} X \wedge \partial_{\phi} X \right\|} =\int{\Omega}{\sqrt{g}}\\
            &=\int{\Omega}{\sqrt{ \rho^{2}(\partial_{\theta}\rho)^{2} + \rho^{4}\sin^{2}(\phi) + \rho^{2}(\partial_{\phi}\rho)^{2}\sin^{2}(\phi)}d\theta d\phi}
    \end{array}
    :label: msarea

The volume of the space enclosed within the shape is easier to express:

.. math::
    \mathcal{V}(\rho)
    = \int{\Omega}{\int_{0}^{\rho(\theta,\phi)} r^{2}\sin(\phi) dr d\theta d\phi}
    = \frac{1}{3}\int{\Omega}{\rho^{3} \sin(\phi) d\theta d\phi}
    :label: msvolume

Derivatives
~~~~~~~~~~~

In order to use a Newton based interior point optimization algorithm, one must be able to evaluate the derivatives of :math:`\mathcal{A}` and :math:`\mathcal{V}` with respect to :math:`rho`.
Concerning the area, we have the following result:

.. math::
   \forall v\in C^{1}(\Omega) \ , \ \langle d\mathcal{A}(\rho),v\rangle
   = \int{\Omega}{\frac{1}{2} \frac{ d\bar{g}(\rho)(v)}{\sqrt{g}}d\theta d\phi }

Where :math:`\bar{g}` is the application mapping the :math:`(\theta,\phi) \mapsto g(\theta,\phi)` scalar field to :math:`\rho`.
This leads to the following expression, easy to transpose in a freefem script using:

.. math::
    \begin{array}{r c l}
        \forall v\in C^{1}(\Omega)& &\\
        \langle d\mathcal{A}(\rho),v\rangle &=& \int{\Omega}{ \left(2\rho^{3}\sin^{2}(\phi) + \rho(\partial_{\theta}\rho)^{2} + \rho(\partial_{\phi}\rho)^{2}\sin^{2}(\phi) \right) v} \\
        & & +\int{\Omega}{\ \rho^{2}\partial_{\theta}\rho\partial_{\theta} v\ + \ \rho^{2}\partial_{\phi}\rho\sin^{2}(\phi)\partial_{\phi} v }
    \end{array}
    :label: msdarea

With a similar approach, one can derive an expression for second order derivatives.
However, comporting no specific difficulties, the detailed calculus are tedious, the result is that these derivatives can be written using a :math:`3\times 3` matrix :math:`\mathbf{B}` whose coefficients are expressed in term of :math:`\rho` and its derivatives with respect to :math:`\theta` and :math:`\phi`, such that:

.. math::
   \forall (w,v)\in C^{1}(\Omega)\ ,\ d^{2}\mathcal{A}(\rho)(w,v) = \int{\Omega}
   {
      \left(\begin{array}{c c c} w & \partial_{\theta} w & \partial_{\phi} w \end{array}\right)
      \mathbf{B}
   } \left( \begin{array}{c} v \\ \partial_{\theta} v \\ \partial_{\phi} v \end{array} \right) d\theta d\phi
   :label: msd2area

Deriving the volume function derivatives is again an easier task.
We immediately get the following expressions:

.. math::
   \begin{array}{r c l}
      \forall v\ ,\ \langle d\mathcal{V}(\rho),v\rangle & = & \int{\Omega}{\rho^{2}\sin(\phi)v\ d\theta d\phi} \\
      \forall w,v\ , d^{2}\mathcal{V}(\rho)(w,v) & = & \int{\Omega}{2\rho\sin(\phi)wv\ d\theta d\phi}
   \end{array}
   :label: msdvolume

The problem and its script
~~~~~~~~~~~~~~~~~~~~~~~~~~

The whole code is available in :ref:`IPOPT minimal surface & volume example <exampleIPOPTMinimalSurfaceVolume>`.
We propose to solve the following problem:

.. tip::

   Given a positive function :math:`\rho_{\mathrm{object}}` piecewise continuous, and a scalar :math:`\mathcal{V}_{\mathrm{max}} > \mathcal{V}(\rho_{\mathrm{object}})`, find :math:`\rho_{0}` such that:

   .. math::
      \rho_{0} = \underset{\rho\in C^{1}(\Omega)}{\operatorname{argmin}}\ \mathcal{A}(\rho)\ ,\ \mathrm{s.t.}\ \rho_{0}\geq\rho_{\mathrm{object}} \ \mathrm{and\ } \mathcal{V}(\rho_{0})\leq \mathcal{V}_{\mathrm{max}}

   If :math:`\rho_{\mathrm{object}}` is the spherical parametrization of the surface of a 3-dimensional object (domain) :math:`\mathcal{O}`, it can be interpreted as finding the surface with minimum area enclosing the object with a given maximum volume. If :math:`\mathcal{V}_{\mathrm{max}}` is close to :math:`\mathcal{V}(\rho_{\mathrm{object}})`, so should be :math:`\rho_{0}` and :math:`\rho_{\mathrm{object}}`. With higher values of :math:`\mathcal{V}_{\mathrm{max}}`, :math:`\rho` should be closer to the unconstrained minimum surface surrounding :math:`\mathcal{O}` which is obtained as soon as :math:`\mathcal{V}_{\mathrm{max}} \geq \frac{4}{3}\pi \|\rho_{\mathrm{object}}\|_{\infty}^{3}` (sufficient but not necessary).

   It also could be interesting to solve the same problem with the constraint :math:`\mathcal{V}(\rho_{0})\geq \mathcal{V}_{\mathrm{min}}` which leads to a sphere when :math:`\mathcal{V}_{\mathrm{min}} \geq \frac{1}{6}\pi \mathrm{diam}(\mathcal{O})^{3}` and moves toward the solution of the unconstrained problem as :math:`\mathcal{V}_{\mathrm{min}}` decreases.

   We start by meshing the domain :math:`[0,2\pi]\times\ [0,\pi]`, then a periodic P1 finite elements space is defined.

   .. code-block:: freefem
      :linenos:

      load "msh3";
      load "medit";
      load "ff-Ipopt";

      // Parameters
      int nadapt = 3;
      real alpha = 0.9;
      int np = 30;
      real regtest;
      int shapeswitch = 1;
      real sigma = 2*pi/40.;
      real treshold = 0.1;
      real e = 0.1;
      real r0 = 0.25;
      real rr = 2-r0;
      real E = 1./(e*e);
      real RR = 1./(rr*rr);

      // Mesh
      mesh Th = square(2*np, np, [2*pi*x, pi*y]);

      // Fespace
      fespace Vh(Th, P1, periodic=[[2, y], [4, y]]);
      //Initial shape definition
      //outside of the mesh adaptation loop to initialize with the previous optimial shape found on further iterations
      Vh startshape = 5;

   We create some finite element functions whose underlying arrays will be used to store the values of dual variables associated to all the constraints in order to reinitialize the algorithm with it in the case where we use mesh adaptation. Doing so, the algorithm will almost restart at the accuracy level it reached before mesh adaptation, thus saving many iterations.

   .. code-block:: freefem
      :linenos:

      Vh uz = 1., lz = 1.;
      rreal[int] lm = [1];

   Then, follows the mesh adaptation loop, and a rendering function, :freefem:`Plot3D`, using 3D mesh to display the shape it is passed with :freefem:`medit` (the :freefem:`movemesh23` procedure often crashes when called with ragged shapes).

   .. code-block:: freefem
      :linenos:

      for(int kkk = 0; kkk < nadapt; ++kkk){
         int iter=0;
         func sin2 = square(sin(y));

         // A function which transform Th in 3d mesh (r=rho)
         //a point (theta,phi) of Th becomes ( r(theta,phi)*cos(theta)*sin(phi) , r(theta,phi)*sin(theta)*sin(phi) , r(theta,phi)*cos(phi) )
         //then displays the resulting mesh with medit
         func int Plot3D (real[int] &rho, string cmm, bool ffplot){
            Vh rhoo;
            rhoo[] = rho;
            //mesh sTh = square(np, np/2, [2*pi*x, pi*y]);
            //fespace sVh(sTh, P1);
            //Vh rhoplot = rhoo;
            try{
               mesh3 Sphere = movemesh23(Th, transfo=[rhoo(x,y)*cos(x)*sin(y), rhoo(x,y)*sin(x)*sin(y), rhoo(x,y)*cos(y)]);
               if(ffplot)
                  plot(Sphere);
               else
                  medit(cmm, Sphere);
            }
            catch(...){
               cout << "PLOT ERROR" << endl;
            }
            return 1;
         }
      }

   Here are the functions related to the area computation and its shape derivative, according to equations :eq:`msarea` and :eq:`msdarea`:

   .. code-block:: freefem
      :linenos:

      // Surface computation
      //Maybe is it possible to use movemesh23 to have the surface function less complicated
      //However, it would not simplify the gradient and the hessian
      func real Area (real[int] &X){
         Vh rho;
         rho[] = X;
         Vh rho2 = square(rho);
         Vh rho4 = square(rho2);
         real res = int2d(Th)(sqrt(rho4*sin2 + rho2*square(dx(rho)) + rho2*sin2*square(dy(rho))));
         ++iter;
         if(1)
            plot(rho, value=true, fill=true, cmm="rho(theta,phi) on [0,2pi]x[0,pi] - S="+res, dim=3);
         else
            Plot3D(rho[], "shape_evolution", 1);
         return res;
      }

      func real[int] GradArea (real[int] &X){
         Vh rho, rho2;
         rho[] = X;
         rho2[] = square(X);
         Vh sqrtPsi, alpha;
         {
            Vh dxrho2 = dx(rho)*dx(rho), dyrho2 = dy(rho)*dy(rho);
            sqrtPsi = sqrt(rho2*rho2*sin2 + rho2*dxrho2 + rho2*dyrho2*sin2);
            alpha = 2.*rho2*rho*sin2 + rho*dxrho2 + rho*dyrho2*sin2;
         }
         varf dArea (u, v)
            = int2d(Th)(
                 1./sqrtPsi * (alpha*v + rho2*dx(rho)*dx(v) + rho2*dy(rho)*sin2*dy(v))
            )
            ;

         real[int] grad = dArea(0, Vh);
         return grad;
      }

   The function returning the hessian of the area for a given shape is a bit blurry, thus we won't show here all of equation :eq:`msd2area` coefficients definition, they can be found in the :freefem:`edp` file.

   .. code-block:: freefem
      :linenos:

      matrix hessianA;
      func matrix HessianArea (real[int] &X){
         Vh rho, rho2;
         rho[] = X;
         rho2 = square(rho);
         Vh sqrtPsi, sqrtPsi3, C00, C01, C02, C11, C12, C22, A;
         {
            Vh C0, C1, C2;
            Vh dxrho2 = dx(rho)*dx(rho), dyrho2 = dy(rho)*dy(rho);
            sqrtPsi = sqrt( rho2*rho2*sin2 + rho2*dxrho2 + rho2*dyrho2*sin2);
            sqrtPsi3 = (rho2*rho2*sin2 + rho2*dxrho2 + rho2*dyrho2*sin2)*sqrtPsi;
            C0 = 2*rho2*rho*sin2 + rho*dxrho2 + rho*dyrho2*sin2;
            C1 = rho2*dx(rho);
            C2 = rho2*sin2*dy(rho);
            C00 = square(C0);
            C01 = C0*C1;
            C02 = C0*C2;
            C11 = square(C1);
            C12 = C1*C2;
            C22 = square(C2);
            A = 6.*rho2*sin2 + dxrho2 + dyrho2*sin2;
         }
         varf d2Area (w, v)
            =int2d(Th)(
                 1./sqrtPsi * (
                    A*w*v
                  + 2*rho*dx(rho)*dx(w)*v
                  + 2*rho*dx(rho)*w*dx(v)
                  + 2*rho*dy(rho)*sin2*dy(w)*v
                  + 2*rho*dy(rho)*sin2*w*dy(v)
                  + rho2*dx(w)*dx(v)
                  + rho2*sin2*dy(w)*dy(v)
               )
               + 1./sqrtPsi3 * (
                    C00*w*v
                  + C01*dx(w)*v
                  + C01*w*dx(v)
                  + C02*dy(w)*v
                  + C02*w*dy(v)
                  + C11*dx(w)*dx(v)
                  + C12*dx(w)*dy(v)
                  + C12*dy(w)*dx(v)
                  + C22*dy(w)*dy(v)
               )
            )
            ;
         hessianA = d2Area(Vh, Vh);
         return hessianA;
      }

   And the volume related functions:

   .. code-block:: freefem
      :linenos:

      // Volume computation
      func real Volume (real[int] &X){
         Vh rho;
         rho[] = X;
         Vh rho3 = rho*rho*rho;
         real res = 1./3.*int2d(Th)(rho3*sin(y));
         return res;
      }

      func real[int] GradVolume (real[int] &X){
         Vh rho;
         rho[] = X;
         varf dVolume(u, v) = int2d(Th)(rho*rho*sin(y)*v);
         real[int] grad = dVolume(0, Vh);
         return grad;
      }

      matrix hessianV;
      func matrix HessianVolume(real[int] &X){
         Vh rho;
         rho[] = X;
         varf d2Volume(w, v) = int2d(Th)(2*rho*sin(y)*v*w);
         hessianV = d2Volume(Vh, Vh);
         return hessianV;
      }

   If we want to use the volume as a constraint function we must wrap it and its derivatives in some **FreeFEM** functions returning the appropriate types.
   It is not done in the above functions in cases where one wants to use it as a fitness function.
   The lagrangian hessian also has to be wrapped since the Volume is not linear with respect to :math:`\rho`, it has some non-null second order derivatives.

   .. code-block:: freefem
      :linenos:

      func real[int] ipVolume (real[int] &X){ real[int] vol = [Volume(X)]; return vol; }
      matrix mdV;
      func matrix ipGradVolume (real[int] &X) { real[int,int] dvol(1,Vh.ndof); dvol(0,:) = GradVolume(X); mdV = dvol; return mdV; }
      matrix HLagrangian;
      func matrix ipHessianLag (real[int] &X, real objfact, real[int] &lambda){
         HLagrangian = objfact*HessianArea(X) + lambda[0]*HessianVolume(X);
         return HLagrangian;
      }

   The :freefem:`ipGradVolume` function could pose some troubles during the optimization process because the gradient vector is transformed in a sparse matrix, so any null coefficient will be discarded.
   Here we create the IPOPT structure manually and use the :freefem:`checkindex` named-parameter to avoid bad indexing during copies.
   This gradient is actually dense, there is no reason for some components to be constantly zero:

   .. code-block:: freefem
      :linenos:

      int[int] gvi(Vh.ndof), gvj=0:Vh.ndof-1;
      gvi = 0;

   These two arrays will be passed to IPOPT with :freefem:`structjacc=[gvi,gvj]`.
   The last remaining things are the bound definitions.
   The simple lower bound must be equal to the components of the P1 projection of :math:`\rho_{object}`.
   And we choose :math:`\alpha\in [0,1]` to set :math:`\mathcal{V}_{\mathrm{max}}` to :math:`(1-\alpha) \mathcal{V}(\rho_{object}) + \alpha\frac{4}{3}\pi \|\rho_{\mathrm{object}}\|_{\infty}^{3}`:

   .. code-block:: freefem
      :linenos:

      func disc1 = sqrt(1./(RR+(E-RR)*cos(y)*cos(y)))*(1+0.1*cos(7*x));
      func disc2 = sqrt(1./(RR+(E-RR)*cos(x)*cos(x)*sin2));

      if(1){
         lb = r0;
         for (int q = 0; q < 5; ++q){
            func f = rr*Gaussian(x, y, 2*q*pi/5., pi/3.);
            func g = rr*Gaussian(x, y, 2*q*pi/5.+pi/5., 2.*pi/3.);
            lb = max(max(lb, f), g);
         }
         lb = max(lb, rr*Gaussian(x, y, 2*pi, pi/3));
      }
      lb = max(lb, max(disc1, disc2));
      real Vobj = Volume(lb[]);
      real Vnvc = 4./3.*pi*pow(lb[].linfty,3);

      if(1)
         Plot3D(lb[], "object_inside", 1);
      real[int] clb = 0., cub = [(1-alpha)*Vobj + alpha*Vnvc];

   Calling IPOPT:

   .. code-block:: freefem
      :linenos:

      int res = IPOPT(Area, GradArea, ipHessianLag, ipVolume, ipGradVolume,
         rc[], ub=ub[], lb=lb[], clb=clb, cub=cub, checkindex=1, maxiter=kkk<nadapt-1 ? 40:150,
         warmstart=kkk, lm=lm, uz=uz[], lz=lz[], tol=0.00001, structjacc=[gvi,gvj]);
      cout << "IPOPT: res =" << res << endl ;

      // Plot
      Plot3D(rc[], "Shape_at_"+kkk, 1);
      Plot3D(GradArea(rc[]), "ShapeGradient", 1);

   Finally, before closing the mesh adaptation loop, we have to perform the said adaptation.
   The mesh is adaptated with respect to the :math:`X=(\rho, 0, 0)` (in spherical coordinates) vector field, not directly with respect to :math:`\rho`, otherwise the true curvature of the 3D-shape would not be well taken into account.

   .. code-block:: freefem
      :linenos:

      if (kkk < nadapt-1){
         Th = adaptmesh(Th, rc*cos(x)*sin(y), rc*sin(x)*sin(y), rc*cos(y),
            nbvx=50000, periodic=[[2, y], [4, y]]);
         plot(Th, wait=true);
         startshape = rc;
         uz = uz;
         lz = lz;
      }

   Here are some pictures of the resulting surfaces obtained for decreasing values of :math:`\alpha` (and a slightly more complicated object than two orthogonal discs).
   We return to the enclosed object when :math:`\alpha=0`:

   .. figure:: images/minsurf3D.jpg
      :width: 100%


The nlOpt optimizers
--------------------

The :freefem:`ff-NLopt` package provides a **FreeFEM** interface to the free/open-source library for nonlinear optimization, easing the use of several different free optimization (constrained or not) routines available online along with the PDE solver.
All the algorithms are well documented in `NLopt documentation <https://nlopt.readthedocs.io/en/latest/>`__, therefore no exhaustive information concerning their mathematical specificities will be found here and we will focus on the way they are used in a **FreeFEM** script.
If needing detailed information about these algorithms, visit the website where a description of each of them is given, as well as many bibliographical links.

Most of the gradient based algorithms of NLopt uses a full matrix approximation of the Hessian, so if you’re planning to solve a large scale problem, use the IPOPT optimizer which definitely surpass them.

All the NLopt features are identified that way:

.. code-block:: freefem
   :linenos:

   load "ff-NLopt"
   //define J, u, and maybe grad(J), some constraints etc...
   real min = nloptXXXXXX(J, u, //Unavoidable part
      grad=<name of grad(J)>, //if needed
      lb= //Lower bounds array
      ub= //Upper bounds array
      ... //Some optional arguments:
      //Constraints functions names,
      //Stopping criteria,
      //Algorithm specific parameters,
      //Etc...
   );

:freefem:`XXXXXX` refers to the algorithm tag (not necessarily 6 characters long).
:freefem:`u` is the starting position (a :freefem:`real[int]` type array) which will be overwritten by the algorithm, the value at the end being the found :math:`argmin`.
And as usual, :freefem:`J` is a function taking a :freefem:`real[int]` type array as argument and returning a :freefem:`real`.
:freefem:`grad`, :freefem:`lb` and :freefem:`ub` are "half-optional" arguments, in the sense that they are obligatory for some routines but not all.

The possible optionally named parameters are the following, note that they are not used by all algorithms (some do not support constraints, or a type of constraints, some are gradient-based and others are derivative free, etc…).
One can refer to the table after the parameters description to check which are the named parameters supported by a specific algorithm.
Using an unsupported parameter will not stop the compiler work, seldom breaks runtime, and will just be ignored.
When it is obvious you are missing a routine, you will get a warning message at runtime (for example if you pass a gradient to a derivative free algorithm, or set the population of a non-genetic one, etc…).
In the following description, :math:`n` stands for the dimension of the search space.

**Half-optional parameters :**

-  :freefem:`grad=` The name of the function which computes the gradient of the cost function (prototype should be :freefem:`real[int]` :math:`\rightarrow` :freefem:`real[int]`, both argument and result should have the size :math:`n`).
   This is needed as soon as a gradient-based method is involved, which is ignored if defined in a derivative free context.
-  :freefem:`lb`/:freefem:`ub` = Lower and upper bounds arrays ( :freefem:`real[int]` type) of size :math:`n`.
   Used to define the bounds within which the search variable is allowed to move.
   Needed for some algorithms, optional, or unsupported for others.
-  :freefem:`subOpt` : Only enabled for the Augmented Lagrangian and MLSL methods who need a sub-optimizer in order to work.
   Just pass the tag of the desired local algorithm with a :freefem:`string`.

**Constraints related parameters (optional - unused if not specified):**

-  :freefem:`IConst`/:freefem:`EConst` : Allows to pass the name of a function implementing some inequality (resp. equality) constraints on the search space.
   The function type must be :freefem:`real[int]` :math:`\rightarrow` :freefem:`real[int]` where the size of the returned array is equal to the number of constraints (of the same type - it means that all of the constraints are computed in one vectorial function).
   In order to mix inequality and equality constraints in a same minimization attempt, two vectorial functions have to be defined and passed.
   See example :eq:`varineqex` for more details about how these constraints have to be implemented.
-  :freefem:`gradIConst`/:freefem:`gradEConst` : Use to provide the inequality (resp. equality) constraints gradient.
   These are :freefem:`real[int]` :math:`\rightarrow` :freefem:`real[int,int]` type functions.
   Assuming we have defined a constraint function (either inequality or equality) with :math:`p` constraints, the size of the matrix returned by its associated gradient must be :math:`p\times n` (the :math:`i`-th line of the matrix is the gradient of the :math:`i`-th constraint).
   It is needed in a gradient-based context as soon as an inequality or equality constraint function is passed to the optimizer and ignored in all other cases.
-  :freefem:`tolIConst`/:freefem:`tolEConst` : Tolerance values for each constraint.
   This is an array of size equal to the number of inequality (resp. equality) constraints.
   Default value is set to :math:`10^{-12}` for each constraint of any type.

**Stopping criteria :**

-  :freefem:`stopFuncValue` : Makes the algorithm end when the objective function reaches this :freefem:`real` value.
-  :freefem:`stopRelXTol` : Stops the algorithm when the relative moves in each direction of the search space is smaller than this :freefem:`real` value.
-  :freefem:`stopAbsXTol` : Stops the algorithm when the moves in each direction of the search space is smaller than the corresponding value in this :freefem:`real[int]` array.
-  :freefem:`stopRelFTol` : Stops the algorithm when the relative variation of the objective function is smaller than this :freefem:`real` value.
-  :freefem:`stopAbsFTol` : Stops the algorithm when the variation of the objective function is smaller than this :freefem:`real` value.
-  :freefem:`stopMaxFEval` : Stops the algorithm when the number of fitness evaluations reaches this :freefem:`integer` value.
-  :freefem:`stopTime` : Stops the algorithm when the optimization time in seconds exceeds this :freefem:`real` value.
   This is not a strict maximum: the time may exceed it slightly, depending upon the algorithm and on how slow your function evaluation is.

   Note that when an AUGLAG or MLSL method is used, the meta-algorithm and the sub-algorithm may have different termination criteria.
   Thus, for algorithms of this kind, the following named parameters has been defined (just adding the SO prefix - for Sub-Optimizer) to set the ending condition of the sub-algorithm (the meta one uses the ones above): :freefem:`SOStopFuncValue`, :freefem:`SOStopRelXTol`, and so on… If these are not used, the sub-optimizer will use those of the master routine.

**Other named parameters :**

-  :freefem:`popSize` : :freefem:`integer` used to change the size of the sample for stochastic search methods.
   Default value is a peculiar heuristic to the chosen algorithm.
-  :freefem:`SOPopSize` : Same as above, but when the stochastic search is passed to a meta-algorithm.
-  :freefem:`nGradStored` : The number (:freefem:`integer` type) of gradients to "remember" from previous optimization steps: increasing this increases the memory requirements but may speed convergence.
   It is set to a heuristic value by default.
   If used with AUGLAG or MLSL, it will only affect the given subsidiary algorithm.

The following table sums up the main characteristics of each algorithm, providing the more important information about which features are supported by which algorithm and what are the unavoidable arguments they need.
More details can be found in `NLopt documentation <https://nlopt.readthedocs.io/en/latest/>`__.

.. figure:: images/nlopttab.png
   :height: 22cm

.. tip:: Variational inequality

   Let :math:`\Omega` be a domain of :math:`\mathbb{R}^{2}`, :math:`f_{1}, f_{2}\in L^{2}(\Omega)` and :math:`g_{1}, g_{2} \in L^{2}(\partial\Omega)` four given functions with :math:`g_{1}\leq g_{2}` almost everywhere.

   We define the space:

   .. math::
      V = \left\lbrace (v_{1},v_{2})\in H^{1}(\Omega)^{2} ; v_{1}\vert_{\partial\Omega}=g_{1}, v_{2}\vert_{\partial\Omega}=g_{2}, v_{1}\leq v_{2}\ \mathrm{a.e.}\ \right\rbrace

   as well as the function :math:`J:H^{1}(\Omega)^{2}\longrightarrow \mathbb{R}`:

   .. math::
      J(v_{1},v_{2}) = \displaystyle{\frac{1}{2}\int_{\Omega}\vert\nabla v_{1}\vert^{2} - \int_{\Omega} f_{1}v_{1} + \frac{1}{2}\int_{\Omega}\vert\nabla v_{2}\vert^{2} - \int_{\Omega} f_{2}v_{2}}
      :label: varineqex

   The problem consists in finding (numerically) two functions :math:`(u_{1},u_{2}) = \underset{(v_{1},v_{2})\in V}{\operatorname{argmin}} J(v_{1},v_{2})`.

   This can be interpreted as finding :math:`u_{1}, u_{2}` as close as possible (in a certain sense) to the solutions of the Laplace equation with respectively :math:`f_{1}, f_{2}` second members and :math:`g_{1}, g_{2}` Dirichlet boundary conditions with the :math:`u_{1}\leq u_{2}` almost everywhere constraint.

   Here is the corresponding script to treat this variational inequality problem with one of the NLOpt algorithms.

   .. code-block:: freefem
      :linenos:

      //A brief script to demonstrate how to use the freefemm interfaced nlopt routines
      //The problem consist in solving a simple variational inequality using one of the
      //optimization algorithm of nlopt. We restart the algorithlm a few times after
      //performing some mesh adaptation to get a more precise output

      load "ff-NLopt"

      // Parameters
      int kas = 3; //choose of the algorithm
      int NN = 10;
      func f1 = 1.;
      func f2 = -1.;
      func g1 = 0.;
      func g2 = 0.1;
      int iter = 0;
      int nadapt = 2;
      real starttol = 1e-6;
      real bctol = 6.e-12;

      // Mesh
      mesh Th = square(NN, NN);

      // Fespace
      fespace Vh(Th, P1);
      Vh oldu1, oldu2;

      // Adaptation loop
      for (int al = 0; al < nadapt; ++al){
         varf BVF (v, w) = int2d(Th)(0.5*dx(v)*dx(w) + 0.5*dy(v)*dy(w));
         varf LVF1 (v, w) = int2d(Th)(f1*w);
         varf LVF2 (v, w) = int2d(Th)(f2*w);
         matrix A = BVF(Vh, Vh);
         real[int] b1 = LVF1(0, Vh), b2 = LVF2(0, Vh);

         varf Vbord (v, w) = on(1, 2, 3, 4, v=1);

         Vh In, Bord;
         Bord[] = Vbord(0, Vh, tgv=1);
         In[] = Bord[] ? 0:1;
         Vh gh1 = Bord*g1, gh2 = Bord*g2;

         func real J (real[int] &X){
            Vh u1, u2;
            u1[] = X(0:Vh.ndof-1);
            u2[] = X(Vh.ndof:2*Vh.ndof-1);
            iter++;
            real[int] Au1 = A*u1[], Au2 = A*u2[];
            Au1 -= b1;
            Au2 -= b2;
            real val = u1[]'*Au1 + u2[]'*Au2;
            if (iter%10 == 9)
               plot(u1, u2, nbiso=30, fill=1, dim=3, cmm="adapt level "+al+" - iteration "+iter+" - J = "+val, value=1);
            return val;
         }

         varf dBFV (v, w) = int2d(Th)(dx(v)*dx(w)+dy(v)*dy(w));
         matrix dA = dBFV(Vh, Vh);
         func real[int] dJ (real[int] &X){
            Vh u1, u2;
            u1[] = X(0:Vh.ndof-1);
            u2[] = X(Vh.ndof:2*Vh.ndof-1);

            real[int] grad1 = dA*u1[], grad2 = dA*u2[];
            grad1 -= b1;
            grad2 -= b2;
            real[int] Grad(X.n);
            Grad(0:Vh.ndof-1) = grad1;
            Grad(Vh.ndof:2*Vh.ndof-1) = grad2;
            return Grad;
         }

         func real[int] IneqC (real[int] &X){
            real[int] constraints(Vh.ndof);
            for (int i = 0; i < Vh.ndof; ++i) constraints[i] = X[i] - X[i+Vh.ndof];
            return constraints;
         }

         func real[int,int] dIneqC (real[int] &X){
            real[int, int] dconst(Vh.ndof, 2*Vh.ndof);
            dconst = 0;
            for(int i = 0; i < Vh.ndof; ++i){
               dconst(i, i) = 1.;
               dconst(i, i+Vh.ndof) = -1.;
            }
            return dconst;
         }

         real[int] BordIndex(Th.nbe); //Indexes of border d.f.
         {
            int k = 0;
            for (int i = 0; i < Bord.n; ++i) if (Bord[][i]){ BordIndex[k] = i; ++k; }
         }

         func real[int] BC (real[int] &X){
            real[int] bc(2*Th.nbe);
            for (int i = 0; i < Th.nbe; ++i){
               int I = BordIndex[i];
               bc[i] = X[I] - gh1[][I];
               bc[i+Th.nbe] = X[I+Th.nv] - gh2[][I];
            }
            return bc;
         }

         func real[int, int] dBC(real[int] &X){
            real[int, int] dbc(2*Th.nbe,2*Th.nv);
            dbc = 0.;
            for (int i = 0; i < Th.nbe; ++i){
               int I = BordIndex[i];
               dbc(i, I) = 1.;
               dbc(i+Th.nbe, I+Th.nv) = 1.;
            }
            return dbc;
         }

         real[int] start(2*Vh.ndof), up(2*Vh.ndof), lo(2*Vh.ndof);

         if (al == 0){
            start(0:Vh.ndof-1) = 0.;
            start(Vh.ndof:2*Vh.ndof-1) = 0.01;
         }
         else{
            start(0:Vh.ndof-1) = oldu1[];
            start(Vh.ndof:2*Vh.ndof-1) = oldu2[];
         }

         up = 1000000;
         lo = -1000000;
         for (int i = 0; i < Vh.ndof; ++i){
            if (Bord[][i]){
               up[i] = gh1[][i] + bctol;
               lo[i] = gh1[][i] - bctol;
               up[i+Vh.ndof] = gh2[][i] + bctol;
               lo[i+Vh.ndof] = gh2[][i] - bctol;
            }
         }

         real mini = 1e100;
         if (kas == 1)
            mini = nloptAUGLAG(J, start, grad=dJ, lb=lo,
               ub=up, IConst=IneqC, gradIConst=dIneqC,
               subOpt="LBFGS", stopMaxFEval=10000, stopAbsFTol=starttol);
         else if (kas == 2)
            mini = nloptMMA(J, start, grad=dJ, lb=lo, ub=up, stopMaxFEval=10000, stopAbsFTol=starttol);
         else if (kas == 3)
            mini = nloptAUGLAG(J, start, grad=dJ, IConst=IneqC,
               gradIConst=dIneqC, EConst=BC, gradEConst=dBC,
               subOpt="LBFGS", stopMaxFEval=200, stopRelXTol=1e-2);
         else if (kas == 4)
            mini = nloptSLSQP(J, start, grad=dJ, IConst=IneqC,
               gradIConst=dIneqC, EConst=BC, gradEConst=dBC,
               stopMaxFEval=10000, stopAbsFTol=starttol);
         Vh best1, best2;
         best1[] = start(0:Vh.ndof-1);
         best2[] = start(Vh.ndof:2*Vh.ndof-1);

         Th = adaptmesh(Th, best1, best2);
         oldu1 = best1;
         oldu2 = best2;
      }

Optimization with MPI
---------------------

The only quick way to use the previously presented algorithms on a parallel architecture lies in parallelizing the used cost function (which is in most real life cases, the expensive part of the algorithm).
Somehow, we provide a parallel version of the CMA-ES algorithm.
The parallelization principle is the trivial one of evolving/genetic algorithms: at each iteration the cost function has to be evaluated :math:`N` times without any dependence at all, these :math:`N` calculus are then equally distributed to each process.
Calling the MPI version of CMA-ES is nearly the same as calling its sequential version (a complete example of use can be found in the :ref:`CMAES MPI variational inequality example <exampleCMAESMPIVariationalInequality>`):

.. code-block:: freefem
   :linenos:

   load "mpi-cmaes"
   ... // Define J, u and all here
   real min = cmaesMPI(J, u, stopTolFun=1e-6, stopMaxIter=3000);
   cout << "minimum value is " << min << " for u = " << u << endl;

If the population size is not changed using the :freefem:`popsize` parameter, it will use the heuristic value slightly changed to be equal to the closest greatest multiple of the size of the communicator used by the optimizer.
The **FreeFEM** :freefem:`mpicommworld` is used by default.
The user can specify his own MPI communicator with the named parameter :freefem:`comm=`, see the MPI section of this manual for more information about communicators in **FreeFEM**.
