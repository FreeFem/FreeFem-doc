.. role:: freefem(code)
  :language: freefem

Variational Inequality
======================

We present, a classical example of variational inequality.

Let us denote :math:`\mathcal{C} = \{ u\in H^1_0(\Omega), u \le g \}`

The problem is:

.. math::
   u = arg \min_{u\in \mathcal{C}} J(u) = \frac{1}{2} \int_\Omega \nabla u . \nabla u - \int_\Omega f u

where :math:`f` and :math:`g` are given function.

The solution is a projection on the convex :math:`\mathcal{C}` of :math:`f^\star` for the scalar product :math:`((v,w)) = \int_\Omega \nabla v . \nabla w` of :math:`H^1_0(\Omega)` where :math:`f^\star` is solution of:

.. math::
   (f^\star, v ) = \int_{\Omega}{f v}, \forall v \in H^1_0(`\Omega)

The projection on a convex satisfy clearly :math:`\forall v \in \mathcal{C}, \quad (( u -v , u - \tilde{f} )) \leq 0`, and after expanding, we get the classical inequality:

.. math::
   \forall v \in \mathcal{C}, \quad \int_\Omega \nabla(u -v) \nabla u \leq \int_\Omega (u-v) f

We can also rewrite the problem as a saddle point problem:

Find :math:`\lambda, u` such that:

.. math::
   \max_{\lambda\in L^2(\Omega), \lambda\geq 0} \min_{u\in H^1_0(\Omega)} \mathcal{L}(u,\lambda) = \frac{1}{2} \int_\Omega \nabla u . \nabla u - \int_\Omega f u + \int_{\Omega} \lambda (u-g)^+

where :math:`((u-g)^+ = max(0,u-g)`.

This saddle point problem is equivalent to find :math:`u, \lambda` such that:

.. math::
    \left\{
    \begin{array}{cc}
        \displaystyle \int_\Omega \nabla u . \nabla v + \lambda v^+ \,d\omega= \int_\Omega f u , &\forall v \in H^1_0(\Omega) \cr
        \displaystyle \int_\Omega \mu (u-g)^+ = 0 , & \forall \mu \in L^2(\Omega) , \mu \geq 0, \lambda \geq 0,
    \end{array}\right.

An algorithm to solve the previous problem is:

1. k=0, and choose :math:`\lambda_0` belong :math:`H^{-1}(\Omega)`

2. Loop on :math:`k = 0, .....`

   -  set :math:`\mathcal{I}_{k} = \{ x \in \Omega / \lambda_{k} + c * ( u_{k+1} - g) \leq 0 \}`
   -  :math:`V_{g,k+1} = \{ v\in H^1_0(\Omega) / v = g` on :math:`{I}_{k} \}`,
   -  :math:`V_{0,k+1} = \{ v\in H^1_0(\Omega) / v = 0` on :math:`{I}_{k} \}`,
   -  Find :math:`u_{k+1} \in V_{g,k+1}` and :math:`\lambda_{k+1} \in H^{-1}(\Omega)` such that

      .. math::
         \left\{\begin{array}{cc}
            \displaystyle \int_\Omega \nabla u_{k+1}. \nabla v_{k+1} \,d\omega = \int_\Omega f v_{k+1} , &\forall v_{k+1} \in V_{0,k+1} \cr
            \displaystyle <\lambda_{k+1},v> = \int_\Omega \nabla u_{k+1}. \nabla v - f v \,d\omega &
         \end{array}\right.

      where :math:`<,>` is the duality bracket between :math:`H^{1}_0(\Omega)` and :math:`H^{-1}(\Omega)`, and :math:`c` is a penalty constant (large enough).

You can find all the mathematics about this algorithm in [ITO2003]_ [HINTERMULLER2002]_.

Now how to do that in **FreeFEM**? The full example is:

.. tip:: Variational inequality

   .. code-block:: freefem
      :linenos:

      load "medit"

      // Parameters
      real eps = 1e-5;
      real c = 1000; //penalty parameter of the algoritm
      real tgv = 1e30; //a huge value for exact penalization
      func f = 1; //right hand side function
      func fd = 0; //Dirichlet boundary condition function

      // Mesh
      mesh Th = square(20, 20);

      // Fespace
      fespace Vh(Th, P1);
      int n = Vh.ndof; //number of degree of freedom
      Vh uh, uhp; //u^n+1 and u^n
      Vh Ik; //to define the set where the containt is reached.
      Vh g = 0.05; //discret function g
      Vh lambda = 0;

      // Problem
      varf a (uh, vh)
          = int2d(Th)(
                dx(uh)*dx(vh)
              + dy(uh)*dy(vh)
          )
          - int2d(Th)(
                f*vh
          )
          + on(1, 2, 3, 4, uh=fd)
          ;

      //the mass Matrix construction
      varf vM (uh, vh) = int2d(Th)(uh*vh);

      //two versions of the matrix of the problem
      matrix A = a(Vh, Vh, tgv=tgv, solver=CG); //one changing
      matrix AA = a(Vh, Vh, solver=CG); //one for computing residual

      matrix M = vM(Vh, Vh); //to do a fast computing of L^2 norm : sqrt(u'*(w=M*u))

      real[int] Aiin(n);
      real[int] Aii = A.diag; //get the diagonal of the matrix
      real[int] rhs = a(0, Vh, tgv=tgv);

      // Initialization
      Ik = 0;
      uhp = -tgv;

      // Loop
      for(int iter = 0; iter < 100; ++iter){
          // Update
          real[int] b = rhs; //get a copy of the Right hand side
          real[int] Ak(n); //the complementary of Ik ( !Ik = (Ik-1))
          Ak = 1.; Ak -= Ik[];
          //adding new locking condition on b and on the diagonal if (Ik ==1 )
          b = Ik[] .* g[]; b *= tgv; b -= Ak .* rhs;
          Aiin = Ik[] * tgv; Aiin += Ak .* Aii; //set Aii= tgv i in Ik
          A.diag = Aiin; //set the matrix diagonal
          set(A, solver=CG); //important to change preconditioning for solving

          // Solve
          uh[] = A^-1* b; //solve the problem with more locking condition

          // Residual
          lambda[] = AA * uh[]; //compute the residual (fast with matrix)
          lambda[] += rhs; //remark rhs = -\int f v

          Ik = (lambda + c*( g- uh)) < 0.; //the new locking value

          // Plot
          plot(Ik, wait=true, cmm=" lock set ", value=true, fill=true);
          plot(uh, wait=true, cmm="uh");

          // Error
          //trick to compute L^2 norm of the variation (fast method)
          real[int] diff(n), Mdiff(n);
          diff = uh[] - uhp[];
          Mdiff = M*diff;
          real err = sqrt(Mdiff'*diff);
          cout << "|| u_{k=1} - u_{k} ||_2 = " << err << endl;

          // Stop test
          if(err < eps) break;

          // Update
          uhp[] = uh[];
      }

      // Plot
      medit("uh", Th, uh);

   .. note:: As you can see on this example, some vector, or matrix operator are not implemented so a way is to skip the expression and we use operator :freefem:`+=`,  :freefem:`-=` to merge the result.
