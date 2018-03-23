
## Variational inequality

We present, a classical example of variational inequality.

Let us denote  $\mathcal{C} = \{ u\in H^1_0(\Omega), u \le g \}$

The problem is :
$$
 u = arg \min_{u\in \mathcal{C}}  J(u) = \frac{1}{2} \int_\Omega \nabla u . \nabla u - \int_\Omega f u
$$
where $f$ and $g$ are given function.

The solution is a projection on the convex $\mathcal{C}$ of $f^\star$ for the scalar product $((v,w)) = \int_\Omega \nabla v . \nabla w$ of $H^1_0(\Omega)$ where ${f^\star} $ is solution of $ ((f^\star, v )) = \int_\Omega f v, \forall v \in H^1_0(\Omega)$.

The projection on a convex satisfy clearly $\forall v \in \mathcal{C}, \quad (( u -v ,  u - \tilde{f}  )) \leq 0$, and after expanding, we get the classical inequality

$$\forall v \in \mathcal{C}, \quad \int_\Omega \nabla(u -v) \nabla u  \leq  \int_\Omega (u-v) f .$$

We can also rewrite the problem as a saddle point problem

Find $\lambda, u$ such that:
$$
  \max_{\lambda\in L^2(\Omega), \lambda\geq 0}  \min_{u\in H^1_0(\Omega)}  \mathcal{L}(u,\lambda) = \frac{1}{2} \int_\Omega \nabla u . \nabla u - \int_\Omega f u  + \int_{\Omega} \lambda (u-g)^+
$$
where $((u-g)^+ = max(0,u-g) $

This saddle point problem is equivalent to find $ u, \lambda $ such that:

\begin{equation}
 \left\{
\begin{array}{cc}
\displaystyle \int_\Omega \nabla u . \nabla v + \lambda v^+ \,d\omega= \int_\Omega f u  , &\forall v \in H^1_0(\Omega) \cr
\displaystyle \int_\Omega   \mu (u-g)^+ = 0  , & \forall \mu \in L^2(\Omega) , \mu \geq 0, \lambda \geq 0,
 \end{array}\right.
\end{equation}

An algorithm to solve the previous problem is:

1. k=0, and choose, $\lambda_0$ belong $H^{-1}(\Omega)$

2. Loop on $k = 0, .....$

	* set $\mathcal{I}_{k} = \{ x \in \Omega / \lambda_{k} + c * ( u_{k+1} - g)  \leq 0 \}$
	* $V_{g,k+1} = \{ v\in H^1_0(\Omega) / v = g$ on ${I}_{k} \}$,
	* $V_{0,k+1} = \{ v\in H^1_0(\Omega) / v = 0$ on ${I}_{k} \}$,
	* Find $u_{k+1} \in V_{g,k+1}$ and  $\lambda_{k+1} \in H^{-1}(\Omega)$ such that

		$$
		\left\{\begin{array}{cc}
		\displaystyle  \int_\Omega \nabla u_{k+1}. \nabla v_{k+1}   \,d\omega = \int_\Omega f v_{k+1}  , &\forall v_{k+1} \in V_{0,k+1} \cr
		\displaystyle  <\lambda_{k+1},v>  =  \int_\Omega \nabla u_{k+1}. \nabla v  -  f v \,d\omega &
		 \end{array}\right.
		$$

		where $<,>$ is the duality bracket between $H^{1}_0(\Omega)$ and  $H^{-1}(\Omega)$, and $c$ is a penalty constant (large enough).

You can find all the mathematics about this algorithm in \cite{ItoKunisch} 38 $\codered$.

Now how to do that in FreeFem++. The full example is:

 __Example 9.25__ VI.edp

```freefem
mesh Th=square(20,20);
real eps=1e-5;
fespace Vh(Th,P1); // P1 FE space
int n = Vh.ndof; // number of Degree of freedom
Vh uh,uhp; // solution and previous one
Vh Ik; // to def the set where the containt is reached.
real[int] rhs(n); // to store the right and side of the equation
real c=1000; // the penalty parameter of the algoritm
func f=1; // right hand side function
func fd=0; // Dirichlet boundary condition function
Vh g=0.05; // the discret function g

real[int] Aii(n),Aiin(n); // to store the diagonal of the matrix 2 version

real tgv = 1e30; // a huge value for exact penalization
// of boundary condition
// the variatonal form of the problem:
varf a(uh,vh) = // definition of the problem
    int2d(Th)( dx(uh)*dx(vh) + dy(uh)*dy(vh) ) // bilinear form
  - int2d(Th)( f*vh ) // linear form
  + on(1,2,3,4,uh=fd) ; // boundary condition form


// two version of the matrix of the problem
matrix A=a(Vh,Vh,tgv=tgv,solver=CG); // one changing
matrix AA=a(Vh,Vh,solver:GC); // one for computing residual

 // the mass Matrix construction:
varf vM(uh,vh) = int2d(Th)(uh*vh);
matrix M=vM(Vh,Vh); // to do a fast computing of $L^2$ norm : sqrt( u'*(w=M*u))

Aii=A.diag; // get the diagonal of the matrix (appear in version 1.46-1)

rhs = a(0,Vh,tgv=tgv);
Ik =0;
uhp=-tgv; // previous value is
Vh lambda=0;
for(int iter=0;iter<100;++iter)
{
  real[int] b(n) ; b=rhs; // get a copy of the Right hand side
  real[int] Ak(n); // the complementary of Ik ( !Ik = (Ik-1))
 // Today the operator Ik- 1. is not implement so we do:
  Ak= 1.; Ak  -= Ik[]; // build Ak  = ! Ik
 // adding new locking condition on b and on the diagonal if (Ik ==1 )
  b = Ik[] .* g[];      b *= tgv;     b  -=  Ak .* rhs;
  Aiin = Ik[] *  tgv;      Aiin  +=  Ak  .* Aii; //set Aii= tgv  $ i \in Ik $
  A.diag = Aiin; // set the matrix diagonal  (appear in version 1.46-1)
  set(A,solver=CG); // important to change preconditioning for solving
  uh[] = A^-1* b; // solve the problem with more locking condition
  lambda[] = AA * uh[]; // compute the residual ( fast with matrix)
  lambda[] += rhs; // remark rhs = $-\int f v $

  Ik = ( lambda + c*( g- uh)) < 0.; // the new of locking value

   plot(Ik, wait=1,cmm=" lock set ",value=1,ps="VI-lock.eps",fill=1 );
   plot(uh,wait=1,cmm="uh",ps="VI-uh.eps");
 // trick to compute  $L^2$ norm of the variation (fast method)
      real[int] diff(n),Mdiff(n);
      diff= uh[]-uhp[];
      Mdiff = M*diff;
      real err = sqrt(Mdiff'*diff);
  cout << "  || u_{k=1} - u_{k} ||_2 " << err << endl;
  if(err< eps) break; // stop test
  uhp[]=uh[] ; // set the previous solution
}
savemesh(Th,"mm",[x,y,uh*10]); // for medit plotting
```

Remark, as you can see on this example, some vector , or matrix operator are not implemented
so a way is to skip the expression and we use operator `:::freefem +=`,  `:::freefem -=` to merge
the result.
