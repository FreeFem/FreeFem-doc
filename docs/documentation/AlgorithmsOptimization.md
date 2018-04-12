$$\newcommand{\MyInt}[2]{\displaystyle{\int_{#1} #2 } }$$

## Conjugate Gradient/GMRES

Suppose we want to solve the Euler problem (here $x$ has nothing to do with the reserved variable for the first coordinate in FreeFem++): find $ x\in \R^n $ such that

\begin{equation}
\label{eqn:dJ=0}
\nabla J(x) = \left(\frac{\p J}{\p x_i} (\vec{x})\right) = 0
\end{equation}

where $ J$ is a functional (to minimize for example) from $ \R^n$ to $ \R$.

If the function is convex we can use the conjugate gradient to solve the problem, and we just need the function (named `:::freefem dJ` for example) which compute $\nabla J$, so the parameters are the name of that function with prototype `:::freefem func real[int] dJ(real[int] & xx);`

which compute $\nabla J$, and a vector `:::freefem x` of type (of course the number 20 can be changed) `:::freefem real[int] x(20);` to initialize the process and get the result.

Given an initial value $\vec{x}^{(0)}$, a maximum number $i_{\max}$ of iterations, and an error tolerance $0<\epsilon<1$:

Put $\vec{x}=\vec{x}^{(0)}$ and write

```freefem
NLCG(dJ, x, precon=M, nbiter=imax, eps=epsilon, stop=stopfunc);
```

will give the solution of $\vec{x}$ of $\nabla J(\vec{x})=0$. We can omit parameters `:::freefem precon, nbiter, eps, stop`. Here $M$ is the preconditioner whose default is the identity matrix.

The stopping test is

\[
\| \nabla J(\vec{x})\|_P\le \epsilon\| \nabla J(\vec{x}^{(0)})\|_P
\]

Writing the minus value in `:::freefem eps=`, i.e.,

```freefem
NLCG(dJ, x, precon=M, nbiter=imax, eps=-epsilon);
```

we can use the stopping test :

\[
\| \nabla J(\vec{x})\|_P^2\le \epsilon
\]

The parameters of these three functions are:

* `:::freefem nbiter=` set the number of iteration (by default 100)

* `:::freefem precon=` set the preconditioner function (`:::freefem P` for example) by default it is the identity, remark the prototype is `:::freefem func real[int] P(real[int] &x)`.

* `:::freefem eps=` set the value of the stop test $\varepsilon$ ($=10^{-6}$ by default) if positive then relative test $||\nabla J(x)||_P\leq \varepsilon||\nabla J(x_0)||_P$, otherwise the absolute test is $||\nabla J(x)||_P^2\leq |\varepsilon|$.

* `:::freefem veps=` set and return the value of the stop test, if positive then relative test $||\nabla J(x)||_P\leq \varepsilon||\nabla J(x_0)||_P$, otherwise the absolute test is $||\nabla J(x)||_P^2\leq |\varepsilon|$. The return value is minus the real stop test (remark: it is useful in loop).

* `:::freefem stop=` `:::freefem stopfunc` $\codered$ add your test function to stop before (after version 3.31) $\codered$. The prototype for the function `:::freefem stopfunc` is

	```freefem
	 func bool stopfunc(int iter, real[int] u, real[int] g)
	```

	where `:::freefem u` is the current solution, and `:::freefem g` the current gradient not preconditioned.

__Example algo.edp__

For a given function $b$, let us find the minimizer $u$ of the functional

\begin{eqnarray*}
 J(u) &=& \frac{1}{2}\int_{\Omega} f(|\nabla u|^2) - \int_{\Omega} u b \\
 f(x) &=& ax + x-\ln(1+x), \quad f'(x) = a+\frac{x}{1+x}, \quad f''(x) = \frac{1}{(1+x)^2}
\end{eqnarray*}

under the boundary condition $u=0$ on $\p\Omega$.

```freefem
func real J(real[int] & u)
{
Vh w;w[]=u; // copy array u in the finite element function w
real r=int2d(Th)(0.5*f( dx(w)*dx(w) + dy(w)*dy(w) ) - b*w) ;
cout << "J(u) =" << r << " " << u.min << " " << u.max << endl;
return r;
}
// -----------------------

Vh u=0; // the current value of the solution
Ph alpha; // of store $df(|\nabla u|^2)$
int iter=0;
alpha=df( dx(u)*dx(u) + dy(u)*dy(u) ); // optimization

func real[int] dJ(real[int] & u)
{
int verb=verbosity; verbosity=0;
Vh w;w[]=u; // copy array u in the finite element function w
alpha=df( dx(w)*dx(w) + dy(w)*dy(w) ); // optimization
varf au(uh,vh) = int2d(Th)( alpha*( dx(w)*dx(vh) + dy(w)*dy(vh) ) - b*vh)
+ on(1,2,3,4,uh=0);
u= au(0,Vh);
verbosity=verb;
return u; // warning no return of local array
}
```

We want to construct also a preconditioner $C$ with solving the problem:

find $u_h \in V_{0h}$ such that :
\[
\forall v_h \in V_{0h}, \quad \int_\Omega \alpha \nabla u_h . \nabla v_h = \int_\Omega b v_h
\]
where $ \alpha=f'(|\nabla u|^2)$.

```freefem
varf alap(uh,vh)= int2d(Th)( alpha *( dx(uh)*dx(vh) + dy(uh)*dy(vh) ))
 + on(1,2,3,4,uh=0);

varf amass(uh)= int2d(Th)( uh*vh) + on(1,2,3,4,uh=0);

matrix amass = alap(Vh,Vh,solver=CG);
matrix Alap= alap(Vh,Vh,solver=Cholesky,factorize=1);

// The preconditionner function
func real[int] C(real[int] & u)
{
 real[int] w = u;
 u = Alap^-1*w;
 return u; // No return of local array variable
}
```

To solve the problem, we make 10 iteration of the conjugate gradient,
recompute the preconditioner and restart the conjugate gradient:

```freefem
 verbosity=5;
 int conv=0;
 real eps=1e-6;
 for(int i=0;i<20;i++)
 {
 conv=NLCG(dJ,u[],nbiter=10,precon=C,veps=eps); // \index{veps=}\index{NLCG}
 if (conv) break; // if converge break loop
 alpha=df( dx(u)*dx(u) + dy(u)*dy(u) ); // recompute alpha optimization
 Alap = alap(Vh,Vh,solver=Cholesky,factorize=1);
 cout << " restart with new preconditionner " << conv
 << " eps =" << eps << endl;
 }

 plot (u,wait=1,cmm="solution with NLCG");
```

$/codered$ END OF EXAMPLE (FRANCK)

For a given symmetric positive matrix $A$, consider the quadratic form

\[
J(\vec{x})=\frac{1}{2}\vec{x}^TA\vec{x}-\vec{b}^T\vec{x}
\]

then $J(\vec{x})$ is minimized by the solution $\vec{x}$ of $A\vec{x}=\vec{b}$. In this case, we can use the function `:::freefem LinearCG` $\codered$ LINEAR OR AFFINE ?? (FRANCK)

```freefem
AffineCG(A, x, precon=M, nbiter=imax, eps=±epsilon, stop=stp);
```

If $A$ is not symmetric, we can use GMRES(Generalized Minimum Residual) algorithm by

```freefem
AffineGMRES(A, x, precon=M, nbiter=imax, eps=±epsilon);
```

Also, we can use the non-linear version of GMRES algorithm
(the functional $J$ is just convex)

```freefem
AffineGMRES(dJ, x, precon=M, nbiter=imax, eps=±epsilon);
```

For detail of these algorithms, refer to \cite{Lucquin} $\codered$[Chapter IV, 1.3].

## Algorithms for Unconstrained Optimization

Two algorithms of COOOL a package \cite{coool} $\codered$ are interfaced with the Newton Raphson method (call `:::freefem Newton`) and the `:::freefem BFGS` method. These two ones are directly available in FreeFem (no dynamical link to load). Be careful with these algorithms, because their implementation uses full matrices. We also provide several optimization algorithms from the NLopt library \cite{nlopt} $\codered$ as well as an interface for Hansen's implementation of CMAES (a MPI version of this one is also available). These last algorithms can be found as dynamical links in the `:::freefem example++-load` folder as the `:::freefem ff-NLopt` and `:::freefem CMA_ES` files (`:::freefem CMA_ES_MPI` from the `:::freefem example++-mpi` folder for the mpi version).

## Example of utilization for BFGS or CMAES

__Example algo.edp__

```freefem
real[int] b(10),u(10);
func real J(real[int] & u)
{
real s=0;
for (int i=0;i<u.n;i++)
s +=(i+1)*u[i]*u[i]*0.5 - b[i]*u[i];
cout << "J ="<< s << " u =" << u[0] << " " << u[1] << "...\n" ;
return s;
}

// the grad of J (this is a affine version (the RHS is in )
func real[int] DJ(real[int] &u)
{
for (int i=0;i<u.n;i++)
u[i]=(i+1)*u[i]-b[i];
return u; // return of global variable ok
};

b=1; u=2; // set right hand side and initial gest
BFGS(J,dJ,u,eps=1.e-6,nbiter=20,nbiterline=20);
cout << "BFGS: J(u) = " << J(u) << endl;

```

%\end{example} $\codered$

Using the CMA evolution strategy is almost the same, except that, as it is a derivative free optimizer, the `:::freefem dJ` argument is omitted and there are some other named parameters to control the behaviour of the algorithm. With the same objective function as above, an example of utilization would be (see `:::freefem cmaes-VarIneq.edp` for a complete example):

```freefem
load "ff-cmaes"
... //define J, u and all here
real min = cmaes(J,u,stopTolFun=1e-6,stopMaxIter=3000);
cout << "minimal value is " << min << " for u = " << u << endl;
```

This algorithm works with a normal multivariate distribution in the parameters space and try to adapt its covariance matrix using the information provides by the successive function evaluations (see \cite{hansen} $\codered$ for more details). Thus, some specific parameters can be passed to control the starting distribution, size of the sample generations etc... Named parameters for this are the following :

 * `:::freefem seed=` Seed for random number generator (`:::freefem val` is an integer). No specified value will lead to a clock based seed initialization.

 * `:::freefem initialStdDev=` Value for the standard deviations of the initial covariance matrix ( `:::freefem val` is a real). If the value $\sigma$ is passed, the initial covariance matrix will be set to $\sigma I$. The expected initial distance between initial $X$ and the $argmin$ should be roughly initialStdDev. Default is 0.3.

 * `:::freefem initialStdDevs=` Same as above except that the argument is an array allowing to set a value of the initial standard deviation for each parameter. Entries differing by several orders of magnitude should be avoided (if it can't be, try rescaling the problem).

 * `:::freefem stopTolFun=` Stops the algorithm if function values differences are smaller than the passed one, default is $10^{-12}$.

 * `:::freefem stopTolFunHist=` Stops the algorithm if function value differences of the best values are smaller than the passed one, default is 0 (unused).

 * `:::freefem stopTolX=` Stopping criteria triggered if step sizes in the parameters space are
 smaller than this real value, default is 0.

 * `:::freefem stopTolXFactor=` Stopping criteria triggered when the standard deviation increases more than this value. The default value is $10^{3}$.

 * `:::freefem stopMaxFunEval=` Stops the algorithm when `:::freefem stopMaxFunEval` function evaluations have been done. Set to $900(n+3)^{2}$ by default, where $n$ is the parameters space dimension.

 * `:::freefem stopMaxIter=` Integer stopping the search when `:::freefem stopMaxIter` generations has been sampled. Unused by default.

 * `:::freefem popsize=` Integer value used to change the sample size. The default value is $4+ \lfloor 3\ln (n) \rfloor$, see \cite{hansen} $\codered$ for more details. Increasing the population size usually improves the global search capabilities at the cost of an at most linear reduction of the convergence speed with respect to `:::freefem popsize`.

 * `:::freefem paramFile=` This `:::freefem string` type parameter allows the user to pass all the parameters using an extern file as in Hansen's original code. More parameters related to the CMA-ES algorithm can be changed with this file. A sample of it can be found in the  `:::freefem examples++-load/ffCMAES/` folder under the name `:::freefem initials.par`. Note that the parameters passed to the CMAES function in the FreeFem script will be ignored if an input parameters file is given.

## IPOPT

The  `:::freefem ff-Ipopt` package is an interface for the IPOPT \cite{ipopt} $\codered$ optimizer. IPOPT is a software library for large scale, non-linear, constrained optimization. Detailed informations about it are in \cite{ipopt} $\codered$ and [https://projects.coin-or.org/Ipopt](https://projects.coin-or.org/Ipopt). It implements a primal-dual interior point method along with filter method based line searchs.
IPOPT needs a direct sparse symmetric linear solver. If your version of FreeFem has been compiled with the  `:::freefem --enable-downlad` tag, it will automatically be linked with a sequential version of MUMPS. An alternative to MUMPS would be to download the HSL subroutines (see [http://www.coin-or.org/Ipopt/documentation/node16.html](http://www.coin-or.org/Ipopt/documentation/node16.html)) and place them in the `:::freefem /ipopt/Ipopt-3.10.2/ThirdParty/HSL` directory of the FreeFem++ downloads folder before compiling.

### Short description of the algorithm

In this section, we give a very brief glimpse at the underlying mathematics of IPOPT. For a deeper introduction on interior methods for nonlinear smooth optimization, one may consults \cite{ipintro} $\codered$, or \cite{ipopt} $\codered$ for more IPOPT specific elements. IPOPT is designed to perform optimization for both equality and inequality constrained problems. Though, nonlinear inequalities are rearranged before the beginning of the optimization process in order to restrict the panel of nonlinear constraints to those of the equality kind. Each nonlinear inequality ones are transformed into a pair of simple bound inequality and nonlinear equality constraint by the introduction of as many slack variables as is needed : $c_{i}(x)\leq 0$ becomes $c_{i}(x) + s_{i} = 0$ and $s_{i}\leq 0$, where $s_{i}$ is added to the initial variables of the problems $x_{i}$. Thus, for convenience, we will assume that the minimization problem does not contain any nonlinear inequality constraint. It means that, given a function $f:\mathbb{R}^{n}\mapsto\mathbb{R}$, we want to find :

\begin{equation}\label{minimproblem}
 \begin{array} {c}
 x_{0} = \underset{x\in V}{\operatorname{argmin}} f(x) \\
 \mathrm{with}\ V = \left\lbrace x\in\R^{n}\ \vert\ c(x)= 0 \ \text{and}\ x_{l}\leq x\leq x_{u}\right\rbrace
 \end{array}
\end{equation}

Where $c:\R^{n}\rightarrow\R^{m}$ and $x_{l},x_{u}\in\R^{n}$ and inequalities hold componentwise.
The $f$ function as well as the constraints $c$ should be twice-continuously differentiable.

As a barrier method, interior points algorithms try to find a Karush-Kuhn-Tucker point for (\ref{minimproblem}) $\codered$ by solving a sequence of problems, unconstrained with respect to the inequality constraints, of the form :

\begin{equation}\label{barrier}
\mathrm{for\ a\ given\ }\mu > 0,\ \mathrm{find}\ x_{\mu} = \underset{x\in\R^{n}\ \vert\ c(x)=0}{\operatorname{argmin}}\ B(x,\mu)
\end{equation}

Where $\mu$ is a positive real number and $B(x,\mu) = f(x) - \displaystyle{\mu\sum_{i=1}^{n} \ln (x_{u,i}-x_{i})} - \displaystyle{\mu\sum_{i=1}^{m} \ln(x_{i}-x_{l,i})}$.

The remaining equality constraints are handled with the usual Lagrange multipliers method. If the sequence of barrier parameters $\mu$ converge to 0, intuition suggests that the sequence of minimizers of (\ref{barrier}) $\codered$ converge to a local constrained minimizer of (\ref{minimproblem}) $\codered$. For a given $\mu$, (\ref{barrier}) $\codered$ is solved by finding $(x_{\mu},\lambda_{\mu})\in\R^{n}\times\R^{m}$ such that :

\begin{equation}\label{muproblem}
 \nabla B(x_{\mu},\mu) + \displaystyle{\sum_{i=1}^{m}\lambda_{\mu,i}\nabla c_{i}(x_{\mu})}
 = \nabla B(x_{\mu},\mu) + J_{c}(x_{\mu})^{T}\lambda_{\mu} = 0\quad \text{and}\quad c(x_{\mu}) = 0
\end{equation}

The derivations for $\nabla B$ only holds for the $x$ variables, so that :

\[
\nabla B(x,\mu) = \nabla f(x) + \left(\begin{matrix}\mu/(x_{u,1}-x_{1}) \\ \vdots \\ \mu/(x_{u,n}-x_{n})\end{matrix}\right) - \left(\begin{matrix}\mu/(x_{1}-x_{l,1}) \\ \vdots \\ \mu/(x_{n}-x_{l,n})\end{matrix}\right)
\]

If we respectively call $z_{u}(x,\mu) = \left(\mu/(x_{u,1}-x_{1}),\dots, \mu/(x_{u,n}-x_{n})\right)$ and $z_{l}(x,\mu)$ the other vector appearing in the above equation, then the optimum $(x_{\mu},\lambda_{\mu})$
satisfies :

\begin{equation}\label{muproblemlambda}
\nabla f(x_{\mu}) + J_{c}(x_{\mu})^{T}\lambda_{\mu}+ z_{u}(x_{\mu},\mu) - z_{l}(x_{\mu},\mu) = 0 \quad \text{and} \quad c(x_{\mu}) = 0
\end{equation}

In this equation, the $z_l$ and $z_u$ vectors seems to play the role of Lagrange multipliers for the simple bounds inequalities, and indeed, when $\mu\rightarrow 0$, they converge toward some suitable Lagrange multipliers for the KKT conditions, provided some technical assumptions are fulfilled (see \cite{ipintro} $\codered$).

Equation \ref{muproblemlambda} is solved by performing a Newton method in order to find a solution of (\ref{muproblem}) for each of the decreasing values of $\mu$. Some order 2 conditions are also taken into account to avoid convergence to local maximizer, see \cite{ipintro} $\codered$ for precision about them. In the most classical IP algorithms, the Newton method is directly applied to (\ref{muproblem}). This is in most case inefficient due to frequent computation of infeasible points. These difficulties are avoided in Primal-Dual interior points methods where (\ref{muproblem}) is transformed into an extended system where $z_u$ and $z_l$ are treated as unknowns and the barrier problems are finding $(x,\lambda,z_u,z_l)\in\R^n\times\R^m\times\R^n\times\R^n$ such that :

\begin{equation}\label{PrimalDualIPBarrierProblem}
\left\lbrace\begin{array}{rcl}
\nabla f(x) + J_{c}(x)^{T}\lambda+ z_{u} - z_{l} & = & 0 \\
c(x) & = & 0 \\
(X_u - X) z_u - \mu e & = & 0 \\
(X - X_l) z_l - \mu e & = & 0
\end{array}\right.
\end{equation}

Where if $a$ is a vector of $\R^n$, $A$ denotes the diagonal matrix $A=(a_i \delta_{ij})_{1\leq i,j\leq n}$ and $e\in\R^{n} = (1,1,\dots,1)$. Solving this nonlinear system by the Newton methods is known as being the _primal-dual_ interior points method. Here again, more details are available in \cite{ipintro} $\codered$. Most actual implementations introduce features in order to globalize the convergence capability of the method, essentially by adding some line-search steps to the Newton algorithm, or by using trust regions. For the purpose of IPOPT, this is achieved by a _filter line search_ methods, the details of which can be found in \cite{iplinesearch} $\codered$.

More IPOPT specific features or implementation details can be found in \cite{ipopt} $\codered$. We will just retain that IPOPT is a smart Newton method for solving constrained optimization problem, with global convergence capabilities due to a robust line search method (in the sense that the algorithm will convergence no matter the initializer). Due to the underlying Newton method, the optimization process requires expressions of all derivatives up to the order 2 of the fitness function as well as those of the constraints. For problems whose hessian matrices are difficult to compute or lead to high dimensional dense matrices, it is possible to use a BFGS approximation of these objects at the cost of a much slower convergence rate.

### IPOPT in FreeFem++

Calling the IPOPT optimizer in a FreeFem++ script is done with the  `:::freefem IPOPT` function included in the `:::freefem ff-Ipopt` dynamic library. IPOPT is designed to solve constrained minimization problem in the form :

$$
\begin{array} {r l}
 	\mathrm{find}& x_{0} = \underset{x\in\R^{n}}{\operatorname{argmin}} f(x) \\
 	\mathrm{s.t.}&\left\lbrace \begin{array}{l r} \forall i\leq n,\ x_{i}^{\mathrm{lb}}\leq x_{i}\leq x_{i}^{\mathrm{ub}} & \mathrm{\ (simple\ bounds)} \\
	 \forall i\leq m,\ c_{i}^{\mathrm{lb}}\leq c_{i}(x)\leq c_{i}^{\mathrm{ub}} & \mathrm{(constraints\ functions)}
\end{array}.
$$

Where $\mathrm{ub}$ and $\mathrm{lb}$ stand for "upper bound" and "lower bound". If for some $i, 1\leq i\leq m$ we have $c_{i}^{\mathrm{lb}} = c_{i}^{\mathrm{ub}}$, it means that $c_{i}$ is an equality constraint, and an inequality one if $c_{i}^{\mathrm{lb}} < c_{i}^{\mathrm{ub}}$.

There are different ways to pass the fitness function and constraints. The more general one is to define the functions using the keyword  `:::freefem func`. Any returned matrix must be a sparse one (type  `:::freefem matrix`, not a  `:::freefem real[int,int]`) :

```freefem
func real J(real[int] &X) {...} // Fitness Function, returns a scalar
func real[int] gradJ(real[int] &X) {...} // Gradient is a vector

func real[int] C(real[int] &X) {...} // Constraints
func matrix jacC(real[int] &X) {...} // Constraints jacobian
```

!!!Warning
	In the current version of FreeFem++, returning a `:::freefem matrix` object local to a function block leads to undefined results. For each sparse matrix returning function you define, an extern matrix object has to be declared, whose associated function will overwrite and return on each call. Here is an example for `:::freefem jacC` :

	```freefem
	 matrix jacCBuffer; //just declare, no need to define yet
	 func matrix jacC(real[int] &X)
	 {
	 ...//fill jacCBuffer
	 return jacCBuffer;
	 }
	```

!!!Warning
	IPOPT requires the structure of each matrix at the initialization of the algorithm. Some errors may occur if the matrices are not constant and are built with the `:::freefem matrix A = [I,J,C]` syntax, or with an intermediary full matrix (`:::freefem real[int,int]`), because any null coefficient is discarded during the construction of the sparse matrix. It is also the case when making matrices linear combinations, for which any zero coefficient will result in the suppression of the matrix from the combination. Some controls are available to avoid such problems. Check the named parameters descriptions (`:::freefem checkindex`, `:::freefem structhess` and `:::freefem structjac` can help). We strongly advice to use  `:::freefem varf` as much as possible for the matrix forging.

$\codered$ PAS SUR DE LA OU SE TROUVE LA FIN DU 2EME WARNING

The hessian returning function is somewhat different because it has to be the hessian of the lagrangian function :

$$
(x,\sigma_{f},\lambda)\mapsto\sigma_{f}\nabla^{2}f(x)+\displaystyle{\sum_{i=1}^{m}\lambda_{i}\nabla^{2}c_{i}(x)} where \lambda\in\R^{m} and \sigma\in\R
$$

where $\lambda\in\R^{m}$ and $\sigma\in\R$ $\codered$ FIX EQUATION WHERE=... (FRANCK). Your hessian function should then have the following prototype :

```freefem
matrix hessianLBuffer; // Just to keep it in mind
func matrix hessianL(real[int] &X,real sigma,real[int] &lambda) {...}
```

If the constraints functions are all affine, or if there are only simple bounds constraints or no constraint at all, the lagrangian hessian is equal to the fitness function hessian, one can then omit the `:::freefem sigma` and `:::freefem lambda` parameters :

```freefem
matrix hessianJBuffer;
func matrix hessianJ(real[int] &X) {...} // Hessian prototype when constraints are affine
```

When these functions are defined, IPOPT is called this way :

```freefem
real[int] Xi = ... ; //starting point
IPOPT(J,gradJ,hessianL,C,jacC,Xi, /*some named parameters*/ );
```

If the hessian is omitted, the interface will tell IPOPT to use the (L)BFGS approximation (it can also be enabled with a named parameter, see further). Simple bounds or unconstrained problems do not require the constraints part, so the following expressions are valid :

```freefem
IPOPT(J,gradJ,C,jacC,Xi, ... ); // IPOPT with BFGS
IPOPT(J,gradJ,hessianJ,Xi, ... ); // Newton IPOPT without constraints
IPOPT(J,gradJ,Xi, ... ); // BFGS, no constraints
```

Simple bounds are passed using the `:::freefem lb` and `:::freefem ub` named parameters, while constraints bounds are passed with the `:::freefem clb` and `:::freefem cub` ones. Unboundedness in some directions can be achieved by using the $1e^{19}$ and $-1e^{19}$ values that IPOPT recognizes as $+\infty$ and $-\infty$ :

```freefem
real[int] xlb(n),xub(n),clb(m),cub(m);
... //fill the arrays...
IPOPT(J,gradJ,hessianL,C,jacC,Xi,lb=xlb,ub=xub,clb=clb,cub=cub, /*some other named parameters*/ );
```

__P2 fitness function and affine constraints function :__ In the case where the fitness function or constraints function can be expressed respectively in the following forms :

$$
\begin{array}{c c}
	\forall x\in\R^{n},\ f(x) = \frac{1}{2}\left\langle Ax,x \right\rangle + \left\langle b,x\right\rangle & (A,b)\in\mathcal{M}_{n,n}(\R)\times\R^{n} \\
	\mathrm{or} ,\ C(x) = Ax + b & (A,b)\in\mathcal{M}_{n,m}(\R)\times\R^{m}
\end{array}
$$

where $A$ and $b$ are constant, it is possible to directly pass the $(A,b)$ pair instead of defining 3 (or 2) functions. It also indicates to IPOPT that some objects are constant and that they have to be evaluated only once, thus avoiding multiple copies of the same matrix. The syntax is :

```freefem
//Affine constraints with "standard" fitness function
matrix A= ... ; // Linear part of the constraints
real[int] b = ... ; // Constant part of constraints
IPOPT(J,gradJ,hessianJ, [A,b] ,Xi, /*bounds and named params*/);
// [b,A] would work as well... Scatterbrains pampering...
```

Note that if you define the constraints in this way, they doesn't contribute to the hessian, so the hessian should only take one  `:::freefem real[int]` as argument.

```freefem
// Affine constraints and P2 fitness func:
matrix A= ... ; // Bilinear form matrix
real[int] b = ... ; // Linear contribution to f
matrix Ac= ... ; // Linear part of the constraints
real[int] bc= ... ; // Constant part of constraints
IPOPT([A,b], [Ac,bc] ,Xi, /*bounds and named params*/);
```

If both objective and constraints functions are given this way, it automatically activates the IPOPT  `:::freefem mehrotra_algorithm` option (better for linear and quadratric programming according to the documentation). Otherwise, this option can only be set through the option file (see the named parameters section).

A spurious case is the one of defining $f$ in this manner while using standard functions for the constraints :

```freefem
matrix A= ... ; // Bilinear form matrix
real[int] b = ... ; // Linear contribution to f
func real[int] C(real[int] &X) {...} // Constraints
func matrix jacC(real[int] &X) {...} // Constraints jacobian
IPOPT([A,b],C,jacC,Xi, /*bounds and named params*/);
```

Indeed, when passing `:::freefem [A,b]` in order to define $f$, the lagrangian hessian is automatically build has the constant $x \mapsto A$ function, with no way to add possible constraints contributions, leading to incorrect second order derivatives. So, a problem should be defined like that in only two cases :

1. constraints are nonlinear but you want to use the BFGS mode (then add `:::freefem bfgs=1` to the named parameter),
2. constraints are affine, but in this case, why not passing them in the same way?

Here are some other valid definitions of the problem (cases when $f$ is a pure quadratic or linear form, or $C$ a pure linear function, etc...) :

```freefem
// Pure quadratic f - A is a matrix:
IPOPT(A, /*constraints args*/, Xi, /*bounds and named params*/);
// Pure linear f - b i a real[int] :
IPOPT(b, /*constraints args*/, Xi, /*bounds and named params*/);
// linear constraints - Ac is a matrix
IPOPT(/*fitness func args*/, Ac, Xi, /*bounds and named params*/);
```

__Returned Value :__ The `:::freefem IPOPT` function returns an error code of type `:::freefem int`. A zero value is obtained when the algorithm succeeds and positive values reflects the fact that IPOPT encounters minor troubles. Negative values reveals more problematic cases. The associated IPOPT return tags are listed in the table below. The [IPOPT pdf documentation](https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw) provides more accurate description of these return status :

|Success|Failures|
|:----|:----|
|0 `:::freefem Solve_Succeeded`<br>1 `:::freefem Solved_To_Acceptable_Level`<br>2 `:::freefem Infeasible_Problem_Detected`<br>3 `:::freefem Search_Direction_Becomes_Too_Small`<br>4 `:::freefem Diverging_Iterates`<br>5 `:::freefem User_Requested_Stop`<br>6 `:::freefem Feasible_Point_Found`|-1 `:::freefem Maximum_Iterations_Exceeded`<br>-2 `:::freefem Restoration_Failed`<br>-3 `:::freefem Error_In_Step_Computation`<br>-4 `:::freefem Maximum_CpuTime_Exceeded`|

|Problem definition issues|Critical errors|
|:----|:----|
|-10 `:::freefem Not_Enough_Degrees_Of_Freedom`<br>-11 `:::freefem Invalid_Problem_Definition`<br>-12 `:::freefem Invalid_Option`<br>-13 `:::freefem Invalid_Number_Detected`|-100 `:::freefem Unrecoverable_Exception`<br>-101 `:::freefem NonIpopt_Exception_Thrown`<br>-102 `:::freefem Insufficient_Memory`<br>-199 `:::freefem Internal_Error`|

__Named Parameters :__ The available named parameters in this interface are those we thought to be the most subject to variations from one optimization to another, plus a few ones that are interface specific. Though, as one could see at [http://www.coin-or.org/Ipopt/documentation/node59.html](http://www.coin-or.org/Ipopt/documentation/node59.html), there are many parameters that can be changed within IPOPT, affecting the algorithm behaviour. These parameters can still be controlled by placing an option file in the execution directory. Note that [IPOPT's pdf documentation](https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw) may provides more informations than the previously mentioned online version for certain parameters. The in-script available parameters are :


 * `:::freefem lb`, `:::freefem ub` : `:::freefem real[int]` for lower and upper simple bounds upon the search variables, must be of size $n$ (search space dimension). If two components of same index in these arrays are equal then the corresponding search variable is fixed. By default IPOPT will remove any fixed variable from the optimization process and always use the fixed value when calling functions. It can be changed using the `:::freefem fixedvar` parameter.

 * `:::freefem clb`, `:::freefem cub` : `:::freefem real[int]` of size $m$ (number of constraints) for lower and upper constraints bounds. Equality between two components of same index $i$ in `:::freefem clb` and  `:::freefem cub` reflect an equality constraint.

 * `:::freefem structjacc` : To pass the greatest possible structure (indexes of non null coefficients) of the constraints jacobian under the form `:::freefem [I,J]` where `:::freefem I` and `:::freefem J` are two integer arrays. If not defined, the structure of the constraints jacobian, evaluated in `:::freefem Xi`, is used (no issue if the jacobian is constant or always defined with the same `:::freefem varf`, hazardous if it is with triplet array or if a full matrix is involved).

 * `:::freefem structhess` : Same as above but for the hessian function (unused if $f$ is P2 or less and constraints are affine). Here again, keep in mind that it is the hessian of the lagrangian function (which is equal to the hessian of $f$ only if constraints are affine). If no structure is given with this parameter, the lagrangian hessian is evaluated on the starting point, with $\sigma=1$ and $\lambda = (1,1,\dots,1)$ (it is safe if all the constraints and fitness function hessians are constant or build with `:::freefem varf`, and here again less reliable if built with triplet array or full matrix).

 * `:::freefem checkindex` : A `:::freefem bool` that triggers an index dichotomic search when matrices are copied from FreeFem functions to IPOPT arrays. It is used to avoid wrong index matching when some null coefficients are removed from the matrices by FreeFem. It will not solve the problems arising when a too small structure has been given at the initialization of the algorithm. Enabled by default (except in cases where all matrices are obviously constant).

 * `:::freefem warmstart` : If set to `:::freefem true`, the constraints dual variables $\lambda$, and simple bounds dual variables are initialized with the values of the arrays passed to `:::freefem lm`, `:::freefem lz` and `:::freefem uz` named parameters (see below).

 * `:::freefem lm` :  `:::freefem real[int]` of size $m$, which is used to get the final values of the constraints dual variables $\lambda$ and/or initialize them in case of a warm start (the passed array is also updated to the last dual variables values at the end of the algorithm).

 * `:::freefem lz`, `:::freefem uz` : `:::freefem real[int]` of size $n$ to get the final values and/or initialize (in case of warm start) the dual variables associated to simple bounds.

 * `:::freefem tol` : `:::freefem real`, convergence tolerance for the algorithm, the default value is $10^{-8}$.

 * `:::freefem maxiter` : `:::freefem int`, maximum number of iterations with 3000 as default value.

 * `:::freefem maxcputime` : `:::freefem real` value, maximum runtime duration. Default is $10^{6}$ (almost 11 days and a half).

 * `:::freefem bfgs` : `:::freefem bool` enabling or not the (low-storage) BFGS approximation of the lagrangian hessian. It is set to false by default, unless there is no way to compute the hessian with the functions that have been passed to IPOPT.

 * `:::freefem derivativetest` : Used to perform a comparison of the derivatives given to IPOPT with finite differences computation. The possible `:::freefem string` values are : `:::freefem "none"` (default),  `:::freefem "first-order"`, `:::freefem "second-order"` and `:::freefem "only-second-order"`. The associated derivative error tolerance can be changed via the option file. One should not care about any error given by it before having tried, and failed, to perform a first optimization.

 * `:::freefem dth` :  Perturbation parameter for the derivative test computations with finite differences. Set by default to $10^{-8}$.

 * `:::freefem dttol` : Tolerance value for the derivative test error detection (default value unknown yet, maybe $10^{-5}$).

 * `:::freefem optfile` : `:::freefem string` parameter to specify the IPOPT option file name. IPOPT will look for a `:::freefem ipopt.opt` file by default. Options set in the file will overwrite those defined in the FreeFem script.

 * `:::freefem printlevel` : An `:::freefem int` to control IPOPT output print level, set to 5 by default, the possible values are from 0 to 12. A description of the output informations is available in the [PDF documentation](https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw) of IPOPT.

 * `:::freefem fixedvar` : `:::freefem string` for the definition of simple bounds equality constraints treatment : use `:::freefem "make_parameter"` (default value) to simply remove them from the optimization process (the functions will always be evaluated with the fixed value for those variables), `:::freefem "make_constraint"` to treat them as any other constraint or `:::freefem "relax_bounds"` to relax fixing bound constraints.

 * `:::freefem mustrategy` : a `:::freefem string` to choose the update strategy for the barrier parameter $\mu$. The two possible tags are `:::freefem "monotone"`, to use the monotone (Fiacco-McCormick) strategy, or `:::freefem "adaptive"` (default setting).

 * `:::freefem muinit` : `:::freefem real` positive value for the barrier parameter initialization. It is only relevant when `:::freefem mustrategy` has been set to `:::freefem monotone`.

 * `:::freefem pivtol` : `:::freefem real` value to set the pivot tolerance for the linear solver. A smaller number pivots for sparsity, a larger number pivots for stability. The value has to be in the $[0,1]$ interval and is set to $10^{-6}$ by default.

 * `:::freefem brf` : Bounds relax factor : before starting the optimization, the bounds given by the user are relaxed. This option sets the factor for this relaxation. If it is set to zero, then the bounds relaxation is disabled. This `:::freefem real` has to be positive and its default value is $10^{-8}$.

 * `:::freefem objvalue` : An identifier to a `:::freefem real` type variable to get the last value of the objective function (best value in case of succes).

 * `:::freefem mumin` : Minimal value for the barrier parameter $\mu$, a `:::freefem real` with $10^{-11}$ as default value.

 * `:::freefem linesearch` : A boolean which disables the line search when set to `:::freefem false`. The line search is activated by default. When disabled, the method becomes a standard Newton algorithm other the primal-dual system. The global convergence is then no longer assured, meaning that many initializers could lead to diverging iterates. But on the other hand, it can be useful when trying to catch a precise local minimum without having some out of control process making the iterate caught by some other near optimum.

## Some short examples using IPOPT

!!!example
	__Example IpoptVI.edp__

	A very simple example consisting in, given two functions $f $ and $g$ (defined on $\Omega\subset\R^{2}$), minimizing $J(u) = \displaystyle{\frac{1}{2}\int_{\Omega} \vert\nabla u\vert^{2} - \int_{\Omega}fu}\ $, with $u\leq g$ almost everywhere :

	```freefem
	load "ff-Ipopt"; // Load the interface
	int nn=20; // Mesh quality
	mesh Th=square(nn,nn); // Build a square mesh
	fespace Vh(Th,P1); // Finite element space

	func f = 1.; // rhs function
	real r=0.03,s=0.1; // Some parameters for g
	func g = r - r/2*exp(-0.5*(square(x-0.5)+square(y-0.5))/ square(s));
	// g is constant minus a gaussian

	macro Grad(u) [dx(u),dy(u)]//the gradient operator
	varf vP(u,v) = int2d(Th)(Grad(u)'*Grad(v)) - int2d(Th)(f*v);
	```

	Here we build the matrix and second member associated to the functional to minimize once and for all. The  `:::freefem [A,b]` syntax for the fitness function is then used to pass it to IPOPT.

	```freefem
	matrix A = vP(Vh,Vh,solver=CG);
	real[int] b = vP(0,Vh);
	```

	We use simple bounds to impose the boundary condition $u=0$ on $\partial\Omega$, as well as the $u\leq g$ condition.

	```freefem
	Vh lb=-1.e19; // Lower-unbounded in the interior
	Vh ub=g; // Upper-bounded by g in the interior
	varf vGamma(u,v) = on(1,2,3,4,u=1);
	real[int] onGamma=vGamma(0,Vh);
	ub[] = onGamma ? 0. : ub[]; // Enforcing the boundary condition
	lb[] = onGamma ? 0. : lb[];

	Vh u=0; // Starting point
	IPOPT([A,b],u[],lb=lb[],ub=ub[]); // Solve the problem
	plot(u,wait=1);
	```

!!!example
	__Example IpoptVI2.edp__
	Let $\Omega$ be a domain of $\mathbb{R}^{2}$, $f_{1}, f_{2}\in L^{2}(\Omega)$ and $g_{1}, g_{2} \in L^{2}(\partial\Omega)$ four given functions with $g_{1}\leq g_{2}$ almost everywhere.
	We define the space :

	$$
	V = \left\lbrace (v_{1},v_{2})\in H^{1}(\Omega)^{2} ; v_{1}\vert_{\partial\Omega}=g_{1}, v_{2}\vert_{\partial\Omega}=g_{2}, v_{1}\leq v_{2}\ \mathrm{a.e.}\ \right\rbrace
	$$

	as well as the functional $J:H^{1}(\Omega)^{2}\longrightarrow \mathbb{R}$:

	$$
	J(v_{1},v_{2}) = \displaystyle{\frac{1}{2}\int_{\Omega}\vert\nabla v_{1}\vert^{2} - \int_{\Omega} f_{1}v_{1} + \frac{1}{2}\int_{\Omega}\vert\nabla v_{2}\vert^{2} - \int_{\Omega} f_{2}v_{2}}
	$$

	The problem consists in finding (numerically) two functions $(u_{1},u_{2}) = \underset{(v_{1},v_{2})\in V}{\operatorname{argmin}} J(v_{1},v_{2})$.

	```freefem
	load "ff-IpOpt";

	mesh Th=square(10,10);
	fespace Vh(Th,[P1,P1] );
	fespace Wh(Th,[P1] );
	int iter=0;

	func f1 = 10; // Right hand sides
	func f2 = -15;
	func g1 = -0.1; // Boundary conditions functions
	func g2 = 0.1;

	while(++iter) // Mesh adaptation loop
	{
	macro Grad(u) [dx(u),dy(u)] // Gradient macro
	varf vP([u1,u2],[v1,v2]) = int2d(Th)(Grad(u1)'*Grad(v1)+ Grad(u2)'*Grad(v2))
	- int2d(Th)(f1*v1+f2*v2);

	matrix A = vP(Vh,Vh); // Fitness function matrix...
	real[int] b = vP(0,Vh); // And linear form

	int[int] II1=[0],II2=[1]; // Constraints matrix
	matrix C1 = interpolate (Wh,Vh, U2Vc=II1);
	matrix C2 = interpolate (Wh,Vh, U2Vc=II2);
	matrix CC = -1*C1 + C2; // u2 - u1 >0
	Wh cl=0; // Constraints lower bounds (no upper bounds)

	// Boundary conditions
	varf vGamma([u1,u2],[v1,v2]) = on(1,2,3,4,u1=1,u2=1);
	real[int] onGamma=vGamma(0,Vh);
	Vh [ub1,ub2]=[g1,g2];
	Vh [lb1,lb2]=[g1,g2];
	ub1[] = onGamma ? ub1[] : 1e19 ; // Unbounded in interior
	lb1[] = onGamma ? lb1[] : -1e19 ;

	Vh [u1,u2]=[0,0]; // Starting point

	IPOPT([b,A],CC,u1[],lb=lb1[],ub=ub1[],clb=cl[]);

	plot(u1,u2,wait=1,nbiso=60,dim=3);
	if(iter > 1) break;
	Th= adaptmesh(Th,[u1,u2],err=0.004,nbvx=100000);
	}
	```

	<center>

	|Fig. 8.1: Numerical Approximation of the Variational Inequality|
	|:----:|
	|![VarIneqFill](images/VarIneqFill.jpg)|
	|![VarIneqIso](images/VarIneqIso.jpg)|

	</center>

## 3D constrained minimal surface with IPOPT

### Area and volume expressions

This example aimed at numerically solving some constrained minimal surface problems with the IPOPT algorithm. We restrain to $C^{k}$ ($k\geq 1$), closed, spherically parametrizable surfaces, i.e. surfaces $S$ such that :

$$
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
$$

Where the components are expressed in the spherical coordinate system. Let's call $\Omega$ the $[0,2\pi ]\times[0,\pi]$ angular parameters set. In order to exclude self crossing and opened shapes, the following assumptions upon $\rho$ are made :

$$
\rho \geq 0\ \ \mathrm{and}\ \ \forall \phi, \rho(0,\phi) = \rho(2\pi,\phi)
$$

For a given function $\rho$ the first fundamental form (the metric) of the defined surface has the following matrix representation :

\begin{equation}\label{msfff}
G =
\left(
\begin{array}{c c}
	\rho^{2}\sin^{2}(\phi) + (\partial_{\theta}\rho)^{2} &\partial_{\theta}\rho\partial_{\phi}\rho \\
	\partial_{\theta}\rho\partial_{\phi}\rho & \rho^{2} + (\partial_{\phi}\rho)^{2} \\
\end{array}
\right)
\end{equation}

This metric is used to express the area of the surface. Let $g=\det(G)$, then we have :

\begin{equation}\label{msarea}
	\mathcal{A}(\rho)
	= \MyInt{\Omega}{\left\| \partial_{\theta} X \wedge \partial_{\phi} X \right\|}
	=\MyInt{\Omega}{\sqrt{g}}
	=\MyInt{\Omega}{\sqrt{ \rho^{2}(\partial_{\theta}\rho)^{2} + \rho^{4}\sin^{2}(\phi) + \rho^{2}(\partial_{\phi}\rho)^{2}\sin^{2}(\phi)}d\theta d\phi}
\end{equation}

The volume of the space enclosed within the shape is easier to express :

\begin{equation}\label{msvolume}
	\mathcal{V}(\rho)
	= \MyInt{\Omega}{\int_{0}^{\rho(\theta,\phi)} r^{2}\sin(\phi) dr d\theta d\phi}
	= \frac{1}{3}\MyInt{\Omega}{\rho^{3} \sin(\phi) d\theta d\phi}
\end{equation}

### Derivatives

In order to use a newton based interior point optimization algorithm, one must be able to evaluate the derivatives of $\mathcal{A}$ and $\mathcal{V}$ with respect to $rho$. Concerning the area we have the following result :

$$
	\forall v\in C^{1}(\Omega) \ , \ \langle d\mathcal{A}(\rho),v\rangle
	= \MyInt{\Omega}{\frac{1}{2} \frac{ d\bar{g}(\rho)(v)}{\sqrt{g}}d\theta d\phi }
$$

Where $\bar{g}$ is the application mapping the $(\theta,\phi) \mapsto g(\theta,\phi)$ scalar field to $\rho$. This leads to the following expression, easy to transpose in a freefem script using :

\begin{equation}\label{msdarea}
	\begin{array}{r c l}
		\forall v\in C^{1}(\Omega) \ , \ \langle d\mathcal{A}(\rho),v\rangle &
		= & \MyInt{\Omega}{ \left(2\rho^{3}\sin^{2}(\phi) + \rho(\partial_{\theta}\rho)^{2} + \rho(\partial_{\phi}\rho)^{2}\sin^{2}(\phi) \right) v} \\
		& & +\MyInt{\Omega}{\ \rho^{2}\partial_{\theta}\rho\partial_{\theta} v\ + \ \rho^{2}\partial_{\phi}\rho\sin^{2}(\phi)\partial_{\phi} v }
	\end{array}
\end{equation}

With a similar approach, one can derive an expression for second order derivatives. Though comporting no specific difficulties, the detailed calculus are tedious, the result is that
these derivatives can be write using a $3\times 3$ matrix $\mathbf{B}$ whose coefficients are expressed in term of $\rho$ and its derivatives with respect to $\theta$ and $\phi$, such that :

\begin{equation}\label{msd2area}
	\forall (w,v)\in C^{1}(\Omega)\ ,\ d^{2}\mathcal{A}(\rho)(w,v) = \MyInt{\Omega}
	{
		\left(\begin{array}{c c c} w & \partial_{\theta} w & \partial_{\phi} w \end{array}\right)
		\mathbf{B}
	}	\left( \begin{array}{c} v \\ \partial_{\theta} v \\ \partial_{\phi} v \end{array} \right) d\theta d\phi
\end{equation}

Deriving the volume function derivatives is again an easier task. We immediately get the following expressions :

\begin{equation}\label{msdvolume}
	\begin{array}{r c l}
		\forall v\ ,\ \langle d\mathcal{V}(\rho),v\rangle & = & \MyInt{\Omega}{\rho^{2}\sin(\phi)v\ d\theta d\phi} \\
		\forall w,v\ , d^{2}\mathcal{V}(\rho)(w,v) & = & \MyInt{\Omega}{2\rho\sin(\phi)wv\ d\theta d\phi}
	\end{array}
\end{equation}

### The problem and its script :

The whole code is available in  `:::freefem IpoptMinSurfVol.edp`. We propose to solve the following problem :

__Example__

Given a positive function $\rho_{\mathrm{object}}$ piecewise continuous, and a scalar $\mathcal{V}_{\mathrm{max}} > \mathcal{V}(\rho_{\mathrm{object}})$, find $\rho_{0}$ such that :

$$
\rho_{0} = \underset{\rho\in C^{1}(\Omega)}{\operatorname{argmin}}\ \mathcal{A}(\rho)\ ,\ \mathrm{s.t.}\ \rho_{0}\geq\rho_{\mathrm{object}} \ \mathrm{and\ } \mathcal{V}(\rho_{0})\leq \mathcal{V}_{\mathrm{max}}
$$

If $\rho_{\mathrm{object}}$ is the spherical parametrization of the surface of a 3-dimensional object (domain) $\mathcal{O}$, it can be interpreted as finding the surface with minimal area enclosing the object with a given maximal volume. If $\mathcal{V}_{\mathrm{max}}$ is close to $\mathcal{V}(\rho_{\mathrm{object}})$, so should be $\rho_{0}$ and $\rho_{\mathrm{object}}$. With higher values of $\mathcal{V}_{\mathrm{max}}$, $\rho$ should be closer to the unconstrained minimal surface surrounding $\mathcal{O}$ which is obtained as soon as $\mathcal{V}_{\mathrm{max}} \geq \frac{4}{3}\pi \|\rho_{\mathrm{object}}\|_{\infty}^{3}$ (sufficient but not necessary).

It also could be interesting to solve the same problem with the constraint $\mathcal{V}(\rho_{0})\geq \mathcal{V}_{\mathrm{min}}$ which lead to a sphere when $\mathcal{V}_{\mathrm{min}} \geq \frac{1}{6}\pi \mathrm{diam}(\mathcal{O})^{3} $ and move toward the solution of the unconstrained problem as $\mathcal{V}_{\mathrm{min}}$ decreases.

We start by meshing the domain $[0,2\pi ]\times\ [0,\pi ]$, then a periodic P1 finite elements space is defined.

```freefem
load "msh3";
load "medit";
load "ff-Ipopt";

int np=40; // Initial mesh quality parameter
mesh Th = square(2*np,np,[2*pi*x,pi*y]);

fespace Vh(Th,P1,periodic=[[2,y],[4,y]]);
Vh startshape=5; // Initial shape
```

We create some finite element functions whose underlying arrays will be used to store the values of dual variables associated to all the constraints in order to reinitialize the algorithm with it in the case where we use mesh adaptation. Doing so, the algorithm will almost restart at the accuracy level it reached before mesh adaptation, thus saving many iterations.

```freefem
Vh uz=1.,lz=1.; // Simple bounds dual variable
real[int] lm=[1]; // Dual variable for volume constraint
```

Then, follows the mesh adaptation loop, and a rendering function, `:::freefem Plot3D`, using 3D mesh to display the shape it is passed with `:::freefem medit` (the `:::freefem movemesh23` procedure often crashes when called with ragged shapes).

```freefem
int nadapt=1;
for(int kkk=0;kkk<nadapt;++kkk) // Mesh adaptation loop
{

int iter=0;//iterations count
func sin2 = square(sin(y)); // A function that will be often used

func int Plot3D(real[int] &rho,string cmm,bool ffplot) {...} // See the .edp file
```

Here are the functions related to the area computation and its shape derivative, according to equations \ref{msarea} and \ref{msdarea} $\codered$ :

```freefem
func real Area(real[int] &X)
{
 Vh rho;
 rho[] = X;
 Vh rho2 = square(rho);
 Vh rho4 = square(rho2);
 real res = int2d(Th)( sqrt( rho4*sin2
 + rho2*square(dx(rho))
 + rho2*sin2*square(dy(rho)) )
 );
 ++iter;
 plot(rho, ... /*some parameters*/ ... );
 return res;
}

func real[int] GradArea(real[int] &X) // The gradient
{
 Vh rho,rho2;
 rho[] = X;
 rho2[] = square(X);
 Vh sqrtPsi,alpha; // Psi is actually det(G)
 {//some optimizations
 Vh dxrho2 = dx(rho)*dx(rho), dyrho2 = dy(rho)*dy(rho);
 sqrtPsi = sqrt( rho2*rho2*sin2 + rho2*dxrho2 + rho2*dyrho2*sin2 );
 alpha = 2.*rho2*rho*sin2 + rho*dxrho2 + rho*dyrho2*sin2;
 }
 varf dSurface(u,v) =
 int2d(Th)(1./sqrtPsi*(alpha*v+rho2*dx(rho)*dx(v)+rho2*dy(rho)*sin2*dy(v)));
 real[int] grad = dSurface(0,Vh);
 return grad;
}
```

The function returning the hessian of the area for a given shape is a bit blurry, thus we won't show here all of equation \ref{msd2area} $\codered$ coefficients definition, they can be found in the `:::freefem edp` file.

```freefem
matrix hessianA; // The global matrix buffer

func matrix HessianArea(real[int] &X)
{
 Vh rho,rho2;
 rho[] = X;
 rho2 = square(rho);
 Vh sqrtPsi,sqrtPsi3,C00,C01,C02,C11,C12,C22,A;
 {
 ... // Definition of the above functions
 }
 varf d2Area(w,v) =
 int2d(Th)(
 1./sqrtPsi * (A*w*v + 2*rho*dx(rho)*dx(w)*v + 2*rho*dx(rho)*w*dx(v)
 + 2*rho*dy(rho)*sin2*dy(w)*v + 2*rho*dy(rho)*sin2*w*dy(v)
 + rho2*dx(w)*dx(v) + rho2*sin2*dy(w)*dy(v))
 + 1./sqrtPsi3 * (C00*w*v + C01*dx(w)*v + C01*w*dx(v) + C02*dy(w)*v
 + C02*w*dy(v) + C11*dx(w)*dx(v)
 + C12*dx(w)*dy(v) + C12*dy(w)*dx(v) + C22*dy(w)*dy(v))
 );//end of int2d and varf
 hessianA = d2Area(Vh,Vh);
 return hessianA;
}
```

And the volume related functions :

```freefem
func real Volume(real[int] &X)
{
 Vh rho;
 rho[]=X;
 Vh rho3=rho*rho*rho;
 real res = 1./3.*int2d(Th)(rho3*sin(y));
 return res;
}

func real[int] GradVolume(real[int] &X)
{
 Vh rho;
 rho[]=X;
 varf dVolume(u,v) = int2d(Th)(rho*rho*sin(y)*v);
 real[int] grad = dVolume(0,Vh);
 return grad;
}

matrix hessianV; // Buffer
func matrix HessianVolume(real[int] &X)
{
 Vh rho;
 rho[]=X;
 varf d2Volume(w,v) = int2d(Th)(2*rho*sin(y)*v*w);
 hessianV = d2Volume(Vh,Vh);
 return hessianV;
}
```

If we want to use the volume as a constraint function we must wrap it and its derivatives in some FreeFem++ functions returning the appropriate types. It is not done in the above functions in case where one wants to use it as fitness function. The lagrangian hessian also have to be wrapped since the Volume is not linear with respect to $\rho$, it has some non-null second order derivatives.

```freefem
func real[int] ipVolume(real[int] &X) {real[int] vol = [Volume(X)]; return vol;}

matrix mdV; // Buffer
func matrix ipGradVolume(real[int] &X)
{//transforms a vector into a sparse matrix
 real[int,int] dvol(1,Vh.ndof);
 dvol(0,:)=GradVolume(X);
 mdV=dvol;
 return mdV;
}

matrix HLagrangian; // Buffer
func matrix ipHessianLag(real[int] &X,real objfact,real[int] &lambda)
{
 HLagrangian = objfact*HessianArea(X) + lambda[0]*HessianVolume(X);
 return HLagrangian;
}
```

The `:::freefem ipGradVolume` function could bring some troubles during the optimization process because the gradient vector is transformed in a sparse matrix, so any null coefficient will be discarded. We are here obliged to give IPOPT the structure by hand and use the `:::freefem checkindex` named-parameter to avoid bad indexing during copies. This gradient is actually dense, there is no reason for some components to be constantly zero :

```freefem
// Sparse structure of a dense vector
int[int] gvi(Vh.ndof),gvj=0:Vh.ndof-1;
gvi=0; // Only one line
```

These two arrays will be passed to IPOPT with `:::freefem structjacc=[gvi,gvj]`. The last remaining things are the bounds definition. The simple lower bounds must be equal to the components of the P1 projection of $\rho_{object}$. And we choose $\alpha\in [0,1]$ to set $\mathcal{V}_{\mathrm{max}}$ to $(1-\alpha) \mathcal{V}(\rho_{object}) + \alpha\frac{4}{3}\pi \|\rho_{\mathrm{object}}\|_{\infty}^{3}$ :

```freefem
real e=0.1,r0=0.25,rr=2-r0;
real E=1./(e*e),RR=1./(rr*rr);
// An indented disc
func disc1 = sqrt(1./(RR+(E-RR)*cos(y)*cos(y)))*(1+0.1*cos(9*x));
// Almost a standard disc
func disc2 = sqrt(1./(RR+(E-RR)*cos(x)*cos(x)*sin2)) ;
Vh lb = max(disc1, disc2); // Glue the object parts
real Vobj = Volume(lb[]); // Object volume
real Vnvc = 4./3.*pi*pow(lb[].linfty,3); // V for no volume constraint
real alpha=0.1;
Plot3D(lb[],"object_inside",0);
real[int] clb=0.,cub=[(1-alpha)*Vobj + alpha*Vnvc];
```

Calling IPOPT :

```freefem
IPOPT(Area,GradArea,ipHessianLag,
 ipVolume,ipGradVolume,rc[], // Functions and starting point
 ub=ub[],lb=lb[],clb=clb,cub=cub, // Simple bounds and volume bounds
 checkindex=1,structjacc=[gvi,gvj], // For safe matrices copies
 maxiter=kkk<nadapt-1 ? 40:150, // Accurate optim only for last mesh adaptation iteration
 warmstart=kkk,lm=lm,uz=uz[],lz=lz[], // Warmstart handling
 tol=0.00001);

Plot3D(rc[],"Shape_at_"+kkk,0); // Displays current solution
```

At last, before closing the mesh adaptation loop, we have to perform the said adaptation. The mesh is adaptated with respect to the $X=(\rho,0,0)$ (in spherical coordinates) vector field, not directly with respect to $\rho$, otherwise the true curvature of the 3D-shape would not be well taken into account.

```freefem
if(kkk<nadapt-1)
{
 Th = adaptmesh(Th,
 rc*cos(x)*sin(y),//X
 rc*sin(x)*sin(y),//Y
 rc*cos(y),//Z
 nbvx=50000,
 periodic=[[2,y],[4,y]]); // Keeps mesh peridicity
 plot(Th);
 startshape = rc; // Shape interpolation on the new mesh
 uz=uz; // Dual variables interpolation
 lz=lz;
} // end if
} // en of mesh adaptation loop
```

Here are some pictures of the resulting surfaces obtained for decreasing values of $\alpha$ (and a slightly more complicated object than two orthogonal discs). We get back the enclosed object when $\alpha=0$ :

<center>
![minsurf3D](images/minsurf3D.jpg)
</center>


## The nlOpt optimizers

The `:::freefem ff-NLopt` package provides a FreeFem interface to the free/open-source library for nonlinear optimization, thus easing the use of several different free optimization (constrained or not) routines available online along with the PDE solver. All the algorithms are well documented in \cite{nlopt} $\codered$, thus no exhaustive informations concerning their mathematical specificities will be found here and we will focus on the way they are called in a FreeFem script. One needing detailed informations about these algorithms should visit the said cite where a description of each of them is given, as well as many bibliographical links.
Most of the gradient based algorithm of nlOpt uses a full matrix approximation of the hessian, so if you're planing to solve a large scale problem, our advise would be to use the IPOPT optimizer which definitely surpass them. Finally, an examples of use can be found in the `:::freefem examples++-load/` directory under the name `:::freefem VarIneq2.edp`.

All the nlOpt features are called that way :

```freefem
load "ff-NLopt"
... // Define J, u, and maybe grad(J), some constraints etc...
real min = nloptXXXXXX(J,u, // Unavoidable part
grad = <name of grad(J)> , // if needed
lb = // Lower bounds array
ub= // Upper bounds array
... // Some optional arguments :
// Constraints functions names,
// Stopping criterions,
// Algo. specific parameters,
// Etc...
);
```

`:::freefem XXXXXX` refers to the algorithm tag (not necessarily 6 characters long). `:::freefem u` is the starting position (a `:::freefem real[int]` type array) which will be overwritten by the algorithm, the value at the end being the found argmin. And as usual, `:::freefem J` is a function taking a `:::freefem real[int]` type array as argument and returning a `:::freefem real`. `:::freefem grad`, `:::freefem lb` and `:::freefem ub` are "half-optionnal" arguments, in the sense that they are obligatory for some routines but not all.

The possible optional named parameters are the following, note that they are not used by all algorithm s (some does not supports constraints, or a type of constraints, some are gradient-based and other are derivative free, etc...). One can refer to the table after the parameters description to check which are the named parameters supported by a specific algorithm. Using an unsupported parameter will not stop the compiler work and seldom breaks runtime, it will just be ignored. That said, when it is obvious you are missusing a routine, you will get a warning message at runtime (for exemple if you pass a gradient to a derivative free algorithm, or set the population of a non-genetic one, etc...). In the following description, $n$ stands for the dimension of the search space.

__Half-optional parameters :__


* `:::freefem grad=` The name of the function which computes the gradient of the cost function (prototype should be  `:::freefem real[int]` $\rightarrow$  `:::freefem real[int]`, both argument and result should have the size $n$). This is needed as soon s a gradient-based method is involved, ignored if defined in a derivative free context.

* `:::freefem lb`/`:::freefem ub` = Lower and upper bounds arrays ( `:::freefem real[int]` type) of size $n$. Used to define the bounds within which the search variable is allowed to move. Needed for some algorithms, optional or unsupported for others.

* `:::freefem subOpt` : Only enabled for the Augmented Lagrangian and MLSL method who need a sub-optimizer in order to work. Just pass the tag of the desired local algorithm with a  `:::freefem string`.

__Constraints related parameters (optional - unused if not specified):__

* `:::freefem IConst`/`:::freefem EConst` : Allows to pass the name of a function implementing some inequality (resp. equality) constraints on the search space. The function type must be `:::freefem real[int]` $\rightarrow$ `:::freefem real[int]` where the size of the returned array is equal to the number of constraints (of the same type - it means that all the constraints are computed in one vectorial function). In order to mix inequality and equality constraints in a same minimization attempt, two vectorial functions have to be defined and passed. See example \ref{varineqex} $\codered$ for more details about how these constraints have to be implemented.

* `:::freefem gradIConst`/`:::freefem gradEConst` : Use to provide the inequality (resp. equality) constraints gradient. These are  `:::freefem real[int]` $\rightarrow$  `:::freefem real[int,int]` type functions. Assuming we have defined a constraint function (either inequality or equality) with $p$ constraints, the size of the matrix returned by its associated gradient must be $p\times n$ (the $i$-th line of the matrix is the gradient of the $i$-th constraint). It is needed in a gradient-based context as soon as an inequality or equality constraint function is passed to the optimizer and ignored in all other cases.

* `:::freefem tolIConst`/`:::freefem tolEConst` : Tolerance values for each constraint. This is an array of size equal to the number of inequality (resp. equality) constraints. Default value os set to $10^{-12}$ for each constraint of any type.

__Stopping criteria :__

* `:::freefem stopFuncValue` : Makes the algorithm end when the objective function reaches this  `:::freefem real` value.

* `:::freefem stopRelXTol` : Stops the algorithm when the relative moves in each direction of the search space is smaller than this  `:::freefem real` value.

* `:::freefem stopAbsXTol` : Stops the algorithm when the moves in each direction of the search space is smaller than the corresponding value in this  `:::freefem real[int]` array.

* `:::freefem stopRelFTol` : Stops the algorithm when the relative variation of the objective function is smaller than this  `:::freefem real` value.

* `:::freefem stopAbsFTol` : Stops the algorithm when the variation of the objective function is smaller than this  `:::freefem real` value.

* `:::freefem stopMaxFEval` : Stops the algorithm when the number of fitness evaluations reaches this  `:::freefem integer` value.

* `:::freefem stopTime` : Stops the algorithm when the otpimization time in second exceeds this  `:::freefem real` value. This is not a strict maximum: the time may exceed it slightly, depending upon the algorithm and on how slow your function evaluation is.

	Note that when an AUGLAG or MLSL method is used, the meta-algorithm and the sub-algorithm may have different termination criteria. Thus, for algorithms of this kind, the following named parameters has been defined (just adding the SO prefix - for Sub-Optimizer) to set the ending condition of the sub-algorithm (the meta one uses the ones above) :  `:::freefem SOStopFuncValue`, `:::freefem SOStopRelXTol`, and so on... If these ones are not used, the sub-optimizer will use those of the master routine.

 __Other named parameters :__

* `:::freefem popSize` : `:::freefem integer` used to change the size of the sample for stochastic search methods. Default value is a peculiar heuristic to the chosen algorithm.

* `:::freefem SOPopSize` : Same as above, but when the stochastic search is passed to a meta-algorithm.

* `:::freefem nGradStored` : The number (`:::freefem integer` type) of gradients to "remember" from previous optimization steps: increasing this increases the memory requirements but may speed convergence. It is set to a heuristic value by default. If used with AUGLAG or MLSL, it will only affect the given subsidiary algorithm.

The following table sums up the main characteristics of each algorithm, providing the more important information about which features are supported by which algorithm and what are the unavoidable arguments they need. More details can be found in \cite{nlopt} $\codered$.

<center>
![nlopttab](images/nlopttab.png)
</center>

$\codered$

%\begin{example}[VarIneq.edp]
%\label{varineqex}
%Let $\Omega$ be a domain of $\mathbb{R}^{2}$, $f_{1}, f_{2}\in L^{2}(\Omega)$ and $g_{1}, g_{2} \in L^{2}(\partial\Omega)$ four given functions with $g_{1}\leq g_{2}$ almost everywhere.
%We define the space :
%$$V = \left\lbrace (v_{1},v_{2})\in H^{1}(\Omega)^{2} ; v_{1}\vert_{\partial\Omega}=g_{1}, v_{2}\vert_{\partial\Omega}=g_{2}, v_{1}\leq v_{2}\ \mathrm{a.e.}\ \right\rbrace$$
%as well as the functional $J:H^{1}(\Omega)^{2}\longrightarrow \mathbb{R}$:
%$$J(v_{1},v_{2}) = \displaystyle{\frac{1}{2}\int_{\Omega}\vert\nabla v_{1}\vert^{2} - \int_{\Omega} f_{1}v_{1} + \frac{1}{2}\int_{\Omega}\vert\nabla v_{2}\vert^{2} - \int_{\Omega} f_{2}v_{2}}$$
%The problem consists in finding (numerically) two functions $(u_{1},u_{2}) = \underset{(v_{1},v_{2})\in V}{\operatorname{argmin}} J(v_{1},v_{2}) $.
%
%This can be interpreted as finding $u_{1}, u_{2}$ as close as possible (in a certain sense) to the solutions of Laplace equation with respectively $f_{1}, f_{2}$ second members
% and $g_{1}, g_{2}$ Dirichlet boundary conditions with the $u_{1}\leq u_{2}$ almost everywhere constraint.
%\end{example}
%
%Here is the corresponding script to treat this variational inequality problem with one of the NLOpt algorithms. The whole script is in the `:::freefem VarIneq2.edp`.
%```freefem
% load "ff-NLopt" //Loading nlopt
% int NN = 10;
% mesh Th = square(NN,NN);
% func f1=1.; //Datas
% func f2=-1.;
% func g1=0.;
% func g2=0.1;
% int nadapt=3,iter=0;
% real starttol=1e-12,bctol=6.e-12;
% fespace Vh(Th,P1);
%```

%
%Instead of defining the $J$ functional with explicit integral on the domain, we will write it in quadratic/linear form, since the corresponding matrix and vector will not change
%along with the iterations and the FreeFem matrix-vector product is optimized, this method will significantly speed up the approximation of the solution.
%
%```freefem
% varf BVF(v,w) = int2d(Th)(0.5*dx(v)*dx(w) + 0.5*dy(v)*dy(w)); //bilinear
% varf LVF1(v,w) = int2d(Th)(f1*w);//and linear parts
% varf LVF2(v,w) = int2d(Th)(f2*w);
% matrix A = BVF(Vh,Vh);//corresponding matrix
% real[int] b1 = LVF1(0,Vh) , b2 = LVF2(0,Vh); //and vectors
%```

%Then follows a trick that we will use later to pass the boundary conditions as non-linear equality constraints (although these are of course linear ones, but this approach
%demonstrate how to use all the features of NLOpt in one shot).
%```freefem
% varf Vbord(v,w) = on(1,2,3,4,v=1);
% Vh In,Bord;
% Bord[] = Vbord(0,Vh,tgv=1); //Bord=1 on the border
% In[] = Bord[] ? 0:1; //and 0 inside
% Vh gh1=Bord*g1,gh2=Bord*g2; //gh1=g1 on the border and 0 inside
%
% //The functional :
% func real J(real[int] &X)
% {
% Vh u1,u2;
% u1[] = X(0:Vh.ndof-1); //splitting the vector
% u2[] = X(Vh.ndof:2*Vh.ndof-1);
% iter++;
%	 real[int] Au1 = A*u1[], Au2 = A*u2[];
%	 Au1 -= b1;
%	 Au2 -= b2;
%	 real val = u1[]'*Au1 + u2[]'*Au2;
%	 plot(u1,u2,...); //parameters omitted
%	 return val;
% }
%```

%To be able to play with gradient-based algorithm, we must have a gradient ready :
%```freefem
% varf dBFV(v,w) = int2d(Th)(dx(v)*dx(w)+dy(v)*dy(w));
% matrix dA = dBFV(Vh,Vh);
% func real[int] dJ(real[int] &X)
% {
%	 Vh u1,u2;
%	 u1[] = X(0:Vh.ndof-1);
%	 u2[] = X(Vh.ndof:2*Vh.ndof-1);
%	 real[int] grad1 = dA*u1[], grad2 = dA*u2[];
%	 grad1 -= b1;
%	 grad2 -= b2;
%	 real[int] Grad(X.n);
%	 Grad(0:Vh.ndof-1) = grad1;
%	 Grad(Vh.ndof:2*Vh.ndof-1) = grad2;
%	 return Grad;
% }
%```

%As we use P1 Lagrange finite elements, the $u_{1}\leq u_{2} $ constraint can be expressed directly upon the degree of freedom. In NLOpt, one constraint must have the form $f(X)\leq 0$.
%The function is vectorial in order to handle multiple constraints :
%```freefem
% func real[int] IneqC(real[int] &X)
% {
%	 real[int] constraints(Vh.ndof);
%	 for(int i=0;i<Vh.ndof;++i) constraints[i] = X[i] - X[i+Vh.ndof];
%	 return constraints; //u1i - u2i <= 0
% }
%```

%One also has to define a gradient for each constraint if a gradient search is to be used. Constraint gradient must be a function which returns a matrix. The $i$-th line of this matrix is the gradient of
%the $i$-th constraint with respect to the variables of the search space.
%```freefem
% func real[int,int] dIneqC(real[int] &X)
% {
%	 real[int,int] dconst(Vh.ndof,2*Vh.ndof);//so sparse... shame!
%	 dconst=0;
%	 for(int i=0;i<Vh.ndof;++i)
%	 {
%		 dconst(i,i) = 1.;
%		 dconst(i,i+Vh.ndof) = -1.;
%	 }
%	 return dconst;
% }
%```

%As mentioned earlier, we also define a set of equality constraints and their derivatives in order to handle the boundary conditions. One could also use the  `:::freefem lb` and  `:::freefem ub` parameters to impose them, this also would avoid the use of large sparse matrix returned by the gradient function (it can be found in the script file, where we set  `:::freefem lb` to $g-\epsilon$ and  `:::freefem ub` to $g+\epsilon$).
%Here again, we can work directly on the degree of freedom because of the P1 finite element. The script could not work with other finite elements without modifications.
%```freefem
% real[int] BordIndex(Th.nbe); //Indexes of border d.f.
%
% {
%	 int k=0; //filling the index array
%	 for(int i=0;i<Bord.n;++i) if(Bord[][i]) {BordIndex[k]=i; ++k;}
% }
%
% func real[int] BC(real[int] &X)
% {
%	 real[int] bc(2*Th.nbe);
%	 for(int i=0;i<Th.nbe;++i)
%	 {
%		 int I = BordIndex[i];
%		 bc[i] = X[I] - gh1[][I];
%		 bc[i+Th.nbe] = X[I+Th.nv] - gh2[][I];
%	 }
%	 return bc;
% }
% func real[int,int] dBC(real[int] &X) //Gradient of boundary conditions
% {
%	 real[int,int] dbc(2*Th.nbe,2*Th.nv);//even sparser...:(
%	 dbc=0.;
%	 for(int i=0;i<Th.nbe;++i)
%	 {
%		 int I=BordIndex[i];
%		 dbc(i,I) = 1.;
%		 dbc(i+Th.nbe,I+Th.nv) = 1.;
%	 }
%	 return dbc;
% }
%```

%Now we fill an initial vector and start the algorithm :
%```freefem
% real[int] start(2*Vh.ndof);
% start(0:Vh.ndof-1) = 0.; //starting with u1=u2=0 is bad
% start(Vh.ndof:2*Vh.ndof-1) = 0.01; //a small difference is needed
%
% real mini = nloptAUGLAG(J, start, grad=dJ,
% IConst=IneqC, gradIConst=dIneqC,
% EConst=BC, gradEConst=dBC,
% subOpt="LBFGS",
% stopMaxFEval=10000, stopAbsFTol=starttol
% );
%```

%Here we wrap a BFGS algorithm (which only works with unconstraint problems) in the augmented lagrangian method for both equality and inequality constraints. Another legit choice would be to use
%SLSQP, the only gradient-based search of the package supporting these two kinds of constraints (note that the  `:::freefem subOpt` parameter is gone ) :
%```freefem
% real mini = nloptSLSQP(J, start, grad=dJ,
% IConst=IneqC, gradIConst=dIneqC,
% EConst=BC, gradEConst=dBC,
% stopMaxFEval=10000, stopAbsFTol=starttol
% );
%```

%One can also try to pass derivative free algorithms to AUGLAG and discard the gradients to see what happens, or directly one which handles the two kinds of constraints. Just check
%the tag in the table. The true script also has a loop in which the algorithm is restarted after mesh adaptation, nothing was said about that to avoid complicating an already relatively long
%script comment.

## Optimization with MPI

The only quick way to use the previously presented algorithms on a parallel architecture lies in parallelizing the used cost function (which is in most real life case, the expensive part of the algorithm). Somehow, we provide a parallel version of the CMA-ES algorithm. The parallelization principle is the trivial one of evolving/genetic algorithms : at each iteration the cost function has to be evaluated $N$ times without any dependence at all, these $N$ calculus are then equally distributed to each processes. Calling the MPI version of CMA-ES is nearly the same as calling its sequential version (a complete example of use can be found in the `:::freefem cmaes-mpi-VarIneq.edp` file):

```freefem
load "mpi-cmaes"
... // Define J, u and all here
real min = cmaesMPI(J,u,stopTolFun=1e-6,stopMaxIter=3000);
cout << "minimal value is " << min << " for u = " << u << endl;
```

If the population size is not changed using the `:::freefem popsize` parameter, it will use the heuristic value slightly changed to be equal to the closest greater multiple of the size of the communicator used by
the optimizer. The FreeFem `:::freefem mpicommworld` is used by default. The user can specify his own MPI communicator with the named parameter `:::freefem comm=`, see the MPI section of this manual for more informations about communicators in FreeFem++.
