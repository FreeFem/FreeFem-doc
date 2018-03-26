FreeFem++ also solves evolution problems such as the heat equation:

\begin{eqnarray}
	\frac{\p u}{\p t}-\mu\Delta u &=& f & \textrm{ in }\Omega\times ]0,T[\label{eqn::heatequation}\\
	u(\mathbf{x},0) &=& u_0(\mathbf{x}) & \textrm{ in }\Omega\\
	\left(\p u/\p n\right)(\mathbf{x},t) &=& 0 & \textrm{ on }\p\Omega\times ]0,T[
\end{eqnarray}

with a positive viscosity coefficient $\mu$ and homogeneous Neumann boundary conditions.

We solve \eqref{eqn::heatequation} by FEM in space and finite differences in time.

We use the definition of the partial derivative of the solution in the time derivative,

$$
\frac{\p u}{\p t}(x,y,t) = \lim_{\tau \to 0}\frac{u(x,y,t)-u(x,y,t-\tau )}{\tau }
$$

which indicates that $u^m(x,y)=u(x,y,m\tau )$ will satisfy approximatively

$$
\frac{\p u}{\p t}(x,y,m\tau )\simeq \frac{u^m(x,y)-u^{m-1}(x,y)}{\tau }
$$

The time discretization of heat equation \eqref{eqn::heatequation} is as follows, $\forall m=0,\cdots,[T/\tau ]$:

\begin{eqnarray}
	\frac{u^{m+1}-u^{m}}{\tau }-\mu\Delta u^{m+1} &=& f^{m+1} & \textrm{ in }\Omega\\
	u^0(\mathbf{x}) &=& u_0(\mathbf{x}) & \textrm{ in }\Omega\\
	\p u^{m+1}/\p n(\mathbf{x}) &=& 0 & \textrm{ on }\p\Omega
\end{eqnarray}

which is so-called _backward Euler method_ for \eqref{eqn::heatequation}.

To obtain the variational formulation, multiply with the test function $v$ both sides of the equation:

\begin{equation*}
\int_{\Omega}\{u^{m+1}v-\tau \Delta u^{m+1}v\}=\int_{\Omega}\{u^m+\tau f^{m+1}\}v
\end{equation*}

By the divergence theorem, we have

\begin{equation*}
	\int_{\Omega}\{u^{m+1}v+\tau\nabla u^{m+1}\cdot \nabla v\}
	-\int_{\p\Omega} \tau \left( \p u^{m+1}/\p n\right) v
	=\int_{\Omega }\{u^mv+\tau f^{m+1}v\}
\end{equation*}

By the boundary condition $\p u^{m+1}/\p n=0$, it follows that

\begin{equation}
	\int_{\Omega} \{u^{m+1}v+\tau \nabla u^{m+1}\cdot \nabla v\}
	-\int_{\Omega }\{u^mv+\tau f^{m+1}v\}
	=0
	\label{eqn::heatequationBWE}
\end{equation}

Using the identity just above, we can calculate the finite element approximation $u_h^m$ of $u^m$ in a step-by-step manner with respect to $t$.

!!!question "Example"
	We now solve the following example with the exact solution $u(x,y,t)=tx^4$, \Omega = ]0,1[^2.

	\begin{eqnarray*}
		\frac{{\p u}}{{\p t}} - \mu \Delta u &=& x^4 - \mu 12tx^2 & \textrm{ in }\Omega\times ]0,3[\\
		u(x,y,0) &=& 0 & \textrm{ on }\Omega\\
		\left. u \right|_{\p\Omega} &=& t*x^4
	\end{eqnarray*}

	<!--- __ --->

	```freefem
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
	```

	In the last statement, the $L^2$-error $\left(\int_{\Omega}\left| u-tx^4\right|^2\right)^{1/2}$ is calculated at $t=m\tau, \tau =0.1$. At $t=0.1$, the error is 0.000213269. The errors increase with $m$ and 0.00628589 at $t=3$.

	The iteration of the backward Euler \eqref{eqn::heatequationBWE} is made by [`:::freefem for` loop](../reference/Loops/#for).

	!!! note
		The stiffness matrix in the loop is used over and over again. FreeFem++ support reuses of stiffness matrix.

### Mathematical Theory on Time Difference Approximations.

In this section, we show the advantage of implicit schemes. Let $V, H$ be separable Hilbert space and $V$ is dense in $H$. Let $a$ be a continuous bilinear form over $V \times V$ with coercivity and symmetry.

Then $\sqrt{a(v,v)}$ become equivalent to the norm $\| v\|$ of $V$.

__Problem Ev$(f,\Omega)$__: For a given $f\in L^2(0,T;V'),\, u^0\in H$

\begin{eqnarray}
\frac{d}{dt}(u(t),v)+a(u(t),v)&=&( f(t),v)\qquad \forall v\in V,\quad a.e. \, t\in [0,T]\\
u(0)&=&u^0\nonumber
\end{eqnarray}

where $V'$ is the dual space of $V$.

Then, there is an unique solution $u\in L^{\infty}(0,T;H)\cap L^2(0,T;V)$.

Let us denote the time step by $\tau>0$, $N_T=[T/\tau]$. For the discretization, we put $u^n = u(n\tau)$ and consider the time difference for each $\theta\in [0,1]$

\begin{eqnarray}
\label{eqn::t-method}
\frac{1}{\tau}\left( u_h^{n+1}-u_h^n,\phi_i\right)
+a\left( u_h^{n+\theta},\phi_i\right)=\langle f^{n+\theta},\phi_i\rangle\\
i=1,\cdots, m,\quad n=0,\cdots, N_T\nonumber\\
u_h^{n+\theta}=\theta u_h^{n+1}+(1-\theta)u_h^n,\quad
f^{n+\theta}=\theta f^{n+1}+(1-\theta)f^n\nonumber
\end{eqnarray}

Formula \eqref{eqn::t-method} is the _forward Euler scheme_ if $\theta=0$, _Crank-Nicolson scheme_ if $\theta=1/2$, the _backward Euler scheme_ if $\theta=1$.

Unknown vectors $u^n=(u_h^1,\cdots,u_h^M)^T$ in

$$
u_h^n(x)=u^n_1\phi_1(x)+\cdots+u^n_m\phi_m(x),\quad u^n_1,\cdots,u^n_m\in \R
$$

are obtained from solving the matrix

\begin{eqnarray}
\label{eqn::Evolution-1}
(M+\theta\tau A)u^{n+1}=\{M-(1-\theta)\tau A\}u^n
+\tau\left\{\theta f^{n+1}+(1-\theta)f^n\right\}\\
M=(m_{ij}),\quad m_{ij}=(\phi_j,\phi_i),\qquad
A=(a_{ij}),\quad a_{ij}=a(\phi_j,\phi_i)\nonumber
\end{eqnarray}

Refer [TABATA1994](#TABATA1994), pp.70--75 for solvability of \eqref{eqn::Evolution-1}. The stability of \eqref{eqn::Evolution-1} is in [TABATA1994](#TABATA1994), Theorem 2.13:

Let $\{\mathcal{T}_h\}_{h\downarrow 0}$ be regular triangulations (see [Regular Triangulation](../documentation/MeshGeneration/#regular-triangulation-htriangle)). Then there is a number $c_0>0$ independent of $h$ such that,

\begin{eqnarray}
|u_h^n|^2\le
\left\{
\begin{array}{lr}
\frac{1}{\delta}\left\{
|u^0_h|^2+\tau \sum_{k=0}^{n-1}\|f^{k+\theta}\|^2_{V_h'}
\right\}&\theta\in [0,1/2)\\
|u^0_h|^2+\tau \sum_{k=0}^{n-1}\|f^{k+\theta}\|^2_{V_h'}&\theta\in [1/2,1]
\end{array}
\right.
\end{eqnarray}

if the following are satisfied:

1. When $\theta\in [0,1/2)$, then we can take a time step $\tau$ in such a way that

	\begin{eqnarray}
	\tau <\frac{2(1-\delta)}{(1-2\theta)c_0^2}h^2
	\end{eqnarray}

	for arbitrary $\delta\in (0,1)$.

2. When $1/2\leq \theta\leq 1$, we can take $\tau$ arbitrary.

!!!question "Example"
	```freefem
	// Parameters
	real tau = 0.1;
	real theta = 0.;

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
	```

	|<a name="Fig19">Fig. 19</a>: $\max_{x\in\Omega}\vert u_h^n(\theta)-u_{ex}(n\tau)\vert\max_{x\in\Omega}\vert u_{ex}(n\tau)\vert$ at $n=0,1,\cdots,29$|
	|:----:|
	|![TimeDifference](images/EvolutionProblems_TimeDifference.png)|

	We can see in <a href="Fig19">Fig. 19</a> that $u_h^n(\theta)$ become unstable at $\theta=0.4$, and figures are omitted in the case $\theta<0.4$.

### Convection

The hyperbolic equation

\begin{equation}
\label{eqn::conv}
\p_t u +\mathbf{\alpha} \cdot \nabla u=f;\ \textrm{ for a vector-valued function }\mathbf{\alpha}
\end{equation}

appears frequently in scientific problems, for example in the Navier-Stokes equations, in the Convection-Diffusion equation, etc.

In the case of 1-dimensional space, we can easily find the general solution $(x,t)\mapsto u(x,t)=u^0(x-\alpha t)$ of the following equation, if $\alpha$ is constant,

\begin{eqnarray}
	\label{eqn::conv0}
	\p_t u +\alpha\p_x u &=& 0\\
	u(x,0) &=& u^0(x),
\end{eqnarray}

because $\p_t u +\alpha\p_x u=-\alpha\dot{u}^0+a\dot{u}^0=0$, where $\dot{u}^0=du^0(x)/dx$.

Even if $\alpha$ is not constant, the construction works on similar principles. One begins with the ordinary differential equation (with the convention that $\alpha$ is prolonged by zero apart from $(0,L)\times (0,T)$):

$$
\dot{X}(\tau )=+\alpha(X(\tau ),\tau ),\ \tau \in (0,t)\quad X(t)=x
$$

In this equation $\tau$ is the variable and $x,t$ are parameters, and we denote the solution by $X_{x,t}(\tau )$. Then it is noticed that $(x,t)\rightarrow v(X(\tau),\tau)$ in $\tau=t$ satisfies the equation

$$
\p _{t}v+\alpha\p _{x}v=\p _{t}X\dot{v}+a\p _{x}X\dot{v}%
=0
$$

<!--- __ -->

and by the definition $\p _{t}X=\dot{X}=+\alpha$ and $\p_{x}X=\p _{x}x$ in $\tau=t$, because if $\tau =t$ we have $X(\tau )=x$.

The general solution of \eqref{eqn::conv0} is thus the value of the boundary condition in $X_{x, t}(0)$, that is to say $u(x,t)=u^{0}(X_{x,t}(0))$ where $X_{x,t}(0)$ is on the $x$ axis, $u(x,t)=u^{0}(X_{x,t}(0))$ if $X_{x,t}(0)$ is on the axis of $t$.

In higher dimension $\Omega \subset R^{d},~d=2,3$, the equation for the convection is written

$$
\p _{t}u+\mathbf{\alpha}\cdot \nabla u=0\hbox{ in }\Omega \times (0,T)
$$

<!--- __ --->

where $\mathbf{a}(x,t)\in R^{d}$.

FreeFem++ implements the Characteristic-Galerkin method for convection operators. Recall that the equation \eqref{eqn::conv} can be discretized as

$$
\frac{Du}{Dt} = f\;\;\textrm{i.e. }\frac{du}{dt}\left( {X(t),t} \right) = f\left(X( t ),t \right)\textrm{ where }\frac{dX}{dt}( t ) = \mathbf{\alpha}( {X(t),t})
$$

where $D$ is the total derivative operator. So a good scheme is one step of backward convection by the method of Characteristics-Galerkin

\begin{eqnarray}
\frac{1}{{\tau }}\left(u^{m + 1}(x) - u^m(X^m(x))\right) = f^m (x)
\label{eqn::Charac}
\end{eqnarray}

where $X^m (x)$ is an approximation of the solution at $t = m\tau $ of the ordinary differential equation

\[
\frac{d\mathbf{X}}{dt}(t) = \mathbf{\alpha}^m(\mathbf{X}(t)),\, \mathbf{X}((m + 1)\tau ) = x.
\]

where $\mathbf{\alpha}^m(x)=(\alpha_1(x,m\tau ),\alpha_2(x,m\tau))$. Because, by Taylor's expansion, we have

\begin{eqnarray}
u^m(\mathbf{X}(m\tau ))&=&
u^m(\mathbf{X}((m+1)\tau )) -
\tau \sum_{i=1}^d \frac{\p u^m}{\p x_i}(\mathbf{X}((m+1)\tau ))
\frac{\p X_i}{\p t}((m+1)\tau )
+o(\tau )\nonumber\\
&=&u^m(x)-\tau \mathbf{\alpha}^m(x)\cdot \nabla u^m(x)+o(\tau )
\label{eqn::conv1}
\end{eqnarray}

where $X_i(t)$ are the i-th component of $\mathbf{X}(t)$, $u^m(x)=u(x,m\tau )$ and we used the chain rule and $x=\mathbf{X}((m+1)\tau )$. From \eqref{eqn::conv1}, it follows that

\begin{eqnarray}
u^m(X^m(x))=u^m(x)-\tau \mathbf{\alpha}^m(x)\cdot \nabla u^m(x)+o(\tau ).
\end{eqnarray}

Also we apply Taylor's expansion for $t\mapsto u^m(x-\mathbf{\alpha}^m(x)t),\, 0\le t\le \tau $, then

$$
u^m(x-\mathbf{\alpha}\tau )=u^m(x)-\tau \mathbf{\alpha}^m(x)\cdot \nabla u^m(x)+o(\tau ).
$$

Putting

`:::freefem convect` $\left( {\mathbf{\alpha},-\tau ,u^m } \right)\approx u^m \left(x - \mathbf{\alpha}^m\tau \right)$

we can get the approximation

$u^m \left( {X^m( x )} \right) \approx$ `:::freefem convect` $\left( {[a_1^m ,a_2^m],-\tau ,u^m } \right)\;\;\textrm{by }X^m \approx x \mapsto x- \tau [a_1^m(x) ,a_2^m(x)]$

A classical convection problem is that of the "rotating bell" (quoted from [LUCQUIN1998](#LUCQUIN1998), p.16).

Let $\Omega$ be the unit disk centered at 0, with its center rotating with speed $\alpha_1 = y,\, \alpha_2 = -x$. We consider the problem \eqref{eqn::conv} with $f=0$ and the initial condition $u(x,0)=u^0(x)$, that is, from \eqref{eqn::Charac}

$u^{m + 1}(x) = u^m(X^m(x))\approx$ `:::freefem convect`$(\mathbf{\alpha},-\tau ,u^m)$

The exact solution is $u(x, t) = u(\mathbf{X}(t))$ where $\mathbf{X}$ equals $x$ rotated around the origin by an angle $\theta = -t$ (rotate in clockwise). So, if $u^0$ in a 3D perspective looks like a bell, then $u$ will have exactly the same shape, but rotated by the same amount. The program consists in solving the equation until $T = 2\pi$, that is for a full revolution and to compare the final solution with the initial one; they should be equal.

 __Example 9.20__ convect.edp

```freefem
border C(t=0, 2*pi) { x=cos(t);  y=sin(t); }; // the unit circle
mesh Th = buildmesh(C(70)); // triangulates the disk
fespace Vh(Th,P1);
Vh u0 = exp(-10*((x-0.3)^2 +(y-0.3)^2)); // give $u^0$

real dt = 0.17,t=0; // time step
Vh a1 = -y, a2 = x; // rotation velocity
Vh u; // $u^{m+1}$
for (int m=0; m<2*pi/dt ; m++) {
    t += dt;
    u=convect([a1,a2],-dt,u0); // $u^{m+1}=u^m(X^m(x))$
    u0=u; // m++
    plot(u,cmm=" t="+t + ", min=" + u[].min + ", max=" +  u[].max,wait=0);
};
```

!!! note

	The scheme `:::freefem convect` is unconditionally stable, then the bell become lower and lower (the maximum of $u^{37}$ is $0.406$ as shown in Fig. \ref{BellLast} 9.21 $\codered$).

|Fig. 9.20: $u^0=e^{-10((x-0.3)^2 +(y-0.3)^2)}$|
|:----:|
|![BellInit](images/BellInit.png)|

|Fig. 9.21: The bell at $t=6.29$|
|:----:|
|![BellLast](images/BellLast.png)|

### 2D Black-Scholes equation for an European Put option

In mathematical finance, an option on two assets is modeled by a Black-Scholes equations in two space
variables, (see for example Wilmott et al\cite{wilmott} $\codered$ or Achdou et al \cite{achdou} $\codered$).

\begin{eqnarray}
 &&\p _t u + \frac{{\left( {\sigma _1 x } \right)^2 }}{2}\frac{{\p ^2 u}}{{\p x^2 }} + \frac{{\left( {\sigma _2 y } \right)^2 }}{2}\frac{{\p ^2 u}}{{\p y^2 }} \\
 &&{\rm{      }} + \rho x y \frac{{\p ^2 u}}{{\p x \p y }} + rS_1 \frac{{\p u}}{{\p x }} + rS_2 \frac{{\p u}}{{\p y }} - rP = 0 \nonumber
\end{eqnarray}

which is to be integrated in $\left( {0,T} \right) \times \R^ +   \times \R^ +$
subject to, in the case of a put

\begin{eqnarray}
u\left( {x , y ,T} \right) = \left( {K - \max \left( {x ,y } \right)} \right)^ +  .
\end{eqnarray}

Boundary conditions for this problem may not be so easy to device. As in the one dimensional case the PDE contains boundary conditions on the axis $x_1 = 0$ and on the axis $x_2 = 0$, namely two one dimensional Black-Scholes equations driven respectively by the data $u\left( {0, + \infty ,T} \right)$ and $u\left( { + \infty ,0,T} \right)$. These will be automatically accounted for because they are embedded in the PDE. So if we do nothing in the variational form (i.e. if we take a Neumann boundary condition at these two axis in the strong form) there will be no disturbance to these. At infinity in one of the variable, as in 1D, it makes sense to impose $u=0$. We take

\begin{eqnarray}
\sigma _1  = 0.3,\;\;\sigma _2  = 0.3,\;\;\rho  = 0.3,\;\;r = 0.05,\;\;K = 40,\;\;T = 0.5
\end{eqnarray}

An implicit Euler scheme is used and a mesh adaptation is done every 10 time steps. To have an unconditionally stable scheme, the first order terms are treated by the Characteristic Galerkin method, which, roughly, approximates

\begin{eqnarray}
\frac{{\p u}}{{\p t}} + a_1 \frac{{\p u}}{{\p x}} + a_2 \frac{{\p u}}{{\p y}} \approx \frac{1}{{\tau }}\left( {u^{n + 1} \left( x \right) - u^n \left( {x - \mathbf{\alpha}\tau } \right)} \right)
\end{eqnarray}

 __Example 9.21__ BlackSchol.edp

```freefem
// file BlackScholes2D.edp
int m=30,L=80,LL=80, j=100;
real sigx=0.3, sigy=0.3, rho=0.3, r=0.05, K=40, dt=0.01;
mesh th=square(m,m,[L*x,LL*y]);
fespace Vh(th,P1);

Vh u=max(K-max(x,y),0.);
Vh xveloc, yveloc, v,uold;

for (int n=0; n*dt <= 1.0; n++)
{
  if(j>20)  { th = adaptmesh(th,u,verbosity=1,abserror=1,nbjacoby=2,
              err=0.001, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
              splitpbedge=1, maxsubdiv=5,rescaling=1) ;
     j=0;
     xveloc = -x*r+x*sigx^2+x*rho*sigx*sigy/2;
     yveloc = -y*r+y*sigy^2+y*rho*sigx*sigy/2;
     u=u;
    };
  uold=u;
  solve eq1(u,v,init=j,solver=LU) = int2d(th)( u*v*(r+1/dt)
        + dx(u)*dx(v)*(x*sigx)^2/2 + dy(u)*dy(v)*(y*sigy)^2/2
        + (dy(u)*dx(v) + dx(u)*dy(v))*rho*sigx*sigy*x*y/2)
        - int2d(th)( v*convect([xveloc,yveloc],dt,w)/dt) + on(2,3,u=0);

  j=j+1;
};
plot(u,wait=1,value=1);
```

Results are shown on Fig. \ref{blackScholesE} 9.21 $\codered$).

|Fig. 9.22: The adapted triangulation|
|:----:|
|![BSth](images/BSth.png)|

|Fig. 9.23: The level line of the European basquet put option|
|:----:|
|![BSval](images/BSval.png)|

## References

<a name="TABATA1994">[TABATA1994]</a> TABATA, M. Numerical solutions of partial differential equations II. Iwanami Applied Math, 1994.

[LUCQUIN1998] PIRONNEAU, O. et LUCQUIN-DESREUX, B. Introduction to scientific computing. Wiley, 1998.
