# Classification of the equations

**Summary :** _It is usually not easy to determine the type of a system. Yet the approximations and algorithms suited to the problem depend on its type:_

* _Finite Elements compatible (LBB conditions) for elliptic systems_
* _Finite difference on the parabolic variable and a time loop on each elliptic subsystem of parabolic systems; better stability diagrams when the schemes are implicit in time._
* _Upwinding, Petrov-Galerkin, Characteristics-Galerkin, Discontinuous-Galerkin, Finite Volumes for hyperbolic systems plus, possibly, a time loop._

_When the system changes type, then expect difficulties (like shock discontinuities) !_

**Elliptic, parabolic and hyperbolic equations**

A partial differential equation (PDE) is a relation between a function
of several variables and its derivatives.

 $$
 F(\varphi(x),{\p\varphi\over\p
 x_1}(x),\cdots,{\p\varphi\over\p
 x_d}(x),{\p^2\varphi\over\p
 x^2_1}(x),\cdots,{\p^m\varphi\over\p x^m_d}(x)) =
 0\quad\forall x\in\Omega\subset \R^d.
 $$

The range of $x$ over which the equation is taken, here $\Omega$, is called the _domain_ of the PDE.
The highest derivation index, here $m$, is called the _order_. If $F$ and $\varphi$ are vector valued functions, then the PDE is actually a _system_ of PDEs.

Unless indicated otherwise, here by convention _one_ PDE corresponds to one scalar valued $F$ and $\varphi$.
If $F$ is linear with respect to its arguments, then the PDE is said to be _linear_.

The general form of a second order, linear scalar PDE is ${\p^2\varphi\over\p x_i\p x_j}$ and $A:B$ means $\sum^d_{i,j=1} a_{ij} b_{ij}.$

 $$\alpha\varphi + a\cdot\nabla\varphi + B :\nabla(\nabla\varphi) =
 f{\quad\hbox{ in }\quad}\Omega\subset \R^d,
 $$

where $f(x),\alpha(x)\in \R, a(x)\in \R^d, B(x)\in \R^{d\times d}$
are the PDE _coefficients_.
If the coefficients are independent of $x$, the PDE is said to have _constant coefficients_.

To a PDE we associate a quadratic form, by replacing $\varphi$ by $1$, $\p\varphi/\p x_i$ by $z_i$ and $\p^2\varphi/\p x_i\p x_j$ by $z_i z_j$, where $z$ is a vector in $\R^d$ :

$$\alpha + a\cdot z + z^T Bz = f.
$$

If it is the equation of an ellipse (ellipsoid if $d \geq 2$),
the PDE is said to be _elliptic_;
if it is the equation of a parabola or a hyperbola, the PDE is said to
be _parabolic_ or _hyperbolic_. If $A \equiv 0$, the degree is
no longer 2 but 1, and for reasons that will appear more clearly
later, the PDE is still said to be hyperbolic.

These concepts can be generalized to systems, by studying whether or
not the polynomial system $P(z)$ associated with the PDE system has
branches at infinity (ellipsoids have no branches at infinity,
paraboloids have one, and hyperboloids have several).

If the PDE is not linear, it is said to be _non linear_.
Those are said to be locally elliptic, parabolic, or hyperbolic
according to the type of the linearized equation.

For example, for the non linear equation

 $${\p^2\varphi\over\p t^2} - {\p\varphi\over\p
 x}{\p^2\varphi\over\p x^2} = 1,
 $$

we have $d=2, x_1 = t, x_2 = x$ and its linearized form is:

 $${\p^2 u\over\p t^2} - {\p u\over\p x}
 {\p^2\varphi\over\p x^2} - {\p\varphi\over\p
 x}{\p^2 u\over\p x^2} = 0,
 $$

which for the unknown $u$ is locally elliptic if ${\p\varphi\over\p x} < 0$ and locally hyperbolic if ${\p\varphi\over\p x} > 0$.

**Examples**

Laplace's equation is elliptic:

 $$\Delta\varphi \equiv {\p^2\varphi\over\p x^2_1} +
 {\p^2\varphi\over\p x^2_2} + \cdots +
 {\p^2\varphi\over\p x^2_d} = f, \ \ \ \forall x
 \in \Omega\subset \R^d.
 $$

The _heat_ equation is parabolic in $Q = \Omega\times]0,T[\subset \R^{d+1}$~:

 $${\p\varphi\over\p t} - \mu\Delta\varphi = f\quad\forall
 x\in\Omega\subset  \R^d, \quad\forall t\in]0,T[.
 $$

If $\mu >0$,  the _wave_ equation is hyperbolic:

 $${\p^2\varphi\over\p t^2} - \mu\Delta\varphi =
 f{\quad\hbox{~in~}\quad}  Q.
 $$

The _convection diffusion_ equation is parabolic if $\mu \neq 0$
and hyperbolic otherwise:

 $${\p\varphi\over\p t} + a\nabla\varphi -
 \mu\Delta\varphi = f.
 $$

The _biharmonic_ equation is elliptic:

 $$\Delta(\Delta\varphi) = f{\quad\hbox{~in~}\quad}\Omega.
 $$

**Boundary conditions**

A relation between a function and its derivatives is not sufficient to define the function. Additional information on the boundary $\Gamma=\p\Omega$ of
$\Omega$, or on part of $\Gamma$ is necessary.
Such information is called a _boundary condition_.
For example,

 $$\varphi(x) \ \hbox{given},\ \forall x\in \Gamma,
 $$

is called a _Dirichlet boundary condition_.
The _Neumann_ condition is

$\codered$ **check formula below with ~~~**
 $${\p\varphi\over\p n}(x) \ \hbox{given on }\
 \Gamma \hbox{~~~(or~~} n\cdot B\nabla\varphi,\hbox{given on }\
 \Gamma\hbox{ for a general second order PDE)}
$$

where $n$ is the normal at $x\in\Gamma$ directed towards the exterior of $\Omega$ (by definition ${\p\varphi\over\p n}=\nabla\varphi\cdot n$).

Another classical condition, called a _Robin_ (or _Fourier_) condition is written as:
 
$$\varphi(x) + \beta(x) {\p\varphi\over\p n}(x) \
 \hbox{given on}\
 \Gamma.
 $$

Finding a set of boundary conditions  that defines a unique
$\varphi$ is a difficult art.
In general, an elliptic equation is well posed (_i.e._ $\varphi$
is unique) with one Dirichlet, Neumann or Robin conditions on the whole boundary.

Thus, Laplace's equations  is well posed with a Dirichlet or Neumann condition but also with :

$$\varphi \ \hbox{given on}\ \Gamma_1,\quad {\p\varphi\over
 \p n} \
 \hbox{given on}\ \Gamma_2, \quad \Gamma_1\cup\Gamma_2 =
 \Gamma,\quad{\dot{\Gamma_1}\cap\dot{\Gamma_2}} =
 \emptyset.
$$

Parabolic and hyperbolic equations rarely require boundary conditions
on all of  $\Gamma\times]0,T[$. For instance, the heat equation
is well posed with :

$$\varphi \ \hbox{given at}\ t=0 \ \hbox{and Dirichlet or Neumann or mixed conditions on}\
 \p\Omega.
$$

Here $t$ is time so the first condition is called an initial condition. The whole set of conditions are also called Cauchy conditions.

The wave equation  is well posed with :

$$\varphi \ \hbox{and}\ {\p\varphi\over\p t} \
 \hbox{given at}\ t=0
 \ \hbox{and Dirichlet or Neumann or mixed conditions on}\
 \p\Omega.
$$
