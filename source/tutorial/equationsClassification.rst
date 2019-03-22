Classification of partial differential equations
================================================

**Summary :**
*It is usually not easy to determine the type of a system.*
*Yet the approximations and algorithms suited to the problem depend on its type:*

-  *Finite Elements compatible (LBB conditions) for elliptic systems*
-  *Finite difference on the parabolic variable and a time loop on each elliptic subsystem of parabolic systems; better stability diagrams when the schemes are implicit in time.*
-  *Upwinding, Petrov-Galerkin, Characteristics-Galerkin, Discontinuous-Galerkin, Finite Volumes for hyperbolic systems plus, possibly, a time loop.*

*When the system changes type, then expect difficulties (like shock discontinuities) !*

**Elliptic, parabolic and hyperbolic equations**

A partial differential equation (PDE) is a relation between a function of several variables and its derivatives.

.. math::
   F\left(\varphi(x),{\partial\varphi\over\partial
   x_1}(x),\cdots,{\partial\varphi\over\partial
   x_d}(x),{\partial^2\varphi\over\partial
   x^2_1}(x),\cdots,{\partial^m\varphi\over\partial x^m_d}(x)\right) =
   0,\quad\forall x\in\Omega\subset \mathbb{R}^d

The range of :math:`x` over which the equation is taken, here :math:`\Omega`, is called the *domain* of the PDE.
The highest derivation index, here :math:`m`, is called the *order*.
If :math:`F` and :math:`\varphi` are vector valued functions, then the PDE is actually a *system* of PDEs.

Unless indicated otherwise, here by convention *one* PDE corresponds to one scalar valued :math:`F` and :math:`\varphi`.
If :math:`F` is linear with respect to its arguments, then the PDE is said to be *linear*.

The general form of a second order, linear scalar PDE is :math:`{\partial^2\varphi\over\partial x_i\partial x_j}` and :math:`A:B` means :math:`\sum^d_{i,j=1} a_{ij} b_{ij}.`

.. math::
   \alpha\varphi + a\cdot\nabla\varphi + B :\nabla(\nabla\varphi) =
   f{\quad\hbox{ in }\quad}\Omega\subset \mathbb{R}^d,

where :math:`f(x),\alpha(x)\in \mathbb{R}`, :math:`a(x)\in \mathbb{R}^d`, :math:`B(x)\in \mathbb{R}^{d\times d}` are the PDE *coefficients*.
If the coefficients are independent of :math:`x`, the PDE is said to have *constant coefficients*.

To a PDE we associate a quadratic form, by replacing :math:`\varphi` by :math:`1`, :math:`\partial\varphi/\partial x_i` by :math:`z_i` and :math:`\partial^2\varphi/\partial x_i\partial x_j` by :math:`z_i z_j`, where :math:`z` is a vector in :math:`\mathbb{R}^d`:

.. math::
   \alpha + A\cdot z + z^T Bz = f.

If it is the equation of an ellipse (ellipsoid if :math:`d \geq 2`), the PDE is said to be *elliptic*; if it is the equation of a parabola or a hyperbola, the PDE is said to be *parabolic* or *hyperbolic*.

If :math:`A \equiv 0`, the degree is no longer 2 but 1, and for reasons that will appear more clearly later, the PDE is still said to be hyperbolic.

These concepts can be generalized to systems, by studying whether or not the polynomial system :math:`P(z)` associated with the PDE system has branches at infinity (ellipsoids have no branches at infinity, paraboloids have one, and hyperboloids have several).

If the PDE is not linear, it is said to be *non-linear*.
Those are said to be locally elliptic, parabolic, or hyperbolic according to the type of the linearized equation.

For example, for the non-linear equation

.. math::
   {\partial^2\varphi\over\partial t^2} - {\partial\varphi\over\partial x}{\partial^2\varphi\over\partial x^2} = 1

we have :math:`d=2`, :math:`x_1 = t`, :math:`x_2 = x` and its linearized form is:

.. math::
   {\partial^2 u\over\partial t^2} - {\partial u\over\partial x}{\partial^2\varphi\over\partial x^2} - {\partial\varphi\over\partial x}{\partial^2 u\over\partial x^2} = 0

which for the unknown :math:`u` is locally elliptic if :math:`{\partial\varphi\over\partial x} < 0` and locally hyperbolic if :math:`{\partial\varphi\over\partial x} > 0`.

.. tip:: Laplace’s equation is elliptic:

   .. math::
      \Delta\varphi \equiv {\partial^2\varphi\over\partial x^2_1}
          + {\partial^2\varphi\over\partial x^2_2}
          + \cdots
          + {\partial^2\varphi\over\partial x^2_d} = f,\ \forall x
          \in \Omega\subset \mathbb{R}^d

.. tip:: The *heat* equation is parabolic in :math:`Q = \Omega\times]0,T[\subset \mathbb{R}^{d+1}`:

   .. math::
      {\partial\varphi\over\partial t} - \mu\Delta\varphi = f
          \ \forall x\in\Omega\subset \mathbb{R}^d, \ \forall t\in]0,T[

.. tip:: If :math:`\mu>0`, the *wave* equation is hyperbolic:

   .. math::
      {\partial^2\varphi\over\partial t^2} - \mu\Delta\varphi
          = f{\ \hbox{ in }\ } Q.

.. tip:: The *convection diffusion* equation is parabolic if :math:`\mu \neq 0` and hyperbolic otherwise:

   .. math::
      {\partial\varphi\over\partial t}
          + a\nabla\varphi
          - \mu\Delta\varphi
          = f

.. tip:: The *biharmonic* equation is elliptic:

   .. math::
      \Delta(\Delta\varphi) = f{\ \hbox{ in }\ }\Omega.

**Boundary conditions**

A relation between a function and its derivatives is not sufficient to define the function.
Additional information on the boundary :math:`\Gamma=\partial\Omega` of :math:`\Omega`, or on part of :math:`\Gamma` is necessary.
Such information is called a *boundary condition*.

For example:

.. math::
   \varphi(x) \ \hbox{given},\ \forall x\in \Gamma,

is called a *Dirichlet boundary condition*. The *Neumann* condition is

.. math::
   {\partial\varphi\over\partial \mathbf{n}}(x) \ \hbox{given on }\
   \Gamma \hbox{ (or } \mathbf{n}\cdot B\nabla\varphi,\hbox{given on }\
   \Gamma\hbox{ for a general second order PDE)}

where :math:`\mathbf{n}` is the normal at :math:`x\in\Gamma` directed towards the exterior of :math:`\Omega` (by definition :math:`{\partial\varphi\over\partial \mathbf{n}}=\nabla\varphi\cdot \mathbf{n}`).

Another classical condition, called a *Robin* (or *Fourier*) condition is written as:

.. math::
   \varphi(x) + \beta(x) {\partial\varphi\over\partial n}(x) \ \hbox{given on}\ \Gamma.

Finding a set of boundary conditions that defines a unique :math:`\varphi` is a difficult art.

In general, an elliptic equation is well posed (*i.e.* :math:`\varphi` is unique) with one Dirichlet, Neumann or Robin condition on the whole boundary.

Thus, Laplace’s equation is well posed with a Dirichlet or Neumann condition but also with :

.. math::
   \varphi \ \hbox{given on}\ \Gamma_1,\ {\partial\varphi\over\partial n} \ \hbox{given on}\ \Gamma_2, \ \Gamma_1\cup\Gamma_2 =\Gamma,\ {\dot{\Gamma_1}\cap\dot{\Gamma_2}} =\emptyset.

Parabolic and hyperbolic equations rarely require boundary conditions on all of :math:`\Gamma\times]0,T[`.
For instance, the heat equation is well posed with :

.. math::
   \varphi \ \hbox{given at}\ t=0 \ \hbox{and Dirichlet or Neumann or mixed conditions on}\
   \partial\Omega.

Here :math:`t` is time so the first condition is called an initial condition.
The whole set of conditions is also called Cauchy condition.

The wave equation is well posed with :

.. math::
   \varphi \ \hbox{and}\ {\partial\varphi\over\partial t} \
   \hbox{given at}\ t=0
   \ \hbox{and Dirichlet or Neumann or mixed conditions on}\
   \partial\Omega.
