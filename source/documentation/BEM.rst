.. role:: freefem(code)
  :language: freefem

.. _BEM:

The Boundary Element Method
===========================

Introduction to the Boundary Element Method (BEM)
-------------------------------------------------

Model problem
~~~~~~~~~~~~~

The model problem we consider here is the scattering of an incoming acoustic wave :math:`u_\text{inc}` by an obstacle :math:`\Omega`. Thus, we want to solve the following homogeneous Helmholtz equation written in terms of the scattered field :math:`u`:

.. figure:: images/BEM_figdomainbem.png
    :name: BEMfigdomainbem
    :width: 30%

.. math::
  \left \{
  \begin{aligned}
  - \Delta u - k^2 u &= 0 \;\; &\text{in} \;\; &\mathbb{R}^3 \backslash \Omega \\
  u &= - u_\text{inc}  \;\; &\text{on} \;\; &\Gamma\\
  &\text{+ radiation condition}\hspace{-2.8cm}
  \end{aligned}
  \right .

with the Sommerfeld radiation condition at infinity, which states that there can be only outgoing waves at infinity:

.. math::
  \lim_{|\boldsymbol{x}| \rightarrow \infty} |\boldsymbol{x}| \left( \frac{\partial}{\partial |\boldsymbol{x}|} - \imath k \right) u(\boldsymbol{x}) = 0

and where the total field :math:`u_\text{tot} = u_\text{inc} + u`.

If the wavenumber :math:`k` is **constant** in :math:`\mathbb{R}^3 \backslash \Omega`, the boundary element method can be applied. It consists in reformulating the problem in terms of unknowns on the boundary :math:`\Gamma` of :math:`\Omega`.  

First, let us introduce the *Green kernel* :math:`\mathcal{G}`, which for the helmholtz equation in 3D is

.. math::
  :label: eq_p

  \mathcal{G}(\boldsymbol{x}) = \exp(\imath k |\boldsymbol{x}|) / (4 \pi |\boldsymbol{x}|).

Let us also introduce the *Single Layer Potential* :math:`\text{SL}`, which for any :math:`q \in H^{-1/2}(\Gamma)` is defined as

.. math::
  :label: eq_sl

  \text{SL}(q)(\boldsymbol{x}) = \int_{\Gamma} \mathcal{G}(\boldsymbol{x}-\boldsymbol{y}) q(\boldsymbol{y}) d\sigma(\boldsymbol{y}), \quad \forall \boldsymbol{x} \in \mathbb{R}^3 \backslash \Gamma.

An interesting property of :math:`\text{SL}` is that it produces solutions of the PDE at hand in :math:`\mathbb{R}^3 \backslash \Gamma` which satisfy the necessary conditions at infinity (here the Sommerfeld radiation condition for the Helmholtz equation).

Thus, we now need to find a so-called *ansatz* :math:`p \in H^{-1/2}(\Gamma)` such that :math:`\forall \boldsymbol{x} \in \mathbb{R}^3 \backslash \Omega`

.. math::
  :label: eq_pv

  u(\boldsymbol{x}) = \text{SL}(p)(\boldsymbol{x}) = \int_{\Gamma} \mathcal{G}(\boldsymbol{x}-\boldsymbol{y}) p(\boldsymbol{y}) d\sigma(\boldsymbol{y}).

In order to find :math:`p`, we define a variational problem by multiplying :eq:`eq_pv` by a test function `q` and integrating over :math:`\Gamma`:

.. math::
  \int_{\Gamma} u(\boldsymbol{x}) q(\boldsymbol{x}) d\sigma(\boldsymbol{x}) =
  \int_{\Gamma \times \Gamma} \frac{\exp(\imath k |\boldsymbol{x}-\boldsymbol{y}|)}{4 \pi |\boldsymbol{x}-\boldsymbol{y}|} p(\boldsymbol{y}) q(\boldsymbol{x}) d\sigma(\boldsymbol{x,y}) \quad \forall q : \Gamma \rightarrow \mathbb{C}.

Using the Dirichlet boundary condition :math:`u = - u_\text{inc}` on :math:`\Gamma`, we end up with the following variational problem to solve: find :math:`p : \Gamma \rightarrow \mathbb{C}` such that

.. math::
  :label: eq_bem

  \int_{\Gamma \times \Gamma} \frac{\exp(\imath k |\boldsymbol{x}-\boldsymbol{y}|)}{4 \pi |\boldsymbol{x}-\boldsymbol{y}|} p(\boldsymbol{y}) q(\boldsymbol{x}) d\sigma(\boldsymbol{x,y}) = - \int_{\Gamma} u_\text{inc}(\boldsymbol{x}) q(\boldsymbol{x}) d\sigma(\boldsymbol{x}) \quad \forall q : \Gamma \rightarrow \mathbb{C}.

Note that knowing :math:`p` on :math:`\Gamma`, we can indeed compute :math:`u` anywhere using the *potential* formulation :eq:`eq_pv`. Thus, we essentially gained one space dimension, as we only have to solve for :math:`p : \Gamma \rightarrow \mathbb{C}` in :eq:`eq_bem`.

Of course, this inherent benefit of the boundary element method comes with a drawback: after discretization of :eq:`eq_bem`, for example with piecewise linear continuous (P1) functions on :math:`\Gamma`, we end up with a linear system whose matrix is **full**: because :math:`\mathcal{G}(\boldsymbol{x}-\boldsymbol{y})` never vanishes, every interaction coefficient is nonzero. Thus, the matrix :math:`A` of the linear system can be very costly to store (:math:`n^2` coefficients) and invert (factorization in :math:`\mathcal{O}(n^3)`) (:math:`n` is the size of the linear system).  
Moreover, compared to the finite element method, the matrix coefficients are much more expensive to compute because of the double integral and the evaluation of the Green function :math:`\mathcal{G}`. Plus, the choice of the quadrature formulas has to be made with extra care because of the singularity of :math:`\mathcal{G}`.