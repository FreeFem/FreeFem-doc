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

Let us also introduce the *Single Layer Potential* :math:`\operatorname{SL}`, which for any :math:`q \in H^{-1/2}(\Gamma)` is defined as

.. math::
  :label: eq_sl

  \operatorname{SL}(q)(\boldsymbol{x}) = \int_{\Gamma} \mathcal{G}(\boldsymbol{x}-\boldsymbol{y}) q(\boldsymbol{y}) d\sigma(\boldsymbol{y}), \quad \forall \boldsymbol{x} \in \mathbb{R}^3 \backslash \Gamma.

An interesting property of :math:`\text{SL}` is that it produces solutions of the PDE at hand in :math:`\mathbb{R}^3 \backslash \Gamma` which satisfy the necessary conditions at infinity (here the Helmholtz equation and the Sommerfeld radiation condition).

Thus, we now need to find a so-called *ansatz* :math:`p \in H^{-1/2}(\Gamma)` such that :math:`\forall \boldsymbol{x} \in \mathbb{R}^3 \backslash \Omega`

.. math::
  :label: eq_pv

  u(\boldsymbol{x}) = \operatorname{SL}(p)(\boldsymbol{x}) = \int_{\Gamma} \mathcal{G}(\boldsymbol{x}-\boldsymbol{y}) p(\boldsymbol{y}) d\sigma(\boldsymbol{y}).

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

Boundary Integral Operators
-------------------------------------------------

In order to solve our model Dirichlet problem, we have used the **Single Layer Potential** :math:`\operatorname{SL}`:

.. math::
  q \mapsto \operatorname{SL}(q)(\boldsymbol{x}) = \int_{\Gamma} \mathcal{G}(\boldsymbol{x}-\boldsymbol{y}) q(\boldsymbol{y}) d\sigma(\boldsymbol{y}).

Depending on the choice of the boundary integral formulation or boundary condition, the **Double Layer Potential** :math:`\operatorname{DL}` can also be of use:

.. math::
  q \mapsto \operatorname{DL}(q)(\boldsymbol{x}) = \int_{\Gamma} \frac{\partial}{\partial \boldsymbol{n} (\boldsymbol{y})} \mathcal{G}(\boldsymbol{x}-\boldsymbol{y}) q(\boldsymbol{y}) d\sigma(\boldsymbol{y}).

Similarly, we have used the **Single Layer Operator** :math:`\mathcal{SL}` in our variational problem

.. math::
  p, q \mapsto \mathcal{SL}(p,q) = \int_{\Gamma \times \Gamma} p(\boldsymbol{x}) q(\boldsymbol{y}) \mathcal{G}(\boldsymbol{x - y}) d \sigma(\boldsymbol{x,y}).

There are three other building blocks that can be of use in the boundary element method, and depending on the problem and the choice of the formulation a boundary integral method makes use of one or a combination of these building blocks:

the **Double Layer Operator** :math:`\mathcal{DL}`:

.. math::
  p, q \mapsto \mathcal{DL}(p,q) = \int_{\Gamma \times \Gamma} p(\boldsymbol{x}) q(\boldsymbol{y}) \frac{\partial}{\partial \boldsymbol{n} (\boldsymbol{y})} \mathcal{G}(\boldsymbol{x - y}) d \sigma(\boldsymbol{x,y})

the **Transpose Double Layer Operator** :math:`\mathcal{TDL}`:

.. math::
  p, q \mapsto \mathcal{TDL}(p,q) = \int_{\Gamma \times \Gamma} p(\boldsymbol{x}) q(\boldsymbol{y}) \frac{\partial}{\partial \boldsymbol{n} (\boldsymbol{x})} \mathcal{G}(\boldsymbol{x - y}) d \sigma(\boldsymbol{x,y})

the **Hypersingular Operator** :math:`\mathcal{HS}`:

.. math::
  p, q \mapsto \mathcal{HS}(p,q) = \int_{\Gamma \times \Gamma} p(\boldsymbol{x}) q(\boldsymbol{y})  \frac{\partial}{\partial \boldsymbol{n} (\boldsymbol{x})} \frac{\partial}{\partial \boldsymbol{n} (\boldsymbol{y})} \mathcal{G}(\boldsymbol{x - y}) d \sigma(\boldsymbol{x,y})

the BEMTool library
-------------------------------------------------

In order to compute the coefficients of the BEM matrix, **FreeFEM** is interfaced with the boundary element library `BEMTool`_. **BEMTool** is a general purpose header-only C++ library written by Xavier Claeys, which handles

- BEM Potentials and Operators for Laplace, Yukawa, Helmholtz and Maxwell equations
- both in 2D and in 3D
- 1D, 2D and 3D triangulations
- :math:`\mathbb{P}_k`-Lagrange for :math:`k = 0,1,2` and surface :math:`\mathbb{RT}_0`

Although **BEMTool** can compute the BEM matrix coefficients by accurately and efficiently evaluating the boundary integral operator, it is very costly and often prohibitive to compute and store all :math:`n^2` coefficients of the matrix. Thus, we have to rely on a *matrix compression* technique. To do so, **FreeFEM** relies on the **Hierarchical Matrix**, or **H-Matrix** format.

.. _BEMTool: https://github.com/xclaeys/BemTool