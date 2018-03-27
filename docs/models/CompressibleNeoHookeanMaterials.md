_Written by Alex Sadovsky_

$$
\def\bR{{\bf R}}
\def\bP{{\bf P}}
\def\bZ{{\bf Z}}
\def\bC{{\bf C}}
\def\VS{\bR^2}
\def\SVS{\underline V}
\def\SO{{\bf SO}}
\def\Sym{{\bf Sym}}
\def\qi{{\bf i}}
\def\qj{{\bf j}}
\def\qk{{\bf k}}
\def\ec{\hat{\bf e}}
\def\xc{\hat{\bf x}}
\def\bdr{\partial}
\def\PD{\partial_}
\def\strain{\underline \epsilon}
\def\stress{\underline \sigma}
\def\strainrate{\underline \epsilon^.}
\def\stressrate{\underline \sigma^.}
\def\stiff{\; \underline{\underline C}\;}
\def\comply{\underline{\underline \kappa}\;}
\def\Id{{\bf I}}
\def\Div{\nabla \cdot}
\def\Grad{\vec{\nabla}}
\def\rot{\nabla \times}
\def\lap{\triangle}
\def\tr{{\bf tr}\;}
\def\udH{\underline H}
\def\refX{\vec X}
\def\Jac{\overline{J}}
\def\spatx{\vec x}
\def\ani{\overline a}
\def\mat{\left[\begin{array}}
\def\tam{\end{array}\right]}
\def\arr{\left.\begin{array}}
\def\rra{\end{array}\right\}}
\def\arl{\left\{\begin{array}}
\def\lra{\end{array}\right.}
\def\ar{\begin{array}}
\def\ra{\end{array}}
\def\const{\mbox{ const.}}
\def\eps{\; \epsilon}
\def\sig{\; \sigma}
\def\th{\theta}
\def\sgn{\mbox{sgn}}
\def\qed{\; Q.E.D.\\}
\def\ranqe{\end{eqnarray}}
\def\ol{\overline}
\def\ul{\underline}
\def\bB{{\bf B}}
\def\bC{{\bf C}}
\def\bD{{\bf D}}
\def\bE{{\bf E}}
\def\bF{{\bf F}}
\def\bK{{\bf K}}
\def\bP{{\bf P}}
\def\bS{{\bf S}}
\def\bT{{\bf T}}
\def\bsig{{\bf \sigma}}
$$

## Notation

In what follows, the symbols $\vec{u}, \bF, \bB, \bC, \stress$ denote, respectively, the displacement field, the deformation gradient, the left Cauchy-Green strain tensor $\bB = \bF \bF^T$, the right Cauchy-Green strain tensor $\bC =\bF^T \bF$, and the Cauchy stress tensor.

We also introduce the symbols $I_1 := \tr \bC$ and $J := \det\bF$. Use will be made of the identity
\begin{equation}
{\PD{}J \over \PD{}\bC} = J \bC^{-1}
\end{equation}

The symbol $\Id$ denotes the identity tensor. The symbol $\Omega_{0}$ denotes the reference configuration of the body to be deformed. The unit volume in the reference (resp., deformed) configuration is denoted $dV$ (resp., $dV_{0}$); these two are related by
$$
dV = J dV_{0},
$$
which allows an integral over $\Omega$ involving the Cauchy stress $\bT$ to be rewritten as an integral of the Kirchhoff stress $\kappa = J \bT$ over $\Omega_{0}$.

## Recommended References

For an exposition of nonlinear elasticity and of the underlying linear- and tensor algebra, see [OGDEN1984](#OGDEN1984). For an advanced mathematical analysis of the Finite Element Method, see [RAVIART1998](#RAVIART1998).

## A Neo-Hookean Compressible Material

_Constitutive Theory and Tangent Stress Measures_

The strain energy density function is given by
\begin{equation}
W = {\mu \over 2}(I_1 - \tr \Id - 2 \ln J)
\end{equation}
(see [HORGAN2004](#HORGAN2004), formula (12)).

The corresponding 2nd Piola-Kirchoff stress tensor is given by
\begin{equation}
\bS_{n} := {\PD{} W \over \PD{}\bE} (\bF_{n})
=
\mu (\Id - \bC^{-1})
\end{equation}
The Kirchhoff stress, then, is
\begin{equation}
\kappa
= \bF \bS \bF^{T}
= \mu (\bB  - \Id)
\end{equation}
The tangent Kirchhoff stress tensor at $\bF_{n}$ acting on
$
\delta \bF_{n+1}
$ is, consequently,
\begin{equation}
{\PD{} \kappa \over \PD{} \bF} (\bF_{n}) \delta \bF_{n+1}
=
\mu
\left[
\bF_{n} (\delta \bF_{n+1})^T
+
\delta \bF_{n+1} (\bF_{n})^T
\right]
\end{equation}

_The Weak Form of the BVP in the Absence of Body (External) Forces_

The $\Omega_0$ we are considering is an elliptical annulus, whose
boundary consists of two concentric ellipses (each allowed to be a
circle as a special case), with the major axes parallel.  Let $P$ denote the dead stress load (traction) on a portion
$\partial \Omega_0^{t}$ (= the inner ellipse) of the boundary
$\partial \Omega_0$.  On the rest of the boundary, we prescribe zero displacement.

The weak formulation of the boundary value
problem is
$$
\arr{lll}
0
& = &
\int_{\Omega_0}
\kappa[\bF]
\:
:
\:
\left\{
(\Grad \otimes \vec{w}) (\bF)^{-1}
\right\}\\
& - & \int_{\PD{} \Omega_0^{t}} P \cdot \hat{N}_0\\
\rra
$$
{\em
For brevity, in the rest of this section we assume $P = 0$.  The provided
FreeFem++ code, however, does not rely on this assumption and allows
for a general value and direction of $P$.}

Given a Newton approximation $\vec{u}_n$ of the displacement field
$\vec{u}$ satisfying the BVP, we seek the correction $\delta \vec{u}_{n+1}$ to
obtain a better approximation
$$
\vec{u}_{n+1} = \vec{u}_{n} + \delta \vec{u}_{n+1}
$$
by solving the weak formulation
\begin{equation}
\arr{lll}
0
& = &
\int_{\Omega_0}
\kappa[\bF_{n} + \delta \bF_{n+1}]
\:
:
\:
\left\{
(\Grad \otimes \vec{w}) (\bF_{n} + \delta \bF_{n+1})^{-1}
\right\}
- \int_{\PD{} \Omega_0} P \cdot \hat{N}_0
\\
& = &
\int_{\Omega_0}
\left\{
\kappa[\bF_{n}] +
{\PD{} \kappa \over \PD{} \bF}[\bF_{n}]
\delta \bF_{n+1}
\right\}
\:
:
\:
\left\{
(\Grad \otimes \vec{w})
(\bF_{n} + \delta \bF_{n+1})^{-1}
\right\}
\\
& = &
\int_{\Omega_0}
\left\{
\kappa[\bF_{n}] +
{\PD{} \kappa \over \PD{} \bF}[\bF_{n}]
\delta \bF_{n+1}
\right\}
\:
:
\:
\left\{
(\Grad \otimes \vec{w}) (\bF_{n}^{-1} + \bF_{n}^{-2} \delta \bF_{n+1})
\right\}
\\
\\
& = &
\int_{\Omega_0}
\kappa[\bF_{n}]
\:
:
\:
\left\{
(\Grad \otimes \vec{w})
\bF_{n}^{-1}
\right\}\\
&-&
\int_{\Omega_0}
\kappa[\bF_{n}]
\:
:
\:
\left\{
(\Grad \otimes \vec{w})
(\bF_{n}^{-2} \delta \bF_{n+1})
\right\}\\
&+&
\int_{\Omega_0}
\left\{
{\PD{} \kappa \over \PD{} \bF}[\bF_{n}]
\delta \bF_{n+1}
\right\}
\:
:
\:
\left\{
(
\Grad \otimes \vec{w})
\bF_{n}^{-1}
\right\}
\\
\rra
\quad
\mbox{for all test functions $\vec{w}$,}
\end{equation}
where we have taken
$$
\delta \bF_{n+1} = \Grad \otimes \delta \vec{u}_{n+1}
$$

{\bf Note:}  Contrary to standard notational use, the symbol $\delta$
here bears no variational context.  By $\delta$ we mean simply an
increment in the sense of Newton's Method.  The role of a variational virtual displacement here
is played by $\vec{w}$.

## An Approach to Implementation in FreeFem++

The associated file is  `:::freefem examples++-tutorial/nl-elast-neo-Hookean.edp`.

Introducing the code-like notation, where a string in $< >$'s is to be
read as one symbol, the individual components of the tensor
\begin{equation}
<TanK>
 :=
{\PD{} \kappa \over \PD{} \bF}[\bF_{n}]
\delta \bF_{n+1}
\end{equation}
will be implemented as the macros $<TanK11>, <TanK12>, \ldots$.

The individual components of the tensor quantities
$$
\bD_{1} :=
\bF_{n} (\delta \bF_{n+1})^T
+
\delta \bF_{n+1} (\bF_{n})^T,
$$
$$
\bD_{2} :=
\bF_{n}^{-T} \delta \bF_{n+1},
$$
$$
\bD_{3} :=
(\Grad \otimes \vec{w})
\bF_{n}^{-2} \delta \bF_{n+1},
$$
and
$$
\bD_{4} :=
(\Grad \otimes \vec{w})
\bF_{n}^{-1},
$$
will be implemented as the macros
\begin{equation}

\arr{l}
<d1Aux11>, <d1Aux12>, \quad \ldots \quad, <d1Aux22>,\\
<d2Aux11>, <d2Aux12>, \quad \ldots \quad, <d2Aux22>\\
<d3Aux11>, <d3Aux12>, \quad \ldots \quad, <d3Aux22>\\
<d4Aux11>, <d4Aux12>, \quad \ldots \quad, <d4Aux22>\\
\rra,
\end{equation}
respectively.

In the above notation, the tangent Kirchhoff stress term becomes
\begin{equation}
{\PD{} \kappa \over \PD{} \bF} (\bF_{n})
\: \delta \bF_{n+1}
=
\mu
\: \bD_{1}
\end{equation}
while the weak BVP formulation acquires the form
\begin{equation}
\arr{lll}
0 & = &
\int_{\Omega_0}
\kappa[\bF_{n}]
\:
:
\:
\bD_{4}
\\
&-&
\int_{\Omega_0}
\kappa[\bF_{n}]
\:
:
\:
\bD_{3}
\\
&+&
\int_{\Omega_0}
\left\{
{\PD{} \kappa \over \PD{} \bF}[\bF_{n}]
\delta \bF_{n+1}
\right\}
\:
:
\:
\bD_{4}
\\
\rra
\quad
\mbox{for all test functions $\vec{w}$}
\end{equation}

## References

<a name="OGDEN1984">[OGDEN1984]</a> OGDEN, Ray W. Non-linear elastic deformations. 1984.

<a name="RAVIART1998">[RAVIART1998]</a> RAVIART, Pierre-Arnaud, THOMAS, Jean-Marie, CIARLET, Philippe G., et al. Introduction à l'analyse numérique des équations aux dérivées partielles. Paris : Dunod, 1998.

<a name="HORGAN2004">[HORGAN2004]</a> HORGAN, Cornelius O. et SACCOMANDI, Giuseppe. Constitutive models for compressible nonlinearly elastic materials with limiting chain extensibility. Journal of Elasticity, 2004, vol. 77, no 2, p. 123-138.
