$$
\newcommand{\vecttwo}{\left|\begin{array}{c}#1\\#2\end{array}\right.}
$$

# Finite Elements

As stated in [tutorials](../tutorial), FEM approximates all functions $w$ as
$$
w(x,y)\simeq w_0\phi_0(x,y)+w_1\phi_1(x,y)+\cdots+w_{M-1}\phi_{M-1}(x,y)
$$
with finite element basis functions $\phi_k(x,y)$ and numbers $w_k$ ($k=0,\cdots,M-1$). The functions $\phi_k(x,y)$ are constructed from the triangle $T_{i_k}$, and called _shape functions_.

In FreeFem++ the finite element space
$$
V_h=\left\{w\left|\; w_0\phi_0+w_1\phi_1+\cdots+w_{M-1}\phi_{M-1},\, w_i\in \R\right.\right\}
$$
 is easily created by :

```freefem
fespace IDspace(IDmesh,<IDFE>);
```

or with $\ell$ pairs of periodic boundary condition in 2d :

```freefem
fespace IDspace(IDmesh,<IDFE>,
	periodic=[[la$_1$,sa$_1$],[lb$_1$,sb$_1$],
			...
			[la$_k$,sa$_k$],[lb$_k$,sb$_\ell$]]);
```

and in 3d :

```freefem
fespace IDspace(IDmesh,<IDFE>,
	periodic=[[la$_1$,sa$_1$,ta$_1$],[lb$_1$,sb$_1$,tb$_1$],
			...
			[la$_k$,sa$_k$,ta$_k$],[lb$_k$,sb$_\ell$,tb$_\ell$]]);
```

where `:::freefem IDspace` is the name of the space (e.g. `:::freefem Vh`), `:::freefem IDmesh` is the name of the associated mesh and `:::freefem <IDFE>` is an identifier of finite element type.

In 2D we have a pair of periodic boundary condition, if $[la_i, sa_i]$, $[lb_i, sb_i]$ is a pair of `:::freefem int`, and the 2 labels $la_i$ and $lb_i$ refer to 2 pieces of boundary to be in equivalence. If $[la_i, sa_i]$, $[lb_i, sb_i]$ is a pair of `:::freefem real`, then $sa_i$ and $sb_i$ give two common abscissa on the two boundary curve, and two points are identified as one if the two abscissa are equal.

In 2D, we have a pair of periodic boundary condition, if $[la_i, sa_i, ta_i]$, $[lb_i, sb_i, tb_i]$ is a pair of `:::freefem int`, the 2 labels $la_i$ and $lb_i$ define the 2 piece of boundary to be in equivalence. If $[la_i, sa_i, ta_i]$, $[lb_i, sb_i, tb_i]$ is a pair of `:::freefem real`, then $sa_i$, $ta_i$ and $sb_i$, $tb_i$ give two common parameters on the two boundary surface, and two points are identified as one if the two parameters are equal.

!!! note
	The 2D mesh of the two identified borders must be the same, so to be sure, use the parameter `:::freefem fixedborder=true` in `:::freefem buildmesh` command (see \ref{buildmesh fixedborder} $\codered$) like in example `:::freefem periodic2bis.edp` (see \ref{exm:periodic4bis} $\codered$).


As of today, the known types of finite element are:

* `:::freefem [P0, P03d]` piecewise constant discontinuous finite element (2d, 3d), the degrees of freedom are the barycenter element value.
	\begin{equation*}
		\label{eq:P0}
		P0_{h} = \left\{ v \in L^2(\Omega) \left|\; \textrm{for all }K \in \mathcal{T}_{h}\;\;\textrm{there is }\alpha_{K}\in \R : \;\; v_{|K} = \alpha_{K } \right.\right\}
	\end{equation*}

* `:::freefem [P1, P13d]` piecewise linear continuous finite element (2d, 3d), the degrees of freedom are the vertices values.
	\begin{equation*}
		\label{eq:P1}
		P1_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h},\ v_{|K} \in P_{1} \right.\right\}
	\end{equation*}

* `:::freefem [P1dc]` piecewise linear discontinuous finite element
	\begin{equation*}
		\label{eq:P1dc}
		P1dc_{h} = \left\{ v \in L^{2}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{1} \right.\right\}
	\end{equation*}

 	!!! warning
		Due to interpolation problem, the degree of freedom is not the vertices but three vectices move inside with $T(X)= G + .99 (X-G)$ where $G$ is the barycenter.

* `:::freefem [P1b, P1b3d]` piecewise linear continuous finite element plus bubble (2d, 3d)

	**The 2d case:**
	\begin{equation*}
		\label{eq:P1b}
		P1b_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{1} \oplus \mathrm{Span}\{ \lambda^{K}_{0} \lambda^{K}_{1} \lambda^{K}_{2} \} \right.\right\}
	\end{equation*}

	**The 3d case:**
	\begin{equation*}
		\label{eq:P1b-3d}
		P1b_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{1} \oplus \mathrm{Span}\{ \lambda^{K}_{0} \lambda^{K}_{1} \lambda^{K}_{2} \lambda^{K}_{3} \} \right.\right\}
	\end{equation*}
	where $\lambda^{K}_{i}, i=0,..,d$ are the $d+1$ barycentric coordinate functions of the element $K$ (triangle or tetrahedron).

* `:::freefem P1bl,P1bl3d` piecewise linear continuous finite element plus linear bubble (2d, 3d). The bulle is build by spliting the $K$ a barycenter in $d+1$ sub element. (need `:::freefem load "Element_P1bl"`)

* `:::freefem [P2, P23d]` piecewise $P_{2}$ continuous finite element (2d, 3d)
	\begin{equation*}
		P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{2} \right.\right\}
	\end{equation*}
	where $P_{2}$ is the set of polynomials of $\R^{2}$ of degrees $\le 2$.

* `:::freefem [P2b]` piecewise $P_{2} $ continuous finite element plus bubble
	\begin{equation*}
		P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{2} \oplus \mathrm{Span}\{ \lambda^{K}_{0} \lambda^{K}_{1} \lambda^{K}_{2} \} \right.\right\}
	\end{equation*}

* `:::freefem [P2dc]` piecewise $P_{2}$ discontinuous finite element
	\begin{equation*}
		P2dc_{h} = \left\{ v \in L^{2}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{2} \right.\right\}
	\end{equation*}

	!!! warning
		Due to interpolation problem, the degree of freedom is not the six P2 nodes but six nodes move	inside with $T(X)= G + .99 (X-G) $ where $G$ is the barycenter.

 * `:::freefem [P3]` piecewise $P_{3}$ continuous finite element (2d) (need `:::freefem load "Element_P3"`)
	\begin{equation*}
		P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{3} \right.\right\}
	\end{equation*}
	where $P_{3}$ is the set of polynomials of $\R^{2}$ of degrees $\le 3$.

 * `:::freefem [P3dc]` piecewise $P_{3}$ discontinuous finite element (2d) (need `:::freefem load "Element_P3dc"`)
	\begin{equation*}
		P2_{h} = \left\{ v \in L^2(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{3} \right.\right\}
	\end{equation*}
	where $P_{3}$ is the set of polynomials of $\R^{2}$ of degrees $\le 3$.

 * `:::freefem [P4]` piecewise $P_{4}$ continuous finite element (2d) (need `:::freefem load "Element_P4"`)
	\begin{equation*}
		P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h},\ v_{|K} \in P_{4} \right.\right\}
	\end{equation*}
	where $P_{4}$ is the set of polynomials of $\R^{2}$ of degrees $\le 4$.

 * `:::freefem [P4dc]` piecewise $P_{4}$ discontinuous finite element (2d) (need `:::freefem load "Element_P4dc"`)
	\begin{equation*}
		P2_{h} = \left\{ v \in L^2(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{3} \right.\right\}
	\end{equation*}
	where $P_{4}$ is the set of polynomials of $\R^{2}$ of degrees $\le 3$.

* `:::freefem [P0Edge]` piecewise $P_{0}$ discontinuous finite element (2d) contant on each edge of the mesh.

* `:::freefem [P1Edge]` piecewise $P_{1}$ discontinuous finite element (2d) (need `:::freefem load "Element_PkEdge"`) $P_1$ on each edge of the mesh.

* `:::freefem [P2Edge]` piecewise $P_{2}$ discontinuous finite element (2d) (need `:::freefem load "Element_PkEdge"`) $P_2$ on each edge of the mesh.

* `:::freefem [P3Edge]` piecewise $P_{3}$ discontinuous finite element (2d) (need `:::freefem load "Element_PkEdge"`) $P_3$ on each edge of the mesh.

* `:::freefem [P4Edge]` piecewise $P_{4}$ discontinuous finite element (2d) (need `:::freefem load "Element_PkEdge"`) $P_4$ on each edge of the mesh.

* `:::freefem [P5Edge]` piecewise $P_{5}$ discontinuous finite element (2d) (need `:::freefem load "Element_PkEdge"`) $P_5$ on each edge of the mesh.

* `:::freefem [P2Morley]` piecewise $P_{2}$ non conform finite element (2d) (need `:::freefem load "Morley"`)
	\begin{equation*}
		P2_{h} = \left\{ v \in L^2(\Omega) \left|\; \forall K \in \mathcal{T}_{h}, \ v_{|K} \in P_{3},
		\left\{\begin{array}{c}
			v \mbox{ continuous at vertices,}\\
			\p_n{v} \mbox{ continuous at middle of edge,}
		\end{array}\right.
		\right.\right\}
	\end{equation*}
	where
	$P_{2}$ is the set of polynomials of $\R^{2}$ of degrees $\le 2$.

	!!! warning
		To build the interplant of a function $u$ (scalar) for this finite element, we need the function and 2 partial derivatives $(u,u_x, u_y)$, so this vectorial finite element with 3 components $(u,u_x,u_y)$.

		See example `:::freefem bilapMorley.edp` of \verb!examples++-load! $\codered$ for solving BiLaplacien problem :

	```freefem
		load "Morley"
		fespace Vh(Th,P2Morley); //The Morley finite element space
		macro bilaplacien(u,v) (dxx(u)*dxx(v) + dyy(u)*dyy(v) + 2.*dxy(u)*dxy(v)) //
		real f = 1;
		Vh [u, ux, uy], [v, vx, vy];

		solve bilap ([u, ux, uy], [v, vx, vy])
			= int2d(Th)(
				bilaplacien(u,v)
			)
			- int2d(Th)(
				f*v
			)
			+ on(1, 2, 3, 4, u=0, ux=0, uy=0)
			;
	```

 * `:::freefem [HCT]` $P_3$ $C^1$ conform finite element (2d) (need `:::freefem load "Element_HCT"`) one 3 sub triangles.

	Let call $\mathcal{T}^\triangle_{h}$ the sub mesh of $\mathcal{T}_{h}$ where all triangle are split in 3 at the barycenter.
	\begin{equation}
		PHCT_{h} = \left\{ v \in C^1(\Omega) \left|\; \forall K \in \mathcal{T}^\triangle_{h}, \ v_{|K} \in P_{3} \right.\right\}
	\end{equation}
	where
	$P_{3}$ is the set of polynomials of $\R^{2}$ of degrees $\le 3$. The degree of freedom are the value and derivative at vertices and normal derivative a middle edge point of initial meshes,and thank to \cite{HCT} $\codered$.

	!!!warning
		To build the interplant of a function $u$ (scalar) for this finite element, we need the function and 2 partial derivatives $(u,u_x, u_y)$, so this vectorial finite element with 3 components $(u,u_x,u_y)$ like in previous Finite Element.

 * `:::freefem [P2BR]` (need `:::freefem load "BernadiRaugel"`) the Bernadi Raugel Finite Element is a Vectorial element (2d) with 2 components, See Bernardi, C., Raugel, G.: Analysis of some finite elements for the Stokes problem. Math. Comp. 44, 71-79 (1985).$\codered$ It is a 2d coupled FE, with the Polynomial space is $ P_1^2$ + 3 normals bubbles edges function $(P_2)$ and the degre of freedom is 6 values at of the $2$ components at the $3$ vertices and the $3$ flux on the $3$ edges. So the number degrees of freedom is 9.

 * `:::freefem [RT0, RT03d]` Raviart-Thomas finite element of degree $0$.

	**The 2d case:**
	\begin{equation}
	\label{eq:RT0}
	RT0_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in	\mathcal{T}_{h} ,\ \mathbf{v}_{|K}(x,y) =
		\vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{x}{y} \right.\right\}
	\end{equation}
	**The 3d case:**
	\begin{equation}
		\label{eq:RT03d}
		RT0_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in \mathcal{T}_{h},\ \mathbf{v}_{|K}(x,y,z) =
			\vectthree{\alpha^1_{K}}{\alpha^2_{K}}{\alpha^3_{K}} + \beta_{K}\vectthree{x}{y}{z} \right.\right\}
	\end{equation}
	where by writing $\textrm{div }\mathbf{w}=\sum_{i=1}^d\p w_i/\p x_i$ with $\mathbf{w}=(w_i)_{i=1}^d$:
	$$
	H(\textrm{div})=\left\{\mathbf{w}\in L^{2}(\Omega)^d\left|\textrm{div } \mathbf{w}\in L^{2}(\Omega)\right.\right\}
	$$
	and where $\alpha^1_{K}$, $\alpha^2_{K}$, $\alpha^3_{K}$, $\beta_{K}$ are real numbers.

* `:::freefem [RT0Ortho]` Raviart-Thomas Orthogonal, or Nedelec finite element type I of degree $0$ in dimension 2
    \begin{equation}
         RT0Ortho{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{-y}{x}  \right.\right\}
         \label{RT0Ortho}
     \end{equation}      

* `:::freefem [Edge03d]` 3d Nedelec finite element or Edge Element of degree $0$.

     **The 3d case:**
     \begin{equation}
         Edge0_{h} = \left\{ \mathbf{v} \in H(\textrm{Curl}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \mathbf{v}_{|K}(x,y,z) =
         \vectthree{\alpha^1_{K}}{\alpha^2_{K}}{\alpha^3_{K}} + \vectthree{\beta^1_{K}}{\beta^2_{K}}{\beta^3_{K}}\times\vectthree{x}{y}{z}  \right.\right\}
         \label{eq:Edge03d}
     \end{equation}
      where by writing $\textrm{curl}\mathbf{w}=\vectthree{\p w_2/\p x_3-\p w_3/\p x_2}{\p w_3/\p x_1-\p w_1/\p x_3}{\p w_1/\p x_2-\p w_2/\p x_1}$ with
      $ \mathbf{w}=(w_i)_{i=1}^d$,
      $$
      H(\textrm{curl})=\left\{\mathbf{w}\in L^{2}(\Omega)^d\left|
      \textrm{curl } \mathbf{w}\in L^{2}(\Omega)^d
      \right.\right\}
      $$
      and
      $\alpha^1_{K},\alpha^2_{K},\alpha^3_{K},\beta^1_{K},\beta^2_{K},\beta^3_{K}$ are real numbers.

 * `:::freefem [Edge13d]` (need `:::freefem load "Element_Mixte3d"`) 3d Nedelec finite element or Edge Element of degree $1$.

 * `:::freefem [Edge23d]` (need `:::freefem load "Element_Mixte3d"`) 3d Nedelec finite element or Edge Element of degree $2$.

* `:::freefem [P1nc]` piecewise linear element continuous at the middle of edge only in 2D (Crouzeix-Raviart Finite Element 2d).

* `:::freefem [P2pnc]` piecewise quadratic plus a bubble P3 element with the continuity of the 2 moments on each edge (version 3.59) (need `:::freefem load "Element_P2pnc"`

* `:::freefem [RT1]` (need `:::freefem load "Element_Mixte"`)
     \begin{equation}
         RT1_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_1^2,P_0,  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{x}{y}   \right.\right\}
         \label{eq:RT1}
     \end{equation}

* `:::freefem [RT1Ortho]` (need `:::freefem load "Element_Mixte"`, version 3.13, dimension 2 $\codered$)
         \begin{equation}
         RT1_{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h},  \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_1^2,P_0,  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{-y}{x}   \right.\right\}
         \label{eq:RT1Ortho}
     \end{equation}

  * `:::freefem [RT2]` (need `:::freefem load "Element_Mixte"`)
     \begin{equation}
         RT2_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad   \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_2^2, P_1,  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{x}{y}   \right.\right\}
         \label{eq:RT2}
     \end{equation}

* `:::freefem [RT2Ortho]` (need `:::freefem load "Element_Mixte"`, version 3.59, dimension 2 $\codered$)
     \begin{equation}
         RT2_{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h} ,  \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_2^2, P_1,  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{-y}{x}   \right.\right\}
         \label{eq:RT1Ortho}
     \end{equation}

* `:::freefem [BDM1]` (need `:::freefem load "Element_Mixte"`, version 3.13, dimension 2 $\codered$) the Brezzi-Douglas-Marini finite element
     \begin{equation}
         BDM1_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad   \mathbf{v}_{|K} \in P_1^2
         \right.\right\}
         \label{eq:BDM1}
     \end{equation}

* `:::freefem [BDM1Ortho]` (need `:::freefem load "Element_Mixte"`, version 3.13, dimension 2 $\codered$) the Brezzi-Douglas-Marini Orthogonal also call Nedelec of type II , finite element
       \begin{equation}
         BDM1Ortho_{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad   \mathbf{v}_{|K} \in P_1^2
         \right.\right\}
         \label{eq:BDM1Ortho}
     \end{equation}

* `:::freefem [FEQF]` (need `:::freefem load "Element_QF"`, $\codered$ version 3.45, dimension 2 or 3) the finite element to store function at default quadrature points (so the quadrature is `:::freefem qf5pT` in 2d and is `:::freefem qfV5` in 3d). For over quadrature you have the following correspondance finite element, quadrature formula.

	* `:::freefem FEQF1` $\mapsto$ `:::freefem qf1pT`,
	* `:::freefem FEQF2` $\mapsto$ `:::freefem qf2pT`,
	* `:::freefem FEQF5` $\mapsto$ `:::freefem qf5pT`,
	* `:::freefem FEQF7` $\mapsto$ `:::freefem qf7pT`,
	* `:::freefem FEQF9` $\mapsto$ `:::freefem qf9pT`,
 	* `:::freefem FEQF13d` $\mapsto$ `:::freefem qfV1`,
	* `:::freefem FEQF23d` $\mapsto$ `:::freefem qfV2`,
	* `:::freefem FEQF53d` $\mapsto$ `:::freefem qfV5`

   You can use this element element of do optimization to store and reuse function with long formula in non linear process in integral.
