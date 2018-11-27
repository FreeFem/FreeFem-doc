Here mathematical expressions and corresponding __`FreeFem++`__ commands are explained.

## Generalities

 * [$\delta_{ij}$] Kronecker delta ($0$ if $i\neq j$, 1 if $i=j$ for integers $i,j$)
 * [$\forall$] for all
 * [$\exists$] there exists
 * [i.e.] that is
 * [PDE] partial differential equation (with boundary conditions)
 * [$\emptyset$] the empty set
 * [$\N$] the set of integers ($a\in \N\Leftrightarrow$ `:::freefme int a`), `:::freefem int` means `:::cpp long int` inside __`FreeFem++`__
 * [$\R$] the set of real numbers ($a\in \R\Leftrightarrow$ `:::freefem real a`), `:::cpp double` inside __`FreeFem++`__
 * [$\C$] the set of complex numbers ($a\in \C\Leftrightarrow$ `:::freefem complex a`), `:::cpp complex<double>`
 * [$\R^d$] $d$-dimensional Euclidean space

## Sets, Mappings, Matrices, Vectors

Let $E,\, F,\, G$ be three sets and $A$ the subset of $E$.


 * [$\{x\in E|\; P\}$] the subset of $E$ consisting of the elements possessing the property $P$
 * [$E\cup F$] the set of elements belonging to $E$ or $F$
 * [$E\cap F$] the set of elements belonging to $E$ and $F$
 * [$E\setminus A$] the set $\{x\in E|\; x\not\in A\}$
 * [$E+F$] $E\cup F$ with $E\cap F=\emptyset$
 * [$E\times F$] the Cartesian product of $E$ and $F$
 * [$E^n$] the $n$-th power of $E$ ($E^2=E\times E$, $E^n=E\times E^{n-1}$)
 * [$f:\; E\to F$] the mapping form $E$ into $F$, i.e.,
	$E\ni x\mapsto f(x)\in F$
 * [$I_E$ or $I$] the identity mapping in $E$,i.e., $I(x)=x\quad \forall x\in E$
 * [$f\circ g$] for $f:\; F\to G$ and $g:\; E\to F$, $E\ni x\mapsto (f\circ g)(x)=f(g(x))\in G$ (see [Elementary function](/reference/Types/#elementary-functions))
 * [$f|_A$] the restriction of $f:\; E\to F$ to the subset $A$ of $E$
 * [$\{a_k\}$] column vector with components $a_k$
 * [$(a_k)$] row vector with components $a_k$
 * [$(a_{k})^T$] denotes the transpose of a matrix $(a_{k})$, and is $\{a_{k}\}$
 * [$\{a_{ij}\}$] matrix with components $a_{ij}$, and $(a_{ij})^T=(a_{ji})$


## Numbers

For two real numbers $a,b$


 * [\quad]$[a,b]$ is the interval $\{x\in \R|\; a\le x\le b\}$
 * [\quad]$]a,b]$ is the interval $\{x\in \R|\; a< x\le b\}$
 * [\quad]$[a,b[$ is the interval $\{x\in \R|\; a\le x< b\}$
 * [\quad]$]a,b[$ is the interval $\{x\in \R|\; a< x< b\}$


## Differential Calculus

 * [$\p f/\p x$] the partial derivative of $f:\R^d\to \R$ with respect to $x$ (`:::freefem dx(f)`)
 * [$\nabla f$] the gradient of $f:\Omega\to \R$,i.e., $\nabla f=(\p f/\p x,\, \p f/\p y)$
 * [div$\mathbf{f}$ or $\nabla.\mathbf{f}$] the divergence of $\mathbf{f}:\Omega\to \R^d$, i.e., div$\mathbf{f}=\p f_1/\p x+\p f_2/\p y$
 * [$\Delta f$] the Laplacian of $f:\; \Omega\to \R$, i.e., $\Delta f=\p^2f/\p x^2+\p^2 f/\p y^2$


## Meshes

 * [$\Omega$] usually denotes a domain on which PDE is defined
 * [$\Gamma$] denotes the boundary of $\Omega$,i.e., $\Gamma=\p\Omega$ (keyword `:::freefem border`, see [Border](/MeshGeneration/#border))
 * [$\mathcal{T}_h$] the triangulation of $\Omega$, i.e., the set of triangles $T_k$, where $h$ stands for mesh size (keyword `:::freefem mesh`, `:::freefem buildmesh`, see [Mesh Generation](/MeshGeneration/#commands-for-mesh-generation))
 * [$n_t$] the number of triangles in $\mathcal{T}_h$ (get by `:::freefem Th.nt`)
 * [$\Omega_h$] denotes the approximated domain $\Omega_h=\cup_{k=1}^{n_t}T_k$ of $\Omega$. If $\Omega$ is polygonal domain, then it will be $\Omega=\Omega_h$
 * [$\Gamma_h$] the boundary of $\Omega_h$
 * [$n_v$] the number of vertices in $\mathcal{T}_h$ (get by `:::freefem Th.nv`)
 * [$n_{be}$] the number of boundary element in $\mathcal{T}_h$ (get by `:::freefem Th.nbe`)
 * [$|\Omega_h|$] the measure (area or volume) in $\mathcal{T}_h$ (get by `:::freefem Th.measure`)
 * [$|\partial \Omega_h|$] the measure of the border (length or area) in $\mathcal{T}_h$ (get by `:::freefem Th.bordermeasure`)
 * [$h_{min}$] the minimum edge size of $\mathcal{T}_h$ (get by `:::freefem Th.hmin`)
 * [$h_{max}$] the maximum edge size of $\mathcal{T}_h$ (get by `:::freefem Th.hmax`)
 * [[$q^iq^j$]] the segment connecting $q^i$ and $q^j$
 * [$q^{k_1},q^{k_2},q^{k_3}$] the vertices of a triangle $T_k$ with anti-clock direction (get the coordinate of $q^{k_j}$ by (`:::freefem Th[k-1][j-1].x, Th[k-1][j-1].y`)
 * [$I_{\Omega}$] the set $\{i\in \N|\; q^i\not\in \Gamma_h\}$


## Finite Element Spaces

 * [$L^2(\Omega)$] the set $\displaystyle{\left\{w(x,y)\left|\; \int_{\Omega}|w(x,y)|^2\d x\d y<\infty\right.\right\}}$
	\begin{eqnarray*}
	&&\textrm{norm:}\; \| w\|_{0,\Omega}=\left(\int_{\Omega}|w(x,y)|^2\d x\d y\right)^{1/2}\\
	&&\textrm{scalar product:}\; (v,w)=\int_{\Omega}vw
	\end{eqnarray*}
 * [$H^1(\Omega)$] the set $\displaystyle{\left\{w\in L^2(\Omega)\left|\; \int_{\Omega}\left(|\p w/\p x|^2+|\p w/\p y|^2\right)\d x\d y <\infty\right.\right\}}$
	\begin{eqnarray*}
	&&\textrm{norm:}\; \| w\|_{1,\Omega}=\left(\| w\|_{0,\Omega}^2+\|\nabla u\|_{0.\Omega}^2\right)^{1/2}
	\end{eqnarray*}
	<!--- __ --->
 * [$H^m(\Omega)$] the set $\displaystyle{\left\{w\in L^2(\Omega)\left|\; \int_{\Omega}\frac{\p^{|\alpha|} w}{\p x^{\alpha_1}\p y^{\alpha_2}}\in L^2(\Omega)\quad\forall \alpha=(\alpha_1,\alpha_2)\in \N^2,\, |\alpha|=\alpha_1+\alpha_2\right.\right\}}$
	\begin{eqnarray*}
	&&\textrm{scalar product:}\; (v,w)_{1,\Omega}=
	\sum_{|\alpha|\le m}\int_{\Omega} D^{\alpha}v D^{\alpha}w
	\end{eqnarray*}
 * [$H^1_0(\Omega)$] the set $\left\{w\in H^1(\Omega)\left|\; u=0\quad \textrm{on }\Gamma\right.\right\}$ * [$L^2(\Omega)^2$] denotes $L^2(\Omega)\times L^2(\Omega)$, and also $H^1(\Omega)^2=H^1(\Omega)\times H^1(\Omega)$
 * [$V_h$] denotes the finite element space created by `:::freefem fespace Vh(Th, *)` in __`FreeFem++`__ (see [Finite Elements](/FiniteElement/) for `*`)
 * [$\Pi_h f$] the projection of the function $f$ into $V_h$ (`:::freefem func f=x^2*y^3; Vh v = f;}` means $v = Pi_h (f) * [\{v\}]$ for FE-function $v$ in $V_h$ means the column vector $(v_1,\cdots,v_M)^T$ if $v=v_1\phi_1+\cdots+v_M\phi_M$, which is shown by `:::freefem fespace Vh(Th, P2); Vh v; cout << v[] << endl;`
