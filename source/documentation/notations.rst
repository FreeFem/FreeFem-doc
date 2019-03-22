.. role:: freefem(code)
  :language: freefem

.. role:: cpp(code)
  :language: c

Notations
=========

Here mathematical expressions and corresponding **FreeFEM** commands are explained.

Generalities
------------

-  [:math:`\delta_{ij}`] Kronecker delta (:math:`0` if :math:`i\neq j`, 1 if :math:`i=j` for integers :math:`i,j`)
-  [:math:`\forall`] for all
-  [:math:`\exists`] there exists
-  [i.e.] that is
-  [PDE] partial differential equation (with boundary conditions)
-  [:math:`\emptyset`] the empty set
-  [:math:`\mathbb{N}`] the set of integers (:math:`a\in \mathbb{N}\Leftrightarrow` :freefem:`int a`), :freefem:`int` means :cpp:`long int` inside **FreeFEM**
-  [:math:`\mathbb{R}`] the set of real numbers (:math:`a\in \mathbb{R}\Leftrightarrow` :freefem:`real a`), :cpp:`double` inside **FreeFEM**
-  [:math:`\mathbb{C}`] the set of complex numbers (:math:`a\in \mathbb{C}\Leftrightarrow` :freefem:`complex a`), :cpp:`complex<double>`
-  [:math:`\mathbb{R}^d`] :math:`d`-dimensional Euclidean space

Sets, Mappings, Matrices, Vectors
---------------------------------

Let :math:`E,\, F,\, G` be three sets and :math:`A` the subset of :math:`E`.

-  [:math:`\{x\in E|\; P\}`] the subset of :math:`E` consisting of the elements possessing the property :math:`P`
-  [:math:`E\cup F`] the set of elements belonging to :math:`E` or :math:`F`
-  [:math:`E\cap F`] the set of elements belonging to :math:`E` and :math:`F`
-  [:math:`E\setminus A`] the set :math:`\{x\in E|\; x\not\in A\}`
-  [:math:`E+F`] :math:`E\cup F` with :math:`E\cap F=\emptyset`
-  [:math:`E\times F`] the Cartesian product of :math:`E` and :math:`F`
-  [:math:`E^n`] the :math:`n`-th power of :math:`E` (:math:`E^2=E\times E`, :math:`E^n=E\times E^{n-1}`)
-  [:math:`f:\; E\to F`] the mapping form :math:`E` into :math:`F`, i.e., :math:`E\ni x\mapsto f(x)\in F`
-  [:math:`I_E` or :math:`I`] the identity mapping in :math:`E`,i.e., :math:`I(x)=x\quad \forall x\in E`
-  [:math:`f\circ g`] for :math:`f:\; F\to G` and :math:`g:\; E\to F`, :math:`E\ni x\mapsto (f\circ g)(x)=f(g(x))\in G` (see :ref:`Elementary function <typeElementaryFunctions>`)
-  [:math:`f|_A`] the restriction of :math:`f:\; E\to F` to the subset :math:`A` of :math:`E`
-  [:math:`\{a_k\}`] column vector with components :math:`a_k`
-  [:math:`(a_k)`] row vector with components :math:`a_k`
-  [:math:`(a_{k})^T`] denotes the transpose of a matrix :math:`(a_{k})`, and is :math:`\{a_{k}\}`
-  [:math:`\{a_{ij}\}`] matrix with components :math:`a_{ij}`, and :math:`(a_{ij})^T=(a_{ji})`

Numbers
-------

For two real numbers :math:`a,b`

-  :math:`[a,b]` is the interval :math:`\{x\in \mathbb{R}|\; a\le x\le b\}`
-  :math:`]a,b]` is the interval :math:`\{x\in \mathbb{R}|\; a< x\le b\}`
-  :math:`[a,b[` is the interval :math:`\{x\in \mathbb{R}|\; a\le x< b\}`
-  :math:`]a,b[` is the interval :math:`\{x\in \mathbb{R}|\; a< x< b\}`

Differential Calculus
---------------------

-  [:math:`\partial f/\partial x`] the partial derivative of :math:`f:\mathbb{R}^d\to \mathbb{R}` with respect to :math:`x` (:freefem:`dx(f)`)
-  [:math:`\nabla f`] the gradient of :math:`f:\Omega\to \mathbb{R}`,i.e., :math:`\nabla f=(\partial f/\partial x,\, \partial f/\partial y)`
-  [:math:`\text{div}(\mathbf{f})` or :math:`\nabla.\mathbf{f}`] the divergence of :math:`\mathbf{f}:\Omega\to \mathbb{R}^d`, i.e., :math:`\text{div}(\mathbf{f})=\partial f_1/\partial x+\partial f_2/\partial y`
-  [:math:`\Delta f`] the Laplacian of :math:`f:\; \Omega\to \mathbb{R}`, i.e., :math:`\Delta f=\partial^2f/\partial x^2+\partial^2 f/\partial y^2`

Meshes
------

-  [:math:`\Omega`] usually denotes a domain on which PDE is defined
-  [:math:`\Gamma`] denotes the boundary of :math:`\Omega`,i.e., :math:`\Gamma=\partial\Omega` (keyword :freefem:`border`, see :ref:`Border <meshBorder>`)
-  [:math:`\mathcal{T}_h`] the triangulation of :math:`\Omega`, i.e., the set of triangles :math:`T_k`, where :math:`h` stands for mesh size (keyword :freefem:`mesh`, :freefem:`buildmesh`, see :ref:`Mesh Generation <meshGeneration>`)
-  [:math:`n_t`] the number of triangles in :math:`\mathcal{T}_h` (get by :freefem:`Th.nt`)
-  [:math:`\Omega_h`] denotes the approximated domain :math:`\Omega_h=\cup_{k=1}^{n_t}T_k` of :math:`\Omega`.
   If :math:`\Omega` is polygonal domain, then it will be :math:`\Omega=\Omega_h`
-  [:math:`\Gamma_h`] the boundary of :math:`\Omega_h`
-  [:math:`n_v`] the number of vertices in :math:`\mathcal{T}_h` (get by :freefem:`Th.nv`)
-  [:math:`n_{be}`] the number of boundary element in :math:`\mathcal{T}_h` (get by :freefem:`Th.nbe`)
-  [:math:`|\Omega_h|`] the measure (area or volume) in :math:`\mathcal{T}_h` (get by :freefem:`Th.measure`)
-  [:math:`|\partial \Omega_h|`] the measure of the border (length or area) in :math:`\mathcal{T}_h` (get by :freefem:`Th.bordermeasure`)
-  [:math:`h_{min}`] the minimum edge size of :math:`\mathcal{T}_h` (get by :freefem:`Th.hmin`)
-  [:math:`h_{max}`] the maximum edge size of :math:`\mathcal{T}_h` (get by :freefem:`Th.hmax`)
-  [[:math:`q^iq^j`]] the segment connecting :math:`q^i` and :math:`q^j`
-  [:math:`q^{k_1},q^{k_2},q^{k_3}`] the vertices of a triangle :math:`T_k` with anti-clock direction (get the coordinate of :math:`q^{k_j}` by :freefem:`(Th[k-1][j-1].x, Th[k-1][j-1].y)`)
-  [:math:`I_{\Omega}`] the set :math:`\{i\in \mathbb{N}|\; q^i\not\in \Gamma_h\}`

Finite Element Spaces
---------------------

-  [:math:`L^2(\Omega)`] the set :math:`\displaystyle{\left\{w(x,y)\left|\; \int_{\Omega}|w(x,y)|^2\text{d} x\text{d} y<\infty\right.\right\}}`

  .. math::

    \textrm{norm:}\; \| w\|_{0,\Omega}&=\left(\int_{\Omega}|w(x,y)|^2\text{d} x\text{d} y\right)^{1/2}\\
    \textrm{scalar product:}\; (v,w)&=\int_{\Omega}vw

-  [:math:`H^1(\Omega)`] the set :math:`\displaystyle{\left\{w\in L^2(\Omega)\left|\; \int_{\Omega}\left(|\partial w/\partial x|^2+|\partial w/\partial y|^2\right)\text{d} x\text{d} y <\infty\right.\right\}}`

  .. math::

    \textrm{norm:}\; \| w\|_{1,\Omega}=\left(\| w\|_{0,\Omega}^2+\|\nabla u\|_{0.\Omega}^2\right)^{1/2}

-  [:math:`H^m(\Omega)`] the set :math:`\displaystyle{\left\{w\in L^2(\Omega)\left|\; \int_{\Omega}\frac{\partial^{|\alpha|} w}{\partial x^{\alpha_1}\partial y^{\alpha_2}}\in L^2(\Omega)\quad\forall \alpha=(\alpha_1,\alpha_2)\in \mathbb{N}^2,\, |\alpha|=\alpha_1+\alpha_2\right.\right\}}`

  .. math::

    \textrm{scalar product:}\; (v,w)_{1,\Omega}=
    \sum_{|\alpha|\le m}\int_{\Omega} D^{\alpha}v D^{\alpha}w

-  [:math:`H^1_0(\Omega)`] the set :math:`\left\{w\in H^1(\Omega)\left|\; u=0\quad \textrm{on }\Gamma\right.\right\}`

   [:math:`L^2(\Omega)^2`] denotes :math:`L^2(\Omega)\times L^2(\Omega)`, and also :math:`H^1(\Omega)^2=H^1(\Omega)\times H^1(\Omega)`
-  [:math:`V_h`] denotes the finite element space created by :freefem:`fespace Vh(Th, *)` in **FreeFEM** (see :ref:`Finite Elements <finiteElement>` for ``*``)
-  [:math:`\Pi_h f`] the projection of the function :math:`f` into :math:`V_h` (:freefem:`func f=x^2*y^3; Vh v = f;`) means :math:`v = Pi_h (f) * [\{v\}]` for FE-function :math:`v` in :math:`V_h` means the column vector :math:`(v_1,\cdots,v_M)^T` if :math:`v=v_1\phi_1+\cdots+v_M\phi_M`, which is shown by :freefem:`fespace Vh(Th, P2); Vh v; cout << v[] << endl;`
