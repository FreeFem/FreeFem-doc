Consider an elastic plate whose displacement change vertically, which is made up of three plates of different materials, welded on each other.

Let $\Omega_i$, $i=1,2,3$ be the domain occupied by $i$-th material with tension $\mu_i$ (see [Soap film](StaticProblems/#soap-film)).

The computational domain $\Omega$ is the interior of $\overline{\Omega_1}\cup \overline{\Omega_2}\cup \overline{\Omega_3}$. The vertical displacement $u(x,y)$ is obtained from

\begin{eqnarray}
	\label{eqn::transm-1}
	-\mu_i\Delta u &=& f & \textrm{ in }\Omega_i\\
	\mu_i\p_n u|_{\Gamma_{i}} &=& -\mu_j\p_n u|_{\Gamma_{j}} & \textrm{ on }\overline{\Omega_{i}}\cap\overline{\Omega_{j}} \textrm{ if }1\le i< j\le 3
	\label{eqn::transm-2}
\end{eqnarray}

where $\p_n u|_{\Gamma_{i}}$ denotes the value of the normal derivative $\p_n u$ on the boundary $\Gamma_i$ of the domain $\Omega_i$.

By introducing the characteristic function $\chi_i$ of $\Omega_i$, that is,

\begin{equation}
\chi_i(x)=1\ \textrm{ if }x\in\Omega_i;\
\chi_i(x)=0\ \textrm{ if }x\not\in\Omega_i
\end{equation}

we can easily rewrite \eqref{eqn::transm-1} and \eqref{eqn::transm-2} to the weak form. Here we assume that $u=0$ on $\Gamma=\p\Omega$.

Transmission problem: For a given function $f$, find $u$ such that

\begin{equation}
	a(u,v) = \ell(f,v) \textrm{ for all }v\in H^1_0(\Omega)
\end{equation}
\begin{eqnarray}
	a(u,v) &=& \int_{\Omega}\mu \nabla u\cdot \nabla v\nonumber\\
	\ell(f,v) &=& \int_{\Omega}fv\nonumber
\end{eqnarray}

where $\mu=\mu_1\chi_1+\mu_2\chi_2+\mu_3\chi_3$. Here we notice that $\mu$ become the discontinuous function.

With dissipation, and at the thermal equilibrium, the temperature equation is:

$\codered$ missing in the originial doc!

This example explains the definition and manipulation of _region_, i.e. sub-domains of the whole domain. Consider this L-shaped domain with 3 diagonals as internal boundaries, defining 4 sub-domains:

```freefem
// Mesh
border a(t=0, 1){x=t; y=0;};
border b(t=0, 0.5){x=1; y=t;};
border c(t=0, 0.5){x=1-t; y=0.5;};
border d(t=0.5, 1){x=0.5; y=t;};
border e(t=0.5, 1){x=1-t; y=1;};
border f(t=0, 1){x=0; y=1-t;};
border i1(t=0, 0.5){x=t; y=1-t;};
border i2(t=0, 0.5){x=t; y=t;};
border i3(t=0, 0.5){x=1-t; y=t;};
mesh th = buildmesh(a(6) + b(4) + c(4) +d(4) + e(4)
	+ f(6) + i1(6) + i2(6) + i3(6));

// Fespace
fespace Ph(th, P0); //constant discontinuous functions / element
Ph reg=region; //defined the P0 function associated to region number

// Plot
plot(reg, fill=true, wait=true, value=true);

```

|<a name="Fig30">Fig. 30</a>: The function `:::freefem reg`|
|:----:|
|![Region](images/TransmissionProblem_Region.png)|

`:::freefem region` is a keyword of FreeFem++ which is in fact a variable depending of the current position (is not a function today, use `:::freefem Ph reg=region;` to set a function). This variable value returned is the number of the sub-domain of the current position. This number is defined by `:::freefem buildmesh` which scans while building the mesh all its connected component.

So to get the number of a region containing a particular point one does:

```freefem
// Characteristic function
int nupper = reg(0.4, 0.9); //get the region number of point (0.4,0.9)
int nlower = reg(0.9, 0.1); //get the region number of point (0.4,0.1)
cout << "nlower = " <<  nlower << ", nupper = " << nupper<< endl;
Ph nu = 1 + 5*(region==nlower) + 10*(region==nupper);

// Plot
plot(nu, fill=true,wait=true);
```

|<a name"=Fig31">Fig. 31</a>: The function `:::freefem nu`|
|:----:|
|![Nu](images/TransmissionProblem_Nu.png)|

This is particularly useful to define discontinuous functions such as might occur when one part of the domain is copper and the other one is iron, for example.

We this in mind we proceed to solve a Laplace equation with discontinuous coefficients ($\nu$ is 1, 6 and 11 below).

```freefem
// Problem
solve lap (u, v)
	= int2d(th)(
		  nu*(dx(u)*dx(v) + dy(u)*dy(v))
	)
	+ int2d(th)(
		- 1*v
	)
	+ on(a, b, c, d, e, f, u=0)
	;

// Plot
plot(u);
```

|<a name="Fig32">Fig. 32</a>: The isovalue of the solution $u$|
|:----:|
|![U](images/TransmissionProblem_U.png)|
