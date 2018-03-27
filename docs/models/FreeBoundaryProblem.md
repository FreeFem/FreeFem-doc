The domain $\Omega$ is defined with:

```freefem
// Parameters
real L = 10; //length
real hr = 2.1; //left height
real hl = 0.35; //right height
int n = 4;

// Mesh
border a(t=0, L){x=t; y=0;}; //bottom: Gamma_a
border b(t=0, hr){x=L; y=t;}; //right: Gamma_b
border f(t=L, 0){x=t; y=t*(hr-hl)/L+hl;}; //free surface: Gamma_f
border d(t=hl, 0){x=0; y=t;}; // left: Gamma_d
mesh Th = buildmesh(a(10*n) + b(6*n) + f(8*n) + d(3*n));
plot(Th);
```

|<a name="Fig33">Fig. 33</a>: The mesh of the domain $\Omega$|
|:----:|
|![Mesh1](images/FreeBoundary_Mesh1.png)|

The free boundary problem is:

Find $u$ and $\Omega$ such that:
\begin{equation*}
	\left\{
	\begin{array}{rcll}
		-\Delta u &=& 0 & \mbox{ in }\Omega\\
		u &=& y & \mbox{ on }\Gamma_b\\
		\p u \over \p n &=& 0 & \mbox{ on }\Gamma_d\cup\Gamma_a\\
		\p u \over \p n &=& {q \over K} n_x & \mbox{ on }\Gamma_f\\
		u &=& y & \mbox{ on }\Gamma_f
	\end{array}
	\right.
\end{equation*}

We use a fixed point method;

$\Omega^0 = \Omega$

In two step, fist we solve the classical following problem:

\begin{equation*}
	\left\{
	\begin{array}{rcll}
		-\Delta u &=& 0 & \mbox{ in }\Omega^n\\
		u &=& y & \mbox{ on }\Gamma^n_b\\
		\p u \over \p n &=& 0 & \mbox{ on }\Gamma^n_d\cup\Gamma^n_a\\
		u &=& y & \mbox{ on }\Gamma^n_f
	\end{array}
	\right.
\end{equation*}

The variational formulation is:

Find $u$ on $V=H^1(\Omega^n)$, such than $u=y$ on $\Gamma^n_b$ and $\Gamma^n_f$
$$
\int_{\Omega^n}\nabla u \nabla u' = 0,\ \forall u' \in V \mbox{ with } u' =0 \mbox{ on }\Gamma^n_b \cup \Gamma^n_f
$$

And secondly to construct a domain deformation $\mathcal{F}(x,y)=[x,y-v(x,y)]$ where $v$ is solution of the following problem:

\begin{equation*}
	\left\{
	\begin{array}{rcll}
		-\Delta v &=& 0 & \mbox{ in }\Omega^n\\
		v &=& 0 & \mbox{ on }\Gamma^n_a\\
		\p v \over \p n &=& 0 & \mbox{ on }\Gamma^n_b\cup\Gamma^n_d\\
		\p v \over \p n &=& {\p u \over \p n} - {q\over K} n_x & \mbox{ on }\Gamma^n_f
	\end{array}
	\right.
\end{equation*}

The variational formulation is:

Find $v$ on $V$, such than $v=0$ on $\Gamma^n_a$
$$
\int_{\Omega^n} \nabla v \nabla v' = \int_{\Gamma_f^n}({\p u \over \p n} - { q\over K} n_x )v',\ \quad \forall v' \in V \mbox{ with } v' =0 \mbox{ on }\Gamma^n_a
$$

Finally the new domain $\Omega^{n+1} = \mathcal{F}(\Omega^n)$

!!!question "Free boundary"
	The FreeFem++ implementation is:

	```freefem
	// Parameters
	real L = 10; //length
	real hr = 2.1; //left height
	real hl = 0.35; //right height
	int n = 4;

	real q = 0.02; //incoming flow
	real K = 0.5; //permeability

	// Mesh
	border a(t=0, L){x=t; y=0;}; //bottom: Gamma_a
	border b(t=0, hr){x=L; y=t;}; //right: Gamma_b
	border f(t=L, 0){x=t; y=t*(hr-hl)/L+hl;}; //free surface: Gamma_f
	border d(t=hl, 0){x=0; y=t;}; // left: Gamma_d
	mesh Th = buildmesh(a(10*n) + b(6*n) + f(8*n) + d(3*n));
	plot(Th);

	// Fespace
	fespace Vh(Th, P1);
	Vh u, v;
	Vh uu, vv;

	// Problem
	problem Pu (u, uu, solver=CG)
		= int2d(Th)(
			  dx(u)*dx(uu)
			+ dy(u)*dy(uu)
		)
		+ on(b, f, u=y)
		;

	problem Pv (v, vv, solver=CG)
		= int2d(Th)(
			  dx(v)*dx(vv)
			+ dy(v)*dy(vv)
		)
		+ on(a, v=0)
		+ int1d(Th, f)(
			  vv*((q/K)*N.y - (dx(u)*N.x + dy(u)*N.y))
		)
		;

	// Loop
	int j = 0;
	real errv = 1.;
	real erradap = 0.001;
	while (errv > 1e-6){
		// Update
		j++;

		// Solve
		Pu;
		Pv;

		// Plot
		plot(Th, u, v);

		// Error
		errv = int1d(Th, f)(v*v);

		// Movemesh
		real coef = 1.;
		real mintcc = checkmovemesh(Th, [x, y])/5.;
		real mint = checkmovemesh(Th, [x, y-v*coef]);

		if (mint < mintcc || j%10==0){ //mesh too bad => remeshing
			Th = adaptmesh(Th, u, err=erradap);
			mintcc = checkmovemesh(Th, [x, y])/5.;
		}

		while (1){
			real mint = checkmovemesh(Th, [x, y-v*coef]);

			if (mint > mintcc) break;

			cout << "min |T| = " << mint << endl;
			coef /= 1.5;
		}

		Th=movemesh(Th, [x, y-coef*v]);

		// Display
		cout << endl << j << " - errv = " << errv << endl;
	}

	// Plot
	plot(Th);
	plot(u, wait=true);
	```

	|<a name="Fig34">Fig. 34</a>: The final solution on the new domain $\Omega^{72}$|
	|:----:|
	|![Sol](images/FreeBoundary_Sol.png)|

	|<a name="Fig35">Fig. 35</a>: The adapted mesh of the domain $\Omega^{72}$|
	|:----:|
	|![Mesh2](images/FreeBoundary_Mesh2.png)|
