## Useful macro

```ffpp
macro grad(u) [dx(u), dy(u)] //
macro Grad(U) [grad(U#x), grad(U#y)] //

macro div(u, v) (dx(u) + dy(v)) //
macro Div(U) (dx(U#x) + dy(U#y)) //
```

## Laplacian
Solve the academic problem :
$$
\left\{
\begin{array}{rcll}
	\Delta u &=& f & \text{ on }\Omega\\
	u &=& 0 & \text{ on }\Gamma\\
\end{array}
\right.
$$

```ffpp
// Mesh
int nn = 10;
mesh Th = square(nn, nn);

// Fespace
fespace Uh(Th, P1);
Uh u, uh;

// Macro
macro grad(u) [dx(u), dy(u)] //

// Problem
problem Laplacian (u, uh)
	= int2d(Th)(
		  grad(u)' * grad(uh)
	)
	+ int2d(Th)(
		  1. * uh
	)
	+ on(1, 2, 3, 4, u=0)
	;

// Solve
Laplacian;

// Plot
plot(u, cmm="Laplacian");
```

## Laplacian with `varf`
Same problem as [Laplacian](#laplacian) but with varf usage.

```ffpp
// Mesh
int nn = 10;
mesh Th = square(nn, nn);

// Fespace
fespace Uh(Th, P1);
Uh u;

// Macro
macro grad(u) [dx(u), dy(u)] //

// Varf
varf vLaplacian (u, uh)
	= int2d(Th)(
		  grad(u)' * grad(uh)
	)
	+ int2d(Th)(
		  1. * uh
	)
	+ on(1, 2, 3, 4, u=0)
	;

// Build matrix and vectors
matrix<real> A = vLaplacian(Uh, Uh);
real[int] b = vLaplacian(0, Uh);

// Solve
u[] = A^-1 * b;

//Plot
plot(u, cmm="Laplacian");
```

## Linear elasticity

```ffpp
// Parameters
int nn = 10;
int L = 10.;
int H = 1.;

real Rho = 800.;	//Density
real E = 200.e9;	//Young's modulus
real Nu = 0.3;		//Poisson's ratio

real fx = 0.;
real fy = -9.81 * Rho;

real g = 0.;
real UDx = 0.;
real UDy = 0.;

real T = 100.;
real dt = 1.e-1;

// Mesh
mesh Th = square(L*nn, H*nn, [L*x, H*y]);
int[int] GammaD = [2];
int[int] GammaN = [1, 3, 4];

// Fespace
fespace Vh(Th, [P1, P1]);
Vh [Ux, Uy];
Vh [Upx, Upy];
Vh [Vx, Vy];

// Macro
macro div(U) (dx(U#x) + dy(U#y)) //
macro Epsilon(U) [[dx(U#x), 0.5*(dx(U#y) + dy(U#x))], [0.5*(dx(U#y) + dy(U#x)), dy(U#y)]] //

// Problem
real Lambda = E / (2.*(1.+Nu));
real Mu = E*Nu / ((1.+Nu)*(1.-2.*Nu));
problem ElasticiteLineaire ([Ux, Uy], [Vx, Vy])
	= int2d(Th)(
		  (Rho/dt) * [Ux, Uy]' * [Vx, Vy]
		+ Lambda * div(U) * div(V)
		+ 2.*Mu * (Epsilon(U) : Epsilon(V))
	)
	- int2d(Th)(
		  (Rho/dt) * [Upx, Upy]' * [Vx, Vy]
		+ [fx, fy]' * [Vx, Vy]
	)
	- int1d(Th, GammaN)(
		  g * [N.x, N.y]' * [Vx, Vy]
	)
	+ on(GammaD, Ux=UDx, Uy=UDy)
	;

// Initialization
[Ux, Uy] = [0, 0];

// Iterations
int NbIter = T/dt;
for (int i = 0; i < NbIter; i++){
	// Update
	[Upx, Upy] = [Ux, Uy];
	
	// Solve
	ElasticiteLineaire;
	
	// Movemesh
	Th = movemesh(Th, [x+dt*(Ux-Upx), y+dt*(Uy-Upy)]);
	[Ux, Uy] = [Ux, Uy];
	
	// Plot
	plot([Ux, Uy]);
}
```



