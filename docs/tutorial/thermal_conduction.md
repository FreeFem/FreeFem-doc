# Thermal Conduction

**Summary :** _Here we shall learn how to deal with a time dependent parabolic problem. We shall also show how to treat an axisymmetric problem  and show also how to deal with a nonlinear problem_

**How air cools a plate**

We seek the temperature distribution in a plate $(0,Lx)\times(0,Ly)\times(0,Lz)$
of rectangular cross section $\Omega=(0,6)\times(0,1)$; the plate is
surrounded by air at temperature $u_e$ and
initially at temperature $u=u_0+\frac x L u_1$. In the plane perpendicular to the plate
at $z=Lz/2$,  the temperature varies little with
the coordinate $z$; as a first approximation the problem is 2D.

We must solve the temperature equation in $\Omega$ in a time interval (0,T).

\begin{eqnarray}&&
    \p_t u -\n\cdot(\kappa\n u)=0 \hbox{ in } \Omega\times(0,T),
    \cr&&
    u(x,y,0)=u_0+x u_1
    \cr&&
    \kappa\frac{\p u}{\p n} +\alpha(u-u_e)=0\hbox{ on } \Gamma\times(0,T).
\end{eqnarray}

Here the diffusion $\kappa$ will take two values, one below the middle horizontal line and ten times less
above, so as to simulate a thermostat.
The term $\alpha(u-u_e)$ accounts for the loss of temperature by convection in air.  Mathematically
this boundary condition is of Fourier (or Robin, or mixed) type.

The variational formulation is  in $L^2(0,T;H^1(\Omega))$; in loose terms and after applying an implicit Euler
finite difference approximation in time; we shall seek $u^n(x,y)$ satisfying for all $w\in H^1(\Omega)$:

\[
    \int_\Omega(\frac{u^n-u^{n-1}}{\delta t} w + \kappa\n u^n\n w) +\int_\Gamma\alpha(u^n-u_ue)w=0
\]

```freefem
func u0 =10+90*x/6;
func k = 1.8*(y<0.5)+0.2;
real ue = 25, alpha=0.25, T=5, dt=0.1 ;

mesh Th=square(30,5,[6*x,y]);
fespace Vh(Th,P1);
Vh u=u0,v,uold;

problem thermic(u,v)= int2d(Th)(u*v/dt + k*(dx(u) * dx(v) + dy(u) * dy(v)))
                + int1d(Th,1,3)(alpha*u*v)
                - int1d(Th,1,3)(alpha*ue*v)
                - int2d(Th)(uold*v/dt) + on(2,4,u=u0);
ofstream ff("thermic.dat");
for(real t=0;t<T;t+=dt){
    uold=u;  // uold $\equiv  u^{n-1} = u^n \equiv $u
    thermic; // here solve the thermic problem
    ff<<u(3,0.5)<<endl;
    plot(u);
}
```

Notice that we must separate by hand the bilinear part from the linear one.

Notice also that  the way we store the temperature at point (3,0.5) for all times in file `thermic.dat`. Should a one dimensional plot be required, the same procedure can be used.  For instance to print $x\mapsto \frac{\p u}{\p y}(x,0.9)$ one would do

```freefem
for(int i=0;i<20;i++) cout<<dy(u)(6.0*i/20.0,0.9)<<endl;
```

Results are shown on Fig. 3.4.

|Fig. 3.4: Temperature at $T=4.9$|Decay of temperature versus time at $x=3, y=0.5$|
|:----:|:----:|
|![Temperature](images/thermic.svg)|![Temperature decay](images/thermicvst.svg)|

## Axisymmetry: 3D Rod with circular section

Let us now deal with a cylindrical rod instead of a flat plate.  For simplicity we take $\kappa=1$.
In cylindrical coordinates, the Laplace
operator becomes ($r$ is the distance to the axis, $z$ is the distance along the
axis, $\theta$ polar angle in a fixed plane perpendicular to the axis):

$$ \Delta u = {1\over r}\p _r(r\p _r u) + {1\over r^2}\p ^2_{\theta\theta} u
 + \p ^2_{z z}.
$$

Symmetry implies that we loose the dependence with respect to
$\theta$; so the domain $\Omega$ is again a rectangle $]0,R[\times]0,|[$ . We take the convention of
numbering of the edges as in `square()` (1 for the bottom horizontal ...);
the problem is now:

\begin{eqnarray}&&
r\p_t u-\p _r(r\p _r u) - \p _z(r\p _z u) = 0 \hbox{ in } \Omega,
\cr&&
u(t=0) = u_0 + \frac z{L_z} (u_1-u)
\cr&&
u|_{\Gamma_4} = u_0,\quad  u|_{\Gamma_2} = u_1,
\quad \alpha(u-u_e) + {\p u\over \p n} |_{\Gamma_1\cup\Gamma_3} = 0.
\end{eqnarray}

Note that the PDE has been multiplied by $r$.

After discretization in time with an implicit scheme, with time steps `dt`,
in the FreeFem++ syntax $r$ becomes $x$ and $z$ becomes $y$ and the problem is:

```freefem
problem thermaxi(u,v)=int2d(Th)((u*v/dt + dx(u)*dx(v) + dy(u)*dy(v))*x)
                + int1d(Th,3)(alpha*x*u*v) - int1d(Th,3)(alpha*x*ue*v)
                - int2d(Th)(uold*v*x/dt) + on(2,4,u=u0);
```

Notice that the bilinear form degenerates at $x=0$. Still one can prove existence and uniqueness
for $u$ and because of this degeneracy no boundary conditions need to be imposed on $\Gamma_1$.

## A Nonlinear Problem : Radiation

Heat loss through radiation is a loss proportional
to the absolute temperature to the fourth power (Stefan's Law). This adds to the
loss by convection and gives the following boundary condition:

$$\kappa{\p u\over \p n} +\alpha(u-u_e) + c[(u + 273)^4 - (u_e+273)^4] = 0$$

The problem is \x{nonlinear}, and must be solved iteratively. If $m$
denotes the iteration index, a semi-linearization of the radiation condition gives
$${\p u^{m+1}\over \p n} + \alpha(u^{m+1}-u_e)+ c(u^{m+1}-u_e)
(u^m+u_e +546) ((u^m + 273)^2 + (u_e+273)^2)  = 0,
$$

because we have the identity $ a^4 - b^4 = (a-b)(a+b)(a^2+b^2)$.
The iterative process will work with $v=u-u_e$.

```freefem
...
fespace Vh(Th,P1); // Finite element space
real rad=1e-8, uek=ue+273; // Definition of the physical constants
Vh vold,w,v=u0-ue,b;
problem thermradia(v,w)
    = int2d(Th)(v*w/dt + k*(dx(v) * dx(w) + dy(v) * dy(w)))
                + int1d(Th,1,3)(b*v*w)
                - int2d(Th)(vold*w/dt) + on(2,4,v=u0-ue);

for(real t=0;t<T;t+=dt){
    vold=v;
    for(int m=0;m<5;m++){
       b= alpha + rad * (v + 2*uek) * ((v+uek)^2 + uek^2);
       thermradia;
    }
}
vold=v+ue; plot(vold);
```
