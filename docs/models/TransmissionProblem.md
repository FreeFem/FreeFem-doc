
## Transmission Problem

Consider an elastic plate whose displacement change vertically, which is made up of three plates of different materials, welded on each other.
Let $\Omega_i,\, i=1,2,3$ be the domain occupied by $i$-th material with tension $\mu_i$ (see \refSec{Soap Film} $\codered$).
The computational domain $\Omega$ is the interior of $\overline{\Omega_1}\cup \overline{\Omega_2}\cup \overline{\Omega_3}$. The vertical displacement $u(x,y)$ is obtained from

\begin{eqnarray}
-\mu_i\Delta u&=&f~\textrm{in }\Omega_i\\
\mu_i\p_n u|_{\Gamma_{i}}&=&-\mu_j\p_n u|_{\Gamma_{j}}
\quad \textrm{on }\overline{\Omega_{i}}\cap\overline{\Omega_{j}}
\qquad \textrm{if }1\le i< j\le 3
\end{eqnarray}

where $\p_n u|_{\Gamma_{i}}$ denotes the value of the normal derivative $\p_n u$ on the boundary $\Gamma_i$ of the domain $\Omega_i$.

By introducing the characteristic function $\chi_i$ of $\Omega_i$, that is,

\begin{equation}
\chi_i(x)=1\quad\textrm{if }x\in \Omega_i;\qquad
\chi_i(x)=0\quad\textrm{if }x\not\in \Omega_i
\end{equation}

we can easily rewrite (\ref{eqn:transm-1} 9.55 $\codered$) and (\ref{eqn:transm-2} 9.56 $\codered$)
to the weak form. Here we assume that $u=0$ on $\Gamma=\p\Omega$.

problem Transmission: For a given function $f$, find $u$ such that

\begin{eqnarray}
a(u,v)&=&\ell(f,v)\quad \textrm{for all }v\in H^1_0(\Omega)\\
a(u,v)=\int_{\Omega}\mu \nabla u\cdot \nabla v,\quad
\ell(f,v)=\int_{\Omega}fv\nonumber
\end{eqnarray}

where $\mu=\mu_1\chi_1+\mu_2\chi_2+\mu_3\chi_3$. Here we notice that $\mu$ become the discontinuous function.

With dissipation, and at the thermal equilibrium, the temperature equation is: $\codered$

This example explains the definition and manipulation of _region_, i.e. subdomains of the whole domain.
Consider this L-shaped domain with 3 diagonals as internal boundaries, defining 4 subdomains:

```freefem
// example using region keyword
// construct a mesh with 4 regions (sub-domains)
border a(t=0,1){x=t;y=0;};
border b(t=0,0.5){x=1;y=t;};
border c(t=0,0.5){x=1-t;y=0.5;};
border d(t=0.5,1){x=0.5;y=t;};
border e(t=0.5,1){x=1-t;y=1;};
border f(t=0,1){x=0;y=1-t;};
// internal boundary
border i1(t=0,0.5){x=t;y=1-t;};
border i2(t=0,0.5){x=t;y=t;};
border i3(t=0,0.5){x=1-t;y=t;};

mesh th = buildmesh (a(6) + b(4) + c(4) +d(4) + e(4) +
    f(6)+i1(6)+i2(6)+i3(6));
fespace Ph(th,P0); // constant discontinuous functions / element
fespace Vh(th,P1); // $P_1$ continuous functions / element

Ph reg=region; // defined the $P_0$ function associated to region number
plot(reg,fill=1,wait=1,value=1);
```

|Fig. 9.30: The function `:::freefem reg`|
|:----:|
|![region](images/region.png)|

|Fig. 9.31: The function `:::freefem nu`|
|:----:|
|![region_nu](images/region_nu.png)|

`:::freefem region` is a keyword of FreeFem++ which is in fact a variable depending of the current position (is not a function today, use `:::freefem Ph reg=region;` to set a function). This variable value returned is the number of the subdomain of the current position. This number is defined by `:::freefem buildmesh` which scans while building the mesh all its connected component.

So to get the number of a region containing a particular point one does:

```freefem
int nupper=reg(0.4,0.9); // get the region number of point (0.4,0.9)
int nlower=reg(0.9,0.1); // get the region number of point (0.4,0.1)
cout << " nlower " <<  nlower << ", nupper = " << nupper<< endl;
// defined the characteristics functions of upper and lower region
Ph nu=1+5*(region==nlower) + 10*(region==nupper);
plot(nu,fill=1,wait=1);
```

This is particularly useful to define discontinuous functions such as might occur
when one part of the domain is copper and the other one is iron, for example.

We this in mind we proceed to solve a Laplace equation with discontinuous coefficients
($\nu$ is 1, 6 and 11 below).

```freefem
Ph nu=1+5*(region==nlower) + 10*(region==nupper);
plot(nu,fill=1,wait=1);
problem lap(u,v) =   int2d(th)( nu*( dx(u)*dx(v)*dy(u)*dy(v) ))
                   + int2d(-1*v) + on(a,b,c,d,e,f,u=0);
plot(u);
```

|Fig. 9.32: The isovalue of the solution $u$|
|:----:|
|![region_u](images/region_u.png)|
