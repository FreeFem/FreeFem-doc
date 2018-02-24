# The System of Stokes for Fluids

In the case of a flow invariant with respect to the third coordinate
(two-dimensional flow), flows at low Reynolds number (for instance
micro-organisms) satisfy,

$\codered$ **check vector**

\begin{eqnarray*}&&
    \vec -\Delta u + \n p =0
    \cr&&
    \n\cdot \vec u =0
\end{eqnarray*}

where $\vec u=(u_1,u_2)$ is the fluid velocity and $p$ its pressure.

The driven cavity is a standard test. It is a box full of liquid with its lid moving horizontally
at speed one.  The pressure and the velocity must be discretized in compatible fintie
element spaces for the LBB conditions to be satisfied:

\[
    \sup_{p\in P_h}\frac{(\vec u,\n p)}{|p|}\geq \beta|\vec u|~~~\forall \vec u\in U_h
\]

```freefem
//file stokes.edp
int n=3;
mesh Th=square(10*n,10*n);
fespace Uh(Th,P1b); Uh u,v,uu,vv;
fespace Ph(Th,P1);  Ph p,pp;

solve stokes([u,v,p],[uu,vv,pp]) =
    int2d(Th)(dx(u)*dx(uu)+dy(u)*dy(uu) + dx(v)*dx(vv)+ dy(v)*dy(vv)
              + dx(p)*uu + dy(p)*vv + pp*(dx(u)+dy(v))
              -\bf 1e-10*p*pp\tt)
            + on(1,2,4,u=0,v=0) + on(3,u=1,v=0);
plot([u,v],p,wait=1);
```

Remark, we add a stabilization term $\bf{-10e-10*p*pp}$ to fixe the constant part of the pressure.

Results are shown on fig. 3.8


|Fig. 3.8: Solution of Stokes' equations for the driven cavity problem, showing the velocity field and the pressure level lines.|
|:----|
|![Stokes](images/stokes.svg)|
