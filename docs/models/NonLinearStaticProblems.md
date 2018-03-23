
## Nonlinear Static Problems

Here we propose to solve the following non-linear academic problem of minimization
of a functional $$J(u) = \int_\Omega \frac{1}{2} f(|\nabla u|^2) - u*b $$
where $u$ is function of $H^1_0(\Omega)$
and $f$ defined by
$$
f(x) = a*x + x-ln(1+x), \quad f'(x) = a+\frac{x}{1+x}, \quad f''(x) =  \frac{1}{(1+x)^2}
$$

### Newton-Raphson algorithm

Now, we solve the Euler problem $ \nabla J (u) = 0$
with Newton-Raphson algorithm, that is,
$$
u^{n+1} = u^n - ( \nabla^2 J (u^{n}))^{-1}*\nabla J(u^n)
$$

First we introduce the two variational form `:::freefem vdJ` and `:::freefem vhJ` to
compute respectively $ \nabla J$ and $ \nabla^2 J$

```freefem
// method of Newton-Raphson to solve dJ(u)=0; $\codered$
// $$ u^{n+1} = u^n - (\frac{\p dJ}{\p u_i})^{-1}*dJ(u^n) $$
// ---------------------------------------------
  Ph dalpha ; //to store  $2 f''( |\nabla u|^2) $  optimisation


 // the variational form of evaluate dJ = $ \nabla J$
 // --------------------------------------
 // dJ =  f'()*( dx(u)*dx(vh) + dy(u)*dy(vh)
  varf vdJ(uh,vh) =  int2d(Th)( alpha*( dx(u)*dx(vh) + dy(u)*dy(vh) ) - b*vh)
  + on(1,2,3,4, uh=0);


 // the variational form of evaluate ddJ   $= \nabla^2 J$
 // hJ(uh,vh) =    f'()*( dx(uh)*dx(vh) + dy(uh)*dy(vh)
 // + 2*f''()( dx(u)*dx(uh) + dy(u)*dy(uh) ) * (dx(u)*dx(vh) + dy(u)*dy(vh))
  varf vhJ(uh,vh) = int2d(Th)( alpha*( dx(uh)*dx(vh) + dy(uh)*dy(vh) )
   +  dalpha*( dx(u)*dx(vh) + dy(u)*dy(vh)  )*( dx(u)*dx(uh) + dy(u)*dy(uh) ) )
   + on(1,2,3,4, uh=0);

 // the Newton algorithm
  Vh v,w;
  u=0;
  for (int i=0;i<100;i++)
   {
    alpha =     df( dx(u)*dx(u) + dy(u)*dy(u) ) ; // optimization
    dalpha = 2*ddf( dx(u)*dx(u) + dy(u)*dy(u) ) ; // optimization
    v[]= vdJ(0,Vh); // $ v = \nabla J(u) $
    real res= v[]'*v[]; // the dot product
    cout << i <<  " residu^2 = " <<  res  << endl;
    if( res< 1e-12) break;
    matrix H= vhJ(Vh,Vh,factorize=1,solver=LU); //
    w[]=H^-1*v[];
    u[] -= w[];
   }
   plot (u,wait=1,cmm="solution with Newton-Raphson");
```

Remark: This example is in `:::freefem Newton.edp` file of `:::freefem examples++-tutorial` $\codered$ directory.
