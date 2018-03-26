
## Fluid/Structures Coupled Problem

This problem involves the Lamé system of elasticity
and the Stokes system for viscous fluids with velocity $\vec u$ and pressure $p$:

\begin{eqnarray*}
-\Delta \vec u +\vec\nabla p = 0, \,
\nabla\cdot \vec u = 0,\;\hbox{in}\;\Omega,\,
\vec u=\vec u_\Gamma\;\hbox{on}\;\Gamma=\p\Omega
\end{eqnarray*}

where $u_\Gamma$ is the velocity of the boundaries. The
force that the fluid applies to the boundaries is the normal stress

$$
\vec h =(\nabla\vec u +\nabla\vec u^T)\vec n -p\vec n
$$

Elastic solids subject to forces deform: a point in the solid, at (x,y) goes to (X,Y) after. When the displacement vector $\vec v=(v_1,v_2) = (X-x, Y-y)$ is small, Hooke's law relates the stress tensor $\sigma$ inside the solid to the deformation tensor $\epsilon$:

$$
\sigma_{ij} = \lambda \delta_{ij} \nabla.\vec v + 2\mu\epsilon_{ij},
\,
\epsilon_{ij} = {1\over 2}({\p v_i\over\p x_j} +
{\p v_j\over\p x_i} )
$$

where $\delta$ is the Kronecker symbol and where $\lambda, \mu$ are two constants describing the material mechanical properties in terms of the modulus of elasticity, and Young's modulus.

The equations of elasticity are naturally written in variational form for the displacement vector $v(x)\in V$ as

$$
\int_\Omega [2\mu\epsilon_{ij}(\vec v)\epsilon_{ij}(\vec w)
+\lambda \epsilon_{ii}(v)\epsilon_{jj}(\vec w)]
=\int_\Omega \vec g\cdot \vec w +\int_\Gamma \vec h\cdot \vec w,%\`{u}
\forall \vec w\in V
$$

The data are the gravity force $\vec g$ and the boundary stress $\vec h$.

 __Example 9.29__ fluidStruct.edp

In our example the Lamé system and the Stokes system are coupled by a common boundary on which the fluid stress creates a displacement of the boundary and hence changes the shape of the domain where the Stokes problem is integrated. The geometry is that of a vertical driven cavity with an elastic lid. The lid is a beam with weight so it will be deformed by its own weight and by the normal stress due to the fluid reaction. The cavity is the $10 \times 10$ square and the lid is a rectangle of height $l=2$.

A beam sits on a box full of fluid rotating because the left vertical side has velocity one. The beam is bent by its own weight, but the pressure of the fluid modifies the bending.

The bending displacement of the beam is given by (uu,vv) whose solution is given as follows.

```freefem
// Fluid-structure interaction for a weighting beam sitting on a
// square cavity filled with a fluid.

int bottombeam = 2; // label of bottombeam
border a(t=2,0)  { x=0; y=t ;label=1;}; // left beam
border b(t=0,10) { x=t; y=0 ;label=bottombeam;}; // bottom of beam
border c(t=0,2)  { x=10; y=t ;label=1;}; // rigth beam
border d(t=0,10) { x=10-t; y=2; label=3;}; // top beam
real E = 21.5;
real sigma = 0.29;
real mu = E/(2*(1+sigma));
real lambda = E*sigma/((1+sigma)*(1-2*sigma));
real gravity = -0.05;
mesh th = buildmesh( b(20)+c(5)+d(20)+a(5));
fespace Vh(th,P1);
Vh uu,w,vv,s,fluidforce=0;
cout << "lambda,mu,gravity ="<<lambda<< " " << mu << " " << gravity << endl;
// deformation of a beam under its own weight
solve bb([uu,vv],[w,s])  =
    int2d(th)(
                  lambda*div(w,s)*div(uu,vv)
                  +2.*mu*( epsilon(w,s)'*epsilon(uu,vv) )
             )
  + int2d(th) (-gravity*s)
  + on(1,uu=0,vv=0)
  + fluidforce[];
 ;

 plot([uu,vv],wait=1);
 mesh th1 = movemesh(th, [x+uu, y+vv]);
 plot(th1,wait=1);
```

Then Stokes equation for fluids ast low speed are solved in the box below the beam,
but the beam has deformed the box (see border h):

```freefem
//Stokes on square b,e,f,g driven cavite on left side g
border e(t=0,10) { x=t; y=-10; label= 1; }; // bottom
border f(t=0,10) { x=10; y=-10+t ; label= 1; }; // right
border g(t=0,10) { x=0; y=-t ;label= 2;}; // left
border h(t=0,10) { x=t; y=vv(t,0)*( t>=0.001 )*(t <= 9.999);
                    label=3;}; // top of cavity deformed

mesh sh = buildmesh(h(-20)+f(10)+e(10)+g(10));
plot(sh,wait=1);
```

 We use the Uzawa conjugate gradient to solve the Stokes problem like in example \refSec{Uzawa} $\codered$

```freefem
fespace Xh(sh,P2),Mh(sh,P1);
Xh u1,u2,v1,v2;
Mh p,q,ppp;

varf bx(u1,q) = int2d(sh)( -(dx(u1)*q));

varf by(u1,q) = int2d(sh)( -(dy(u1)*q));

varf Lap(u1,u2)= int2d(sh)(  dx(u1)*dx(u2) + dy(u1)*dy(u2) )
                    +  on(2,u1=1) +  on(1,3,u1=0)  ;

Xh bc1; bc1[] = Lap(0,Xh);
Xh brhs;

matrix A= Lap(Xh,Xh,solver=CG);
matrix Bx= bx(Xh,Mh);
matrix By= by(Xh,Mh);
Xh bcx=0,bcy=1;

func real[int] divup(real[int] & pp)
{
  int verb=verbosity;
   verbosity=0;
   brhs[]  = Bx'*pp; brhs[] += bc1[] .*bcx[];
   u1[] = A^-1*brhs[];
   brhs[]  = By'*pp; brhs[] += bc1[] .*bcy[];
   u2[] = A^-1*brhs[];
   ppp[] =   Bx*u1[];
   ppp[] +=  By*u2[];
   verbosity=verb;
   return ppp[] ;
};
```

do a loop on the two problems

```freefem
for(step=0;step<2;++step)
 {
   p=0;q=0;u1=0;v1=0;

   LinearCG(divup,p[],eps=1.e-3,nbiter=50);
   divup(p[]);
```

Now the beam will feel the stress constraint from the fluid:

```freefem
  Vh sigma11,sigma22,sigma12;
  Vh uu1=uu,vv1=vv;

  sigma11([x+uu,y+vv]) = (2*dx(u1)-p);
  sigma22([x+uu,y+vv]) = (2*dy(u2)-p);
  sigma12([x+uu,y+vv]) = (dx(u1)+dy(u2));
```

which comes as a boundary condition to the PDE of the beam:

```freefem
  solve bbst([uu,vv],[w,s],init=i)  =
     int2d(th)(
                  lambda*div(w,s)*div(uu,vv)
                  +2.*mu*( epsilon(w,s)'*epsilon(uu,vv) )
              )
  + int2d(th) (-gravity*s)
  + int1d(th,bottombeam)( -coef*(   sigma11*N.x*w + sigma22*N.y*s
                                   + sigma12*(N.y*w+N.x*s) )  )
  + on(1,uu=0,vv=0);
  plot([uu,vv],wait=1);
  real err = sqrt(int2d(th)( (uu-uu1)^2 + (vv-vv1)^2 ));
  cout <<  " Erreur L2 = " << err << "----------\n";
```

Notice that the matrix generated by bbst is reused (see `:::freefem init=i`).
Finally we deform the beam

```freefem
 th1 = movemesh(th, [x+0.2*uu, y+0.2*vv]);
 plot(th1,wait=1);
 } // end of loop
```

|Fig. 9.29: Fluid velocity and pressure (left) $\codered$ and displacement vector (center) $\codered$ of the structure and displaced geometry (right) $\codered$ in the fluid-structure interaction of a soft side and a driven cavity|
|:----:|
|![fluidstruct1](images/fluidstruct1.png)|
|![fluidstruct2](images/fluidstruct2.png)|
|![fluidstruct3](images/fluidstruct3.png)|
