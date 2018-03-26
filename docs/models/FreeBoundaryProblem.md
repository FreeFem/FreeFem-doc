
## Free Boundary Problem

The domain $\Omega$ is defined with:

```freefem
real L=10; //longueur du domaine
real h=2.1; // hauteur du bord gauche
real h1=0.35; // hauteur du bord droite

// maillage d'un tapeze
border a(t=0,L){x=t;y=0;}; // bottom:  $\Gamma_a$ \hfill
border b(t=0,h1){x=L;y=t;}; // right:  $\Gamma_b$ \hfill
border f(t=L,0){x=t;y=t*(h1-h)/L+h;}; // free surface:  $\Gamma_f$ \hfill
border d(t=h,0){x=0;y=t;}; // left:  $\Gamma_d$ \hfill

int n=4;
mesh Th=buildmesh (a(10*n)+b(6*n)+f(8*n)+d(3*n));
plot(Th,ps="dTh.eps");
```

|Fig. 9.33: The mesh of the domain $\Omega$|
|:----:|
|![dTh](images/dTh.png)|

The free boundary problem is:
Find $u$ and $\Omega$ such that:

$$
\left\{\begin{array}{cl}
\displaystyle - \Delta u = 0  & \mbox{in } \Omega\\
\displaystyle u = y         &\mbox{on } \Gamma_b \\
\displaystyle      {\p u  \over \p n} = 0   &\mbox{on } \Gamma_d \cup \Gamma_a \\
\displaystyle    {\p u  \over \p n} = { q\over K} n_x
          \mbox{and} {u = y}  &\mbox{on} \Gamma_ f
\end{array}\right.
$$

We use a fixed point method;
$\Omega^0 = \Omega$

in two step, fist we solve the classical following problem:

$$
\left\{\begin{array}{rll}
\displaystyle - \Delta u &= 0  & \mbox{in } \Omega^n\\
\displaystyle u &= y         &\mbox{on } \Gamma^n_b \\
\displaystyle      {\p u  \over \p n} &= 0   &\mbox{on } \Gamma^n_d \cup \Gamma^n_a\\
\displaystyle u &= y        &\mbox{on} \Gamma^n_ f
\end{array}\right.
$$

The variational formulation is:

find $u$ on $V=H^1(\Omega^n)$, such than $u=y$ on $\Gamma^n_b$ and $\Gamma^n_f$

$$
 \int_{\Omega^n}  \nabla u \nabla u' = 0,  \quad \forall u' \in V  \mbox{ with }  u' =0 \mbox{ on }
\Gamma^n_b \cup \Gamma^n_f
$$

and secondly to construct a domain deformation $\mathcal{F}(x,y)=[x,y-v(x,y)]$
where $v$ is solution of the following problem:

 $$
\left\{\begin{array}{rll}
\displaystyle - \Delta v &= 0  & \mbox{in } \Omega^n\\
\displaystyle v  &= 0         &\mbox{on } \Gamma^n_a \\
\displaystyle      {\p v \over \p n} &= 0   &\mbox{on } \Gamma^n_b \cup \Gamma^n_d \\
\displaystyle    {\p v  \over \p n}  &=  \displaystyle {\p u  \over \p n} - { q\over K} n_x
            &\mbox{on } \Gamma^n_ f
\end{array}\right. $$

The variational formulation is:

find $v$ on $V$, such than $v=0$ on $\Gamma^n_a$

$$
 \int_{\Omega^n}  \nabla v \nabla v' = \int_{\Gamma_f^n}  ({\p u  \over \p n} - { q\over K} n_x )v',  \quad \forall v' \in V  \mbox{ with }  v' =0 \mbox{ on }
\Gamma^n_a
$$

Finally the new domain $\Omega^{n+1} = \mathcal{F}(\Omega^n)$

 __Example 9.30__ freeboundary.edp

The FreeFem++ implementation is:

```freefem
real q=0.02; //flux entrant
real K=0.5; //permeabilit\'{e}

fespace Vh(Th,P1);
int j=0;

Vh u,v,uu,vv;

problem Pu(u,uu,solver=CG) = int2d(Th)( dx(u)*dx(uu)+dy(u)*dy(uu))
  + on(b,f,u=y) ;

problem Pv(v,vv,solver=CG) = int2d(Th)( dx(v)*dx(vv)+dy(v)*dy(vv))
  +  on (a, v=0) + int1d(Th,f)(vv*((q/K)*N.y- (dx(u)*N.x+dy(u)*N.y)));


real errv=1;
real erradap=0.001;
verbosity=1;
while(errv>1e-6)
{
  j++;
  Pu;
  Pv;
  plot(Th,u,v ,wait=0);
  errv=int1d(Th,f)(v*v);
   real coef=1;

//
  real mintcc = checkmovemesh(Th,[x,y])/5.;
  real mint = checkmovemesh(Th,[x,y-v*coef]);

  if (mint<mintcc ||  j%10==0) { // mesh to bad => remeshing
    Th=adaptmesh(Th,u,err=erradap ) ;
    mintcc = checkmovemesh(Th,[x,y])/5.;
  }

  while (1)
  {
    real mint = checkmovemesh(Th,[x,y-v*coef]);

    if (mint>mintcc) break;

    cout << " min |T]  " << mint << endl;
    coef /= 1.5;
  }

  Th=movemesh(Th,[x,y-coef*v]); // calcul de la deformation
  cout << "\n\n"<<j <<"------------ errv = " << errv << "\n\n";

}
plot(Th,ps="d_Thf.eps");
plot(u,wait=1,ps="d_u.eps");
```

|Fig. 9.34: The final solution on the new domain $\Omega^{72}$|
|:----:|
|![d_u](images/d_u.png)|

|Fig. 9.35: The adapted mesh of the domain $\Omega^{72}$|
|:----:|
|![d_Thf](images/d_Thf.png)|
