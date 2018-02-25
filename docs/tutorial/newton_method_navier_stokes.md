# Newton Method for the Steady Navier-Stokes equations

The problem is find the velocity field $\bm{u}=(u_i)_{i=1}^d$ and the pressure $p$ of a Flow  satisfying in the domain $\Omega \subset  \mathbb{R}^d (d=2,3)$:

\begin{eqnarray*}
  (\bm{u}\cdot\nabla) \bm{u}-\nu \Delta \bm{u}+\nabla p&=&0,\\ \nabla\cdot \bm{u}&=&0
\end{eqnarray*}

where $\nu$ is the viscosity of the fluid, $\nabla = (\p_i )_{i=1}^d $, the dot product is $\cdot$, and  $\Delta = \nabla\cdot\nabla$
with the some boundary conditions ( $\bm{u}$ is  given on $\Gamma$)

The weak form is 
find $\bm{u}, p $ such than for $\forall \bm{v}$ (zero on $\Gamma$), and $\forall  q$ 

\begin{equation}
  \int_\Omega  ((\bm{u}\cdot\nabla) \bm{u} ). \bm v + \nu \nabla \bm{u}:\nabla \bm{v} 
  - p \nabla\cdot \bm{v} - q \nabla\cdot \bm{u} = 0 
\end{equation}

The Newton Algorithm to solve nonlinear Problem is

Find  $u\in  V$ such that $F(u)=0$ where 
$ F : V \mapsto V $. 

1. choose $u_0\in \R^n $ , ;
2. for ( $i =0$; $i$ < niter; $i = i+1$) 
	1. solve $DF(u_i) w_i =  F(u_i)$;
	2. $u_{i+1} = u_i - w_i$;  

break  $|| w_i|| < \varepsilon$.


Where $DF(u)$ is the differential of $F$ at point  $u$, this is a linear application such that:

$
  F(u+\delta) = F(u) + DF(u) \delta + o(\delta) 
$

For Navier Stokes, $F$ and $DF$ are  :

\begin{eqnarray*}
F(\bm{u},p) =  \int_\Omega  &&((\bm{u}\cdot\nabla) \bm{u} ). \bm v + \nu \nabla \bm{u}:\nabla \bm{v} 
  - p \nabla\cdot \bm{v} - q \nabla\cdot \bm{u}\\
DF(\bm{u},p)(\bm{\delta u} ,\delta p)  =  \int_\Omega  &&((\bm{\delta u}\cdot\nabla) \bm{u} ). \bm v + ((\bm{u}\cdot\nabla) \bm{\delta u} ). \bm v \\
 &+& \nu \nabla \bm{\delta u}:\nabla \bm{v}  - \delta p \nabla\cdot \bm{v} - q \nabla\cdot \bm{\delta u}
\end{eqnarray*}

So the Newton algorithm become   
```freefem
\\ file NSNewton.edp 
for (n = 0; n < 15; n++) {
	solve Oseen([du1,du2,dp],[v1,v2,q]) =
          int2d(Th) (  nu*(Grad(du1,du2)'*Grad(v1,v2) )
                      + UgradV(du1,du2, u1, u2)'*[v1,v2]
                      + UgradV( u1, u2,du1,du2)'*[v1,v2]
                      - div(du1,du2)*q - div(v1,v2)*dp 
                      - 1e-8*dp*q // stabilization term 
                     )
        - int2d(Th) (  nu*(Grad(u1,u2)'*Grad(v1,v2) )
                      + UgradV(u1,u2, u1, u2)'*[v1,v2]
                      - div(u1,u2)*q - div(v1,v2)*p 
                     )
        + on(1,du1=0,du2=0) ;
      u1[] -= du1[];  u2[] -= du2[]; p[]  -= dp[];
      err= du1[].linfty + du2[].linfty + dp[].linfty;        
      if(err < eps) break;    
      if( n>3 && err > 10.) break; // Blowup ????
}
```


With the  operator: 

```freefem
macro Grad(u1,u2) [ dx(u1),dy(u1) , dx(u2),dy(u2) ]
macro UgradV(u1,u2,v1,v2) [ [u1,u2]'*[dx(v1),dy(v1)] ,
                            [u1,u2]'*[dx(v2),dy(v2)] ]
macro div(u1,u2)  (dx(u1)+dy(u2))
```

We build a computation mesh the exterior of a 2d cylinder. 

```freefem
real R = 5,L=15;
border cc(t=0,2*pi){ x=cos(t)/2;y=sin(t)/2;label=1;}
border ce(t=pi/2,3*pi/2) { x=cos(t)*R;y=sin(t)*R;label=1;}
border beb(tt=0,1) { real t=tt^1.2; x= t*L; y= -R; label = 1;}
border beu(tt=1,0) { real t=tt^1.2; x= t*L; y= R; label = 1;}
border beo(t=-R,R) {  x= L; y= t; label = 0;}
border bei(t=-R/4,R/4) {  x= L/2; y= t; label = 0;}
mesh Th=buildmesh(cc(-50)+ce(30)+beb(20)+beu(20)+beo(10)+bei(10));
plot(Th);

// Bounding box for the plot
func bb=[[-1,-2],[4,2]];

// FE Space Taylor Hood
fespace Xh(Th,P2); // For volicity 
fespace Mh(Th,P1); // For pressure 
Xh u1,u2,v1,v2,du1,du2,u1p,u2p;
Mh p,q,dp,pp;

// Initial guess with B.C. 
u1 = ( x^2+y^2) > 2;
u2=0;
```

Finally we use  trick to make continuation on the viscosity $\nu$, because the Newton method blowup 
owe start with the final viscosity $\nu$

```freefem
// Physical parameter
real nu= 1./50, nufinal=1/200. ,cnu=0.5;

// Stop test for Newton
real eps=1e-6;

verbosity=0;
while(1)  //  Loop on viscosity
{   int n;
    real err=0; // err on Newton algo ... 
    
      //... put the new the Newton  algo here
        
	if(err < eps)                                                               
	 { // converge decrease $\nu$ (more difficult)                                                                                                                           
	   plot([u1,u2],p,wait=1,cmm=" rey = " + 1./nu , coef=0.3,bb=bb);                                                           
	   if( nu == nufinal) break;                                                                                                
	   if( n < 4) cnu=cnu^1.5; // Fast converge => change faster                                                                
	   nu = max(nufinal, nu* cnu); // New vicosity                                                                              
	   u1p=u1;  u2p=u2;  pp=p; // Save correct solution ...                                                                    
	 }                                                                                                                          
	 else 
	 {   //  blowup increase $\nu$  (more simple)                                                                         
	   assert(cnu< 0.95); // The method  finally  blowup                                                                                               
	   nu = nu/cnu; //  get previous value of viscosity                                                                         
	   cnu= cnu^(1./1.5); // No conv. => change lower                                                                           
	   nu = nu* cnu; // New viscosity                                                                                           
	   cout << " restart nu = " << nu << " Rey= "<< 1./nu << "  (cnu = " << cnu << " ) \n";                                     
	   // Restore a correct solution ..                                                                                           
	   u1=u1p;                                                                                                                  
	   u2=u2p;                                                                                                                  
	   p=pp;                                                                                                                    
	 }     
}        
```

|Fig. 3.10: Mesh and the velocity and pressure at Reynolds 200|
|:----|
|![NSNewtonTh](images/NSNewtonTh.jpg)|
|![NSNewtonUP](images/NSNewtonUP.jpg)|

