\subsection{The System of elasticity}

\paragraph{Elasticity}

Solid objects deform under the action of applied forces:
a point in the solid, originally at $(x,y,z)$ will come to $(X,Y,Z)$
after some time; the vector $\mathbf{u}=(u_1,u_2,u_3) = (X-x, Y-y, Z-z)$ is called
the displacement. When the displacement is small and the solid is
elastic, Hooke's law gives a relationship between the stress tensor
$\sigma(u)=(\sigma_{ij}(u) )$ and the strain tensor $\epsilon(u)=\epsilon_{ij}(u)$
$$
\sigma_{ij}(u) = \lambda \delta_{ij} \nabla.\mathbf{u}+ 2\mu\epsilon_{ij}(u),
$$
where the Kronecker symbol $\delta_{ij} = 1$ if $i=j$, $0$ otherwise, with
$$\epsilon_{ij}(u) = {1\over 2}({\p u_i\over\p x_j} +
{\p u_j\over\p x_i} ),
$$
and where $\lambda, \mu$ are two constants that describe the
mechanical properties of the solid, and are themselves related to the
better known constants $E$, Young's modulus, and $\nu$, Poisson's ratio:
$$ \mu = {E\over 2( 1+\nu)}, \quad \lambda = {E\nu\over (1+\nu)(1-2\nu)}.
$$


 \paragraph{Lam\'e's system}

Let us consider a beam with axis $Oz$ and with perpendicular section
$\Omega$. The components along $x$ and $y$ of the strain ${\bf u}(x)$
in a section $\Omega$ subject to forces ${\bf f}$ perpendicular to the
axis are governed by \\
$$
  -\mu \Delta {\bf u} - (\mu+\lambda)  \nabla (\nabla .{\bf u})={\bf f}~~\hbox{in}~~\Omega,
$$
where $\lambda ,\mu  $ are the Lam\'{e} coefficients introduced above.

Remark, we do not used this equation because the associated  variationnal
form does not give the right boundary condition, we simply use
$$
  - div( \sigma ) = \mathbf{f}  \quad  \mbox{in} \Omega
$$
where the corresponding variationnal form is:
$$
 \int_{\Omega} \sigma(u) : \epsilon(\mathbf{v})\;dx - \int_{\Omega}  \mathbf{v} f \;dx =0;
$$
where $:$  denote the tensor scalar product,   i.e. $ a: b = \sum_{i,j}  a_{ij}b_{ij}$.

So the variationnal form can be written as :
$$
 \int_{\Omega} \lambda \nabla.u   \nabla.v  + 2 \mu \epsilon(\mathbf{u}):\epsilon(\mathbf{v}) \; dx - \int_{\Omega}  \mathbf{v} f  \;dx  =0;
$$
\paragraph{Example}  Consider  elastic plate with the undeformed rectangle shape
$[0,20]\times [-1,1]$.
The body force is the gravity force $\vec f$ and the
boundary force $\vec g$ is zero on lower, upper and right sides.
The left vertical sides of the beam is fixed.
The boundary conditions are
\begin{eqnarray*}
     \sigma . {\bf n}  &=& g = 0    ~~\hbox{on}~~\Gamma_1, \Gamma_4, \Gamma_3, \\
      {\bf u} &=& \mathbf{0} ~~\hbox{on}~~\Gamma_2
 \end{eqnarray*}
Here ${\bf u}=(u,v) $ has two components.\bigskip

The above two equations are strongly coupled by their mixed
derivatives, and thus any iterative solution on each of the
components is risky. One should rather use FreeFem++'s system
approach and write:

\begin{example}[lame.edp]\label{lame.edp}
\bFF
// file lame.edp
mesh Th=square(10,10,[20*x,2*y-1]);
fespace Vh(Th,P2);
Vh u,v,uu,vv;
real sqrt2=sqrt(2.);
macro epsilon(u1,u2)  [dx(u1),dy(u2),(dy(u1)+dx(u2))/sqrt2] // EOM  \index{macro!with parameter}
//the sqrt2 is because we want: epsilon(u1,u2)'* epsilon(v1,v2) $==  \epsilon(\bm{u}): \epsilon(\bm{v})$
macro div(u,v) ( dx(u)+dy(v) ) // EOM


real E = 21e5, nu = 0.28, mu= E/(2*(1+nu));
real lambda = E*nu/((1+nu)*(1-2*nu)), f = -1; //

solve lame([u,v],[uu,vv])= int2d(Th)(
        lambda*div(u,v)*div(uu,vv)
        +2.*mu*( epsilon(u,v)'*epsilon(uu,vv) ) )	
        - int2d(Th)(f*vv)
        + on(4,u=0,v=0);
real coef=100;
plot([u,v],wait=1,ps="lamevect.eps",coef=coef);

mesh th1 = movemesh(Th, [x+u*coef, y+v*coef]);
plot(th1,wait=1,ps="lamedeform.eps");
real dxmin  = u[].min;
real dymin  = v[].min;

cout << " - dep.  max   x = "<< dxmin<< " y=" << dymin << endl;
cout << "   dep.  (20,0)  = " << u(20,0) << " " << v(20,0) << endl;
\eFF
\end{example}

The numerical results are shown on figure \ref{figlame} and the output is:
\bFF
 -- square mesh : nb vertices  =121 ,  nb triangles = 200 ,  nb boundary edges 40
 -- Solve :           min -0.00174137  max 0.00174105
          min -0.0263154  max 1.47016e-29
 - dep.  max   x = -0.00174137 y=-0.0263154
   dep.  (20,0)  = -1.8096e-07 -0.0263154
times: compile 0.010219s, execution 1.5827s
\eFF



\begin{figure}[hbtp]
\begin{center}
\includegraphics[width=15cm]{lamevect}\\
\includegraphics[width=15cm]{lamedeform}

\caption{\label{figlame} Solution of Lam\'e's equations for elasticity for a 2D beam deflected by its
own weight and clamped by its left vertical side; result are shown with a amplification factor equal to  100.
{\em Remark: the size of the arrow  is automatically bound, but the color gives the real length}}
\end{center}
\end{figure}

\subsection{The System of Stokes for Fluids}

In the case of a flow invariant with respect to the third coordinate
(two-dimensional flow), flows at low Reynolds number (for instance
micro-organisms) satisfy,
\begin{eqnarray*}&&
    \vec -\Delta u + \n p =0
    \cr&&
    \n\cdot \vec u =0
\end{eqnarray*}
where $\vec u=(u_1,u_2)$ is the fluid velocity and $p$ its pressure.
\\
The driven cavity is a standard test. It is a box full of liquid with its lid moving horizontally
at speed one.  The pressure and the velocity must be discretized in compatible fintie
element spaces for the LBB conditions to be satisfied:
\[
    \sup_{p\in P_h}\frac{(\vec u,\n p)}{|p|}\geq \beta|\vec u|~~~\forall \vec u\in U_h
\]
\bFF
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
\eFF

Remark, we add a stabilization term {\bf{-10e-10*p*pp}} to fixe the constant part of the pressure.
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=8cm]{stokes}

\caption{\label{figstokes} Solution of Stokes' equations for the driven cavity problem, showing the
velocity field and the pressure level lines.}
\end{center}
\end{figure}

Results are shown on figure \ref{figstokes}

\subsection{A Projection Algorithm for the Navier-Stokes equations }
\paragraph{Summary}\emph{Fluid flows require good algorithms and good triangultions. We show
here an example of a complex algorithm and or first example of \x{mesh adaptation}}.
\\\\
An incompressible viscous fluid satisfies:
$$ \p _t u + u\cdot\nabla u + \nabla p - \nu\Delta u = 0,\quad  \nabla\cdot u=0
\quad  \hbox{ in } \Omega\times ]0,T[,
$$
$$ u|_{t=0} = u^0,\quad  u|_\Gamma = u_\Gamma.
$$
A possible algorithm, proposed by Chorin, is
$$ {1\over \delta t}[u^{m+1} - u^moX^m] + \nabla p^m -\nu\Delta u^m= 0,\quad  u|_\Gamma = u_\Gamma, \nu \p_n u|_{\Gamma_{out}}=0 
  $$
$$ -\Delta p^{m+1} = -\nabla\cdot  u^moX^m, \quad  \p _n p^{m+1} = 0 \mbox{ on } \Gamma , p^{m+1} = 0 \mbox{ on } \Gamma_{out}
$$
where $uoX(x) = u(x-u(x)\delta t)$ since $\p _t u + u\cdot\nabla
u $ is approximated by the method of characteristics, as in the previous section.
\\\\
We use the s Chorin's algorithm with free boundary condition  at outlet (i.e.  $p=0,\nu \p_n u = 0$),  to compute a correction, q,
to the pressure. 
\[
    -\Delta q= \n\cdot\vec u ,  \quad q=0 \mbox{ on } \Gamma_{out}
\]
and define
\[
    u^{m+1}=\tilde u + P \n q\delta t,~~~p^{m+1}=p^m-q
\]
where $\tilde u$ is the $(u^{m+1},v^{m+1})$ of Chorin's algorithm,
and  where $P$ is the  $L^2$ projection with mass lumping ( a sparse matrix). 

\paragraph{The backward facing step}

The geometry is that of a channel with a backward facing step so that
the inflow section is smaller than the outflow section. This geometry
produces a fluid recirculation zone that must be captured correctly.

This can only be done if the triangulation is sufficiently fine, or
well adapted to the flow.

Remark (FH), The are a technical difficulty is the example, is the output B.C., here 
we put $p=0$ and $ \nu \p_n u = 0$.  

\begin{example}[NSprojection.edp]\index{adaptmesh}\index{mesh!adaptation}
\bFF
// file NSprojection.edp
// Version july 2014, 
// FH. Change B.C on u on outpout , more simple .. 
// ............
verbosity=0;
border a0(t=1,0){ x=-2;      y=t;      label=1;}
border a1(t=-2,0){ x=t;    y=0;        label=2;}
border a2(t=0,-0.5){ x=0;      y=t;       label=2;}
border a3(t=0,1){ x=18*t^1.2;  y=-0.5;       label=2;}
border a4(t=-0.5,1){ x=18;     y=t;   label=3;}
border a5(t=1,0){ x=-2+20*t; y=1;        label=4;}
int n=1;
mesh Th= buildmesh(a0(3*n)+a1(20*n)+a2(10*n)+a3(150*n)+a4(5*n)+a5(100*n));

plot(Th);
fespace Vh(Th,P1);
real nu = 0.0025, dt = 0.2; // Reynolds=200

Vh w,u = 0, v =0, p = 0, q=0;

real epsv = 1e-6, epsu = 1e-6, epsp = 1e-6;// Eps CG ..

 // def of Matrix dtMx and dtMy
matrix dtM1x,dtM1y;
 
macro  BuildMat()
 { /* for memory managenemt */
   varf vM(unused,v) = int2d(Th)(v) ;
   varf vdx(u,v) = int2d(Th)(v*dx(u)*dt) ;
   varf vdy(u,v) = int2d(Th)(v*dy(u)*dt) ;

   real[int] Mlump = vM(0,Vh); 
   real[int] one(Vh.ndof); one = 1;  
   real[int] M1 =  one ./ Mlump; 
   matrix dM1 = M1;
   matrix Mdx = vdx(Vh,Vh);
   matrix Mdy = vdy(Vh,Vh);
   dtM1x = dM1*Mdx;
   dtM1y = dM1*Mdy; 
 }// EOF \hfilll
 
BuildMat

real err=1, outflux=1;
for(int n=0;n<300;n++)
 {	
  Vh uold = u,  vold = v, pold=p;
  
  solve pb4u(u,w,init=n,solver=CG,eps=epsu)
        =int2d(Th)(u*w/dt +nu*(dx(u)*dx(w)+dy(u)*dy(w)))
        -int2d(Th)((convect([uold,vold],-dt,uold)/dt-dx(p))*w)
        + on(1,u = 4*y*(1-y)) + on(2,4,u = 0) ;// Neuman on $Gamma_3$
       
  plot(u);

  solve pb4v(v,w,init=n,solver=CG,eps=epsv)
        = int2d(Th)(v*w/dt +nu*(dx(v)*dx(w)+dy(v)*dy(w)))
        -int2d(Th)((convect([uold,vold],-dt,vold)/dt-dy(p))*w)
        +on(1,2,3,4,v = 0);

 solve pb4p(q,w,solver=CG,init=n,eps=epsp) = int2d(Th)(dx(q)*dx(w)+dy(q)*dy(w))
    - int2d(Th)((dx(u)+ dy(v))*w/dt)+ on(3,q=0);

 // to have absolute epsilon in CG algorithm. 
  epsv = -abs(epsv);
  epsu = -abs(epsu);
  epsp = -abs(epsp);

  p = pold-q;
  u[] += dtM1x*q[];
  v[] += dtM1y*q[];
 
  
  if(n%50==49){
    Th = adaptmesh(Th,[u,v],q,err=0.04,nbvx=100000);
    plot(Th, wait=true);
    BuildMat // rebuild mat.   
 }
 
  err = sqrt(int2d(Th)(square(u-uold)+square(v-vold))/Th.area) ;
  outflux = int1d(Th)( [u,v]'*[N.x,N.y]) ;
  cout << " iter " << n << " Err L2 = " << err 
       <<  " flux sortant = "<< outflux << endl; 
  if(err < 1e-3) break;
}
assert(abs(outflux)< 2e-3); // verifaction ... 
plot(p,wait=1,ps="NSprojP.eps");
plot(u,wait=1,ps="NSprojU.eps");
\eFF
\end{example}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=15cm]{NSprojTh}\\
\includegraphics[width=15cm]{NSprojP}\\
\includegraphics[width=15cm]{NSprojU}

\caption{\label{figNSproj} Rannacher's projection algorithm: result on an adapted mesh (top) showing
the pressure (middle) and the horizontal velocity $u$ at Reynolds 400.}
\end{center}
\end{figure}

We show in figure \ref{figNSproj} the numerical results obtained for a Reynolds
number of 400 where mesh adaptation is done after 50 iterations on the first mesh.

\subsection{Newton Method for the Steady Navier-Stokes equations}

The problem is find the velocity field $\bm{u}=(u_i)_{i=1}^d$ and the pressure $p$ of a Flow 
satisfying in the domain $\Omega \subset  \mathbb{R}^d (d=2,3)$:
\Blue{
  \begin{eqnarray*}
    (\bm{u}\cdot\nabla) \bm{u}-\nu \Delta \bm{u}+\nabla p&=&0,\\ \nabla\cdot \bm{u}&=&0
  \end{eqnarray*}
}
where $\nu$ is the viscosity of the fluid, $\nabla = (\p_i )_{i=1}^d $, the dot product is $\cdot$, and  $\Delta = \nabla\cdot\nabla$
with the some boundary conditions ( $\bm{u}$ is  given on $\Gamma$)

\bigskip

The weak form is 
find $\bm{u}, p $ such than for $\forall \bm{v}$ (zero on $\Gamma$), and $\forall  q$ 
\Blue{
  \begin{equation}
    \int_\Omega  ((\bm{u}\cdot\nabla) \bm{u} ). \bm v + \nu \nabla \bm{u}:\nabla \bm{v} 
    - p \nabla\cdot \bm{v} - q \nabla\cdot \bm{u} = 0 
  \end{equation}
}


The Newton Algorithm to solve nonlinear Problem is

Find  $u\in  V$ such that $F(u)=0$ where 
$ F : V \mapsto V $. 

\begin{enumerate} \index{Newton Algorithm}
\item choose $u_0\in \R^n $ , ;
\item for ( $i =0$; $i$ < niter; $i = i+1$) 

\begin{enumerate}
\item solve $DF(u_i) w_i =  F(u_i)$;
\item $u_{i+1} = u_i - w_i$;  
\end{enumerate}
break  $|| w_i|| < \varepsilon$.
\end{enumerate}

Where $DF(u)$ is the differential of $F$ at point  $u$, this is a linear application such that:
{$
  F(u+\delta) = F(u) + DF(u) \delta + o(\delta) 
$}


For Navier Stokes, $F$ and $DF$ are  :
\Blue{\small
\begin{eqnarray*}F(\bm{u},p) =  \int_\Omega  &&((\bm{u}\cdot\nabla) \bm{u} ). \bm v + \nu \nabla \bm{u}:\nabla \bm{v} 
  - p \nabla\cdot \bm{v} - q \nabla\cdot \bm{u}\\
DF(\bm{u},p)(\bm{\delta u} ,\delta p)  =  \int_\Omega  &&((\bm{\delta u}\cdot\nabla) \bm{u} ). \bm v + ((\bm{u}\cdot\nabla) \bm{\delta u} ). \bm v \\
 &+& \nu \nabla \bm{\delta u}:\nabla \bm{v}  - \delta p \nabla\cdot \bm{v} - q \nabla\cdot \bm{\delta u}
\end{eqnarray*}}



So the Newton algorithm become   
\begin{example}[NSNewton.edp] 
\bFF
...   
	for( n=0;n< 15;n++)
	{ solve Oseen([du1,du2,dp],[v1,v2,q]) =
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
      if( n>3 && err > 10.) break; //  blowup ???? 
    }
\eFF


With the  operator: 
\bFF
macro Grad(u1,u2) [ dx(u1),dy(u1) , dx(u2),dy(u2) ]// 
macro UgradV(u1,u2,v1,v2) [ [u1,u2]'*[dx(v1),dy(v1)] ,
                            [u1,u2]'*[dx(v2),dy(v2)] ]// 
macro div(u1,u2)  (dx(u1)+dy(u2))//
\eFF

We build a computation mesh the exterior of a 2d cylinder. 

\bFF
real R = 5,L=15;
border cc(t=0,2*pi){ x=cos(t)/2;y=sin(t)/2;label=1;}
border ce(t=pi/2,3*pi/2) { x=cos(t)*R;y=sin(t)*R;label=1;}
border beb(tt=0,1) { real t=tt^1.2; x= t*L; y= -R; label = 1;}
border beu(tt=1,0) { real t=tt^1.2; x= t*L; y= R; label = 1;}
border beo(t=-R,R) {  x= L; y= t; label = 0;}
border bei(t=-R/4,R/4) {  x= L/2; y= t; label = 0;}
mesh Th=buildmesh(cc(-50)+ce(30)+beb(20)+beu(20)+beo(10)+bei(10));
plot(Th);

// bounding box for the plot
func bb=[[-1,-2],[4,2]];

/  FE Space Taylor Hood
fespace Xh(Th,P2);// for volicity 
fespace Mh(Th,P1);// for pressure 
Xh u1,u2,v1,v2,du1,du2,u1p,u2p;
Mh p,q,dp,pp;

// intial guess with B.C. 
u1 = ( x^2+y^2) > 2;
u2=0;
\eFF

Finally we use  trick to make continuation on the viscosity $\nu$, because the Newton method blowup 
owe start with the final viscosity $\nu$
\bFF
//  Physical parameter
real nu= 1./50, nufinal=1/200. ,cnu=0.5;

// stop test for Newton
real eps=1e-6;

verbosity=0;
while(1)  //  Loop on viscosity
{   int n;
    real err=0; // err on Newton algo ... 
    
      ... put the new the Newton  algo here
        
	if(err < eps)                                                               
	 { // converge decrease $\nu$ (more difficult)                                                                                                                           
	   plot([u1,u2],p,wait=1,cmm=" rey = " + 1./nu , coef=0.3,bb=bb);                                                           
	   if( nu == nufinal) break;                                                                                                
	   if( n < 4) cnu=cnu^1.5; // fast converge => change faster                                                                
	   nu = max(nufinal, nu* cnu); // new vicosity                                                                              
	   u1p=u1;  u2p=u2;  pp=p; //  save correct solution ...                                                                    
	 }                                                                                                                          
	 else 
	 {   //  blowup increase $\nu$  (more simple)                                                                         
	   assert(cnu< 0.95); //  the method  finally  blowup                                                                                               
	   nu = nu/cnu; //  get previous value of viscosity                                                                         
	   cnu= cnu^(1./1.5); // no conv. => change lower                                                                           
	   nu = nu* cnu;  // new viscosity                                                                                           
	   cout << " restart nu = " << nu << " Rey= "<< 1./nu << "  (cnu = " << cnu << " ) \n";                                     
	   // restore a correct solution ..                                                                                           
	   u1=u1p;                                                                                                                  
	   u2=u2p;                                                                                                                  
	   p=pp;                                                                                                                    
	 }     
}        
\eFF
\end{example}

\begin{figure}[htbp]
\begin{center}
\includegraphics[width=8cm]{NSNewtonTh}
\includegraphics[width=8cm]{NSNewtonUP}
\caption{\label{NSNewton} Mesh and the velocity and pressure at Reynolds  $ 200$ }
\end{center}
\end{figure}


%\paragraph{Gag's or Bug with Stokes Equation, with Boundary condition}
%
%
%Denote: \Blue{$\Gamma$} the boundary, \Blue{$\bm{n}$} the unit exterior normal,  
%\Blue{$\varepsilon(\bm{u}) =  ( ^t \nabla \bm{u}+ \nabla \bm{u})/2 $}, \Blue{$I_d$} the Identity matrix.
%
%\begin{enumerate}
%\item With BC, the incompressibility imply
%\Red{$ \int_\Gamma \bm{u}.\bm{n} =0 $} at discrete level. 
%\item 
% The force in the fluid on surface $\Sigma$ is \Blue{$ \int_\Sigma (\nu\; \varepsilon(\bm{u})  -  p I_d) \bm{n}  $}. The Weak form is  
%$\forall \bm{v}, q$ :
%\[  \int_\Omega  \nu \nabla \bm{u}:\nabla \bm{v} 
%  - p \nabla\cdot \bm{v} - q \nabla\cdot \bm{u} = {\int_\Gamma \;{}^t\bm{n}(\nu \nabla \bm{u} -  p I_d )   \bm{v} }
%  \]
%    Or with a more physical formulation (with the true strain $\varepsilon$):
%\[  \int_\Omega  \nu\; \varepsilon( \bm{u}):\varepsilon( \bm{v}) 
%  - p \nabla\cdot \bm{v} - q \nabla\cdot \bm{u} = {\int_\Gamma \;{}^t\bm{n}(\nu \varepsilon (\bm{u}) -  p I_d )   \bm{v} }
%  \]  
%\end{enumerate}
%
%\begin{figure}[htbp]
%\begin{center}
%\includegraphics[width=8cm]{NSnewtonU}
%\includegraphics[width=8cm]{NSnewtonE}
%
%\caption{\label{figconvect} The rotated hill after one revolution, left with Characteristics-Galerkin,
%on the right with Discontinuous $P_1$ Galerkin FEM.}
%\end{center}
%\end{figure}
%



\subsection{A Large Fluid Problem}
A friend of one of us in Auroville-India was building a ramp to access an air conditioned room. As I was visiting the construction site he told me that he expected to cool air escaping by the door to the room to slide down the ramp and refrigerate the feet of the coming visitors.  I told him "no way" and decided to check numerically.  The results are on the front page of this book.
\\
The fluid velocity and pressure are solution of the Navier-Stokes equations with varying density function of the temperature.
\\
The geometry is trapezoidal with prescribed inflow made of cool air at the bottom and warm air above and so are the initial conditions; there is free outflow, slip velocity at the top (artificial) boundary and no-slip at the bottom.  However the Navier-Stokes cum temperature equations have a RANS $k-\epsilon$ model and a Boussinesq approximation for the buoyancy. This comes to
\begin{eqnarray}&&
\p_t\theta+u\n\theta-\n\cdot(\kappa_T^m\n\theta)=0
\cr&&
\p_t u +u\n u -\n\cdot(\mu_T\n u) +\n p+ e(\theta-\theta_0)\vec e_2,~~\n\cdot u=0
\cr&&
\mu_T = c_\mu\frac{k^2}\epsilon,~~\kappa_T=\kappa\mu_T
\cr&&
\p_t k + u\n k + \epsilon  -\n\cdot(\mu_T\n k)  = \frac{\mu_T}2|\n u+\n u^T|^2
\cr&&
\p_t\epsilon+u\n\epsilon + c_2\frac{\epsilon^2} k -\frac{c_\epsilon}{c_\mu}\n\cdot (\mu_T\n\epsilon)= \frac{c_1}2  k|\n u+\n u^T|^2=0
\end{eqnarray}
We use a time discretization which preserves positivity and uses the method of characteristics ($X^m(x)\approx  x-u^m(x)\delta t$)
\begin{eqnarray}&&
\frac 1{\delta t}(\theta^{m+1}-\theta^m \circ X^m)-\n\cdot(\kappa_T^m\n\theta^{m+1})=0
\cr&&
\frac1{\delta t}(u^{m+1}-u^m \circ X^m) -\n\cdot(\mu_T^m\n u^{m+1}) +\n p^{m+1}+ e(\theta^{m+1}-\theta_0)\vec e_2
,~~\n\cdot u^{m+1}=0
\cr&&
\frac1{\delta t}(k^{m+1}-k^m \circ X^m) + k^{m+1}\frac{\epsilon^m}{k^m}  -\n\cdot(\mu_T^m\n k^{m+1})  = \frac{\mu_T^m}2|\n u^m+{\n u^m}^T|^2
\cr&&
\frac1{\delta t}(\epsilon^{m+1}-\epsilon^m \circ X^m) + c_2\epsilon^{m+1}\frac{\epsilon^m} {k^m} -\frac{c_\epsilon}{c_\mu}\n\dot(\mu_T^m\n\epsilon^{m+1})= \frac{c_1}2  k^m|\n u^m+{\n u^m}^T|^2
\cr&&
\mu_T ^{m+1}= c_\mu\frac{{k^{m+1}}^2}{\epsilon^{m+1}},~~\kappa_T^{m+1}=\kappa\mu_T^{m+1}
\end{eqnarray}
In variational form and with appropriated boundary conditions the problem is:

\bFF
real L=6;
border aa(t=0,1){x=t; y=0 ;}
border bb(t=0,14){x=1+t; y= - 0.1*t ;}
border cc(t=-1.4,L){x=15; y=t ;}
border dd(t=15,0){x= t ; y = L;}
border ee(t=L,0.5){ x=0; y=t ;}
border ff(t=0.5,0){ x=0; y=t ;}
int n=8;
mesh Th=buildmesh(aa(n)+bb(9*n) + cc(4*n) + dd(10*n)+ee(6*n) + ff(n));
real s0=clock();

fespace Vh2(Th,P1b); // velocity space
fespace Vh(Th,P1); // pressure space
fespace V0h(Th,P0); // for gradients
Vh2 u2,v2,up1=0,up2=0;
Vh2 u1,v1;
Vh  u1x=0,u1y,u2x,u2y, vv;

real reylnods=500;
//cout << " Enter the reynolds number :"; cin >> reylnods;
assert(reylnods>1 && reylnods < 100000);
up1=0;
up2=0;
func g=(x)*(1-x)*4;  // inflow
Vh p=0,q, temp1,temp=35, k=0.001,k1,ep=0.0001,ep1;
V0h muT=1,prodk,prode, kappa=0.25e-4, stress;
real alpha=0, eee=9.81/303, c1m = 1.3/0.09 ;
real  nu=1, numu=nu/sqrt( 0.09), nuep=pow(nu,1.5)/4.1;
int i=0,iter=0;
real dt=0;
problem TEMPER(temp,q) = // temperature equation
    int2d(Th)(
             alpha*temp*q + kappa * ( dx(temp)*dx(q) + dy(temp)*dy(q) ))
//   + int1d(Th,aa,bb)(temp*q* 0.1)
  + int2d(Th) ( -alpha*convect([up1,up2],-dt,temp1)*q )
   + on(ff,temp=25)
  + on(aa,bb,temp=35) ;

problem kine(k,q)=  // get the kinetic turbulent energy
    int2d(Th)(
             (ep1/k1+alpha)*k*q + muT * ( dx(k)*dx(q) + dy(k)*dy(q) ))
//   + int1d(Th,aa,bb)(temp*q*0.1)
  + int2d(Th) ( prodk*q-alpha*convect([up1,up2],-dt,k1)*q )
   + on(ff,k=0.0001)  + on(aa,bb,k=numu*stress) ;

 problem viscturb(ep,q)= // get the rate of turbulent viscous energy
    int2d(Th)(
             (1.92*ep1/k1+alpha)*ep*q + c1m*muT * ( dx(ep)*dx(q) + dy(ep)*dy(q) ))
//   + int1d(Th,aa,bb)(temp*q*0.1)
  + int2d(Th) ( prode*q-alpha*convect([up1,up2],-dt,ep1)*q )
   + on(ff,ep= 0.0001) + on(aa,bb,ep=nuep*pow(stress,1.5)) ;


 solve NS ([u1,u2,p],[v1,v2,q]) = // Navier-Stokes k-epsilon and Boussinesq
    int2d(Th)(
             alpha*( u1*v1 + u2*v2)
            + muT * (dx(u1)*dx(v1)+dy(u1)*dy(v1)+dx(u2)*dx(v2)+dy(u2)*dy(v2))
 //           ( 2*dx(u1)*dx(v1) + 2*dy(u2)*dy(v2)+(dy(u1)+dx(u2))*(dy(v1)+dx(v2)))
            + p*q*(0.000001)
            - p*dx(v1) - p*dy(v2)
            - dx(u1)*q - dy(u2)*q
           )
  + int1d(Th,aa,bb,dd)(u1*v1* 0.1)
  + int2d(Th) (eee*(temp-35)*v1 -alpha*convect([up1,up2],-dt,up1)*v1
                             -alpha*convect([up1,up2],-dt,up2)*v2 )
   + on(ff,u1=3,u2=0)
  + on(ee,u1=0,u2=0)
  + on(aa,dd,u2=0)
  + on(bb,u2= -up1*N.x/N.y)
  + on(cc,u2=0) ;
 plot(coef=0.2,cmm=" [u1,u2] et p  ",p,[u1,u2],ps="StokesP2P1.eps",value=1,wait=1);
{
  real[int] xx(21),yy(21),pp(21);
  for (int i=0;i<21;i++)
   {
     yy[i]=i/20.;
     xx[i]=u1(0.5,i/20.);
     pp[i]=p(i/20.,0.999);
    }
      cout << " " << yy << endl;
//     plot([xx,yy],wait=1,cmm="u1 x=0.5 cup");
//     plot([yy,pp],wait=1,cmm="pressure y=0.999 cup");
}

dt = 0.05;
int nbiter = 3;
real coefdt = 0.25^(1./nbiter);
real coefcut = 0.25^(1./nbiter) , cut=0.01;
real tol=0.5,coeftol = 0.5^(1./nbiter);
nu=1./reylnods;

for (iter=1;iter<=nbiter;iter++)
{
 cout << " dt = " << dt << " ------------------------ " << endl;
  alpha=1/dt;
  for (i=0;i<=500;i++)
   {
     up1=u1;
     up2=u2;
     temp1=max(temp,25);
     temp1=min(temp1,35);
     k1=k; ep1=ep;
     muT=0.09*k*k/ep;
      NS; plot([u1,u2],wait=1); // Solves Navier-Stokes
     prode =0.126*k*(pow(2*dx(u1),2)+pow(2*dy(u2),2)+2*pow(dx(u2)+dy(u1),2))/2;
     prodk= prode*k/ep*0.09/0.126;
     kappa=muT/0.41;
     stress=abs(dy(u1));
     kine; plot(k,wait=1);
     viscturb; plot(ep,wait=1);
     TEMPER; // solves temperature equation
     if ( !(i % 5)){
         plot(temp,value=1,fill=true,ps="temp_"+iter+"_"+i+".ps");
         plot(coef=0.2,cmm=" [u1,u2] et p  ",p,[u1,u2],ps="plotNS_"+iter+"_"+i+".ps");
     }
     cout << "CPU " << clock()-s0 << "s " << endl;
   }

  if (iter>= nbiter) break;
   Th=adaptmesh(Th,[dx(u1),dy(u1),dx(u1),dy(u2)],splitpbedge=1,
               abserror=0,cutoff=cut,err=tol, inquire=0,ratio=1.5,hmin=1./1000);
 plot(Th,ps="ThNS.eps");
  dt = dt*coefdt;
  tol = tol *coeftol;
  cut = cut *coefcut;
}
cout << "CPU " <<clock()-s0 << "s " << endl;
\eFF


\subsection{An Example with Complex Numbers}

 \index{FE function!complex}\index{complex}\index{Helmholtz}
In a microwave oven heat comes from molecular excitation by an
electromagnetic field. For a plane monochromatic wave, amplitude is
given by Helmholtz's equation:
$$ \beta v + \Delta v = 0.
$$
We consider a rectangular oven where the wave is emitted by part of
the upper wall. So the boundary of the domain is made up of a part
$\Gamma_1$ where $v=0$ and of another part $\Gamma_2=[c,d]$ where for
instance $v=\sin(\pi{y-c\over c-d})$.

Within an object to be cooked, denoted by $B$, the heat source is
proportional to $v^2$.
At equilibrium, one has

$$-\Delta\theta = v^2 I_B, \quad \theta_\Gamma = 0
$$
where $I_B$ is $1$ in the object and $0$ elsewhere.
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=5cm]{rmuonde}
\includegraphics[width=5cm]{imuonde}
\includegraphics[width=5cm]{tempmuonde}
\caption{\label{figmuonde}A microwave oven: real (left) and imaginary (middle) parts
 of wave  and temperature (right).}
\end{center}
\end{figure}

Results are shown on figure \ref{figmuonde}

In the program below $\beta = 1/(1-I/2)$ in the air and $2/(1-I/2)$
in the object ($i=\sqrt{-1}$):

\begin{example}[muwave.edp]
\bFF

// file muwave.edp
real a=20, b=20, c=15, d=8, e=2, l=12, f=2, g=2;
border a0(t=0,1) {x=a*t; y=0;label=1;}
border a1(t=1,2) {x=a; y= b*(t-1);label=1;}
border a2(t=2,3) { x=a*(3-t);y=b;label=1;}
border a3(t=3,4){x=0;y=b-(b-c)*(t-3);label=1;}
border a4(t=4,5){x=0;y=c-(c-d)*(t-4);label=2;}
border a5(t=5,6){ x=0; y= d*(6-t);label=1;}

border b0(t=0,1) {x=a-f+e*(t-1);y=g; label=3;}
border b1(t=1,4) {x=a-f; y=g+l*(t-1)/3; label=3;}
border b2(t=4,5) {x=a-f-e*(t-4); y=l+g; label=3;}
border b3(t=5,8) {x=a-e-f; y= l+g-l*(t-5)/3; label=3;}
int n=2;
mesh Th = buildmesh(a0(10*n)+a1(10*n)+a2(10*n)+a3(10*n)
        +a4(10*n)+a5(10*n)+b0(5*n)+b1(10*n)+b2(5*n)+b3(10*n));
plot(Th,wait=1);
fespace Vh(Th,P1);
real meat =  Th(a-f-e/2,g+l/2).region, air= Th(0.01,0.01).region;
Vh R=(region-air)/(meat-air);

Vh<complex> v,w;
solve muwave(v,w) = int2d(Th)(v*w*(1+R)
                -(dx(v)*dx(w)+dy(v)*dy(w))*(1-0.5i))
   + on(1,v=0) + on(2, v=sin(pi*(y-c)/(c-d)));
Vh vr=real(v), vi=imag(v);
plot(vr,wait=1,ps="rmuonde.ps", fill=true);
plot(vi,wait=1,ps="imuonde.ps", fill=true);

fespace Uh(Th,P1); Uh u,uu, ff=1e5*(vr^2 + vi^2)*R;

solve temperature(u,uu)= int2d(Th)(dx(u)* dx(uu)+ dy(u)* dy(uu))
     - int2d(Th)(ff*uu) + on(1,2,u=0);
plot(u,wait=1,ps="tempmuonde.ps", fill=true);
\eFF
\end{example}

\subsection{Optimal Control} Thanks to the function `BFGS} it
is possible to solve complex nonlinear optimization problem within
FreeFem++. For example consider the following inverse
problem
\[
    \min_{b,c,d\in R}J=\int_E(u-u_d)^2~:~
    -\nabla(\kappa(b,c,d)\cdot\nabla u)=0,~~u|_\Gamma=u_\Gamma
\]
where the desired state $u_d$, the boundary data $u_\Gamma$ and the
observation set $E\subset\Omega$ are all given.  Furthermore let us
assume that
\[
    \kappa(x)=1+bI_B(x)+cI_C(x)+dI_D(x)~~~\forall x\in\Omega
\]
where $B,C,D$ are separated subsets of $\Omega$.
\\\\
To solve this problem by the quasi-Newton BFGS method we need the
derivatives of $J$ with respect to $b,c,d$.  We self explanatory
notations, if $\delta b,\delta c,\delta d$ are variations of
$b,c,d$ we have
\[
    \delta J\approx 2\int_E(u-u_d)\delta u,~~
    -\nabla(\kappa\cdot\nabla\delta u)\approx\nabla(\delta\kappa\cdot\nabla
    u)
    ~~\delta u|_\Gamma=0
\]
Obviously $J'_b$ is equal to $\delta J$ when $\delta b=1,\delta
c=0,\delta d=0$, and so on for $J'_c$ and $J'_d$.
\\\\
All this is implemented in the following program
 \bFF
// file optimcontrol.edp
border aa(t=0, 2*pi) {    x = 5*cos(t);    y = 5*sin(t);  };
border bb(t=0, 2*pi) {    x = cos(t);    y = sin(t);  };
border cc(t=0, 2*pi) {    x = -3+cos(t);    y = sin(t);  };
border dd(t=0, 2*pi) {    x = cos(t);    y = -3+sin(t);  };
mesh th = buildmesh(aa(70)+bb(35)+cc(35)+dd(35));
fespace Vh(th,P1);
Vh Ib=((x^2+y^2)<1.0001),
   Ic=(((x+3)^2+ y^2)<1.0001),
   Id=((x^2+(y+3)^2)<1.0001),
   Ie=(((x-1)^2+ y^2)<=4),
   ud,u,uh,du;
real[int] z(3);
problem A(u,uh) =int2d(th)((1+z[0]*Ib+z[1]*Ic+z[2]*Id)*(dx(u)*dx(uh)
                    +dy(u)*dy(uh))) + on(aa,u=x^3-y^3);
z[0]=2; z[1]=3; z[2]=4;
A; ud=u;
ofstream f("J.txt");
func real J(real[int] & Z)
{
    for (int i=0;i<z.n;i++)z[i]=Z[i];
    A; real s= int2d(th)(Ie*(u-ud)^2);
    f<<s<<"   "; return s;
}

real[int] dz(3), dJdz(3);

problem B(du,uh)
  =int2d(th)((1+z[0]*Ib+z[1]*Ic+z[2]*Id)*(dx(du)*dx(uh)+dy(du)*dy(uh)))
  +int2d(th)((dz[0]*Ib+dz[1]*Ic+dz[2]*Id)*(dx(u)*dx(uh)+dy(u)*dy(uh)))
  +on(aa,du=0);

func real[int] DJ(real[int] &Z)
    {
      for(int i=0;i<z.n;i++)
        { for(int j=0;j<dz.n;j++) dz[j]=0;
          dz[i]=1; B;
          dJdz[i]= 2*int2d(th)(Ie*(u-ud)*du);
      }
     return dJdz;
 }

 real[int] Z(3);
 for(int j=0;j<z.n;j++) Z[j]=1;
 BFGS(J,DJ,Z,eps=1.e-6,nbiter=15,nbiterline=20);
 cout << "BFGS: J(z) = " << J(Z) <<  endl;
 for(int j=0;j<z.n;j++) cout<<z[j]<<endl;
 plot(ud,value=1,ps="u.eps");
\eFF
In this example the sets $B,C,D,E$ are circles of boundaries $bb,cc,dd,ee$ are the domain
$\Omega$ is the circle of boundary $aa$.  The desired state $u_d$ is the solution
of the PDE for $b=2,c=3,d=4$.  The unknowns are packed into array $z$.  Notice that it is
necessary to recopy $Z$ into $z$ because one is a local variable while the other one is global.
The program found $b=2.00125,c=3.00109,d=4.00551$.
Figure \ref{figcontrol} shows $u$ at convergence
and the successive function evaluations of $J$.
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=5cm]{u-bfgs}
\includegraphics[width=5cm]{J-bfgs}
\caption{\label{figcontrol}On the left the level lines of $u$. On the right the successive evaluations of
$J$ by BFGS (5 values above 500 have been removed for readability)}
\end{center}
\end{figure}
Note that an \emph{adjoint state} could have been used. Define $p$ by
\[
-\nabla\cdot(\kappa\nabla p)=2I_E(u-u_d),~~~p|_\Gamma=0
\]
Consequently
\begin{eqnarray}&&
\delta J = -\int_{\Omega} (\nabla\cdot(\kappa\nabla p))\delta u
\cr&&
= \int_\Omega(\kappa\nabla p\cdot\nabla\delta u)
=-\int_\Omega(\delta\kappa\nabla p\cdot\nabla u)
\end{eqnarray}
Then the derivatives are found by setting $\delta b=1, \delta c=\delta d=0$ and so on:
\[
    J'_b=-\int_B \nabla p\cdot\nabla u,~~J'_c=-\int_C \nabla p\cdot\nabla u,~~
    J'_d=-\int_D \nabla p\cdot\nabla u
\]
\paragraph{Remark} As BFGS stores an $M\times M$ matrix where $M$ is the number of
unknowns, it is dangerously expensive to use this method when the unknown $x$ is a
Finite Element Function.  One should use another optimizer such as
the NonLinear Conjugate Gradient `NLCG} (also a key word
of FreeFem++).  See the file algo.edp in the examples folder.

\subsection{A Flow with Shocks}
Compressible Euler  equations should be discretized with Finite Volumes or FEM with flux up-winding scheme but these are not implemented in FreeFem++.  Nevertheless acceptable results can be obtained with the method of characteristics
provided that the mean values $\bar f=\frac12(f^++f^-)$ are used at shocks in the scheme, and finally  mesh adaptation \index{adaptmesh}\index{mesh!adaptation}.%
\begin{eqnarray}\label{euler}&&
    \p_t\rho+\bar u\n\rho + \bar\rho\n\cdot u=0
    \cr&&
   \bar\rho( \p_t u+\frac{\overline{\rho u}}{\bar\rho}\n u +\n p=0
    \cr&&
    \p_t p + \bar u\n p +(\gamma-1)\bar p\n\cdot u =0
\end{eqnarray}
%
One possibility is to couple $u,p$ and then update $\rho$, i.e.
%
\begin{eqnarray}\label{eulalgo}&&
    \frac 1{(\gamma-1)\delta t\bar p^m} (p^{m+1}-p^m \circ X^m) + \n\cdot u^{m+1} =0
    \cr&&
    \frac{\bar\rho^m}{\delta t}(u^{m+1}-u^m \circ {\tilde X}^m ) +\n p^{m+1}=0
    \cr&&
    \rho^{m+1} = \rho^m \circ X^m +
        \frac{\bar\rho^m}{(\gamma-1)\bar p^m}(p^{m+1}-p^m \circ X^m)
\end{eqnarray}
A numerical result is given on Figure \ref{figvfive} and the FreeFem++ script is

%
\begin{figure}[htbp]
\begin{center}
\includegraphics[width=10cm]{mach2r}
\caption{ \label{figvfive} Pressure for a Euler flow around a disk at Mach 2 computed by (\ref{eulalgo}) }
\end{center}
\end{figure}
%
\bFF
verbosity=1;
int anew=1;
real x0=0.5,y0=0, rr=0.2;
border ccc(t=0,2){x=2-t;y=1;};
border ddd(t=0,1){x=0;y=1-t;};
border aaa1(t=0,x0-rr){x=t;y=0;};
border cercle(t=pi,0){ x=x0+rr*cos(t);y=y0+rr*sin(t);}
border aaa2(t=x0+rr,2){x=t;y=0;};
border bbb(t=0,1){x=2;y=t;};

int m=5; mesh Th;
if(anew) Th = buildmesh (ccc(5*m) +ddd(3*m) + aaa1(2*m) + cercle(5*m)
              + aaa2(5*m) + bbb(2*m) );
      else Th = readmesh("Th_circle.mesh"); plot(Th,wait=0);

real dt=0.01, u0=2, err0=0.00625, pena=2;
fespace Wh(Th,P1);
fespace Vh(Th,P1);
Wh u,v,u1,v1,uh,vh;
Vh r,rh,r1;
macro dn(u) (N.x*dx(u)+N.y*dy(u) ) //  def the normal derivative

if(anew){ u1= u0; v1= 0; r1 = 1;}
else {
    ifstream g("u.txt");g>>u1[];
    ifstream gg("v.txt");gg>>v1[];
    ifstream ggg("r.txt");ggg>>r1[];
    plot(u1,ps="eta.eps", value=1,wait=1);
    err0=err0/10; dt = dt/10;
}

problem  eul(u,v,r,uh,vh,rh)
   = int2d(Th)(  (u*uh+v*vh+r*rh)/dt
                  + ((dx(r)*uh+ dy(r)*vh) - (dx(rh)*u + dy(rh)*v))
               )
 + int2d(Th)(-(rh*convect([u1,v1],-dt,r1) + uh*convect([u1,v1],-dt,u1)
                + vh*convect([u1,v1],-dt,v1))/dt)
  +int1d(Th,6)(rh*u)   // +int1d(Th,1)(rh*v)
 + on(2,r=0) + on(2,u=u0) + on(2,v=0);

int j=80;
for(int k=0;k<3;k++)
{
    if(k==20){ err0=err0/10; dt = dt/10; j=5;}
    for(int i=0;i<j;i++){
       eul; u1=u; v1=v; r1=abs(r);
        cout<<"k="<<k<<"  E="<<int2d(Th)(u^2+v^2+r)<<endl;
        plot(r,wait=0,value=1);
}
Th = adaptmesh (Th,r, nbvx=40000,err=err0,
      abserror=1,nbjacoby=2, omega=1.8,ratio=1.8, nbsmooth=3,
      splitpbedge=1, maxsubdiv=5,rescaling=1) ;
 plot(Th,wait=0);
 u=u;v=v;r=r;

savemesh(Th,"Th_circle.mesh");
ofstream f("u.txt");f<<u[];
ofstream ff("v.txt");ff<<v[];
ofstream fff("r.txt");fff<<r[];
r1 = sqrt(u*u+v*v);
plot(r1,ps="mach.eps", value=1);
r1=r;
}
\eFF

 \subsection{Classification of the equations}
\paragraph{Summary}\emph{
It is usually not easy to determine the type of a system.  Yet the approximations
and algorithms suited to the problem depend on its type:
\begin{itemize}
\item Finite Elements compatible (LBB conditions) for elliptic systems
\item Finite difference on the parabolic variable and a time loop on each
elliptic subsystem of parabolic systems; better stability diagrams when the schemes are implicit in time.
\item Upwinding, Petrov-Galerkin, Characteristics-Galerkin, Discontinuous-Galerkin, Finite Volumes
for hyperbolic systems plus, possibly, a time loop.
\end{itemize}
When the system changes type, then expect difficulties (like shock discontinuities)!}

 \paragraph{Elliptic, parabolic and hyperbolic equations}

A partial differential equation (PDE) is a relation between a function
of several variables and its derivatives.
 $$
 F(\varphi(x),{\p\varphi\over\p
 x_1}(x),\cdots,{\p\varphi\over\p
 x_d}(x),{\p^2\varphi\over\p
 x^2_1}(x),\cdots,{\p^m\varphi\over\p x^m_d}(x)) =
 0\quad\forall x\in\Omega\subset \Rel^d.
 $$
 The range of $x$ over which the equation is taken, here $\Omega$, is called
 the \emph{domain} of the PDE.
The highest derivation index, here $m$, is called the {\it
 order}. If $F$ and $\varphi$ are vector valued functions, then the
 PDE is actually a \emph{system} of PDEs.
 \\
Unless indicated otherwise, here by convention \emph{one} PDE
 corresponds to one scalar valued $F$ and $\varphi$.
If $F$ is linear with respect to its arguments, then the PDE is said
 to be \emph{linear}.
 \\
The general form of a second order, linear scalar PDE is
${\p^2\varphi\over\p x_i\p x_j}$ and $A:B$ means
$\sum^d_{i,j=1} a_{ij} b_{ij}.$
 $$\alpha\varphi + a\cdot\nabla\varphi + B :\nabla(\nabla\varphi) =
 f{\quad\hbox{ in }\quad}\Omega\subset \Rel^d,
 $$
where $f(x),\alpha(x)\in \Rel, a(x)\in \Rel^d, B(x)\in \Rel^{d\times d}$
are the PDE \emph{coefficients}.
If the coefficients are independent of $x$, the PDE is said to have
\emph{constant coefficients}.

To a PDE we associate a quadratic form, by replacing
$\varphi$ by $1$,
$\p\varphi/\p x_i$ by $z_i$ and
$\p^2\varphi/\p x_i\p x_j$ by $z_i z_j$, where $z$
is a vector in $\Rel^d$~:
$$\alpha + a\cdot z + z^T Bz = f.
$$
If it is the equation of an ellipse (ellipsoid if $d \geq 2$),
the PDE is said to be {\it elliptic};
if it is the equation of a parabola or a hyperbola, the PDE is said to
be {\it parabolic} or {\it hyperbolic}. If $A \equiv 0$, the degree is
no longer 2 but 1, and for reasons that will appear more clearly
later, the PDE is still said to be hyperbolic.
\\\\
These concepts can be generalized to systems, by studying whether or
not the polynomial system $P(z)$ associated with the PDE system has
branches at infinity (ellipsoids have no branches at infinity,
paraboloids have one, and hyperboloids have several).
\\
If the PDE is not linear, it is said to be \emph{non linear}.
Those are said to be locally elliptic, parabolic, or hyperbolic
according to the type of the linearized equation.
\\
For example, for the non linear equation
 $${\p^2\varphi\over\p t^2} - {\p\varphi\over\p
 x}{\p^2\varphi\over\p x^2} = 1,
 $$
we have $d=2, x_1 = t, x_2 = x$ and its linearized form is:
 $${\p^2 u\over\p t^2} - {\p u\over\p x}
 {\p^2\varphi\over\p x^2} - {\p\varphi\over\p
 x}{\p^2 u\over\p x^2} = 0,
 $$
which for the unknown $u$ is locally elliptic if
${\p\varphi\over\p x} < 0$  and locally hyperbolic if
${\p\varphi\over\p x} > 0$.

\paragraph{Examples}

\noindent{Laplace's} equation is elliptic:
 $$\Delta\varphi \equiv {\p^2\varphi\over\p x^2_1} +
 {\p^2\varphi\over\p x^2_2} + \cdots +
 {\p^2\varphi\over\p x^2_d} = f, \ \ \ \forall x
 \in \Omega\subset \Rel^d.
 $$
The \emph{heat} equation is parabolic in
 $Q = \Omega\times]0,T[\subset \Rel^{d+1}$~:
 $${\p\varphi\over\p t} - \mu\Delta\varphi = f\quad\forall
 x\in\Omega\subset  \Rel^d, \quad\forall t\in]0,T[.
 $$
If $\mu >0$,  the \emph{wave} equation is hyperbolic:
 $${\p^2\varphi\over\p t^2} - \mu\Delta\varphi =
 f{\quad\hbox{~in~}\quad}  Q.
 $$
The \emph{convection diffusion} equation is parabolic if $\mu \neq 0$
and hyperbolic otherwise:
 $${\p\varphi\over\p t} + a\nabla\varphi -
 \mu\Delta\varphi = f.
 $$
The \emph{biharmonic} equation is elliptic:
 $$\Delta(\Delta\varphi) = f{\quad\hbox{~in~}\quad}\Omega.
 $$

\paragraph{Boundary conditions}

A relation between a function and its derivatives is not sufficient to define the function.
Additional information on the boundary $\Gamma=\p\Omega$ of
$\Omega$, or on part of $\Gamma$ is necessary.\\
Such information is called a \emph{boundary condition}.
For example,
 $$\varphi(x) \ \hbox{given},\ \forall x\in \Gamma,
 $$
is called a \emph{\x{Dirichlet} boundary condition}.
The \emph{\x{Neumann}} condition is
 $${\p\varphi\over\p n}(x) \ \hbox{given on }\
 \Gamma \hbox{~~~(or~~} n\cdot B\nabla\varphi,\hbox{given on }\
 \Gamma\hbox{ for a general second order PDE)}
$$
where $n$ is the normal at $x\in\Gamma$
directed towards the exterior of $\Omega$ (by definition
${\p\varphi\over\p n}=\nabla\varphi\cdot n$).

Another classical condition, called a \emph{Robin} (or \emph{Fourier})
condition is written as:
 $$\varphi(x) + \beta(x) {\p\varphi\over\p n}(x) \
 \hbox{given on}\
 \Gamma.
 $$
Finding a set of boundary conditions  that defines a unique
$\varphi$ is a difficult art.\\
In general, an elliptic equation is well posed (\emph{i.e.} $\varphi$
is unique) with one Dirichlet, Neumann or Robin conditions on the whole boundary.
\\
Thus, Laplace's equations  is well posed with
a Dirichlet or Neumann condition but also with
 $$\varphi \ \hbox{given on}\ \Gamma_1,\quad {\p\varphi\over
 \p n} \
 \hbox{given on}\ \Gamma_2, \quad \Gamma_1\cup\Gamma_2 =
 \Gamma,\quad{\dot{\Gamma_1}\cap\dot{\Gamma_2}} =
 \emptyset.$$
Parabolic and hyperbolic equations rarely require boundary conditions
on all of  $\Gamma\times]0,T[$. For instance, the heat equation
is well posed with
 $$\varphi \ \hbox{given at}\ t=0 \ \hbox{and Dirichlet or Neumann or mixed conditions on}\
 \p\Omega.$$
Here $t$ is time so the first condition is called an \x{initial condition}.  The whole set of conditions
are also called \x{Cauchy} conditions.
\\
The wave equation  is well posed with
$$\varphi \ \hbox{and}\ {\p\varphi\over\p t} \
 \hbox{given at}\ t=0
 \ \hbox{and Dirichlet or Neumann or mixed conditions on}\
 \p\Omega.
$$
\subsection{Optimize Time depend schema for Heat equation}


 First, it is possible to define variational forms, and use this forms to build matrix 
  and vector to make very fast script (4 times faster here) (see the example  \hrefexample{chapt3}{Heat.edp}  directory).
  
  For example solve the Thermal Conduction problem of section \ref{ss Thermal Conduction}
  We must solve the temperature equation in $\Omega$ in a time interval (0,T).
{\begin{eqnarray}&&
    \p_t u -\nabla\cdot(\kappa\nabla u)=0 \hbox{ in } \Omega\times(0,T),
    \cr&&
    u(x,y,0)=u_0+x u_1
    \cr&&
   u = 30 \hbox{ on } \Gamma_{24}\times(0,T) , \quad  \kappa\frac{\p u}{\p n} +\alpha(u-u_e)=0\hbox{ on } \Gamma\times(0,T).
\end{eqnarray}}
  
 The variational formulation is  in {$L^2(0,T;H^1(\Omega))$};
  we shall seek $u^n$ satisfying
{\[
\forall w \in V_{0}; \qquad   \int_\Omega \frac{u^n-u^{n-1}}{\delta t} w + \kappa\nabla u^n\nabla w) +\int_\Gamma\alpha(u^n-u_{ue})w=0
\]}
where {$ V_0 = \{w\in H^1(\Omega)/ w_{|\Gamma_{24}}=0\}$}.

So the to code the method with the matrices $A=(A_{ij})$, $M=(M_{ij})$, and  the vectors 
$ u^n, b^n, b',b", b_{cl}$ 
( notation if $w$ is a vector then $w_i$ is a component of the vector).
\def\tgv{\Red{{\frac{1}{\varepsilon}}}}
{$$ u^n = A^{-1} b^n, \quad
  \quad b' = b_0 + M u^{n-1}, 
  \quad b"=  \tgv \; b_{cl} , 
  \quad  b^n_i = \left\{
  \begin{array}{cl}   b"_i  & \mbox{if }\ i \in \Gamma_{24} \\
                       b'_i & \mbox{else } %\not\in \Gamma_{24}
                        \end{array}\right.
                       \label{eq tgv}  $$}
Where with $ \tgv = \Red{\mathtt{tgv}} = \Red{10^{30}}$ :
\begin{eqnarray*}
 A_{ij} &=& \left\{\begin{array}{cl}   \tgv  & \mbox{if } i  \in \Gamma_{24}, \mbox{and}\quad  j=i \\  
\displaystyle 
 {\int_{\Omega} w_j w_i / dt + k (\nabla w_j. \nabla w_i ) + \int_{\Gamma_{13}} \alpha w_j w_i} & \mbox{else } % i  \not\in \Gamma_{24}, \mbox{or}\quad  j\ne i 
 \end{array}\right.  \\ 
 M_{ij} &=& \left\{\begin{array}{cl}   \tgv & \mbox{if } i  \in \Gamma_{24}, \mbox{and}\quad  j=i  \\  
\displaystyle 
  n{\int_{\Omega} w_j w_i / dt}
 & \mbox{else  } %i  \not\in \Gamma_{24}, \mbox{or}  j\ne i 
   \end{array}\right. \\ 
 b_{0,i} &=&  n{\int_{\Gamma_{13}} \alpha u_{ue} w_i } \\
 b_{cl} &=& u^{0}  \quad \mbox{the initial data} 
\end{eqnarray*}


The Fast version script:

\bFF
 ...
Vh u0=fu0,u=u0; 
\eFF

Create three variational formulation, and build the matrices $A$,$M$.\label{matrix-varf}
\bFF
varf vthermic (u,v)= int2d(Th)(u*v/dt 
           + k*(dx(u) * dx(v) + dy(u) * dy(v)))  
  +  int1d(Th,1,3)(alpha*u*v)  + on(2,4,u=1); 
varf vthermic0(u,v) =   int1d(Th,1,3)(alpha*ue*v);
varf vMass (u,v)= int2d(Th)( u*v/dt)  + on(2,4,u=1);

real tgv = 1e30;
matrix A= vthermic(Vh,Vh,tgv=tgv,solver=CG);
matrix M= vMass(Vh,Vh);
\eFF


Now, to build the right hand size we need 4 vectors.

\bFF
real[int]  b0 = vthermic0(0,Vh);//constant part of  RHS 
real[int]  bcn = vthermic(0,Vh); //tgv on Dirichlet part  
// we have for the node $i$ : $i\in \Gamma_{24}  \quad \Leftrightarrow \quad bcn[i] \ne 0 $
real[int]  bcl=tgv*u0[]; //  the Dirichlet B.C. part 

// The fast loop ... 
for(real t=0;t<T;t+=dt){
  real[int] b = b0 ; // for the  RHS
  b += M*u[]; //add the the time dependent part
  b = bcn ? bcl : b;//do $\forall i$:  b[i] =  bcn[i] ? bcl[i] : b[i]  ;    
  u[] = A^-1*b;   //  Solve linear problem  
  plot(u);
}
\eFF 

\subsection{Tutorial to write a transient Stokes solver in matrix form}

Consider the following script to solve a time dependent Stokes problem in a cavity

\bFF
mesh Th=square(10,10);
fespace Vh(Th,P2), Qh(Th,P1);
Vh u,v,uu,vv, uold=0,vold=0;
Qh p,pp;
real nu=0.1, T=1., dt = 0.1;
int m, M= T/dt;
problem stokes(u, v, p, uu, vv, pp)=
  int2d(Th)(  (u*uu+v*vv)/dt + nu*(dx(u)*dx(uu) + dy(u)*dy(uu)
            	 + dx(v)*dx(vv) + dy(v)*dy(vv)) 
                 - p*pp*1.e-6   - p*(dx(uu) +dy(vv))- pp*(dx(u)+dy(v))
           	  ) - int2d(Th)((uold*uu+vold*vv)/dt)
  		 		+ on(1,2,4,u=0,v=0) + on(3,u=1,v=0);

for(m=0;m<M;m++){
	stokes; uold=u; vold=v;
}
plot(p,[u,v],value=true, wait=true, cmm="t="+m*dt);
\eFF

Every iteration is in fact of the form $ A[u,v,p] = B[uold,vold,pold] + b$ where $A,B$ are matrices and $b$ is a vector containing the boundary conditions.
The $A,B,b$ are constructed by

\bFF
fespace Xh(Th,[P2,P2,P1]);
varf aa ([u,v,p],[uu,vv,pp])
  = int2d(Th)(  (u*uu+v*vv)/dt + nu*(dx(u)*dx(uu) + dy(u)*dy(uu)
            	 + dx(v)*dx(vv) + dy(v)*dy(vv)) 
                 - p*pp*1.e-6  - p*(dx(uu) +dy(vv))- pp*(dx(u)+dy(v))
           ) + on(1,2,4,u=0,v=0) + on(3,u=1,v=0);

varf bb ([uold,vold,pold],[uu,vv,pp]) 
 = int2d(Th)((uold*uu+vold*vv)/dt)
;//  + on(1,2,4,uold=0,vold=0) + on(3,uold=0,vold=0);

varf bcl ([uold,vold,pold],[uu,vv,pp]) = on(1,2,4,uold=0,vold=0) + on(3,uold=1,vold=0);

matrix A= aa(Xh,Xh,solver=UMFPACK); 
matrix B= bb(Xh,Xh);
real[int] b=bcl(0,Xh); 
\eFF
Note that the boundary conditions are not specified in $bb$. Removing the comment ``//" would cause the compiler to multiply the diagonal terms corresponding to a Dirichlet degree of freedom by a very large term (tgv); if so $b$ would not be needed, on the condition that $uold=1$ on boundary 3 initially. Note also that b has a tgv on the Dirichlet nodes, by construction, and so does A.

The loop will them be 
\bFF
real[int] sol(Xh.ndof), aux(Xh.ndof);
for(m=0;m<M;m++){
    aux=B*sol;  aux+=b;
    sol=A^-1*aux;
}
\eFF
There is yet a difficulty with the initialization of `sol} and with the solution from `sol}.  For this we need a temporary vector in $X_h$ and here is a solution
\bFF
Xh [w1,w2,wp]=[uold,vold,pp];  
sol=w1[]; // cause also the copy of w2 and wp
for(m=0;m<M;m++){
    aux=B*sol;  aux+=b;
    sol=A^-1*aux;
}
w1[]=sol;  u=w1; v= w2; p=wp;
plot(p,[u,v],value=true, wait=true, cmm="t="+m*dt);
\eFF
The freefem team agrees that the line `sol=w1[];} is mysterious as it copies also w2 and wp into sol. Structured data such as vectors of $X_h$ here cannot be written component by component. Hence `w1=u} is not allowed.

