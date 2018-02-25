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

