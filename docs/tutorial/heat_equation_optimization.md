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

