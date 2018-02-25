# Optimize Time depend schema for Heat equation


First, it is possible to define variational forms, and use this forms to build matrix and vector to make very fast script (4 times faster here) (see the example Heat.edp).
  
For example solve the Thermal Conduction problem of section 3.4
We must solve the temperature equation in $\Omega$ in a time interval (0,T).

\begin{eqnarray}&&
    \p_t u -\nabla\cdot(\kappa\nabla u)=0 \hbox{ in } \Omega\times(0,T),
    \cr&&
    u(x,y,0)=u_0+x u_1
    \cr&&
   u = 30 \hbox{ on } \Gamma_{24}\times(0,T) , \quad  \kappa\frac{\p u}{\p n} +\alpha(u-u_e)=0\hbox{ on } \Gamma\times(0,T).
\end{eqnarray}
  
The variational formulation is  in $L^2(0,T;H^1(\Omega))$; we shall seek $u^n$ satisfying

\[
\forall w \in V_{0}; \qquad   \int_\Omega \frac{u^n-u^{n-1}}{\delta t} w + \kappa\nabla u^n\nabla w) +\int_\Gamma\alpha(u^n-u_{ue})w=0
\]

where $ V_0 = \{w\in H^1(\Omega)/ w_{|\Gamma_{24}}=0\}$.

So the to code the method with the matrices $A=(A_{ij})$, $M=(M_{ij})$, and  the vectors 
$ u^n, b^n, b',b", b_{cl}$ (notation if $w$ is a vector then $w_i$ is a component of the vector).

\def\tgv{\Red{{\frac{1}{\varepsilon}}}}

$$ u^n = A^{-1} b^n, \quad
  \quad b' = b_0 + M u^{n-1}, 
  \quad b"=  \tgv \; b_{cl} , 
  \quad  b^n_i = \left\{
  \begin{array}{cl}   b"_i  & \mbox{if }\ i \in \Gamma_{24} \\
                       b'_i & \mbox{else } %\not\in \Gamma_{24}
                        \end{array}\right.
                       \label{eq tgv}
$$

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

```freefem
 ...
Vh u0=fu0,u=u0; 
```

Create three variational formulation, and build the matrices $A$,$M$.\label{matrix-varf}

```freefem
varf vthermic (u,v)= int2d(Th)(u*v/dt 
           + k*(dx(u) * dx(v) + dy(u) * dy(v)))  
  +  int1d(Th,1,3)(alpha*u*v)  + on(2,4,u=1); 
varf vthermic0(u,v) =   int1d(Th,1,3)(alpha*ue*v);
varf vMass (u,v)= int2d(Th)( u*v/dt)  + on(2,4,u=1);

real tgv = 1e30;
matrix A= vthermic(Vh,Vh,tgv=tgv,solver=CG);
matrix M= vMass(Vh,Vh);
```


Now, to build the right hand size we need 4 vectors.

```freefem
real[int]  b0 = vthermic0(0,Vh); // Constant part of  RHS 
real[int]  bcn = vthermic(0,Vh); // tgv on Dirichlet part  
// we have for the node $i$ : $i\in \Gamma_{24}  \quad \Leftrightarrow \quad bcn[i] \ne 0 $
real[int]  bcl=tgv*u0[]; // The Dirichlet B.C. part 

// The fast loop ... 
for(real t=0;t<T;t+=dt){
  real[int] b = b0 ; // For the  RHS
  b += M*u[]; // Add the the time dependent part
  b = bcn ? bcl : b; // Do $\forall i$:  b[i] =  bcn[i] ? bcl[i] : b[i]  ;    
  u[] = A^-1*b; // Solve linear problem  
  plot(u);
}
```

