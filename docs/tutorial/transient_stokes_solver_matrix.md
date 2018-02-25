# Tutorial to write a transient Stokes solver in matrix form

Consider the following script to solve a time dependent Stokes problem in a cavity

```freefem
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
```

Every iteration is in fact of the form $ A[u,v,p] = B[uold,vold,pold] + b$ where $A,B$ are matrices and $b$ is a vector containing the boundary conditions.
The $A,B,b$ are constructed by

```freefem
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
```

Note that the boundary conditions are not specified in $bb$. Removing the comment "//" would cause the compiler to multiply the diagonal terms corresponding to a Dirichlet degree of freedom by a very large term (tgv); if so $b$ would not be needed, on the condition that $uold=1$ on boundary 3 initially. Note also that b has a tgv on the Dirichlet nodes, by construction, and so does A.

The loop will them be 
```freefem
real[int] sol(Xh.ndof), aux(Xh.ndof);
for(m=0;m<M;m++){
    aux=B*sol;  aux+=b;
    sol=A^-1*aux;
}
```

There is yet a difficulty with the initialization of `sol` and with the solution from `sol`.  For this we need a temporary vector in $X_h$ and here is a solution

```freefem
Xh [w1,w2,wp]=[uold,vold,pp];  
sol=w1[]; // Cause also the copy of w2 and wp
for(m=0;m<M;m++){
    aux=B*sol;  aux+=b;
    sol=A^-1*aux;
}
w1[]=sol;  u=w1; v= w2; p=wp;
plot(p,[u,v],value=true, wait=true, cmm="t="+m*dt);
```

The freefem team agrees that the line `sol=w1[];` is mysterious as it copies also w2 and wp into sol. Structured data such as vectors of $X_h$ here cannot be written component by component. Hence `w1=u` is not allowed.

