## Poisson Equation

```freefem
border aaa(t=0,1){x=t;y=0;};
border bbb(t=0,0.5){x=1;y=t;};
border ccc(t=0,0.5){x=1-t;y=0.5;};
border ddd(t=0.5,1){x=0.5;y=t;};
border eee(t=0.5,1){x=1-t;y=1;};
border fff(t=0,1){x=0;y=1-t;};
mesh Th = buildmesh (aaa(6) + bbb(4) + ccc(4) +ddd(4) + eee(4) + fff(6));
fespace Vh(Th,P1); //to change P1 in P2 to make P2 finite element.
Vh u=0,v;
func f= 1;
func g= 0;
int i=0;
real error=0.1, coef= 0.1^(1./5.);
problem Probem1(u,v,solver=CG,eps=-1.0e-6) =
    int2d(Th)(  dx(u)*dx(v) + dy(u)*dy(v)) 
  + int2d(Th) ( v*f ) 
  + on(aaa,bbb,ccc,ddd,eee,fff,u=g)  ;
  
for (i=0;i< 10;i++)
{   
  real d = clock();
  Probem1; //solves the problem 
  plot(u,Th,wait=1);
  Th=adaptmesh(Th,u,inquire=1,err=error);
  error = error * coef;
} ;
```

Solution on adapted mesh and associated mesh   |  |
:-------------------------:|:-------------------------:
![poisson Associated mesh](images/poisson_associated_mesh.jpg)  |  ![poisson adapted mesh](images/poisson_adapted_mesh.jpg)

## Stoke Equation on Cube

```freefem
load "msh3" load "medit" // dynamics load tools for 3d.
int nn=8;
mesh Th2=square(nn,nn);
fespace Vh2(Th2,P2);  Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1;
mesh3 Th=buildlayers(Th2,nn,zbound=[zmin,zmax],reffacemid=rmid, 
  reffaceup=rup, reffacelow=rdown);
  
medit("c10x10x10",Th); // see the 3d mesh with medit software
fespace VVh(Th,[P2,P2,P2,P1]);

macro Grad(u) [dx(u),dy(u),dz(u)] // EOM
macro div(u1,u2,u3) (dx(u1)+dy(u2)+dz(u3))  // EOM

VVh [u1,u2,u3,p];
VVh [v1,v2,v3,q];
  
solve vStokes([u1,u2,u3,p],[v1,v2,v3,q]) = 
  int3d(Th,qforder=3)( Grad(u1)'*Grad(v1) +  Grad(u2)'*Grad(v2) +  Grad(u3)'*Grad(v3)
                  - div(u1,u2,u3)*q - div(v1,v2,v3)*p + 1e-10*q*p ) 
  + on(2,u1=1.,u2=0,u3=0) + on(1,u1=0,u2=0,u3=0) ;
 plot(p,wait=1, nbiso=5);  // A 3d plot of iso  pressure. in progress... march 2009
 //  to see the 10 cut plan in 2d 
for(int i=1;i<10;i++)
{
 real yy=i/10.; // Compute yy.
 // do 3d -> 2d interpolation.
 ux= u1(x,yy,y); uz= u3(x,yy,y);  p2= p(x,yy,y);
 plot([ux,uz],p2,cmm=" cut y = "+yy,wait= 1);
}
```

Solution on cup plan $y=0.5$ and $10*10*10$ associated mesh  |
:-------------------------:|
![Stokes 3D](images/Stokes3d.jpg)  |
![Stokes 3d mesh](images/Stokes3d-Th.jpg)  |