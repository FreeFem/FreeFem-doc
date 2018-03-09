$\codered$
This script does not run as expected

# A Large Fluid Problem

A friend of one of us in Auroville-India was building a ramp to access an air conditioned room. As I was visiting the construction site he told me that he expected to cool air escaping by the door to the room to slide down the ramp and refrigerate the feet of the coming visitors.  I told him "no way" and decided to check numerically.  The results are on the front page of this book.

The fluid velocity and pressure are solution of the Navier-Stokes equations with varying density function of the temperature.

The geometry is trapezoidal with prescribed inflow made of cool air at the bottom and warm air above and so are the initial conditions; there is free outflow, slip velocity at the top (artificial) boundary and no-slip at the bottom.  However the Navier-Stokes cum temperature equations have a RANS $k-\epsilon$ model and a Boussinesq approximation for the buoyancy. This comes to :

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

In variational form and with appropriated boundary conditions the problem is :

$\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$ $\codered$
```freefem
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

fespace Vh2(Th,P1b); // Velocity space
fespace Vh(Th,P1); // Pressure space
fespace V0h(Th,P0); // For gradients
Vh2 u2,v2,up1=0,up2=0;
Vh2 u1,v1;
Vh  u1x=0,u1y,u2x,u2y, vv;

real reylnods=500;
//cout << " Enter the reynolds number :"; cin >> reylnods;
assert(reylnods>1 && reylnods < 100000);
up1=0;
up2=0;
func g=(x)*(1-x)*4;  // Inflow
Vh p=0,q, temp1,temp=35, k=0.001,k1,ep=0.0001,ep1;
V0h muT=1,prodk,prode, kappa=0.25e-4, stress;
real alpha=0, eee=9.81/303, c1m = 1.3/0.09 ;
real  nu=1, numu=nu/sqrt( 0.09), nuep=pow(nu,1.5)/4.1;
int i=0,iter=0;
real dt=0;
problem TEMPER(temp,q) = // Temperature equation
    int2d(Th)(
             alpha*temp*q + kappa * ( dx(temp)*dx(q) + dy(temp)*dy(q) ))
//   + int1d(Th,aa,bb)(temp*q* 0.1)
  + int2d(Th) ( -alpha*convect([up1,up2],-dt,temp1)*q )
   + on(ff,temp=25)
  + on(aa,bb,temp=35) ;

problem kine(k,q)= // Get the kinetic turbulent energy
    int2d(Th)(
             (ep1/k1+alpha)*k*q + muT * ( dx(k)*dx(q) + dy(k)*dy(q) ))
//   + int1d(Th,aa,bb)(temp*q*0.1)
  + int2d(Th) ( prodk*q-alpha*convect([up1,up2],-dt,k1)*q )
   + on(ff,k=0.0001)  + on(aa,bb,k=numu*stress) ;

 problem viscturb(ep,q)= // Get the rate of turbulent viscous energy
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
```
