MPI Parallel version

A first attempt of parallelization of FreeFem++ is made here with __`:::freefem mpi`__. An extended interface with MPI has been added to FreeFem++ version 3.5,  (see the MPI documentation for the functionality of the language at [http://www.mpi-forum.org/docs/mpi21-report.pdf](http://www.mpi-forum.org/docs/mpi21-report.pdf).


## MPI Keywords

The following keywords and concepts are used:

* `:::freefem mpiComm` to defined a _communication world_
* `:::freefem mpiGroup` to defined a group of _processors_ in the communication world
* `:::freefem mpiRequest` to defined a request to wait for the end of the communication


## MPI Constants

* `:::freefem mpisize` The total number of  _processes_,
* `:::freefem mpirank` the id-number of my current process in ${0,..., mpisize-1}$,
* `:::freefem mpiUndefined` The \verb!MPI_Undefined! constant, $\codered$ CORRIGER LES VERB!
* `:::freefem mpiAnySource` The \verb!MPI_ANY_SOURCE! constant,
* `:::freefem mpiCommWorld` The \verb!MPI_COMM_WORLD! constant,
* [ ... ] and all the keywords of `:::freefem MPI_Op` for the _reduce_ operator:
	`:::freefem mpiMAX`, `:::freefem mpiMIN`, `:::freefem mpiSUM`, `:::freefem mpiPROD`, `:::freefem mpiLAND`, `:::freefem mpiLOR`, `:::freefem mpiLXOR`, `:::freefem mpiBAND`, `:::freefem mpiBXOR`.

## MPI Constructor

```freefem
int[int] proc1=[1,2,3],proc2=[0,4];
mpiGroup grp(procs); // set `:::freefem MPI\_Group` to proc 1,2,3 in `:::freefem MPI\_COMM\_WORLD`
mpiGroup grp1(comm,proc1); // set `:::freefem MPI\_Group` to proc 1,2,3 in comm
mpiGroup grp2(grp,proc2); // set `:::freefem MPI\_Group` to  grp union proc1

mpiComm  comm=mpiCommWorld; //  set a `:::freefem MPI\_Comm` to `:::freefem MPI\_COMM\_WORLD`
mpiComm  ncomm(mpiCommWorld,grp); // set the `:::freefem MPI\_Comm` form grp
// `:::freefem MPI\_COMM\_WORLD`
mpiComm  ncomm(comm,color,key); // \verb!MPI_Comm_split(MPI_Comm comm,!
// \verb! int color, int key, MPI_Comm *ncomm)!
mpiComm  nicomm(processor(local_comm,local_leader),
                processor(peer_comm,peer_leader),tag);
// build \verb! MPI_INTERCOMM_CREATE(local_comm, local_leader, peer_comm,!
// \verb! remote_leader, tag, &nicomm)!
mpiComm  ncomm(intercomm,hight) ; // build using
// \verb!MPI_Intercomm_merge( intercomm, high, &ncomm)!
mpiRequest rq; // defined  an \verb!MPI_Request!
mpiRequest[int] arq(10); // defined an array of 10 \verb!MPI_Request!
```

## MPI Functions

```freefem
mpiSize(comm) ; // return the size of comm (int)
mpiRank(comm) ; // return the rank  in comm (int)

processor(i) // return  processor  i with no Resquest in \verb!MPI_COMM_WORLD!
processor(mpiAnySource) // return processor `:::freefem any source`
// with no Resquest in \verb!MPI_COMM_WORLD!
processor(i,comm) // return processor i with  no Resquest in comm
processor(comm,i) // return processor i with  no Resquest in comm
processor(i,rq,comm) // return  processor i with Resquest rq  in comm
processor(i,rq) // return  processor i with Resquest rq in
// \verb!MPI_COMM_WORLD!
processorblock(i) // return  processor i  in \verb!MPI_COMM_WORLD!
// in block mode for synchronously communication
processorblock(mpiAnySource) // return processor  `:::freefem any source`
// in \verb!MPI_COMM_WORLD!  in block mode for synchronously communication
processorblock(i,comm) // return processor i in in comm  in block mode

mpiBarrier(comm) ; // do a  `:::freefem MPI\_Barrier` on communicator `:::freefem comm`,
mpiWait(rq); // wait on of Request,
mpiWaitAll(arq); // wait add of Request array,
mpiWtime() ; //  return MPIWtime in second (real),
mpiWtick() ; // return MPIWTick in second (real),
```

where a `:::freefem processor` is just  a integer rank, pointer to a \verb!MPI_comm! and pointer to a \verb!MPI_Request!, and `:::freefem processorblock` with a special \verb!MPI_Request!.

## MPI Communicator operator

```freefem
int status;// to get the MPI status of send / recv
processor(10) << a << b; // send a,b asynchronously to the process 1,
processor(10) >> a >> b;//receive a,b synchronously from the process 10,
broadcast(processor(10,comm),a); // broadcast from processor
// of com  to other comm processor
status=Send( processor(10,comm) , a);// send synchronously

//to the process 10 the data a
status=Recv( processor(10,comm) , a);// receive synchronously
// from the process 10 the data a;
status=Isend( processor(10,comm) , a);// send asynchronously to
// the process 10 , the data a without request
status=Isend( processor(10,rq,comm) , a) ; // send asynchronously to to
// the process 10, the data a with request
status=Irecv( processor(10,rq) , a) ; // receive synchronously from
// the process 10, the data a;
 status=Irecv( processor(10) , a) ; // Error
// Error asynchronously without request.
broadcast(processor(comm,a)); //  Broadcast to all process of `:::freefem comm`
```

where the data type of `:::freefem a` can be of type of `:::freefem int`,`:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem double[int]`, `:::freefem complex[int]`, `:::freefem int[int,int]`, `:::freefem double[int,int]`, `:::freefem complex[int,int]`, `:::freefem mesh`, `:::freefem mesh3`, `:::freefem mesh[int]`, `:::freefem mesh3[int]`, `:::freefem matrix`, `:::freefem matrix<complex>`

% because in this case  the communication are multiple (header + data).
$\codered$

```freefem
processor(10,rq) << a ; // send asynchronously to the process 10
// the data a  with request
processor(10,rq) >> a ; // receive asynchronously  from the process 10
//the data a with request
```

If \verb!a,b! are arrays or full matrices of int, real, or complex, we can use the following MPI functions:

```freefem
mpiAlltoall(a,b[,comm]);
mpiAllgather(a,b[,comm]);
mpiGather(a,b,processor(..) );
mpiScatter(a,b,processor(..));
mpiReduce(a,b,processor(..),mpiMAX);
mpiAllReduce(a,b,comm, mpiMAX);
mpiReduceScatter(a,b,comm, mpiMAX);
```

See the `:::freefem examples++-mpi/essai.edp` $\codered$ to test of all this functionality and thank you to Guy-Antoine Atenekeng Kahou for his help to code this interface.

## Schwarz example in parallel
This example is a rewritting of example `:::freefem schwarz-overlap` in section \ref{schwarz-overlap} $\codered$.

```freefem
[examples++-mpi] Hecht%lamboot
LAM 6.5.9/MPI 2 C++/ROMIO - Indiana University
[examples++-mpi] hecht% mpirun -np 2 FreeFem++-mpi schwarz-c.edp
```

```freefem
//  a new coding version c, methode de schwarz in parallele
// with 2 proc.
//  -------------------------------
// F.Hecht december 2003
// ----------------------------------
//  to test the broadcast instruction
//  and array of mesh
//  add add the stop test
//  ---------------------------------

if ( mpisize != 2 ) {
    cout << " sorry, number of processors !=2 " << endl;
     exit(1);}
verbosity=3;
int interior = 2;
int exterior = 1;
border a(t=1,2){x=t;y=0;label=exterior;};
border b(t=0,1){x=2;y=t;label=exterior;};
border c(t=2,0){x=t ;y=1;label=exterior;};
border d(t=1,0){x = 1-t; y = t;label=interior;};
border e(t=0, pi/2){ x= cos(t); y = sin(t);label=interior;};
border e1(t=pi/2, 2*pi){ x= cos(t); y = sin(t);label=exterior;};
int n=4;
mesh[int]  Th(mpisize);
if (mpirank == 0)
 Th[0] = buildmesh( a(5*n) + b(5*n) + c(10*n) + d(5*n));
else
 Th[1] = buildmesh ( e(5*n) + e1(25*n) );

broadcast(processor(0),Th[0]);
broadcast(processor(1),Th[1]);

fespace Vh(Th[mpirank],P1);
fespace Vhother(Th[1-mpirank],P1);

Vh u=0,v;
Vhother U=0;
int i=0;

problem pb(u,v,init=i,solver=Cholesky) =
    int2d(Th[mpirank])( dx(u)*dx(v)+dy(u)*dy(v) )
  - int2d(Th[mpirank])( v)
  + on(interior,u = U)  +  on(exterior,u= 0 ) ;

for ( i=0 ;i< 20; i++)
{
  cout << mpirank << " looP " << i << endl;
   pb;
   //  send u  to the other proc, receive in U
   processor(1-mpirank) << u[];   processor(1-mpirank) >> U[];
   real err0,err1;
   err0 = int1d(Th[mpirank],interior)(square(U-u)) ;
   // send err0  to the other proc, receive in err1
   processor(1-mpirank)<<err0;   processor(1-mpirank)>>err1;
   real err= sqrt(err0+err1);
   cout <<" err = " << err << " err0 = " << err0
         << ", err1 = " << err1 << endl;
   if(err<1e-3) break;
};
if (mpirank==0)
    plot(u,U,ps="uU.eps");
```

### True parallel Schwarz example

This is a explanation of the two script  `:::freefem examples++-mpi/MPIGMRES[2]D.edp` $\codered$, a Schwarz parallel with a complexity almost independent of the number of process (with a coarse grid preconditioner).

Thank you to F. Nataf.

To solve the following Poisson problem
on domain $\Omega$ with boundary $\Gamma$ in $L^2(\Omega)$ :

$$
 -\Delta u = f,  \mbox{ in } \Omega,\mbox{ and } u= g \mbox{ on } \Gamma,
$$

where $f$ and $g$ are two given functions of $L^2(\Omega)$ and of $H^{\frac12}(\Gamma)$,

Lets introduce $(\pi_i)_{i=1,.., N_p}$ a regular partition of the unity of $\Omega$, q-e-d:

$$
\pi_i \in \mathcal{C}^0(\Omega) : \quad \pi_i\ge 0 \mbox{ and } \sum_{i=1}^{N_p} \pi_i =1 .
$$

Denote $\Omega_i$ the sub domain which is the support of $\pi_i$ function and also denote $\Gamma_i$ the boundary of $\Omega_i$.

The parallel Schwarz method is $\codered$
Let $\ell=0$ the iterator and a initial guest $u^0$ respecting the boundary condition (i.e. $u^0_{|\Gamma} = g$).

\begin{eqnarray}
\label{eq:lapl} \forall i = 1 .., N_p: & \displaystyle   -\Delta u_i^\ell = f, \mbox{ in }	 \Omega_i ,&\mbox{ and } u_i^\ell= u^\ell \mbox{ on }  \Gamma_i \setminus \Gamma,\; u_i^\ell=g \mbox{ on }  \Gamma_i \cap  \Gamma     \\
\label{eq:pu1}&u^{\ell+1} = \sum_{i=1}^{N_p} \pi_i u_i^\ell &
\end{eqnarray}

After discretization with the Lagrange finite element method, with a compatible mesh ${\mathcal{T}_h}_i$ of $\Omega_i$, i. e., the exist a global mesh ${\mathcal{T}_h}$ such that  ${\mathcal{T}_h}_i$ is include in ${\mathcal{T}_h}$.

Let us denote:

* ${V_h}_i$ the finite element space corresponding to domain $\Omega_i$,
* ${\mathcal{N}_h}_i$ is the set of the degree of freedom $\sigma_i^k$,
* ${\mathcal{N}^{\Gamma_i}_{hi}}$ is the set of the degree of freedom of ${V_h}_i$ on the boundary $\Gamma_i$ of $\Omega_i$,
* $\sigma_i^k({v_h})$ is the value the degree of freedom $k$,
* ${V_{0h}}_i= \{ {v_h} \in {V_h}_i :\forall k \in {\mathcal{N}^{\Gamma_i}_{hi}}, \quad \sigma_i^k({v_h})=0 \}$,
* The conditional expression $a\;?\;b:c$ is defined like in `:::freefem C` of `:::freefem C++` language by

	$$
	a?b: c \equiv
	\left\{
	\begin{array}{l}
	\mbox{if $a$ is true then return $b$}\\
	\mbox{else return $c$}\\
	\end{array}
	\right..
	$$

!!! note
	We never use finite element space associated to the full domain $\Omega$ because it to expensive.

We have to defined to operator to build the previous algorithm:

We denote ${u_h^{\ell}}_{|i}$ the restriction of $u_h^\ell$ on ${V_h}_i$, so the discrete problem on $\Omega_i$ of problem (10.1 \ref{eq:lapl} $\codered$) is find ${u_h^{\ell}}_{i}\in {V_h}_i$ such that:

$\codered$ COMPILATION ERROR
%% FFCS: this does not compile because of \uh^ and {\mathcal{N}^\Gamma_h}i_. Just wait
%% for a working version.
%%%\begin{equation}
%%%\forall {v_h}_i\in V_{0i}: \int_{\Omega_i} \nabla {v_h}_i . \nabla \uh^{\ell}_{i} = \int_{\Omega_i} f  {v_h}_i ,\quad \forall  k \in {\mathcal{N}^{\Gamma_i}_{hi}}\;:\;  \sigma_i^k(\uh^\ell_i) =  (k\in \Gamma) \; ? \; g_i^k : \sigma_i^k(\uh^{\ell}_{|i})
%%%\end{equation}

where $g_i^k$ is the value of $g$ associated to the degree of freedom $k\in {\mathcal{N}^{\Gamma_i}_{hi}}$.

In FreeFem++, it can be written has with `:::freefem U` is the vector corresponding to ${u_h^{\ell}}_{|i}$ and the vector `:::freefem U1` is the vector corresponding to ${u_h^{\ell}}_{i}$ is the solution of:

```freefem
real[int] U1(Ui.n);
real[int] b= onG .* U;
b  = onG ? b : Bi ;
U1 = Ai^-1*b;
```

where $\mathtt{onG}[i] =(i \in \Gamma_i\setminus\Gamma) ? 1 : 0$, and $\mathtt{Bi}$ the right of side of the problem, are defined by

```freefem
fespace Whi(Thi,P2); // def of the Finite element space.
varf vPb(U,V)= int3d(Thi)(grad(U)'*grad(V)) + int3d(Thi)(F*V) +on(1,U=g) + on(10,U=G);
varf vPbon(U,V)=on(10,U=1)+on(1,U=0);
matrix Ai = vPb(Whi,Whi,solver=sparsesolver);
real[int] onG = vPbon(0,Whi);
real[int] Bi=vPb(0,Whi);
```

where the FreeFem++ label of $\Gamma$ is 1 and the label of $\Gamma_i\setminus \Gamma$ is $10$.

To build the transfer/update part corresponding to (10.2 \ref{eq:pu1} $\codered$) equation on process $i$, let us call `:::freefem njpart` the number the neighborhood of domain of $\Omega_i$ (i.e: $\pi_j$ is none $0$ of $\Omega_i$), we store in an array `:::freefem jpart` of size `:::freefem njpart` all this  neighborhood.
Let us introduce two array of matrix, `:::freefem Smj[j]` to defined the vector to send from $i$ to $j$ a neighborhood process, and the matrix $rMj[j]$ to after to reduce owith neighborhood $j$ domain.

So the tranfert and update part compute $v_i= \pi_i u_i + \sum_{j\in J_i} \pi_j u_j$ and can be write the
FreeFem++ function Update:

```freefem
func bool Update(real[int] &ui, real[int] &vi)
{ int n= jpart.n;
  for(int j=0;j<njpart;++j)  Usend[j][]=sMj[j]*ui;
  mpiRequest[int] rq(n*2);
  for (int j=0;j<n;++j) Irecv(processor(jpart[j],comm,rq[j  ]), Ri[j][]);
  for (int j=0;j<n;++j) Isend(processor(jpart[j],comm,rq[j+n]), Si[j][]);
  for (int j=0;j<n*2;++j) int k= mpiWaitAny(rq);
  // apply the unity local partition .
   vi = Pii*ui;// set to $ \pi_i u_i$
   for(int j=0;j<njpart;++j)  vi += rMj[j]*Vrecv[j][]; // add $\pi_j  u_j$
   return true; }
```

where the buffer are defined by:

```freefem
InitU(njpart,Whij,Thij,aThij,Usend) // defined the send buffer
InitU(njpart,Whij,Thij,aThij,Vrecv) // defined the revc buffer
```

with the following macro definition:

```freefem
macro InitU(n,Vh,Th,aTh,U) Vh[int] U(n); for(int j=0;j<n;++j) {Th=aTh[j];  U[j]=0;}
```

__ First `:::freefem gmres` algorithm:__  you can easily accelerate the fixe point algorithm by using a parallel `:::freefem GMRES` algorithm after the introduction the following affine $\mathcal{A}_i$ operator
sub domain $\Omega_i$.

```freefem
func real[int] DJ0(real[int]& U) {
 real[int] V(U.n), b= onG .* U;
  b = onG ? b : Bi ;
  V = Ai^-1*b;
  Update(V,U);
  V -= U; return V; }
```

Where the  parallel `:::freefem MPIGMRES` or `:::freefem MPICG`  algorithm is just a simple way to solve in parallel the following $A_i x_i = b_i, i = 1, .., N_p$ by just changing the dot product by reduce the local dot product of all process with the following MPI code:

```freefem
template<class R> R ReduceSum1(R s,MPI_Comm * comm)
{   R r=0;
    MPI_Allreduce( &s, &r, 1 ,MPI_TYPE<R>::TYPE(),   MPI_SUM,  *comm );
    return r; }
```

This is done in `:::freefem MPIGC` dynamics library tool.

__ Second `:::freefem gmres` algorithm:__ Use scharwz algorithm as a preconditioner of basic GMRES
method to solving the parallel problem.

```freefem
func real[int] DJ(real[int]& U) // the original problem
{
  ++kiter;
  real[int] V(U.n);
   V =  Ai*U;
  V = onGi ? 0.: V;  // remove boundary term ...
  return V;
}

func real[int] PDJ(real[int]& U) // the preconditioner
{
  real[int] V(U.n);
  real[int] b= onG ? 0. : U;
  V =  Ai^-1*b;
  Update(V,U);
  return U;
}
```

__ Third  `:::freefem gmres` algorithm:__ Add a coarse solver to the previous algorithm

First build a coarse grid on processor 0, and the

```freefem
matrix AC,Rci,Pci;//
if(mpiRank(comm)==0)
  AC = vPbC(VhC,VhC,solver=sparsesolver);// the corase problem

Pci = interpolate(Whi,VhC); // the projection on coarse grid.
Rci = Pci'*Pii; // the Restiction on Process $i$  grid with the partition  $\pi_i$

func bool CoarseSolve(real[int]& V,real[int]& U,mpiComm& comm)
{
   // solvibg the coarse probleme
   real[int] Uc(Rci.n),Bc(Uc.n);
   Uc= Rci*U;
   mpiReduce(Uc,Bc,processor(0,comm),mpiSUM);
   if(mpiRank(comm)==0)
      Uc = AC^-1*Bc;
	  broadcast(processor(0,comm),Uc);
   V = Pci*Uc;
}
```

The New preconditionner

```freefem
func real[int] PDJC(real[int]& U) //
{ // Precon  C1= Precon //, C2  precon Coarse
// Idea : F. Nataf.
  //  0 ~  (I C1A)(I-C2A) => I ~  - C1AC2A +C1A +C2A
  //  New Prec P= C1+C2 - C1AC2   = C1(I- A C2) +C2
  // (  C1(I- A C2) +C2 ) Uo
  //   V =  - C2*Uo
  // ....
  real[int] V(U.n);
  CoarseSolve(V,U,comm);
  V = -V; //  -C2*Uo
  U  += Ai*V; // U =  (I-A C2) Uo
  real[int] b= onG ? 0. :  U;
  U =  Ai^-1*b;	// ( C1( I -A C2) Uo
  V = U -V; //
  Update(V,U);
  return U;
}
```

The code to the 4 algorithms:

```freefem
real epss=1e-6;
int rgmres=0;
if(gmres==1)
  {
   rgmres=MPIAffineGMRES(DJ0,u[],veps=epss,nbiter=300,comm=comm,
                         dimKrylov=100,verbosity=ipart ? 0: 50);
   real[int] b= onG .* u[];
   b  = onG ? b : Bi ;
   v[] = Ai^-1*b;
   Update(v[],u[]);
  }
else if(gmres==2)
  rgmres= MPILinearGMRES(DJ,precon=PDJ,u[],Bi,veps=epss,nbiter=300,comm=comm
                        ,dimKrylov=100,verbosity=ipart ? 0: 50);
else if(gmres==3)
   rgmres= MPILinearGMRES(DJ,precon=PDJC,u[],Bi,veps=epss,nbiter=300,comm=comm,
                          dimKrylov=100,verbosity=ipart ? 0: 50);
else // algo Shwarz for demo ...
   for(int iter=0;iter <10; ++iter)
     ....
```

We  have all ingredient to solve in parallel if we have et the partitions of the unity.
To build this partition we do:
the initial step on process $1$ to build a coarse mesh, ${\mathcal{T}_h}^*$ of the full domain, and build the partition $\pi$  function constant equal to $i$ on each sub domain $\mathcal{O}_i, i =1 ,.., N_p$, of the grid with the `:::freefem Metis` graph partitioner \cite{metis} $\codered$ and on each process $i$ in $1..,N_p$ do

1. Broadcast from process $1$, the mesh ${\mathcal{T}_h}^*$ (call `:::freefem Thii` in FreeFem++ script), and $\pi$ function,

2. remark that the characteristic function  $\mathrm{1\!\!I}_{\mathcal{O}_i}$ of domain $\mathcal{O}_i$, is defined by $(\pi=i)?1:0$,

3. Let us call $\Pi^2_P$ (resp. $\Pi^2_V$) the $L^2$ on $P_h^*$ the space of the constant finite element function per element on ${\mathcal{T}_h}^*$ (resp. $V_h^*$ the space of the affine continuous finite element per element on ${\mathcal{T}_h}^*$). and build in parallel the  $\pi_i$ and $\Omega_i$, such that $\mathcal{O}_i\ \subset \Omega_i$ where $\mathcal{O}_i= supp ((\Pi^2_V \Pi^2_C)^m \mathrm{1\!\!I}_{O_i})$, and $m$ is a the overlaps size on the coarse mesh (generally one), (this is done in function \verb!AddLayers(Thii,suppii[],nlayer,phii[]);! We choose a function $\pi^*_i = (\Pi^2_1 \Pi^2_0)^m \mathrm{1\!\!I}_{\mathcal{O}_i}$
 so the partition of the unity is simply defined by

\begin{equation}
  \pi_i = \frac{\pi_i^*}{\sum_{j=1}^{N_p} \pi_j^*}
\end{equation}

The set $J_i$ of neighborhood of the domain $\Omega_i$, and the local version on $V_{hi}$ can be defined the array `:::freefem jpart` and `:::freefem njpart` with:

\def\piff,#1{$\pi^*_#1$}
$\codered$

```freefem
Vhi pii=\piff,i ; Vhi[int] pij(npij); // local partition of 1 = pii + $\sum_j$ pij[j]
int[int] jpart(npart); int njpart=0;
Vhi sumphi =  \piff,i  ;
for (int i=0;i<npart;++i)
    if(i != ipart ) {
       if(int3d(Thi)( \piff,j)>0) {
         pij[njpart]=\piff,j;
         sumphi[] += pij[njpart][];
         jpart[njpart++]=i;}}}
       pii[]=pii[] ./ sumphi[];
      for (int j=0;j<njpart;++j) pij[j][] = pij[j][] ./ sumphi[];
      jpart.resize(njpart);
 ```

 4. We call ${\mathcal{T}_h}^*_{ij}$ the sub mesh part of ${\mathcal{T}_h}_i$ where $\pi_j$ are none zero. and thanks to the function `:::freefem trunc` to build this array,

```freefem
for(int jp=0;jp<njpart;++jp)
	aThij[jp] = trunc(Thi,pij[jp]>1e-10,label=10);
```

$\codered$

%\item we exchange all sub mesh ${\mathcal{T}_h}_{ij}$ of the process $i$ to to all process $j$
%and we call the mesh on process $i$ the mesh ${\mathcal{T}_h}_{ji}$, this part is the most tricky
%part because we send / receive large  complex data form all to all, with no real order.
%The only way do  that in MPI is to
%make asynchronous  Isend / Irecv with a wait at end otherwise we break all MPI buffer.
%This imply the implementation of this  send/recv/wait
% of meshes (not so simple code), but after this we can write
%

```freefem
%     macro  ISendRecvAny(comm,jpart,Si,Ri) {
%       int n= jpart.n,nn=n+n;
%       mpiRequest[int] rq(nn);
%       for (int j=0;j<n;++j)
%         Irecv(processor(jpart[j],comm,rq[j]),Ri[j]);
%       for (int j=0;j<n;++j)
%         Isend(processor(jpart[j],comm,rq[n+j]),Si[j]);
%       for (int j=0;j<nn;++j)
%         mpiWaitAny(rq); }//EndofMacro
%
%     ISendRecvAny(comm,jpart,aThij,Thij)
```

5. At this step we have all on the coarse mesh, so we can build the fine final mesh by splitting all meshes : \verb!Thi, Thij[j],Thij[j]! with FreeFem++ `:::freefem trunc` mesh function which do restriction and slipping.

6. The construction of the send/recv matrices \verb!sMj! and \verb!rMj! : can done with this code:

```freefem
mesh3 Thij=Thi; // variable meshes
fespace Whij(Thij,Pk);// variable fespace ..
matrix Pii; Whi wpii=pii; Pii = wpii[]; // Diagonal matrix corresponding $\times \pi_i$
matrix[int] sMj(njpart), rMj(njpart); // M send/rend case..
 for(int jp=0;jp<njpart;++jp)
    { int j=jpart[jp];
      Thij = aThij[jp];//change mesh to change Whij,Whij
      matrix I = interpolate(Whij,Whi); // Whij <- Whi
      sMj[jp] = I*Pii;  // Whi -> s Whij
      rMj[jp] = interpolate(Whij,Whi,t=1); }} // Whij -> Whi
```

To buil a not too bad application, I have add code tout change variable from parametre value with the following code

```freefem
include "getARGV.idp"
verbosity=getARGV("-vv",0);
int vdebug=getARGV("-d",1);
int ksplit=getARGV("-k",10);
int nloc = getARGV("-n",25);
string sff=getARGV("-p,","");
int gmres=getARGV("-gmres",3);
bool dplot=getARGV("-dp",0);
int nC = getARGV("-N" ,max(nloc/10,4));
```

$\codered$
%%\include{docFFGUI}

And small include to make graphic in parallel of distribute solution of vector $u$ on mesh $T_h$ with the following interface:

```freefem
include "MPIplot.idp"
func bool plotMPIall(mesh &Th,real[int] & u,string cm)
{ PLOTMPIALL(mesh,Pk, Th, u,{cmm=cm,nbiso=20,fill=1,dim=3,value=1}); return 1;}
```

!!! note
	The `:::freefem cmm=cm,  ... =1` in the macro argument is a way to quote macro argument so the argument is `:::freefem cmm=cm, ... =1`.
