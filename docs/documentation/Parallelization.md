A first attempt of parallelization of FreeFem++ is made here with __`:::freefem mpi`__. An extended interface with MPI has been added to FreeFem++ version 3.5, (see the [MPI documentation](http://mpi-forum.org/docs/) for the functionality of the language).

## MPI

### MPI Keywords

The following keywords and concepts are used:

* `:::freefem mpiComm` to defined a _communication world_
* `:::freefem mpiGroup` to defined a group of _processors_ in the communication world
* `:::freefem mpiRequest` to defined a request to wait for the end of the communication


### MPI Constants

* `:::freefem mpisize` The total number of _processes_,
* `:::freefem mpirank` the id-number of my current process in `{0, ..., mpisize-1}`,
* `:::freefem mpiUndefined` The `:::cpp MPI_Undefined` constant,
* `:::freefem mpiAnySource` The `:::cpp MPI_ANY_SOURCE` constant,
* `:::freefem mpiCommWorld` The `:::cpp MPI_COMM_WORLD` constant,
* [ ... ] and all the keywords of `:::freefem MPI_Op` for the _reduce_ operator:
	`:::freefem mpiMAX`, `:::freefem mpiMIN`, `:::freefem mpiSUM`, `:::freefem mpiPROD`, `:::freefem mpiLAND`, `:::freefem mpiLOR`, `:::freefem mpiLXOR`, `:::freefem mpiBAND`, `:::freefem mpiBXOR`.

### MPI Constructor

```freefem
// Parameters
int[int] proc1 = [1, 2], proc2 = [0, 3];
int color = 1;
int key = 1;

// MPI ranks
cout << "MPI rank = " << mpirank << endl;

// MPI
mpiComm comm(mpiCommWorld, 0, 0); //set a MPI_Comm to MPI_COMM_WORLD

mpiGroup grp(proc1); //set MPI_Group to proc 1,2 in MPI_COMM_WORLD
mpiGroup grp1(comm, proc1); //set MPI_Group to proc 1,2 in comm
mpiGroup grp2(grp, proc2); //set MPI_Group to grp union proc2

mpiComm ncomm1(mpiCommWorld, grp); //set the MPI_Comm form grp

// MPI_COMM_WORLD
mpiComm ncomm2(comm, color, key); //MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *ncomm)
mpiComm nicomm(processor(localcomm, localleader),
				processor(peercomm, peerleader), tag);
//build MPI_INTERCOMM_CREATE(local_comm, local_leader, peer_comm, remote_leader, tag, &nicomm)

mpiComm ncomm3(intercomm, hight); //build using MPI_Intercomm_merge(intercomm, high, &ncomm)
mpiRequest rq; //defined an MPI_Request
mpiRequest[int] arq(10); //defined an array of 10 MPI_Request
```
$\codered$ script bug

### MPI Functions

```freefem
mpiComm Comm(mpiCommWorld, 0, 0);

int MPICommSize = mpiSize(Comm);
int MPIRank = mpiRank(Comm);

if (MPIRank == 0) cout << "MPI Comm size = " << MPICommSize << endl;
cout << "MPI rank in Comm = " << mpiRank(Comm) << endl;

mpiRequest Req;
mpiRequest[int] ReqArray(10);

for (int i = 0; i < MPICommSize; i++){
	 //return processor i with no Resquest in MPI_COMM_WORLD
	processor(i);
	//return processor any source with no Resquest in MPI_COMM_WORLD
	processor(mpiAnySource);
	//return processor i with no Resquest in Comm
	processor(i, Comm);
	//return processor i with no Resquest in Comm
	processor(Comm, i);
	//return processor i with Resquest rq in Comm
	/* processor(i, Req, Comm);
	//return processor i with Resquest rq in MPI_COMM_WORLD
	processor(i, Req); */
	//return processor i in MPI_COMM_WORLD in block mode for synchronously communication
	processorblock(i);
	//return processor any source in MPI_COMM_WORLD in block mode for synchronously communication
	processorblock(mpiAnySource);
	//return processor i in in Comm in block mode
	processorblock(i, Comm);
}

mpiBarrier(Comm); //do a MPI_Barrier on communicator Comm
mpiWaitAny(ReqArray); //wait add of Request array,
mpiWait(Req); //wait on a Request
real t = mpiWtime(); //return MPIWtime in second
real tick = mpiWtick(); //return MPIWTick in second
```

where a `:::freefem processor` is just a integer rank, pointer to a `:::cpp MPI_comm` and pointer to a `:::cpp MPI_Request`, and `:::freefem processorblock` with a special `:::cpp MPI_Request`.

### MPI Communicator operator

```freefem
int status; //to get the MPI status of send / recv
real a, b;

mpiComm comm(mpiCommWorld, 0, 0);
mpiRequest req;

//send a,b asynchronously to the process 1
processor(1) << a << b;
//receive a,b synchronously from the process 10
processor(10) >> a >> b;

//broadcast from processor of comm to other comm processor
broadcast(processor(10, comm), a);
//send synchronously to the process 10 the data a
status = Send(processor(10, comm), a);
//receive synchronously from the process 10 the data a
status = Recv(processor(10, comm), a);

//send asynchronously to the process 10 the data a without request
status = Isend(processor(10, comm), a);
//send asynchronously to the process 10 the data a with request
status = Isend(processor(10, req, comm), a);
//receive asynchronously from the process 10 the data a
status = Irecv(processor(10, req), a);
//Error asynchronously without request.
status = Irecv(processor(10), a);
//Broadcast to all process comm
broadcast(processor(comm, a));
```
$\codered$ script bug

where the data type of `:::freefem a` can be of type of `:::freefem int`,`:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem real[int]`, `:::freefem complex[int]`, `:::freefem int[int,int]`, `:::freefem double[int,int]`, `:::freefem complex[int,int]`, `:::freefem mesh`, `:::freefem mesh3`, `:::freefem mesh[int]`, `:::freefem mesh3[int]`, `:::freefem matrix`, `:::freefem matrix<complex>`

```freefem
//send asynchronously to the process 10 the data a with request
processor(10, req) << a ;
//receive asynchronously from the process 10 the data a with request
processor(10, req) >> a ;
```

If `:::freefem a, b` are arrays or full matrices of `:::freefem int`, `:::freefem real`, or `:::freefem complex`, we can use the following MPI functions:

```freefem
mpiAlltoall(a, b, [comm]);
mpiAllgather(a, b, [comm]);
mpiGather(a, b, processor(...));
mpiScatter(a, b, processor(...));
mpiReduce(a, b, processor(...), mpiMAX);
mpiAllReduce(a, b, comm, mpiMAX);
mpiReduceScatter(a, b, comm, mpiMAX);
```
$\codered$ mpiReduceScatter is commented in parallelempi.cpp

See the `:::freefem examples++-mpi/essai.edp` $\codered$ to test of all this functionality and thank you to Guy-Antoine Atenekeng Kahou for his help to code this interface.

### Schwarz example in parallel
This example is a rewritting of example [Schwarz overlapping](../models/DomainDecomposition/#schwarz-overlapping).

```bash
ff-mpirun -np 2 SchwarzParallel.edp
# OR
mpirun -np 2 FreeFem++-mpi SchwarzParallel.edp
```

```freefem
if (mpisize != 2){
	cout << " sorry, number of processors !=2 " << endl;
 	exit(1);
}

// Parameters
verbosity = 0;
int interior = 2;
int exterior = 1;
int n = 4;

// Mesh
border a(t=1, 2){x=t; y=0; label=exterior;}
border b(t=0, 1){x=2; y=t; label=exterior;}
border c(t=2, 0){x=t; y=1; label=exterior;}
border d(t=1, 0){x=1-t; y=t; label=interior;}
border e(t=0, pi/2){x=cos(t); y=sin(t); label=interior;}
border e1(t=pi/2, 2*pi){x=cos(t); y=sin(t); label=exterior;}
mesh[int] Th(mpisize);
if (mpirank == 0)
	Th[0] = buildmesh(a(5*n) + b(5*n) + c(10*n) + d(5*n));
else
	Th[1] = buildmesh(e(5*n) + e1(25*n));

broadcast(processor(0), Th[0]);
broadcast(processor(1), Th[1]);

// Fespace
fespace Vh(Th[mpirank], P1);
Vh u = 0, v;

fespace Vhother(Th[1-mpirank], P1);
Vhother U = 0;

//Problem
int i = 0;
problem pb (u, v, init=i, solver=Cholesky)
	= int2d(Th[mpirank])(
		  dx(u)*dx(v)
		+ dy(u)*dy(v)
	)
	- int2d(Th[mpirank])(
		  v
	)
	+ on(interior, u=U)
	+ on(exterior, u= 0 )
	;

// Loop
for (i = 0; i < 20; i++){
	cout << mpirank << " - Loop " << i << endl;

	// Solve
	pb;
	//send u to the other proc, receive in U
	processor(1-mpirank) << u[]; processor(1-mpirank) >> U[];

	// Error
	real err0, err1;
	err0 = int1d(Th[mpirank],interior)(square(U - u));
	// send err0 to the other proc, receive in err1
	processor(1-mpirank) << err0; processor(1-mpirank) >> err1;
	real err = sqrt(err0 + err1);
	cout << " err = " << err << " - err0 = " << err0 << " - err1 = " << err1 << endl;
	if (err < 1e-3) break;
}
if (mpirank == 0)
	plot(u, U);
```
$\codered$ script bug

#### True parallel Schwarz example

$\codered$ check this part

_Thank you to F. Nataf_

This is a explanation of the two examples [MPI-GMRES 2D](../examples/#mpi-gmres-2d) and [MPI-GMRES 3D](../examples/#mpi-gmres-3d), a Schwarz parallel with a complexity almost independent of the number of process (with a coarse grid preconditioner).

To solve the following Poisson problem on domain $\Omega$ with boundary $\Gamma$ in $L^2(\Omega)$ :

\begin{eqnarray}
	-\Delta u &=& f & \mbox{ in } \Omega\\
	u &=& g & \mbox{ on } \Gamma
\end{eqnarray}

where $f$ and $g$ are two given functions of $L^2(\Omega)$ and of $H^{\frac12}(\Gamma)$,

Lets introduce $(\pi_i)_{i=1,.., N_p}$ a regular partition of the unity of $\Omega$, q-e-d:
<!--- __ --->
$$
\pi_i \in \mathcal{C}^0(\Omega) : \quad \pi_i\ge 0 \mbox{ and } \sum_{i=1}^{N_p} \pi_i =1 .
$$

Denote $\Omega_i$ the sub domain which is the support of $\pi_i$ function and also denote $\Gamma_i$ the boundary of $\Omega_i$.

The parallel Schwarz method is:

Let $\ell=0$ the iterator and a initial guest $u^0$ respecting the boundary condition (i.e. $u^0_{|\Gamma} = g$).

\begin{eqnarray}
	\forall i = 1 .., N_p:&\nonumber\\
	\displaystyle -\Delta u_i^\ell &=& f &\mbox{ in } \Omega_i\label{eq:lapl}\\
	u_i^\ell &=& u^\ell & \mbox{ on } \Gamma_i \setminus \Gamma\\
	u_i^\ell &=& g & \mbox{ on } \Gamma_i \cap \Gamma
\end{eqnarray}

\begin{equation}
\label{eq:pu1}
u^{\ell+1} = \sum_{i=1}^{N_p} \pi_i u_i^\ell
\end{equation}

After discretization with the Lagrange finite element method, with a compatible mesh ${\mathcal{T}_h}_i$ of $\Omega_i$, i. e., the exist a global mesh ${\mathcal{T}_h}$ such that ${\mathcal{T}_h}_i$ is include in ${\mathcal{T}_h}$.

Let us denote:

* ${V_h}_i$ the finite element space corresponding to domain $\Omega_i$,
* ${\mathcal{N}_h}_i$ is the set of the degree of freedom $\sigma_i^k$,
* ${\mathcal{N}^{\Gamma_i}_{hi}}$ is the set of the degree of freedom of ${V_h}_i$ on the boundary $\Gamma_i$ of $\Omega_i$,
* $\sigma_i^k({v_h})$ is the value the degree of freedom $k$,
* ${V_{0h}}_i= \{ {v_h} \in {V_h}_i :\forall k \in {\mathcal{N}^{\Gamma_i}_{hi}}, \quad \sigma_i^k({v_h})=0 \}$,
* The conditional expression $a\;?\;b:c$ is defined like in `:::c C` of `:::cpp C++` language by

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
	We never use finite element space associated to the full domain $\Omega$ because it is too expensive.

We have to defined to operator to build the previous algorithm:

We denote ${u_h^{\ell}}_{|i}$ the restriction of $u_h^\ell$ on ${V_h}_i$, so the discrete problem on $\Omega_i$ of problem \eqref{eq:lapl} is find ${u_h^{\ell}}_{i}\in {V_h}_i$ such that:

\begin{equation}
\forall {v_h}_i\in V_{0i}:
\int_{\Omega_i} \nabla {v_h}_i \cdot \nabla {u_h}^{\ell}_{i}
= \int_{\Omega_i} f {v_h}_i ,\quad \forall k \in {\mathcal{N}^{\Gamma_i}_{hi}}\;:\; \sigma_i^k({u_h}^\ell_i) = (k\in \Gamma) \; ? \; g_i^k : \sigma_i^k({u_h}^{\ell}_{|i})
\end{equation}

where $g_i^k$ is the value of $g$ associated to the degree of freedom $k\in {\mathcal{N}^{\Gamma_i}_{hi}}$.

In FreeFem++, it can be written has with `:::freefem U` is the vector corresponding to ${u_h^{\ell}}_{|i}$ and the vector `:::freefem U1` is the vector corresponding to ${u_h^{\ell}}_{i}$ is the solution of:

```freefem
real[int] U1(Ui.n);
real[int] b = onG .* U;
b = onG ? b : Bi ;
U1 = Ai^-1*b;
```

where $\mathtt{onG}[i] =(i \in \Gamma_i\setminus\Gamma) ? 1 : 0$, and $\mathtt{Bi}$ the right of side of the problem, are defined by

```freefem
// Fespace
fespace Whi(Thi, P2);

// Problem
varf vPb (U, V)
	= int3d(Thi)(
		  grad(U)'*grad(V)
	)
	+ int3d(Thi)(
		  F*V
	)
	+ on(1, U=g)
	+ on(10, U=G)
	;

varf vPbon (U, V) = on(10, U=1) + on(1, U=0);

matrix Ai = vPb (Whi, Whi, solver=sparsesolver);
real[int] onG = vPbon(0, Whi);
real[int] Bi=vPb(0, Whi);
```

where the FreeFem++ label of $\Gamma$ is 1 and the label of $\Gamma_i\setminus \Gamma$ is $10$.

To build the transfer/update part corresponding to \eqref{eq:pu1} equation on process $i$, let us call `:::freefem njpart` the number the neighborhood of domain of $\Omega_i$ (i.e: $\pi_j$ is none $0$ of $\Omega_i$), we store in an array `:::freefem jpart` of size `:::freefem njpart` all this neighborhood.

Let us introduce two array of matrix, `:::freefem Smj[j]` to defined the vector to send from $i$ to $j$ a neighborhood process, and the matrix $rMj[j]$ to after to reduce owith neighborhood $j$ domain.

So the tranfert and update part compute $v_i= \pi_i u_i + \sum_{j\in J_i} \pi_j u_j$ and can be write the FreeFem++ function Update:

```freefem
func bool Update (real[int] &ui, real[int] &vi){
	int n = jpart.n;
	for (int j = 0; j < njpart; ++j) Usend[j][] = sMj[j]*ui;
	mpiRequest[int] rq(n*2);
	for (int j = 0; j < n; ++j) Irecv(processor(jpart[j], comm,rq[j]), Ri[j][]);
	for (int j = 0; j < n; ++j) Isend(processor(jpart[j], comm, rq[j+n]), Si[j][]);
	for (int j = 0; j < n*2; ++j) int k = mpiWaitAny(rq);
	// apply the unity local partition
	vi = Pii*ui; //set to pi_i u_i
	for (int j = 0; j < njpart; ++j) vi += rMj[j]*Vrecv[j][]; //add pi_j u_j
	return true;
}
```

where the buffer are defined by:

```freefem
InitU(njpart, Whij, Thij, aThij, Usend) //defined the send buffer
InitU(njpart, Whij, Thij, aThij, Vrecv) //defined the revc buffer
```

with the following macro definition:

```freefem
macro InitU(n, Vh, Th, aTh, U) Vh[int] U(n); for (int j = 0; j < n; ++j){Th = aTh[j]; U[j] = 0;}
```

__ First GMRES algorithm:__ you can easily accelerate the fixed point algorithm by using a parallel GMRES algorithm after the introduction the following affine $\mathcal{A}_i$ operator sub domain $\Omega_i$.

```freefem
func real[int] DJ0 (real[int]& U){
	real[int] V(U.n), b = onG .* U;
	b = onG ? b : Bi ;
	V = Ai^-1*b;
	Update(V, U);
	V -= U;
	return V;
}
```

Where the parallel `:::freefem MPIGMRES` or `:::freefem MPICG` algorithm is just a simple way to solve in parallel the following $A_i x_i = b_i, i = 1, .., N_p$ by just changing the dot product by reduce the local dot product of all process with the following MPI code:

```cpp
template<class R> R ReduceSum1(R s, MPI_Comm *comm){
	R r = 0;
	MPI_Allreduce(&s, &r, 1, MPI_TYPE<R>::TYPE(), MPI_SUM, *comm );
	return r;
}
```

This is done in `:::freefem MPIGC` dynamics library tool.

__ Second GMRES algorithm:__ Use scharwz algorithm as a preconditioner of basic GMRES method to solving the parallel problem.

```freefem
func real[int] DJ (real[int]& U){ //the original problem
	++kiter;
	real[int] V(U.n);
	V = Ai*U;
	V = onGi ? 0.: V; //remove boundary term
	return V;
}

func real[int] PDJ (real[int]& U){ //the preconditioner
	real[int] V(U.n);
	real[int] b = onG ? 0. : U;
	V = Ai^-1*b;
	Update(V, U);
	return U;
}
```

__ Third GMRES algorithm:__ Add a coarse solver to the previous algorithm

First build a coarse grid on processor 0, and the

```freefem
matrix AC, Rci, Pci;
if (mpiRank(comm) == 0)
	AC = vPbC(VhC, VhC, solver=sparsesolver); //the coarse problem

Pci = interpolate(Whi, VhC); //the projection on coarse grid
Rci = Pci'*Pii; //the restriction on Process i grid with the partition pi_i

func bool CoarseSolve (real[int]& V, real[int]& U, mpiComm& comm){
	// solving the coarse problem
	real[int] Uc(Rci.n), Bc(Uc.n);
	Uc = Rci*U;
	mpiReduce(Uc, Bc, processor(0, comm), mpiSUM);
	if (mpiRank(comm) == 0)
	Uc = AC^-1*Bc;
	broadcast(processor(0, comm), Uc);
	V = Pci*Uc;
}
```

The New preconditionner

```freefem
func real[int] PDJC (real[int]& U){
	// Idea: F. Nataf.
	// 0 ~ (I C1A)(I-C2A) => I ~ - C1AC2A +C1A +C2A
	// New Prec P= C1+C2 - C1AC2 = C1(I- A C2) +C2
	// ( C1(I- A C2) +C2 ) Uo
	// V = - C2*Uo
	// ....
	real[int] V(U.n);
	CoarseSolve(V, U, comm);
	V = -V; //-C2*Uo
	U += Ai*V; //U = (I-A C2) Uo
	real[int] b = onG ? 0. : U;
	U = Ai^-1*b; //C1( I -A C2) Uo
	V = U - V;
	Update(V, U);
	return U;
}
```

The code of the 4 algorithms:

```freefem
real epss = 1e-6;
int rgmres = 0;
if (gmres == 1){
	rgmres = MPIAffineGMRES(DJ0, u[], veps=epss, nbiter=300,
		comm=comm, dimKrylov=100, verbosity=ipart?0: 50);
	real[int] b = onG .* u[];
	b = onG ? b : Bi ;
	v[] = Ai^-1*b;
	Update(v[], u[]);
}
else if (gmres == 2)
	rgmres = MPILinearGMRES(DJ, precon=PDJ, u[], Bi, veps=epss,
		nbiter=300, comm=comm, dimKrylov=100, verbosity=ipart?0: 50);
else if (gmres == 3)
	rgmres = MPILinearGMRES(DJ, precon=PDJC, u[], Bi, veps=epss,
		nbiter=300, comm=comm, dimKrylov=100, verbosity=ipart?0: 50);
else //algo Shwarz for demo
	for(int iter = 0; iter < 10; ++iter)
		...
```

We have all ingredient to solve in parallel if we have et the partitions of the unity. To build this partition we do:

The initial step on process $1$ to build a coarse mesh, ${\mathcal{T}_h}^*$ of the full domain, and build the partition $\pi$ function constant equal to $i$ on each sub domain $\mathcal{O}_i, i =1 ,.., N_p$, of the grid with the `:::freefem metis` graph partitioner [KARYPIS1995](#KARYPIS1995) and on each process $i$ in $1..,N_p$ do

<!--- ** --->

1. Broadcast from process $1$, the mesh ${\mathcal{T}_h}^*$ (call `:::freefem Thii` in FreeFem++ script), and $\pi$ function,
<!--- *** --->

2. remark that the characteristic function $\mathrm{1\!\!I}_{\mathcal{O}_i}$ of domain $\mathcal{O}_i$, is defined by $(\pi=i)?1:0$,

3. Let us call $\Pi^2_P$ (resp. $\Pi^2_V$) the $L^2$ on $P_h^*$ the space of the constant finite element function per element on ${\mathcal{T}_h}^*$ (resp. $V_h^*$ the space of the affine continuous finite element per element on ${\mathcal{T}_h}^*$) and build in parallel the $\pi_i$ and $\Omega_i$, such that $\mathcal{O}_i\ \subset \Omega_i$ where $\mathcal{O}_i= supp ((\Pi^2_V \Pi^2_C)^m \mathrm{1\!\!I}_{O_i})$, and $m$ is a the overlaps size on the coarse mesh (generally one), (this is done in function `:::freefem AddLayers(Thii,suppii[],nlayer,phii[]);` We choose a function $\pi^*_i = (\Pi^2_1 \Pi^2_0)^m \mathrm{1\!\!I}_{\mathcal{O}_i}$ so the partition of the unity is simply defined by

	<!--- ** --->

	\begin{equation}
	\pi_i = \frac{\pi_i^*}{\sum_{j=1}^{N_p} \pi_j^*}
	\end{equation}

	The set $J_i$ of neighborhood of the domain $\Omega_i$, and the local version on $V_{hi}$ can be defined the array `:::freefem jpart` and `:::freefem njpart` with:

	```freefem
	Vhi pii = piistar;
	Vhi[int] pij(npij); //local partition of 1 = pii + sum_j pij[j]
	int[int] jpart(npart);
	int njpart = 0;
	Vhi sumphi = piistar;
	for (int i = 0; i < npart; ++i)
		if (i != ipart){
			if (int3d(Thi)(pijstar,j) > 0){
				pij[njpart] = pijstar;
				sumphi[] += pij[njpart][];
				jpart[njpart++] = i;
			}
		}
	pii[] = pii[] ./ sumphi[];
	for (int j = 0; j < njpart; ++j)
		pij[j][] = pij[j][] ./ sumphi[];
	jpart.resize(njpart);
 	```

 4. We call ${\mathcal{T}_h}^*_{ij}$ the sub mesh part of ${\mathcal{T}_h}_i$ where $\pi_j$ are none zero. and thanks to the function `:::freefem trunc` to build this array,
	<!--- ** --->

	```freefem
	for(int jp = 0; jp < njpart; ++jp)
		aThij[jp] = trunc(Thi, pij[jp] > 1e-10, label=10);
	```

5. At this step we have all on the coarse mesh, so we can build the fine final mesh by splitting all meshes : `:::freefem Thi, Thij[j], Thij[j]` with FreeFem++ `:::freefem trunc` mesh function which do restriction and slipping.

6. The construction of the send/recv matrices `:::freefem sMj` and `:::freefemrMj`: can done with this code:

	```freefem
	mesh3 Thij = Thi;
	fespace Whij(Thij, Pk);
	matrix Pii; Whi wpii = pii; Pii = wpii[]; //Diagonal matrix corresponding X pi_i
	matrix[int] sMj(njpart), rMj(njpart); //M send/recive case
	 for (int jp = 0; jp < njpart; ++jp){
		 int j = jpart[jp];
		Thij = aThij[jp]; //change mesh to change Whij, Whij
		matrix I = interpolate(Whij, Whi); //Whij <- Whi
		sMj[jp] = I*Pii; //Whi -> s Whij
		rMj[jp] = interpolate(Whij, Whi, t=1); //Whij -> Whi
	}
	```

To buil a not too bad application, all variables come from parameters value with the following code

```freefem
include "getARGV.idp"
verbosity = getARGV("-vv", 0);
int vdebug = getARGV("-d", 1);
int ksplit = getARGV("-k", 10);
int nloc = getARGV("-n", 25);
string sff = getARGV("-p, ", "");
int gmres = getARGV("-gmres", 3);
bool dplot = getARGV("-dp", 0);
int nC = getARGV("-N", max(nloc/10, 4));
```

And small include to make graphic in parallel of distribute solution of vector $u$ on mesh $T_h$ with the following interface:

```freefem
include "MPIplot.idp"
func bool plotMPIall(mesh &Th, real[int] &u, string cm){
	PLOTMPIALL(mesh, Pk, Th, u, {cmm=cm, nbiso=20, fill=1, dim=3, value=1});
	return 1;
}
```

!!! note
	The `:::freefem cmm=cm, ... ` in the macro argument is a way to quote macro argument so the argument is `:::freefem cmm=cm, ...`.

<!--- The part upper needs to be reviewed --->
<!--- The parts lower have been reviewed --->

## Parallel sparse solvers

Parallel sparse solvers use several processors to solve linear systems of equation. Like sequential, parallel linear solvers can be direct or iterative. In __`FreeFem++`__ both are available.

### Using parallel sparse solvers in __`FreeFem++`__

We recall that the `:::freefem solver` parameters are defined in the following commands: `:::freefem solve`, `:::freefem problem`, `:::freefem set` (setting parameter of a matrix) and in the construction of the matrix corresponding to a bilinear form. In these commands, the parameter `:::freefem solver` must be set to `:::freefem sparsesolver` for parallel sparse solver. We have added specify parameters to these command lines for parallel sparse solvers. These are

* `:::freefem lparams` : vector of integer parameters (`l` is for the `C++` type `:::cpp long`)
* `:::freefem dparams` : vector of real parameters
* `:::freefem sparams` : string parameters
* `:::freefem datafilename` : name of the file which contains solver parameters

The following four parameters are only for direct solvers and are vectors. These parameters allow the user to preprocess the matrix (see the section on [sparse direct solver](#sparse-direct-solver) for more information).

* `:::freefem permr` : row permutation (integer vector)
* `:::freefem permc` : column permutation or inverse row permutation (integer vector)
* `:::freefem scaler` : row scaling (real vector)
* `:::freefem scalec` : column scaling (real vector)

There are two possibilities to control solver parameters. The first method defines parameters with `:::freefem lparams`, `:::freefem dparams` and `:::freefem sparams` in `.edp` file.

The second one reads the solver parameters from a data file. The name of this file is specified by `:::freefem datafilename`. If `:::freefem lparams`, `:::freefem dparams`, `:::freefem sparams` or `:::freefem datafilename` is not provided by the user, the solver's default values are used.

To use parallel solver in __`FreeFem++`__, we need to load the dynamic library corresponding to this solver. For example to use [MUMPS](http://mumps.enseeiht.fr/) solver as parallel solver in __`FreeFem++`__, write in the `.edp` file `:::freefem load "MUMPS_FreeFem"`.

If the libraries are not loaded, the default sparse solver will be loaded (default sparse solver is `:::freefem UMFPACK`). The [table 1](#Tab1) gives this new value for the different libraries.

<table>
	<thead>
		<tr>
			<th colspan="3"><a name="Tab1">Table 1</a>: Default sparse solver for real and complex arithmetics when we load a parallel sparse solver library</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td rowspan="2" style="vertical-align: middle !important; font-weight: bold">Libraries</td>
			<td colspan="2" align="center" style="font-weight: bold">Default sparse solver</td>
		</tr>
		<tr>
			<td align="center" style="font-weight: bold">real</td>
			<td align="center" style="font-weight: bold">complex</td>
		</tr>
		<tr>
			<td>MUMPS_FreeFem</td>
			<td align="center">mumps</td>
			<td align="center">mumps</td>
		</tr>
		<tr>
			<td>real_SuperLU_DIST_FreeFem</td>
			<td align="center">SuperLU_DIST</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complex_SuperLU_DIST_FreeFem</td>
			<td align="center">previous solver</td>
			<td align="center">SuperLU_DIST</td>
		</tr>
		<tr>
			<td>real_pastix_FreeFem</td>
			<td align="center">PaStiX</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complex_pastix_FreeFem</td>
			<td align="center">previous solver</td>
			<td align="center">PaStiX</td>
		</tr>
		<tr>
			<td>hips_FreeFem</td>
			<td align="center">hips</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>hypre_FreeFem</td>
			<td align="center">hypre</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>parms_FreeFem</td>
			<td align="center">parms</td>
			<td align="center">previous solver</td>
		</tr>
	</tbody>
</table>

We also add functions (see [Table 2](#Tab2)) with no parameter to change the default sparse solver in the `.edp` file. To use these functions, we need to load the library corresponding to the solver. An example of using different parallel sparse solvers for the same problem is given in [Direct solvers example](../examples/#direct-solvers).

<table>
	<thead>
		<tr>
			<th colspan="3"><a name="Tab2">Table 2</a>: Functions that allow to change the default sparse solver for real and complex arithmetics and the result of these functions</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td rowspan="2" style="vertical-align: bottom !important; font-weight: bold">Function</td>
			<td colspan="2" align="center" style="font-weight: bold">default sparse solver</td>
		</tr>
		<tr>
			<td align="center" style="font-weight: bold">real</td>
			<td align="center" style="font-weight: bold">complex</td>
		</tr>
		<tr>
			<td>defaulttoMUMPS()</td>
			<td align="center">mumps</td>
			<td align="center">mumps</td>
		</tr>
		<tr>
			<td>realdefaulttoSuperLUdist()</td>
			<td align="center">SuperLU_DIST</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complexdefaulttoSuperLUdist()</td>
			<td align="center">previous solver</td>
			<td align="center">SuperLU_DIST</td>
		</tr>
		<tr>
			<td>realdefaultopastix()</td>
			<td align="center">pastix</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complexdefaulttopastix()</td>
			<td align="center">previous solver</td>
			<td align="center">pastix</td>
		</tr>
		<tr>
			<td>defaulttohips()</td>
			<td align="center">hips</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>defaulttohypre()</td>
			<td align="center">hypre</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>defaulttoparms()</td>
			<td align="center">parms</td>
			<td align="center">previous solver</td>
		</tr>
	</tbody>
</table>

!!!example "Test direct solvers"
	```freefem
	load "MUMPS_FreeFem"
	//default solver: real-> MUMPS, complex -> MUMPS
	load "real_SuperLU_DIST_FreeFem"
	//default solver: real-> SuperLU_DIST, complex -> MUMPS
	load "real_pastix_FreeFem"
	//default solver: real-> pastix, complex -> MUMPS

	// Solving with pastix
	{
		matrix A =
			[[1, 2, 2, 1, 1],
			[ 2, 12, 0, 10 , 10],
			[ 2, 0, 1, 0, 2],
			[ 1, 10, 0, 22, 0.],
			[ 1, 10, 2, 0., 22]];

		real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
		b = A*xx;
		cout << "b = " << b << endl;
		cout << "xx = " << xx << endl;

		set(A, solver=sparsesolver, datafilename="ffpastix_iparm_dparm.txt");
		cout << "solve" << endl;
		x = A^-1*b;
		cout << "b = " << b << endl;
		cout << "x = " << endl;
		cout << x << endl;
		di = xx - x;
		if (mpirank == 0){
			cout << "x-xx = " << endl;
			cout << "Linf = " << di.linfty << ", L2 = " << di.l2 << endl;
		}
	}

	// Solving with SuperLU_DIST
	realdefaulttoSuperLUdist();
	//default solver: real-> SuperLU_DIST, complex -> MUMPS
	{
		matrix A =
			[[1, 2, 2, 1, 1],
			[ 2, 12, 0, 10 , 10],
			[ 2, 0, 1, 0, 2],
			[ 1, 10, 0, 22, 0.],
			[ 1, 10, 2, 0., 22]];

		real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
		b = A*xx;
		cout << "b = " << b << endl;
		cout << "xx = " << xx << endl;

		set(A, solver=sparsesolver, datafilename="ffsuperlu_dist_fileparam.txt");
		cout << "solve" << endl;
		x = A^-1*b;
		cout << "b = " << b << endl;
		cout << "x = " << endl;
		cout << x << endl;
		di = xx - x;
		if (mpirank == 0){
			cout << "x-xx = " << endl;
			cout << "Linf = " << di.linfty << ", L2 = " << di.l2 << endl;
		}
	}

	// Solving with MUMPS
	defaulttoMUMPS();
	//default solver: real-> MUMPS, complex -> MUMPS
	{
		matrix A =
			[[1, 2, 2, 1, 1],
			[ 2, 12, 0, 10 , 10],
			[ 2, 0, 1, 0, 2],
			[ 1, 10, 0, 22, 0.],
			[ 1, 10, 2, 0., 22]];

		real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
		b = A*xx;
		cout << "b = " << b << endl;
		cout << "xx = " << xx << endl;

		set(A, solver=sparsesolver, datafilename="ffmumps_fileparam.txt");
		cout << "solving solution" << endl;
		x = A^-1*b;
		cout << "b = " << b << endl;
		cout << "x = " << endl;
		cout << x << endl;
		di = xx - x;
		if (mpirank == 0){
			cout << "x-xx = " << endl;
			cout << "Linf = " << di.linfty << ", L2 " << di.l2 << endl;
		}
	}
	```

### Sparse direct solver

In this section, we present the sparse direct solvers interfaced with __`FreeFem++`__.

#### MUMPS solver

MUltifrontal Massively Parallel Solver ([MUMPS](http://mumps.enseeiht.fr/)) is an open-source library.

This package solves linear system of the form $A \: x = b$ where $A$ is a square sparse matrix with a direct method. The square matrix considered in MUMPS can be either unsymmetric, symmetric positive definite or general symmetric.

The method implemented in MUMPS is a direct method based on a multifrontal approach. It constructs a direct factorization $A \:= \: L\:U$, $A\: = \: L^t \: D \: L$ depending of the symmetry of the matrix $A$.

MUMPS uses the following libraries : [BLAS](http://www.netlib.org/blas/), [BLACS](http://www.netlib.org/blacs/) and [ScaLAPACK](http://www.netlib.org/scalapack/).

!!!warning
	MUMPS does not solve linear system with a rectangular matrix.

__MUMPS parameters:__

There are four input parameters in [MUMPS](http://mumps.enseeiht.fr/index.php?page=doc). Two integers `:::cpp SYM` and `:::cpp PAR`, a vector of integer of size 40 `:::cpp INCTL` and a vector of real of size 15 `:::cpp CNTL`.

The first parameter gives the type of the matrix: 0 for unsymmetric matrix, 1 for symmetric positive matrix and 2 for general symmetric.

The second parameter defined if the host processor work during the factorization and solves steps : `:::cpp PAR=1` host processor working and `:::cpp PAR=0` host processor not working.

The parameter `:::cpp INCTL` and `:::cpp CNTL` is the control parameter of MUMPS. The vectors `:::cpp ICNTL` and `:::cpp CNTL` in MUMPS becomes with index 1 like vector in `Fortran`. For more details see the [MUMPS user's guide](http://mumps.enseeiht.fr/index.php?page=doc).

We describe now some elements of the main parameters of `:::cpp ICNTL` for MUMPS.

* __Input matrix parameter__
	The input matrix is controlled by parameters `ICNTL(5)` and `ICNTL(18)`. The matrix format (resp. matrix pattern and matrix entries) are controlled by `INCTL(5)` (resp. `INCTL(18)`).

	The different values of `ICNTL(5)` are 0 for assembled format and 1 for element format. In the current release of __`FreeFem++`__, we consider that FE matrix or matrix is storage in assembled format. Therefore, `INCTL(5)` is treated as 0 value.

	The main option for `ICNTL(18)`: `INCLTL(18)=0` centrally on the host processor, `ICNTL(18)=3` distributed the input matrix pattern and the entries (recommended option for distributed matrix by developer of MUMPS). For other values of `ICNTL(18)` see the [MUMPS user's guide](http://mumps.enseeiht.fr/index.php?page=doc). These values can be used also in __`FreeFem++`__.

	The default option implemented in __`FreeFem++`__ are `ICNTL(5)=0` and `ICNTL(18)=0`.

* __Preprocessing parameter__
	The preprocessed matrix $A_{p}$ that will be effectively factored is defined by
	$$
	A_{p} = P \: D_r \: A \: Q_c \ D_c P^t
	$$
	where $P$ is the permutation matrix, $Q_c$ is the column permutation, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling.

	The ordering strategy to obtain $P$ is controlled by parameter `ICNTL(7)`. The permutation of zero free diagonal $Q_c$ is controlled by parameter `ICNTL(6)`. The row and column scaling is controlled by parameter `ICNTL(18)`. These option are connected and also strongly related with `ICNTL(12)` (see the [MUMPS user's guide](http://mumps.enseeiht.fr/index.php?page=doc) for more details).

	The parameters `:::freefem permr`, `:::freefem scaler`, and `:::freefem scalec` in __`FreeFem++`__ allow to give permutation matrix($P$), row scaling ($D_r$) and column scaling ($D_c$) of the user respectively.

__Calling MUMPS in `FreeFem++`__

To call MUMPS in __`FreeFem++`__, we need to load the dynamic library `MUMPS_freefem.dylib` (MacOSX), `MUMPS_freefem.so` (Unix) or `MUMPS_freefem.dll` (Windows).

This is done in typing `:::freefem load "MUMPS_FreeFem"` in the `.edp` file. We give now the two methods to give the option of MUMPS solver in __`FreeFem++`__.

* __Solver parameters is defined in .edp file:__
	In this method, we need to give the parameters `:::freefem lparams` and `:::freefem dparams`. These parameters are defined for MUMPS by :

	`:::freefem lparams[0] = SYM`,<br>
	`:::freefem lparams[1] = PAR`,<br>
	$\forall i$ = 1,...,40, `:::freefem lparams[i+1] = ICNTL(i)`<br>
	$\forall i$ = 1,...,15, `:::freefem dparams[i-1] = CNTL(i)`

* __Reading solver parameters on a file:__

	The structure of data file for MUMPS in __`FreeFem++`__ is : first line parameter `SYM` and second line parameter `PAR` and in the following line the different value of vectors `ICNTL` and `CNTL`. An example of this parameter file is given in `:::freefem ffmumpsfileparam.txt`.

	```freefem
	0 /* SYM :: 0 for non symmetric matrix, 1 for symmetric definite positive matrix and 2 general symmetric matrix*/
	1 /* PAR :: 0 host not working during factorization and solves steps, 1 host working during factorization and solves steps*/
	-1 /* ICNTL(1) :: output stream for error message */
	-1 /* ICNTL(2) :: output for diagnostic printing, statics and warning message */
	-1 /* ICNTL(3) :: for global information */
	0 /* ICNTL(4) :: Level of printing for error, warning and diagnostic message */
	0 /* ICNTL(5) :: matrix format : 0 assembled format, 1 elemental format. */
	7 /* ICNTL(6) :: control option for permuting and/or scaling the matrix in analysis phase */
	3 /* ICNTL(7) :: pivot order strategy : AMD, AMF, metis, pord scotch*/
	77 /* ICNTL(8) :: Row and Column scaling strategy */
	1 /* ICNTL(9) :: 0 solve Ax = b, 1 solve the transposed system A^t x = b : parameter is not considered in the current release of freefem++*/
	0 /* ICNTL(10) :: number of steps of iterative refinement */
	0 /* ICNTL(11) :: statics related to linear system depending on ICNTL(9) */
	1 /* ICNTL(12) :: constrained ordering strategy for general symmetric matrix */
	0 /* ICNTL(13) :: method to control splitting of the root frontal matrix */
	20 /* ICNTL(14) :: percentage increase in the estimated working space (default 20\%)*/
	0 /* ICNTL(15) :: not used in this release of MUMPS */
	0 /* ICNTL(16) :: not used in this release of MUMPS */
	0 /* ICNTL(17) :: not used in this release of MUMPS */
	3 /* ICNTL(18) :: method for given : matrix pattern and matrix entries : */
	0 /* ICNTL(19) :: method to return the Schur complement matrix */
	0 /* ICNTL(20) :: right hand side form ( 0 dense form, 1 sparse form) : parameter will be set to 0 for FreeFem++ */
	0 /* ICNTL(21) :: 0, 1 kept distributed solution : parameter is not considered in the current release of FreeFem++ */
	0 /* ICNTL(22) :: controls the in-core/out-of-core (OOC) facility */
	0 /* ICNTL(23) :: maximum size of the working memory in Megabyte than MUMPS can allocate per working processor */
	0 /* ICNTL(24) :: control the detection of null pivot */
	0 /* ICNTL(25) :: control the computation of a null space basis */
	0 /* ICNTL(26) :: This parameter is only significant with Schur option (ICNTL(19) not zero). : parameter is not considered in the current release of FreeFem++ */
	-8 /* ICNTL(27) (Experimental parameter subject to change in next release of MUMPS) :: control the blocking factor for multiple righthand side during the solution phase : parameter is not considered in the current release of FreeFem++ */
	0 /* ICNTL(28) :: not used in this release of MUMPS*/
	0 /* ICNTL(29) :: not used in this release of MUMPS*/
	0 /* ICNTL(30) :: not used in this release of MUMPS*/
	0 /* ICNTL(31) :: not used in this release of MUMPS*/
	0 /* ICNTL(32) :: not used in this release of MUMPS*/
	0 /* ICNTL(33) :: not used in this release of MUMPS*/
	0 /* ICNTL(34) :: not used in this release of MUMPS*/
	0 /* ICNTL(35) :: not used in this release of MUMPS*/
	0 /* ICNTL(36) :: not used in this release of MUMPS*/
	0 /* ICNTL(37) :: not used in this release of MUMPS*/
	0 /* ICNTL(38) :: not used in this release of MUMPS*/
	1 /* ICNTL(39) :: not used in this release of MUMPS*/
	0 /* ICNTL(40) :: not used in this release of MUMPS*/
	0.01 /* CNTL(1) :: relative threshold for numerical pivoting */
	1e-8 /* CNTL(2) :: stopping criteria for iterative refinement */
	-1 /* CNTL(3) :: threshold for null pivot detection */
	-1 /* CNTL(4) :: determine the threshold for partial pivoting */
	0.0 /* CNTL(5) :: fixation for null pivots */
	0 /* CNTL(6) :: not used in this release of MUMPS */
	0 /* CNTL(7) :: not used in this release of MUMPS */
	0 /* CNTL(8) :: not used in this release of MUMPS */
	0 /* CNTL(9) :: not used in this release of MUMPS */
	0 /* CNTL(10) :: not used in this release of MUMPS */
	0 /* CNTL(11) :: not used in this release of MUMPS */
	0 /* CNTL(12) :: not used in this release of MUMPS */
	0 /* CNTL(13) :: not used in this release of MUMPS */
	0 /* CNTL(14) :: not used in this release of MUMPS */
	0 /* CNTL(15) :: not used in this release of MUMPS */
	```

If no solver parameter is given, we used default option of MUMPS solver.

!!!example "MUMPS example"
	A simple example of calling MUMPS in __`FreeFem++`__ with this two methods is given in the [Test solver MUMPS example](../examples/#solver-mumps).

#### SuperLU distributed solver

The package [SuperLU_DIST](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) solves linear systems using LU factorization. It is a free scientific library under BSD license.

This library provides functions to handle square or rectangular matrix in real and complex arithmetics. The method implemented in SuperLU_DIST is a supernodal method. New release of this package includes a parallel symbolic factorization. This scientific library is written in C and MPI for communications.

__SuperLU_DIST parameters:__

We describe now some parameters of SuperLU_DIST. The SuperLU_DIST library use a 2D-logical process group. This process grid is specified by $nprow$ (process row) and $npcol$ (process column) such that $N_{p} = nprow \: npcol$ where $N_{p}$ is the number of all process allocated for SuperLU_DIST.

The input matrix parameters is controlled by "matrix= " in `:::freefem sparams` for internal parameter or in the third line of parameters file. The different value are

* `:::freefem matrix=assembled` global matrix are available on all process
* `:::freefem matrix=distributedglobal` The global matrix is distributed among all the process
* `:::freefem matrix=distributed` The input matrix is distributed (not yet implemented)

The option arguments of SuperLU_DIST are described in the section Users-callable routine of the [SuperLU users' guide](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ug.pdf).

The parameter `Fact` and `TRANS` are specified in __`FreeFem++`__ interfaces to SuperLU_DIST during the different steps. For this reason, the value given by the user for this option is not considered.

The factorization LU is calculated in SuperLU_DIST on the matrix $A_p$.
$$
A_{p} = P_{c} \: P_r \: D_r \: A \: D_{c} \: P_{c}^{t}
$$
where $P_c$ and $P_r$ is the row and column permutation matrix respectively, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling.

The option argument `RowPerm` (resp. `ColPerm`) control the row (resp. column) permutation matrix. $D_r$ and $D_c$ is controlled by the parameter `DiagScale`.

The parameter `:::freefem permr`, `:::freefem permc`, `:::freefem scaler`, and `:::freefem scalec` in __`FreeFem++`__ is provided to give row permutation, column permutation, row scaling and column scaling of the user respectively.

The other parameters for LU factorization are `ParSymFact` and `ReplaceTinyPivot`. The parallel symbolic factorization works only on a power of two processes and need the `ParMetis` ordering. The default option argument of SuperLU_DIST are given in the file `ffsuperlu_dist_fileparam.txt`.

__Calling SuperLU_DIST in __`FreeFem++`____

To call SuperLU_DIST in __`FreeFem++`__, we need to load the library dynamic correspond to interface. This done by the following line `:::freefem load "real_superlu _DIST_FreeFem"` (resp. `:::freefem load "complex_superlu_DIST_FreeFem"`) for real (resp. complex) arithmetics in the file `.edp`.

__Solver parameters is defined in `.edp` file:__

To call SuperLU_DIST with internal parameter, we used the parameters `sparams`. The value of parameters of SuperLU_DIST in `sparams` are defined by :

* `nprow=1`,
* `npcol=1`,
* `matrix= distributedgloba`,
* `Fact= DOFACT`,
* `Equil=NO`,
* `ParSymbFact=NO`,
* `ColPerm= MMD_AT_PLUS_A`,
* `RowPerm= LargeDiag`,
* `DiagPivotThresh=1.0`,
* `IterRefine=DOUBLE`,
* `Trans=NOTRANS`,
* `ReplaceTinyPivot=NO`,
* `SolveInitialized=NO`,
* `PrintStat=NO`,
* `DiagScale=NOEQUIL`

This value correspond to the parameter in the file `ffsuperlu_dist_fileparam.txt`. If one parameter is not specified by the user, we take the default value of SuperLU_DIST.

__Reading solver parameters on a file:__
The structure of data file for SuperLU_DIST in __`FreeFem++`__ is given in the file `ffsuperlu_dist_fileparam.txt` (default value of the __`FreeFem++`__ interface).

```freefem
1 /* nprow : integer value */
1 /* npcol : integer value */
distributedglobal /* matrix input : assembled, distributedglobal, distributed */
DOFACT /* Fact : DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED */
NO /* Equil : NO, YES */
NO /* ParSymbFact : NO, YES */
MMD_AT_PLUS_A /* ColPerm : NATURAL, MMD_AT_PLUS_A, MMD_ATA, METIS_AT_PLUS_A, PARMETIS, MY_PERMC */
LargeDiag /* RowPerm : NOROWPERM, LargeDiag, MY_PERMR */
1.0 /* DiagPivotThresh : real value */
DOUBLE /* IterRefine : NOREFINE, SINGLE, DOUBLE, EXTRA */
NOTRANS /* Trans : NOTRANS, TRANS, CONJ*/
NO /* ReplaceTinyPivot : NO, YES*/
NO /* SolveInitialized : NO, YES*/
NO /* RefineInitialized : NO, YES*/
NO /* PrintStat : NO, YES*/
NOEQUIL /* DiagScale : NOEQUIL, ROW, COL, BOTH*/
```

If no solver parameter is given, we used default option of SuperLU_DIST solver.

!!!example
	A simple example of calling SuperLU_DIST in __`FreeFem++`__ with this two methods is given in the [Solver superLU_DIST example](../examples/#solver-superlu_dist).

#### PaStiX solver

[PaStiX](http://pastix.gforge.inria.fr/files/README-txt.html) (Parallel Sparse matrix package) is a free scientific library under CECILL-C license. This package solves sparse linear system with a direct and block ILU(k) iterative methods. This solver can be applied to a real or complex matrix with a symmetric pattern.

__PaStiX parameters:__

The input `:::freefem matrix` parameter of __`FreeFem++`__ depend on PaStiX interface. `:::freefem matrix = assembled` for non distributed matrix. It is the same parameter for SuperLU_DIST.

There are four parameters in PaStiX : `iparm`, `dparm`, `perm` and `invp`. These parameters are respectively the integer parameters (vector of size 64), real parameters (vector of size 64), permutation matrix and inverse permutation matrix respectively. `iparm` and `dparm` vectors are described in [PaStiX RefCard](https://gforge.inria.fr/docman/?group_id=186&view=listfile&dirid=246).

The parameters `:::freefem permr` and `:::freefem permc` in __`FreeFem++`__ are provided to give permutation matrix and inverse permutation matrix of the user respectively.

__Solver parameters defined in `.edp` file:__

To call PaStiX in __`FreeFem++`__ in this case, we need to specify the parameters `:::freefem lparams` and `:::freefem dparams`. These parameters are defined by :

$\forall i$ = 0,... ,63, `:::freefem lparams[i] = iparm[i]`.

$\forall i$ = 0,... ,63, `:::freefem dparams[i] = dparm[i]`.


__Reading solver parameters on a file:__

The structure of data file for PaStiX parameters in __`FreeFem++`__ is : first line structure parameters of the matrix and in the following line the value of vectors `iparm` and `dparm` in this order.

```freefem
assembled /* matrix input :: assembled, distributed global and distributed */
iparm[0]
iparm[1]
...
...
iparm[63]
dparm[0]
dparm[1]
...
...
dparm[63]
```

An example of this file parameter is given in `ffpastix_iparm_dparm.txt` with a description of these parameters. This file is obtained with the example file `iparm.txt` and `dparm.txt` including in the PaStiX package.

If no solver parameter is given, we use the default option of PaStiX solver.

!!!example
	A simple example of calling PaStiX in __`FreeFem++`__ with this two methods is given in the [Solver PaStiX example](../examples/#solver-pastix).

In [Table 3](#Tab3), we recall the different matrix considering in the different direct solvers.

<table>
	<thead>
		<tr>
			<th colspan="7"><a name="Tab3">Table 3</a>: Type of matrix used by the different direct sparse solver</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td rowspan="2" style="vertical-align: bottom !important; font-weight: bold">direct solver</td>
			<td colspan="3" align="center" style="font-weight: bold">square matrix</td>
			<td colspan="3" align="center" style="font-weight: bold">rectangular matrix</td>
		</tr>
		<tr>
			<td style="font-weight: bold">sym</td>
			<td align="center" style="font-weight: bold">sym pattern</td>
			<td align="center" style="font-weight: bold">unsym</td>
			<td align="center" style="font-weight: bold">sym</td>
			<td align="center" style="font-weight: bold">sym pattern</td>
			<td align="center" style="font-weight: bold">unsym</td>
		</tr>
		<tr>
			<td>SuperLU_DIST</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
		</tr>
		<tr>
			<td>MUMPS</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">no</td>
			<td align="center">no</td>
			<td align="center">no</td>
		</tr>
		<tr>
			<td>pastix</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">no</td>
			<td align="center">no</td>
			<td align="center">no</td>
			<td align="center">no</td>
		</tr>
	</tbody>
</table>


### Parallel sparse iterative solver

Concerning iterative solvers, we have chosen [pARMS](http://www-users.cs.umn.edu/~saad/software/pARMS/), [HIPS](http://hips.gforge.inria.fr/) and [Hypre](https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods).

Each software implements a different type of parallel preconditioner.

So, pARMS implements algebraic domain decomposition preconditioner type such as additive Schwartz [CAI1989](#CAI1989) and interface method; while HIPS implement hierarchical incomplete factorization and finally HYPRE implements multilevel preconditioner are AMG(Algebraic MultiGrid) and parallel approximated inverse.

To use one of these programs in __`FreeFem++`__, you have to install it independently of __`FreeFem++`__. It is also necessary to install the MPI communication library which is essential for communication between the processors and, in some cases, software partitioning graphs like [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) or [Scotch](http://www.labri.fr/perso/pelegrin/scotch/).

All this preconditioners are used with Krylov subspace methods accelerators.

Krylov subspace methods are iterative methods which consist in finding a solution $x$ of linear system $Ax=b$ inside the affine space $x_0+K_m$ by imposing that $b-Ax \bot \mathcal{L}_m$, where $K_m$ is Krylov subspace of dimension $m$ defined by $K_m=\{r_0, Ar_0, A^2r_0,...,A^{m-1}r_0\}$ and $\mathcal{L}_m$ is another subspace of dimension $m$ which depends on type of Krylov subspace. For example in GMRES, $\mathcal{L}_m=AK_m$.

We realized an interface which is easy to use, so that the call of these different softwares in __`FreeFem++`__ is done in the same way. You just have to load the solver and then specify the parameters to apply to the specific solvers. In the rest of this chapter, when we talk about Krylov subspace methods we mean one among GMRES, CG and BICGSTAB.

#### pARMS solver

[pARMS](http://www-users.cs.umn.edu/~saad/software/pARMS/) (parallel Algebraic Multilevel Solver) is a software developed by Youssef Saad and al at University of Minnesota.

This software is specialized in the resolution of large sparse non symmetric linear systems of equation. Solvers developed in pARMS are of type "Krylov's subspace".

It consists of variants of GMRES like FGMRES (Flexible GMRES), DGMRES (Deflated GMRES) [SAAD2003](#SAAD2003) and BICGSTAB. pARMS also implements parallel preconditioner like RAS (Restricted Additive Schwarz) [CAI1989](#CAI1989) and Schur Complement type preconditioner.

All these parallel preconditioners are based on the principle of domain decomposition. Thus, the matrix $A$ is partitioned into sub matrices $A_i$($i=1,...,p$) where p represents the number of partitions one needs. The union of $A_i$ forms the original matrix. The solution of the overall system is obtained by solving the local systems on $A_i$ (see [SMITH1996](#SMITH1996)). Therefore, a distinction is made between iterations on $A$ and the local iterations on $A_i$.

To solve the local problem on $A_i$ there are several preconditioners as __ilut__ (Incomplete LU with threshold), __iluk__ (Incomplete LU with level of fill in) and __ARMS__ (Algebraic Recursive Multilevel Solver).

!!!example "Default parameters"
	```freefem
	load "parms_FreeFem" //Tell FreeFem that you will use pARMS

	// Mesh
	border C(t=0, 2*pi){x=cos(t); y=sin(t); label=1;}
	mesh Th = buildmesh (C(50));

	// Fespace
	fespace Vh(Th, P2);
	Vh u, v;

	// Function
	func f= x*y;

	// Problem
	problem Poisson (u, v, solver=sparsesolver)
		= int2d(Th)(
			  dx(u)*dx(v)
			+ dy(u)*dy(v)
		)
		+ int2d(Th)(
			- f*v
		)
		+ on(1, u=0)
		;

	// Solve
	real cpu = clock();
	Poisson;
	cout << " CPU time = " << clock()-cpu << endl;

	// Plot
	plot(u);
	```

	In line 1, the pARMS dynamic library is loaded with interface __`FreeFem++`__. After this, in line 15 we specify that the bilinear form will be solved by the last sparse linear solver load in memory which, in this case, is pARMS.

	The parameters used in pARMS in this case are the default one since the user does not have to provide any parameter.

	!!!note
		In order to see the plot of a parallel script, run the command `FreeFem++-mpi -glut ffglut script.edp`

Here are some default parameters:

* `solver=FGMRES`,
* `Krylov dimension=30`,
* `Maximum of Krylov=1000`,
* `Tolerance for convergence=$1e-08$`(see book [SAAD2003](#SAAD2003) to understand all this parameters),
* `preconditionner=Restricted Additif Schwarz` [CAI1989](#CAI1989),
* `Inner Krylov dimension=5`,
* `Maximum of inner Krylov dimension=5`,
* `Inner preconditionner=ILUK`.

To specify the parameters to apply to the solver, the user can either give an integer vector for __integer parameters__ and real vectors for __real parameters__ or provide a __file__ which contains those parameters.

!!!example "User specifies parameters inside two vectors"
	Lets us consider Navier-Stokes example. In this example we solve linear systems coming from discretization of Navier-Stokes equations with pARMS. Parameters of solver is specified by user.

	```freefem
	load "parms_FreeFem"

	// Parameters
	real nu = 1.;
	int[int] iparm(16);
	real[int] dparm(6);
	for (int ii = 0; ii < 16; ii++)
		iparm[ii] = -1;
	for (int ii = 0; ii < 6; ii++)
		dparm[ii] = -1.0;
	iparm[0]=0;

	// Mesh
	mesh Th = square(10, 10);
	int[int] wall = [1, 3];
	int inlet = 4;

	// Fespace
	fespace Vh(Th, [P2, P2, P1]);

	// Function
	func uc = 1.;

	varf Stokes ([u, v, p], [ush, vsh, psh], solver=sparsesolver)
		= int2d(Th)(
			  nu*(
				  dx(u)*dx(ush)
				+ dy(u)*dy(ush)
				+ dx(v)*dx(vsh)
				+ dy(v)*dy(vsh)
			)
			- p*psh*(1.e-6)
			- p*(dx(ush) + dy(vsh))
			- (dx(u) + dy(v))*psh
		)
		+ on(wall, wall, u=0., v=0.)
		+ on(inlet, u=uc, v=0)
		;

	matrix AA = Stokes(Vh, Vh);
	set(AA, solver=sparsesolver, lparams=iparm, dparams=dparm); //set pARMS as linear solver
	real[int] bb = Stokes(0, Vh);
	real[int] sol(AA.n);
	sol = AA^-1 * bb;
	```

	We need two vectors to specify the parameters of the linear solver. In line 5-6 of the example, we have declared these vectors(`:::freefem int[int] iparm(16); real[int] dparm(6);`). In line 7-10 we have initialized these vectors by negative values.

	We do this because all parameters values in pARMS are positive and if you do not change the negative values of one entry of this vector, the default value will be set.

	In [table 4](#Tab4) and [table 5](#Tab5), we have the meaning of different entries of these vectors.

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab4">Table 4</a>: Meaning of `:::freefem lparams` corresponding variables</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Entries of `:::freefem iparm`</td>
				<td style="font-weight: bold">Significations of each entries</td>
			</tr>
			<tr>
				<td>`iparm[0]`</td>
				<td>Krylov subspace methods.<br>Different values for this parameters are specify on [table 7](#Tab6)</td>
			</tr>
			<tr>
				<td>`iparm[1]`</td>
				<td>Preconditionner.<br>Different preconditionners for this parameters are specify on [table 7](#Tab7)</td>
			</tr>
			<tr>
				<td>`iparm[2]`</td>
				<td>Krylov subspace dimension in outer iteration: default value 30</td>
			</tr>
			<tr>
				<td>`iparm[3]`</td>
				<td>Maximum of iterations in outer iteration: default value 1000</td>
			</tr>
			<tr>
				<td>`iparm[4]`</td>
				<td>Number of level in arms when used.</td>
			</tr>
			<tr>
				<td>`iparm[5]`</td>
				<td>Krylov subspace dimension in inner iteration: default value 3</td>
			</tr>
			<tr>
				<td>`iparm[6]`</td>
				<td>Maximum of iterations in inner iteration: default value 3</td>
			</tr>
			<tr>
				<td>`iparm[7]`</td>
				<td>Symmetric(=1 for symmetric) or unsymmetric matrix:<br>default value 0(unsymmetric matrix)</td>
			</tr>
			<tr>
				<td>`iparm[8]`</td>
				<td>Overlap size between different subdomain: default value 0(no overlap)</td>
			</tr>
			<tr>
				<td>`iparm[9]`</td>
				<td>Scale the input matrix or not: Default value 1 (Matrix should be scaled)</td>
			</tr>
			<tr>
				<td>`iparm[10]`</td>
				<td>Block size in arms when used: default value 20</td>
			</tr>
			<tr>
				<td>`iparm[11]`</td>
				<td>lfil0 (ilut, iluk, and arms) : default value 20</td>
			</tr>
			<tr>
				<td>`iparm[12]`</td>
				<td>lfil for Schur complement const : default value 20</td>
			</tr>
			<tr>
				<td>`iparm[13]`</td>
				<td>lfil for Schur complement const : default value 20</td>
			</tr>
			<tr>
				<td>`iparm[14]`</td>
				<td>Multicoloring or not in ILU when used : default value 1</td>
			</tr>
			<tr>
				<td>`iparm[15]`</td>
				<td>Inner iteration : default value 0</td>
			</tr>
			<tr>
				<td>`iparm[16]`</td>
				<td>Print message when solving:default 0 (no message print).<br>0: no message is print,<br>1: Convergence informations like number of iteration and residual ,<br>2: Timing for a different step like preconditioner<br>3 : Print all informations.</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab5">Table 5</a>: Significations of `:::freefem dparams` corresponding variables</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Entries of `:::freefem dparm`</td>
				<td style="font-weight: bold">Significations of each entries</td>
			</tr>
			<tr>
				<td>`dparm[0]`</td>
				<td>precision for outer iteration : default value 1e-08</td>
			</tr>
			<tr>
				<td>`dparm[1]`</td>
				<td>precision for inner iteration: default value 1e-2</td>
			</tr>
			<tr>
				<td>`dparm[2]`</td>
				<td>tolerance used for diagonal domain: : default value 0.1</td>
			</tr>
			<tr>
				<td>`dparm[3]`</td>
				<td>drop tolerance droptol0 (ilut, iluk, and arms) : default value 1e-2</td>
			</tr>
			<tr>
				<td>`dparm[4]`</td>
				<td>droptol for Schur complement const: default value 1e-2</td>
			</tr>
			<tr>
				<td>`dparm[5]`</td>
				<td>droptol for Schur complement const: default value 1e-2</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab6">Table 6</a>: Krylov Solvers in pARMS</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Values of `iparm[0]`</td>
				<td style="font-weight: bold">Krylov subspace methods</td>
			</tr>
			<tr>
				<td>0</td>
				<td>FGMRES (Flexible GMRES)</td>
			</tr>
			<tr>
				<td>1</td>
				<td>DGMRES (Deflated GMRES)</td>
			</tr>
			<tr>
				<td>2</td>
				<td>BICGSTAB</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab7">Table 7</a>: Preconditionners in pARMS</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Values of `iparm[1]`</td>
				<td style="font-weight: bold">Preconditionners type</td>
			</tr>
			<tr>
				<td>0</td>
				<td>additive Schwartz preconditioner with ilu0 as local preconditioner</td>
			</tr>
			<tr>
				<td>1</td>
				<td>additive Schwartz preconditioner with iluk as local preconditioner</td>
			</tr>
			<tr>
				<td>2</td>
				<td>additive Schwartz preconditioner with ilut as local preconditioner</td>
			</tr>
			<tr>
				<td>3</td>
				<td>additive Schwartz preconditioner with arms as local preconditioner</td>
			</tr>
			<tr>
				<td>4</td>
				<td>Left Schur complement preconditioner with ilu0 as local preconditioner</td>
			</tr>
			<tr>
				<td>5</td>
				<td>Left Schur complement preconditioner with ilut as local preconditioner</td>
			</tr>
			<tr>
				<td>6</td>
				<td>Left Schur complement preconditioner with iluk as local preconditioner</td>
			</tr>
			<tr>
				<td>7</td>
				<td>Left Schur complement preconditioner with arms as local preconditioner</td>
			</tr>
			<tr>
				<td>8</td>
				<td>Right Schur complement preconditioner with ilu0 as local preconditioner</td>
			</tr>
			<tr>
				<td>9</td>
				<td>Right Schur complement preconditioner with ilut as local preconditioner</td>
			</tr>
			<tr>
				<td>10</td>
				<td>Right Schur complement preconditioner with iluk as local preconditioner</td>
			</tr>
			<tr>
				<td>11</td>
				<td>Right Schur complement preconditioner with arms as local preconditioner</td>
			</tr>
			<tr>
				<td>12</td>
				<td>sch_gilu0, Schur complement preconditioner with global ilu0</td>
			</tr>
			<tr>
				<td>13</td>
				<td>SchurSymmetric GS preconditioner</td>
			</tr>
		</tbody>
	</table>

	We run this example on a cluster paradent of Grid5000 and report results in [table 8](#Tab8).

	<table>
		<thead>
			<tr>
				<th colspan="5"><a name="Tab8">Table 8</a>: Convergence and time for solving linear system</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td colspan="2" align="center" style="font-weight: bold">n=471281</td>
				<td colspan="2" align="center" style="font-weight: bold">nnz=$13\times10^6$</td>
				<td colspan="2" align="center" style="font-weight: bold">Te=571,29</td>
			</tr>
			<tr>
				<td rowspan="2" align="center" style="vertical-align: bottom !important">np</td>
				<td colspan="2" align="center">add(iluk)</td>
				<td colspan="2" align="center">schur(iluk)</td>
			</tr>
			<tr>
				<td align="center">`nit`</td>
				<td align="center">`time`</td>
				<td align="center">`nit`</td>
				<td align="center">`time`</td>
			</tr>
			<tr>
				<td align="center">4</td>
				<td align="center">230</td>
				<td align="center">637.57</td>
				<td align="center">21</td>
				<td align="center">557.8</td>
			</tr>
			<tr>
				<td align="center">8</td>
				<td align="center">240</td>
				<td align="center">364.12</td>
				<td align="center">22</td>
				<td align="center">302.25</td>
			</tr>
			<tr>
				<td align="center">16</td>
				<td align="center">247</td>
				<td align="center">212.07</td>
				<td align="center">24</td>
				<td align="center">167.5</td>
			</tr>
			<tr>
				<td align="center">32</td>
				<td align="center">261</td>
				<td align="center">111.16</td>
				<td align="center">25</td>
				<td align="center">81.5</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Ta9">Table 9</a>: Legend of [table 8](#Tab8)</th>
			</tr>
		</thead>
		<tbody>
		<tr>
			<td>`n`</td>
			<td>matrix size</td>
		</tr>
		<tr>
			<td>`nnz`</td>
			<td>number of non null entries inside matrix</td>
		</tr>
		<tr>
			<td>`nit`</td>
			<td>number of iteration for convergence</td>
		</tr>
		<tr>
			<td>`time`</td>
			<td>Time for convergence</td>
		</tr>
		<tr>
			<td>`Te`</td>
			<td>Time for constructing finite element matrix</td>
		</tr>
		<tr>
			<td>`np`</td>
			<td>number of processor</td>
		</tr>
	</table>

	In this example, we fix the matrix size (in term of finite element, we fix the mesh) and increase the number of processors used to solve the linear system. We saw that, when the number of processors increases, the time for solving the linear equation decreases, even if the number of iteration increases. This proves that, using pARMS as solver of linear systems coming from discretization of partial differential equation in __`FreeFem++`__ can decrease drastically the total time of simulation.

#### Interfacing with HIPS

[HIPS](http://hips.gforge.inria.fr/) (_Hierarchical Iterative Parallel Solver_) is a scientific library that provides an efficient parallel iterative solver for very large sparse linear systems. HIPS is available as free software under the CeCILL-C licence.

HIPS implements two solver classes which are the iteratives class (GMRES, PCG) and the Direct class. Concerning preconditionners, HIPS implements a type of multilevel ILU. For further informations on those preconditionners see the [HIPS documentation](http://hips.gforge.inria.fr/doc/hips_user.pdf).

!!!example "Laplacian 3D solved with HIPS"
	Let us consider the 3D Laplacian example inside __`FreeFem++`__ package where after discretization we want to solve the linear equation with HIPS.

	The following example is a Laplacian 3D using Hips as linear solver. We first load Hips solver at line 2. From line 7 to 18 we specify the parameters for the Hips solver and in line 82 we set these parameters in the linear solver.

	In [Table 10](#Tab10) results of running on Cluster Paradent of Grid5000 are reported. We can see in this running example the efficiency of parallelism.

	```freefem
	load "msh3"
	load "hips_FreeFem" //load Hips library

	// Parameters
	int nn = 10;
	real zmin = 0, zmax = 1;
	int[int] iparm(14);
	real[int] dparm(6);
	for (int iii = 0; iii < 14; iii++) iparm[iii] = -1;
	for (int iii = 0; iii < 6; iii++) dparm[iii] = -1;
	iparm[0] = 0; //use iterative solver
	iparm[1] = 1; //PCG as Krylov method
	iparm[4] = 0; //Matrix are symmetric
	iparm[5] = 1; //Pattern are also symmetric
	iparm[9] = 1; //Scale matrix
	dparm[0] = 1e-13; //Tolerance to convergence
	dparm[1] = 5e-4; //Threshold in ILUT
	dparm[2] = 5e-4; //Threshold for Schur preconditionner

	// Functions
	func ue = 2*x*x + 3*y*y + 4*z*z + 5*x*y + 6*x*z + 1;
	func uex = 4*x + 5*y + 6*z;
	func uey = 6*y + 5*x;
	func uez = 8*z + 6*x;
	func f = -18.;

	// Mesh
	mesh Th2 = square(nn, nn);

	int[int] rup = [0,2], rdown=[0, 1];
	int[int] rmid=[1, 1, 2, 1, 3, 1, 4, 1];

	mesh3 Th=buildlayers(Th2, nn,
		zbound=[zmin, zmax],
		reffacemid=rmid,
		reffaceup = rup,
		reffacelow = rdown);

	// Fespace
	fespace Vh2(Th2, P2);
	Vh2 ux, uz, p2;

	fespace Vh(Th, P2);
	Vh uhe = ue;
	cout << "uhe min = " << uhe[].min << ", max = " << uhe[].max << endl;
	Vh u, v;
	Vh F;

	// Macro
	macro Grad3(u) [dx(u), dy(u), dz(u)] //

	// Problem
	varf va (u, v)
		= int3d(Th)(
			  Grad3(v)' * Grad3(u)
		)
		+ int2d(Th, 2)(
			  u*v
		)
		- int3d(Th)(
			  f*v
		)
		- int2d(Th, 2)(
			  ue*v
			+ (uex*N.x + uey*N.y + uez*N.z)*v
		)
		+ on(1, u=ue);

	varf l (unused, v) = int3d(Th)(f*v);

	real cpu=clock();
	matrix Aa = va(Vh, Vh);

	F[] = va(0, Vh);

	if (mpirank == 0){
		cout << "Size of A = " << Aa.n << endl;
		cout << "Non zero coefficients = " << Aa.nbcoef << endl;
		cout << "CPU TIME FOR FORMING MATRIX = " << clock()-cpu << endl;
	}

	set(Aa, solver=sparsesolver, dparams=dparm, lparams=iparm); //Set hips as linear solver

	// Solve
	u[] = Aa^-1*F[];

	// Plot
	plot(u);
	```

	<table>
		<thead>
			<tr>
				<th colspan="3"><a name="Tab10">Table 10</a>: Legend of this table are give in [table 9](#Tab9)</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td align='center'>$n=4 \times 10^6$</td>
				<td align='center'>$nnz=118 \times 10^6$</td>
				<td align='center'>$Te=221.34$</td>
			</tr>
			<tr>
				<td align='center'>`np`</td>
				<td align='center'>`nit`</td>
				<td align='center'>`time`</td>
			</tr>
			<tr>
				<td align='center'>8</td>
				<td align='center'>190</td>
				<td align='center'>120.34</td>
			</tr>
			<tr>
				<td align='center'>16</td>
				<td align='center'>189</td>
				<td align='center'>61.08</td>
			</tr>
			<tr>
				<td align='center'>32</td>
				<td align='center'>186</td>
				<td align='center'>31.70</td>
			</tr>
			<tr>
				<td align='center'>64</td>
				<td align='center'>183</td>
				<td align='center'>23.44</td>
			</tr>
		</tbody>
	</table>

!!!tips
	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab11">Table 11</a>: Significations of `:::freefem lparams` corresponding to HIPS interface</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Entries of `iparm`</td>
				<td style="font-weight: bold">Significations of each entries</td>
			</tr>
			<tr>
				<td>`iparm[0]`</td>
				<td>Strategy use for solving (Iterative=0 or Hybrid=1 or Direct=2). Defaults values are : Iterative</td>
			</tr>
			<tr>
				<td>`iparm[1]`</td>
				<td>Krylov methods. If iparm[0]=0, give type of Krylov methods: 0 for GMRES, 1 for PCG</td>
			</tr>
			<tr>
				<td>`iparm[2]`</td>
				<td>Maximum of iterations in outer iteration: default value 1000</td>
			</tr>
			<tr>
				<td>`iparm[3]`</td>
				<td>Krylov subspace dimension in outer iteration: default value 40</td>
			</tr>
			<tr>
				<td>`iparm[4]`</td>
				<td>Symmetric(=0 for symmetric) and 1 for unsymmetricmatrix: default value 1 (unsymmetric matrix)</td>
			</tr>
			<tr>
				<td>`iparm[5]`</td>
				<td>Pattern of matrix are symmetric or not: default value 0</td>
			</tr>
			<tr>
				<td>`iparm[6]`</td>
				<td>Partition type of input matrix: default value 0</td>
			</tr>
			<tr>
				<td>`iparm[7]`</td>
				<td>Number of level that use the HIPS locally consistentfill-in: Default value 2</td>
			</tr>
			<tr>
				<td>`iparm[8]`</td>
				<td>Numbering in indices array will start at 0 or 1: Default value 0</td>
			</tr>
			<tr>
				<td>`iparm[9]`</td>
				<td>Scale matrix. Default value 1</td>
			</tr>
			<tr>
				<td>`iparm[10]`</td>
				<td>Reordering use inside subdomains for reducingfill-in: Only use for iterative. Default value 1</td>
			</tr>
			<tr>
				<td>`iparm[11]`</td>
				<td>Number of unknowns per node in the matrix non-zeropattern graph: Default value 1</td>
			</tr>
			<tr>
				<td>`iparm[12]`</td>
				<td>This value is used to set the number of time the normalization is applied to the matrix: Default 2.</td>
			</tr>
			<tr>
				<td>`iparm[13]`</td>
				<td>Level of informations printed during solving: Default 5.</td>
			</tr>
			<tr>
				<td>`iparm[14]`</td>
				<td>HIPS_DOMSIZE Subdomain size</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab12">Table 12</a>: Significations of `:::freefem dparams` corresponding to HIPS interface</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td>`dparm[0]`</td>
				<td>HIPS_PREC: Relative residual norm: Default=1e-9</td>
			</tr>
			<tr>
				<td>`dparm[1]`</td>
				<td>HIPS_DROPTOL0: Numerical threshold in ILUT for interior domain (important : set 0.0 in HYBRID: Default=0.005)</td>
			</tr>
			<tr>
				<td>`dparm[2]`</td>
				<td>HIPS_DROPTOL1 : Numerical threshold in ILUT for Schur preconditioner: Default=0.005</td>
			</tr>
			<tr>
				<td>`dparm[3]`</td>
				<td>HIPS_DROPTOLE : Numerical threshold for coupling between the interior level and Schur: Default 0.005</td>
			</tr>
			<tr>
				<td>`dparm[4]`</td>
				<td>HIPS_AMALG : Numerical threshold for coupling between the interior level and Schur: Default=0.005</td>
			</tr>
			<tr>
				<td>`dparm[5]`</td>
				<td>HIPS_DROPSCHUR : Numerical threshold for coupling between the interior level and Schur: Default=0.005</td>
			</tr>
		</tbody>
	</table>

#### Interfacing with HYPRE

[Hypre](https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) (High Level Preconditioner) is a suite of parallel preconditioner developed at Lawrence Livermore National Lab.

There are two main classes of preconditioners developed in HYPRE: AMG (Algebraic MultiGrid) and Parasails (Parallel Sparse Approximate Inverse).

Now, suppose we want to solve $Ax=b$.

At the heart of AMG there is a series of progressively coarser (smaller) representations of the matrix $A$. Given an approximation $\hat{x}$ to the solution $x$, consider solving the residual equation $Ae=r$ to find the error $e$, where $r=b-A\hat{x}$. A fundamental principle of AMG is that it is an algebraically smooth error. To reduce the algebraically smooth errors further, they need to be represented by a smaller defect equation (coarse grid residual equation) $A_ce_c=r_c$, which is cheaper to solve. After solving this coarse equation, the solution is then interpolated in fine grid represented here by matrix $A$. The quality of AMG depends on the choice of coarsening and interpolating operators.

The _sparse approximate inverse_ approximates the inverse of a matrix $A$ by a sparse matrix $M$. A technical idea to construct matrix $M$ is to minimize the Frobenuis norm of the residual matrix $I-MA$. For more details on this preconditioner technics see [CHOW1997](#CHOW1997).

HYPRE implement three Krylov subspace solvers: GMRES, PCG and BiCGStab.

!!!example "Laplacian 3D solved with HYPRE"
	Let us consider again the 3D Laplacian example inside FreeFem++ package where after discretization we want to solve the linear equation with Hypre. The following example is a Laplacian 3D using Hypre as linear solver. This is the same example as Hips one, so we just show here the lines where we set some Hypre parameters.

	We first load the Hypre solver at line 2. From line 6 to 18 we specifies the parameters to set to Hypre solver and in line 22 we set parameters to Hypre solver.

	It should be noted that the meaning of the entries of these vectors is different from those of Hips. In the case of HYPRE, the meaning of differents entries of vectors `:::freefem iparm` and `:::freefem dparm` are given in [tables 13](#Tab13) to [17](#Tab17).

	In [Table 18](#Tab18) the results of running on Cluster Paradent of Grid5000 are reported. We can see in this running example the efficiency of parallelism, in particular when AMG are use as preconditioner.


	```freefem
	load "msh3"
	load "hipre_FreeFem" //Load Hipre librairy

	// Parameters
	int nn = 10;
	int[int] iparm(20);
	real[int] dparm(6);
	for (int iii = 0; iii < 20; iii++) iparm[iii] = -1;
	for (int iii = 0; iii < 6; iii++) dparm[iii] = -1;
	iparm[0] = 2; //PCG as krylov method
	iparm[1] = 0; //AMG as preconditionner 2: if ParaSails
	iparm[7] = 7; //Interpolation
	iparm[9] = 6; //AMG Coarsen type
	iparm[10] = 1; //Measure type
	iparm[16] = 2; //Additive schwarz as smoother
	dparm[0] = 1e-13; //Tolerance to convergence
	dparm[1] = 5e-4; //Threshold
	dparm[2] = 5e-4; //Truncation factor

	...

	set(Aa, solver=sparsesolver, dparams=dparm, lparams=iparm);
	```

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab13">Table 13</a>: Definitions of common entries of `:::freefem iparms` and `:::freefem dparms` vectors for every preconditioner in HYPRE</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td>`iparms[0]`</td>
				<td>Solver identification: 0: BiCGStab, 1: GMRES, 2: PCG. Default=1</td>
			</tr>
			<tr>
				<td>`iparms[1]`</td>
				<td>Preconditioner identification: 0: BOOMER AMG, 1: PILUT, 2: Parasails, 3: Schwartz Default=0</td>
			</tr>
			<tr>
				<td>`iparms[2]`</td>
				<td>Maximum of iteration: Default=1000</td>
			</tr>
			<tr>
				<td>`iparms[3]`</td>
				<td>Krylov subspace dim: Default= 40</td>
			</tr>
			<tr>
				<td>`iparms[4]`</td>
				<td>Solver print info level: Default=2</td>
			</tr>
			<tr>
				<td>`iparms[5]`</td>
				<td>Solver log: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[6]`</td>
				<td>Solver stopping criteria only for BiCGStab : Default=1</td>
			</tr>
			<tr>
				<td>`dparms[0]`</td>
				<td>Tolerance for convergence: Default=$1.0e-11$</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab14">Table 14</a>: Definitions of other entries of `:::freefem iparms` and `:::freefem dparms` if preconditioner is BOOMER AMG</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td>`iparms[7]`</td>
				<td>AMG interpolation type: Default=6</td>
			</tr>
			<tr>
				<td>`iparms[8]`</td>
				<td>Specifies the use of GSMG - geometrically smooth coarsening and interpolation: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[9]`</td>
				<td>AMG coarsen type: Default=6</td>
			</tr>
			<tr>
				<td>`iparms[10]`</td>
				<td>Defines whether local or global measures are used: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[11]`</td>
				<td>AMG cycle type: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[12]`</td>
				<td>AMG Smoother type: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[13]`</td>
				<td>AMG number of levels for smoothers: Default=3</td>
			</tr>
			<tr>
				<td>`iparms[14]`</td>
				<td>AMG number of sweeps for smoothers: Default=2</td>
			</tr>
			<tr>
				<td>`iparms[15]`</td>
				<td>Maximum number of multigrid levels: Default=25</td>
			</tr>
			<tr>
				<td>`iparms[16]`</td>
				<td>Defines which variant of the Schwartz method isused:<br>
				0: hybrid multiplicative Schwartz method (no overlap across processor boundaries)<br>
				1: hybrid additive Schwartz method (no overlap across processor boundaries)<br>
				2: additive Schwartz method<br>
				3: hybrid multiplicative Schwartz method (with overlap across processor boundaries)<br>
				Default=1
				</td>
			</tr>
			<tr>
				<td>`iparms[17]`</td>
				<td>Size of the system of PDEs: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[18]`</td>
				<td>Overlap for the Schwarz method: Default=1</td>
			</tr>
			<tr>
				<td>Type of domain used for the Schwarz method<br>
				<td>`iparms[19]`</td>
				0: each point is a domain<br>
				1: each node is a domain (only of interest in "systems" AMG)<br>
				2: each domain is generated by agglomeration (default)</td>
			</tr>
			<tr>
				<td>`dparms[1]`</td>
				<td>AMG strength threshold: Default=0.25</td>
			</tr>
			<tr>
				<td>`dparms[2]`</td>
				<td>Truncation factor for the interpolation: Default=1e-2</td>
			</tr>
			<tr>
				<td>`dparms[3]`</td>
				<td>Sets a parameter to modify the definition of strength for diagonal dominant portions of the matrix: Default=0.9</td>
			</tr>
			<tr>
				<td>`dparms[3]`</td>
				<td>Defines a smoothing parameter for the additive Schwartz method. Default=1</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab15">Table 15</a>: Definitions of other entries of `:::freefem iparms` and `:::freefem dparms` if preconditioner is PILUT</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td>`iparms[7]`</td>
				<td>Row size in Parallel ILUT: Default=1000</td>
			</tr>
			<tr>
				<td>`iparms[8]`</td>
				<td>Set maximum number of iterations: Default=30</td>
			</tr>
			<tr>
				<td>`dparms[1]`</td>
				<td>Drop tolerance in Parallel ILUT: Default=$1e-5$</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab16">Table 16</a>: Definitions of other entries of `:::freefem iparms` and `:::freefem dparms` if preconditioner is ParaSails</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td>`iparms[7]`</td>
				<td>Number of levels in Parallel Sparse Approximate inverse: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[8]`</td>
				<td>Symmetric parameter for the ParaSails preconditioner:<br>
				0: nonsymmetric and/or indefinite problem, and nonsymmetric preconditioner<br>
				1: SPD problem, and SPD (factored) preconditioner<br>
				2: nonsymmetric, definite problem, and SPD (factored) preconditioner<br>
				Default=0</td>
			</tr>
			<tr>
				<td>`dparms[1]`</td>
				<td>Filters parameters. The filter parameter is used to drop small nonzeros in the preconditioner, to reduce the cost of applying the preconditioner: Default=0.1</td>
			</tr>
			<tr>
				<td>`dparms[2]`</td>
				<td>Threshold parameter: Default=0.1</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab17">Table 17</a>: Definitions of other entries of `:::freefem iparms` and `:::freefem dparms` if preconditionner is Schwartz</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td>`iparms[7]`</td>
				<td>Defines which variant of the Schwartz method isused:<br>
				0: hybrid multiplicative Schwartz method (no overlap across processor boundaries)<br>
				1: hybrid additive Schwartz method (no overlap across processor boundaries)<br>
				2: additive Schwartz method<br>
				3: hybrid multiplicative Schwartz method (with overlap across processor boundaries)<br>
				Default=1</td>
			</tr>
			<tr>
				<td>`iparms[8]`</td>
				<td>Overlap for the Schwartz method: Default=1</td>
			</tr>
			<tr>
				<td>`iparms[9]`</td>
				<td>Type of domain used for the Schwartz method<br>
				0: each point is a domain<br>
				1: each node is a domain (only of interest in "systems" AMG)<br>
				2: each domain is generated by agglomeration (default)</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="3"><a name="Tab18">Table 18</a>: Convergence and time for solving linear system</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td align='center'>n = $4\times10^6$</td>
				<td align='center'>nnz = $13\times10^6$</td>
				<td align='center'>$Te = 571,29$</td>
			</tr>
			<tr>
				<td rowspan="2" align='center'>np</td>
				<td colspan="2" align='center'>AMG</td>
			</tr>
			<tr>
				<td align='center'>`nit`</td>
				<td align='center'>`time`</td>
			</tr>
			<tr>
				<td align='center'>8</td>
				<td align='center'>6</td>
				<td align='center'>1491.83</td>
			</tr>
			<tr>
				<td align='center'>16</td>
				<td align='center'>5</td>
				<td align='center'>708.49</td>
			</tr>
			<tr>
				<td align='center'>32</td>
				<td align='center'>4</td>
				<td align='center'>296.22</td>
			</tr>
			<tr>
				<td align='center'>64</td>
				<td align='center'>4</td>
				<td align='center'>145.64</td>
			</tr>
		</tbody>
	</table>

#### Conclusion

With the different runs presented here, we wanted to illustrate the gain in time when we increase the number of processors used for the simulations. We saw that in every case the time for the construction of the finite element matrix is constant. This is normal because until now this phase is sequential in FreeFem++. In contrast, phases for solving the linear system are parallel. We saw on several examples presented here that when we increase the number of processors, in general we decrease the time used for solving the linear systems. But this is not true in every case. In several case, when we increase the number of processors the time to convergence also increases. There are two main reasons for this. First, the increase of processors can lead to the increase of volume of exchanged data across processors consequently increasing the time for solving the linear systems.

Furthermore, in decomposition domain type preconditioners, the number of processors generally corresponds to the number of sub domains. In subdomain methods, generally when we increase the number of subdomains we decrease convergence quality of the preconditioner. This can increase the time used for solving linear equations.

To end this, we should note that good use of the preconditioners interfaced in __`FreeFem++`__ is empiric, because it is difficult to know what is a good preconditioner for some type of problems. Although, the efficiency of preconditioners sometimes depends on how its parameters are set. For this reason we advise the user to pay attention to the meaning of the parameters in the user guide of the iterative solvers interfaced in __`FreeFem++`__.

### Domain decomposition

In the previous section, we saw that the phases to construct a matrix are sequential. One strategy to construct the matrix in parallel is to divide geometrically the domain into subdomains. In every subdomain we construct a local submatrix and after that we assemble every submatrix to form the global matrix.

We can use this technique to solve PDE directly in domain $\Omega$. In this case, in every subdomains you have to define artificial boundary conditions to form consistent equations in every subdomains. After this, you solve equation in every subdomains and define a strategy to obtain the global solution.

In terms of parallel programming for __`FreeFem++`__, with MPI, this means that the user must be able to divide processors avaible for computation into subgroups of processors and also must be able to realize different type of communications in __`FreeFem++`__ script. Here is a wrapper of some MPI functions.

#### Communicators and groups

__Groups__

`:::freefem mpiGroup grpe(mpiGroup gp, KN_<long>)`: Create MPI_Group from existing group `:::freefem gp` by
given vector.

__Communicators__

Communicators is an abstract MPI object which allows MPI user to communicate across group of processors. Communicators can be Intra-communicators(involves a single group) or Inter-communicators (involves two groups). When we not specify type of communicator it will be Intra-communicators

__mpiComm cc(mpiComm comm, mpiGroup gp):__ Creates a new communicator.

`:::freefem comm` communicator(handle), `:::freefem gp` group which is a subset of the group of `:::freefem comm` (handle). Return new communicator

__mpiComm cc(mpiGroup gp)__: Same as previous constructor but default `:::freefem comm` here is `MPI_COMM_WORLD`.

__mpiComm cc(mpiComm comm, int color, int key):__ Creates new communicators based on `:::freefem colors` and `:::freefem key`. This constructor is based on MPI\_Comm\_split routine of MPI.

__mpiComm cc(MPIrank p, int key):__ Same constructor than the last one.

Here `:::freefem colors` and `:::freefem comm` is defined in `:::freefem MPIrank`. This constructor is based on `MPI_Comm_split` routine of MPI.

!!!example "Split communicator"
	```freefem
	mpiComm comm(mpiCommWorld, 0, 0);
	int color = mpiRank(comm)%2;
	mpiComm ccc(processor(color, comm), 0);
	mpiComm qpp(comm, 0, 0);
	mpiComm cp(ccc, color, 0);
	```

__mpiComm cc(mpiComm comm, int high):__ Creates an intracommunicator from an intercommunicator. `:::freefem comm` intercommunicator, `:::freefem high`.

Used to order the groups within `:::freefem comm` (logical) when creating the new communicator. This constructor is based on `MPI_Intercomm_merge` routine of MPI.

__mpiComm cc(MPIrank p1, MPIrank p2, int tag):__ This constructor creates an intercommuncator from two intracommunicators. `:::freefem p1` defined local (intra)communicator and rank in `local_comm` of leader (often 0) while `:::freefem p2` defined remote communicator and rank in `peer_comm` of remote leader (often 0). `:::freefem tag` Message tag to use in constructing intercommunicator. This constructor is based on `MPI_Intercomm_create`.

!!!example "Merge"
	```freefem
	mpiComm comm, cc;
	int color = mpiRank(comm)%2;
	int rk = mpiRank(comm);
	int size = mpiSize(comm);
	cout << "Color values: " << color << endl;
	mpiComm ccc(processor((rk<size/2), comm), rk);
	mpiComm cp(cc, color, 0);
	int rleader;
	if (rk == 0){ rleader = size/2; }
	else if (rk == size/2){ rleader = 0; }
	else{ rleader = 3; }
	mpiComm qqp(processor(0, ccc), processor(rleader, comm), 12345);
	int aaa = mpiSize(ccc);
	cout << "Number of processor: " << aaa << endl;
	```

#### Process

In __`FreeFem++`__ we wrap MPI process by function call `:::freefem processor` which create internal __`FreeFem++`__ object call `:::freefem MPIrank`. This mean that do not use `:::freefem MPIrank` in __`FreeFem++`__ script.

`:::freefem processor(int rk)`: Keep process rank inside object `:::freefem MPIrank`. Rank is inside `MPI_COMM_WORLD`.

`:::freefem processor(int rk, mpiComm cc)` and `:::freefem processor(mpiComm cc, int rk)` process rank inside communicator cc.

`:::freefem processor(int rk, mpiComm cc)` and `:::freefem processor(mpiComm cc, int rk)` process rank inside communicator cc.

`:::freefem processorblock(int rk)`: This function is exactlly the same than `:::freefem processor(int rk)` but is use in case of blocking communication.

`:::freefem processorblock(int rk, mpiComm cc)`: This function is exactly the same as `:::freefem processor(int rk, mpiComm cc)` but uses a synchronization point.

#### Points to Points communicators

In __`FreeFem++`__ you can call MPI points to points communications functions.

`:::freefem Send(processor(int rk, mpiComm cc), Data D)` : Blocking send of `:::freefem Data D` to processor of `:::freefem rank rk` inside communicator `:::freefem cc`. Note that `:::freefem Data D` can be: `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem real[int]`, `:::freefem complex[int]`, `:::freefem Mesh`, `:::freefem Mesh3`, `:::freefem Matrix`.

`:::freefem Recv(processor(int rk, mpiComm cc), Data D)`: Receive `:::freefem Data D` from process of rank `:::freefem rk` in communicator `:::freefem cc`. Note that `:::freefem Data D` can be: `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem real[int]`, `:::freefem complex[int]`, `:::freefem Mesh`, `:::freefem Mesh3`, `:::freefem Matrix` and should be the same type than corresponding send.

`:::freefem Isend(processor(int rk, mpiComm cc), Data D)` : Non blocking send of `:::freefem Data D` to processor of `:::freefem rank rk` inside communicator `:::freefem cc`. Note that `:::freefem Data D` can be: `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem real[int]`, `:::freefem complex[int]`, `:::freefem mesh`, `:::freefem mesh3`, `:::freefem matrix`.

`:::freefem Recv(processor(int rk, mpiComm cc), Data D)`: Receive corresponding to send.

#### Global operations

In __`FreeFem++`__ you can call MPI global communication functions.

`:::freefem broadcast(processor(int rk, mpiComm cc), Data D)`: Process `:::freefem rk` Broadcast `:::freefem Data D` to all process inside `:::freefem communicator cc`. Note that `:::freefem Data D` can be: `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem real[int]`, `:::freefem complex[int]`, `:::freefem Mesh`, `:::freefem Mesh3`, `:::freefem Matrix`.

`:::freefem broadcast(processor(int rk), Data D)`: Process `:::freefem rk` Broadcast `:::freefem Data D` to all process inside `MPI_COMM_WORLD`. Note that `:::freefem Data D` can be: `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem int[int]`, `:::freefem real[int]`, `:::freefem complex[int]`, `:::freefem Mesh`, `:::freefem Mesh3`, `:::freefem Matrix`.

`:::freefem mpiAlltoall(Data a, Data b)`: Sends `:::freefem data a` from all to all processes. Receive buffer is `:::freefem Data b`. This is done inside communicator `MPI_COMM_WORLD`.

`:::freefem mpiAlltoall(Data a, Data b, mpiComm cc)`: Sends `:::freefem data a` from all to all processes. Receive buffer is `:::freefem Data b`. This is done inside communicator `cc`.

`:::freefem mpiGather(Data a, Data b, processor(mpiComm, int rk)`: Gathers together values `:::freefem Data a` from a group of processes. Process of rank `:::freefem rk` get data on communicator `:::freefem rk`. This function is like `MPI_Gather`.

`:::freefem mpiAllgather(Data a, Data b)`: Gathers `:::freefem Data a` from all processes and distribute it to all in `:::freefem Data b`. This is done inside communicator `MPI_COMM_WORLD`. This function is like `MPI_Allgather`.

`:::freefem mpiAllgather(Data a, Data b, mpiComm cc)`: Gathers `:::freefem Data a` from all processes and distribute it to all in `:::freefem Data b`. This is done inside `:::freefem communicator cc`. This function is like `MPI_Allgather`.

`:::freefem mpiScatter(Data a,Data b,processor(int rk, mpiComm cc))`: Sends `:::freefem Data a` from one process whith rank `:::freefem rk` to all other processes in group represented by communicator `:::freefem mpiComm cc`.

`:::freefem mpiReduce(Data a, Data b, processor(int rk, mpiComm cc), MPI_Op op)` Reduces values `:::freefem Data a` on all processes to a single value `:::freefem Data b` on process of rank `:::freefem rk` and communicator `:::freefem cc`.

Operation use in reduce is: `:::freefem MPI_Op op` which can be: `:::freefem mpiMAX`, `:::freefem mpiMIN`, `:::freefem mpiSUM`, `:::freefem mpiPROD`, `:::freefem mpiLAND`, `:::freefem mpiLOR`, `:::freefem mpiLXOR`, `:::freefem mpiBAND`, `:::freefem mpiBXOR`, `:::freefem mpiMAXLOC`, `:::freefem mpiMINLOC`.

Note that, for all global operations, only `:::freefem int[int]` and `:::freefem real[int]` are data type take in account in __`FreeFem++`__.

### HPDDM solvers

Real valued problems (diffusion, heat, elasticity and Stokes) and complex valued problems (Maxwell and Helmholtz) are given in both 2D and 3D. We detail here the 3D elasticity problem and the 3D time-dependent heat problem.

!!!example "Elasticity 3D"
	A three dimensional elasticity problem is defined. The solver is a domain decomposition method. Domain decomposition methods are a natural framework for parallel computers. The scripts run on multicores computers (from 2 to tens of thousands of cores). Recall that like in any MPI code the number of MPI processes, `:::freefem mpisize`, is given in the command line via the option `:::freefem -np`. We focus on the script `:::freefem Elasticity3D.edp` but the other scripts have the same structure. The command line to run the example on four processes with `:::freefem ffglut` visualization is: `:::freefem ff-mpirun -np 4 Elasticity3D.edp -glut ffglut`

	```freefem
	load "hpddm" //load HPDDM plugin
	macro partitioner()metis//metis, scotch, or parmetis
	macro dimension()3//2D or 3D
	macro vectorialfe()P1//
	include "macro_ddm.idp" //additional DDM functions

	// Macro
	macro def(i)[i, i#B, i#C] //vector field definition
	macro init(i)[i, i, i] //vector field initialization

	real Sqrt = sqrt(2.0);
	macro epsilon(u) [dx(u), dy(u#B), dz(u#C),
		(dz(u#B) + dy(u#C)) / Sqrt,
		(dz(u) + dx(u#C)) / Sqrt,
		(dy(u) + dx(u#B)) / Sqrt] //
	macro div(u) (dx(u) + dy(u#B) + dz(u#C)) //

	// Parameters
	real f = -9000.0;
	real strain = 100.0;
	real Young = 2.0e11; // steel
	real poisson = 0.35;

	func Pk = [vectorialfe, vectorialfe, vectorialfe];

	string deflation = getARGV("-deflation", "geneo"); //coarse space construction
	int overlap = getARGV("-overlap", 1); //geometric overlap between subdomains
	int fakeInterface = getARGV("-interface", 10); //interface between subdomains
	int s = getARGV("-split", 1); //refinement factor
	int p = getARGV("-hpddm_master_p", 1);

	mpiComm comm;
	bool excluded = splitComm(mpiCommWorld, p, comm, topology = getARGV("-hpddm_master_topology", 0), exclude = (usedARGV("-hpddm_master_exclude") != -1));

	// Display
	if (verbosity > 0 && mpirank == 0){
		cout << " --- " << mpirank << "/" << mpisize;
		cout << " - Elasticity3D.edp - input parameters: refinement factor = " << s << " - overlap = " << overlap << endl;
	}

	// Mesh
	int[int] LL = [2, 3, 2, 1, 2, 2];
	meshN ThBorder, Th = cube(1, 1, 1, [x, y, z]);
	fespace Wh(Th, Pk); //local finite element space

	int[int] arrayIntersection; //ranks of neighboring subdomains
	int[int][int] restrictionIntersection(0); //local-to-neighbors renumbering
	real[int] D; //partition of unity
	{
		meshN ThGlobal = cube(10*getARGV("-global", 5), getARGV("-global", 5), getARGV("-global", 5), [10*x, y, z], label=LL); //global mesh
		build(Th, ThBorder, ThGlobal, fakeInterface, s, overlap, D, arrayIntersection, restrictionIntersection, Wh, Pk, comm, excluded, 3)
	}

	// Problem
	real tmp = 1.0 + poisson;
	real mu = Young / (2.0 * tmp);
	real lambda = Young * poisson / (tmp * (1.0 - 2.0 * poisson));
	real[int] rhs; //local right-hand side
	matrix<real> Mat; //local operator
	{ //local weak form
		meshN ThAugmented = Th + ThBorder;
		varf vPb (def(u), def(v))
			= intN(ThAugmented)(
				  lambda * div(u) * div(v)
				+ 2.0 * mu * (epsilon(u)' * epsilon(v))
			)
			+ intN(ThAugmented)(
				  f * vC
			)
			+ on(1, u=0.0, uB=0.0, uC=0.0)
			;

		fespace WhAugmented(ThAugmented, Pk);
		Mat = vPb(WhAugmented, WhAugmented, tgv=-1);
		real[int] rhsFull = vPb(0, WhAugmented, tgv=-1);
		matrix R = interpolate(Wh, WhAugmented);
		renumbering(Mat, R, rhsFull, rhs);
	}
	ThBorder = cube(1, 1, 1, [x, y, z]);

	dschwarz A(Mat, arrayIntersection, restrictionIntersection, scaling = D);

	set(A, sparams = "-hpddm_schwarz_method ras -hpddm_schwarz_coarse_correction balanced -hpddm_variant right -hpddm_verbosity 1 -hpddm_geneo_nu 10");

	matrix<real> Opt; //local operator with optimized boundary conditions
	dpair ret;
	{
		int solver = getOption("schwarz_method");
		if (solver == 1 || solver == 2 || solver == 4){ //optimized Schwarz methods
			fespace Ph(Th, P0);
			real kZero = getARGV("-kZero", 10.0);
			Ph transmission = 2 * kZero * mu * (2 * mu + lambda) / (lambda + 3 * mu);
			varf vOptimized (def(u), def(v))
				= intN(Th)(
					  lambda * div(u) * div(v)
					+ 2.0 * mu * (epsilon(u)' * epsilon(v))
				)
				+ intN1(Th, fakeInterface)(
					  transmission * (def(u)' * def(v))
				)
				+ on(1, u=0.0, uB=0.0, uC=0.0)
				;
			Opt = vOptimized(Wh, Wh, tgv=-1);
		}
		if (mpisize > 1 && isSetOption("schwarz_coarse_correction")){ //two-level Schwarz methods
			if(excluded)
				attachCoarseOperator(mpiCommWorld, A);
			else {
				varf vPbNoPen (def(u), def(v))
					= intN(Th)(
						  lambda * div(u) * div(v)
						+ 2.0 * mu * (epsilon(u)' * epsilon(v))
					)
					+ on(1, u=0.0, uB=0.0, uC=0.0)
					;
				matrix<real> noPen = vPbNoPen(Wh, Wh, solver=CG);
				if(deflation == "geneo") //standard GenEO, no need for RHS -> deduced from LHS (Neumann matrix)
					attachCoarseOperator(mpiCommWorld, A, A=noPen, ret=ret);
				else if(deflation == "dtn"){
					varf vMass (def(u), def(v)) = intN1(Th, fakeInterface)(u * v);
					matrix<real> massMatrix = vMass(Wh, Wh, solver=CG);
					attachCoarseOperator(mpiCommWorld, A, A=noPen, B=massMatrix, pattern=Opt, ret=ret);
				}
				else if(deflation == "geneo-2") //GenEO-2 for optimized Schwarz methods, need for RHS (LHS is still Neumann matrix)
					attachCoarseOperator(mpiCommWorld, A, A=noPen, B=Opt, pattern=Opt, ret=ret);
			}
		}
	}


	// Solve
	Wh<real> def(u); //local solution

	if(Opt.n > 0) //optimized Schwarz methods
		DDM(A, u[], rhs, excluded=excluded, ret=ret, O=Opt);
	else
		u[] = A^-1 * rhs;

	// Error
	real[int] err(u[].n);
	err = A * u[]; //global matrix-vector product
	err -= rhs;

	// Plot
	plotMPI(Th, u[], "Global solution", Pk, def, real, 3, 1)
	plotMPI(Th, err, "Global residual", Pk, def, real, 3, 1)
	real alpha = 2000.0;
	meshN ThMoved = movemesh3(Th, transfo = [x + alpha*u, y + alpha*uB, z + alpha*uC]);
	u[] = mpirank;
	plotMPI(ThMoved, u[], "Global moved solution", Pk, def, real, 3, 1)
	```

	The macro `:::freefem build` is of particular interest since it handles the data distribution among the `:::freefem mpisize` MPI processes with the following steps:

	* The initial mesh `:::freefem ThGlobal` is partitioned by process 0 into `:::freefem mpisize` submeshes

	* The partition is broadcasted to every process $i$ for 0 < $i$ < `:::freefem mpisize`. From then on, all tasks are parallel.

	* Each process creates the local submesh `:::freefem Th` (if the refinement factor `:::freefem s` defined via the option `:::freefem -split` is larger than 1, each local edge is splitted into $s$ subedges, resulting in each element being split into $s^2$ element in 2D and $s^3$ elements in 3D) so that the collection of these submeshes is an overlapping domain decomposition of a refined mesh. The number of extra layers added to the initial partition is monitored by the option `:::freefem -overlap`.

	* Connectivity structures are created
		* `:::freefem D` is the diagonal of the local partition of unity (see [Distributed vectors in HPDDM](#distributed-vectors-in-hpddm for more details)
		* `:::freefem arrayIntersection` is the list of neighbors of the current subdomain
		* For `:::freefem j` in `:::freefem arrayIntersection`, `:::freefem restrictionIntersection[j]` is the list of the degrees of freedom that belong to the intersection of the current subdomain with its neighbor `:::freefem j`.

	Then, the variational formulation `:::freefem vPb` of a three dimensional elasticity problem is used to assemble a local matrix `:::freefem Mat`. This matrix along with `:::freefem D`, `:::freefem arrayIntersection` and `:::freefem restrictionIntersection` are arguments for the constructor of the distributed matrix `:::freefem A`. This is enough to solve the problem with a one-level additive Schwarz method which can be either ASM or RAS.

	For some problems it is interesting to use optimized interface conditions. When there are many subdomains, it is usually profitable to add a second level to the solver. Options are set in the sequel of the script:

	```freefem
	set(A, sparams="-hpddm_schwarz_method ras -hpddm_schwarz_coarse_correction balanced -hpddm_variant right -hpddm_verbosity 1 -hpddm_geneo_nu 10");
	```

	In the above line, the first option selects the one-level preconditioner `:::freefem ras` (possible choices are `:::freefem ras`, `:::freefem oras`, `:::freefem soras`, `:::freefem asm`, `:::freefem osm` or `:::freefem none`), the second option selects the correction formula for the second level here `:::freefem balanced` (possible options are `:::freefem deflated`, `:::freefem additive` or `:::freefem balanced`), the third option selects right preconditioning, the fourth one is verbosity level of HPDDM (different from the one of __`FreeFem++`__), the fifth one prints all possible options of HPPDM and the last one specifies the number of coarse degrees of freedom per subdomain of the GENEO coarse space. All other options of [HPDDM library](cheatsheet of the HPDDM) can be selected via the __`FreeFem++`__ function `:::freefem set`.

	In the last part of the script, the global linear system is solved by the domain decomposition method defined above.

	```freefem
	// Solve
	Wh<real> def(u); //local solution

	if(Opt.n > 0) //optimized Schwarz methods
		DDM(A, u[], rhs, excluded=excluded, ret=ret, O=Opt);
	else
		u[] = A^-1 * rhs;
	```

### Time dependent problem

!!!example "Heat 3D"
	A three dimensional heat problem

	\[
	\frac{\partial u}{\partial t} - \Delta u = 1,\ \ \ u(0,\cdot) := 0 \text{ in }\Omega\,.
	\]

	is discretized by an implicit Euler scheme. At each time step $n$, we shall seek $u^n(x,y,z)$ satisfying for all $w\in H^1(\Omega)$:

	\[
	\int_\Omega \frac{u^n-u^{n-1}}{\delta t}\,w + \nabla u^n \nabla w = \int_\Omega w ,\ \ \ u^0 := 0 \text{ in }\Omega\,.
	\]

	so that at each time step a linear system

	\[
	(M+dt*K) u^n[] = M*u^{n-1}[] + \delta t*F
	\]

	is solved by a domain decomposition method where $M$ is the mass matrix and $K$ is the rigidity matrix. In order to save computational efforts, the domain decomposition method preconditioner is built only once and then reused for all subsequent solves with matrix $A:=M+dt*K$. The distributed matrix vector product with matrix $M$ is made through the call to the function `:::freefem dmv` using the partition of unity associated to matrix $A$.

	```freefem
	load "hpddm" //load HPDDM plugin
	macro partitioner()metis//metis, scotch, or parmetis
	macro dimension()3//2D or 3D
	include "macro_ddm.idp" //additional DDM functions

	// Macro
	macro def(i)i //scalar field definition
	macro init(i)i //scalar field initialization
	macro grad(u) [dx(u), dy(u), dz(u)] //three-dimensional gradient

	// Parameters
	func Pk = P2; //finite element space

	string deflation = getARGV("-deflation", "geneo"); //coarse space construction
	int overlap = getARGV("-overlap", 1); //geometric overlap between subdomains
	int fakeInterface = getARGV("-interface", 10); //interface between subdomains
	int s = getARGV("-split", 1); //refinement factor
	real dt = getARGV("-dt", 0.01); //time step
	int iMax = getARGV("-iMax", 10); //number of iterations

	mpiComm comm;
	int p = getARGV("-hpddm_master_p", 1);
	bool excluded = splitComm(mpiCommWorld, p, comm, topology = getARGV("-hpddm_master_topology", 0), exclude = (usedARGV("-hpddm_master_exclude") != -1));

	// Display
	if (verbosity > 0 && mpirank == 0){
		cout << " --- " << mpirank << "/" << mpisize;
		cout << " - Heat3D.edp - input parameters: refinement factor = " << s << " - overlap = " << overlap << endl;
	}

	// Mesh
	int[int] LL = [1, 2, 1, 1, 1, 1];
	meshN ThBorder, Th = cube(1, 1, 1, [x, y, z]);
	fespace Wh(Th, Pk); //local finite element space
	int[int] arrayIntersection; //ranks of neighboring subdomains
	int[int][int] restrictionIntersection(0); //local-to-neighbors renumbering
	real[int] D; //partition of unity
	{
		meshN ThGlobal = cube(getARGV("-global", 10), getARGV("-global", 10), getARGV("-global", 10), [x, y, z], label=LL); //global mesh
		build(Th, ThBorder, ThGlobal, fakeInterface, s, overlap, D, arrayIntersection, restrictionIntersection, Wh, Pk, comm, excluded)
	}

	// Problem
	real[int] rhs; // local right-hand side
	matrix<real> Mat; //local operator
	matrix<real> M; //local mass matrix
	{ //local weak form
		meshN ThAugmented = Th + ThBorder;
		varf vPb (u, v)
		 	= intN(ThAugmented)(
				  u * v
				+ dt * (grad(u)' * grad(v))
			)
			+ intN(ThAugmented)(
				  dt * v
			)
			+ on(1, u=0.0)
			;
		fespace WhAugmented(ThAugmented, Pk);
		Mat = vPb(WhAugmented, WhAugmented, tgv=-1);
		real[int] rhsFull = vPb(0, WhAugmented, tgv=-1);
		matrix R = interpolate(Wh, WhAugmented);
		varf vPbM (u, v) = intN(ThAugmented)(u * v);
		M = vPbM(WhAugmented, WhAugmented);
		renumbering(M, R, rhsFull, rhs);
		renumbering(Mat, R, rhsFull, rhs);
	}
	ThBorder = cube(1, 1, 1, [x, y, z]);

	dschwarz A(Mat, arrayIntersection, restrictionIntersection, scaling=D);

	matrix<real> Opt; //local operator with optimized boundary conditions
	dpair ret;
	{
		int solver = getOption("schwarz_method");
		if (solver == 1 || solver == 2 || solver == 4){ //optimized Schwarz methods
			fespace Ph(Th, P0);
			real kZero = getARGV("-kZero", 10.0);
			Ph transmission = kZero;
			varf vOptimized (u, v)
				= intN(Th)(
					  u * v
					+ dt * (grad(u)' * grad(v))
				)
				+ intN1(Th, fakeInterface)(
					  transmission * (u * v)
				)
				+ on(1, u=0.0)
				;
			Opt = vOptimized(Wh, Wh, tgv=-1);
		}
		if (mpisize > 1 && isSetOption("schwarz_coarse_correction")){ //two-level Schwarz methods
			if(excluded)
				attachCoarseOperator(mpiCommWorld, A);
			else {
				varf vPbNoPen (u, v)
					= intN(Th)(
						  u * v
						+ dt * (grad(u)' * grad(v))
					)
					+ on(1, u=0.0)
					;
				matrix<real> noPen = vPbNoPen(Wh, Wh, solver=CG);
				if(deflation == "geneo") //standard GenEO, no need for RHS -> deduced from LHS (Neumann matrix)
					attachCoarseOperator(mpiCommWorld, A, A=noPen, ret = ret);
				else if(deflation == "dtn") {
					varf vMass (def(u), def(v)) = intN1(Th, fakeInterface)(u * v);
					matrix<real> massMatrix = vMass(Wh, Wh, solver=CG);
					attachCoarseOperator(mpiCommWorld, A, A=noPen, B=massMatrix, pattern=Opt, ret=ret);
				}
				else if(deflation == "geneo-2") //GenEO-2 for optimized Schwarz methods, need for RHS (LHS is still Neumann matrix)
					attachCoarseOperator(mpiCommWorld, A, A=noPen, B=Opt, pattern=Opt, ret=ret);
			}
		}
	}

	// Solve
	set(A, sparams="-hpddm_reuse_preconditioner=1");
	Wh<real> def(u) = init(0.0); //local solution
	for (int i = 0; i < iMax; ++i){
		real[int] newRhs(rhs.n);
		dmv(A, M, u[], newRhs); //newRhs = M * u[]
		newRhs += rhs;

		if (Opt.n > 0) //optimized Schwarz methods
			DDM(A, u[], newRhs, excluded=excluded, ret=ret, O=Opt);
		else
			u[] = A^-1 * newRhs;

		plotMPI(Th, u[], "Global solution", Pk, def, real, 3, 0)
	}
	```

### Distributed vectors in HPDDM

We give here some hints on the way vectors are distributed among $np$ processes when using __`FreeFem++`__ interfaced with HPDDM. The set of degrees of freedom ${\mathcal N}$ is decomposed into $np$ overlapping sets $({\mathcal N}_i)_{1\le i\le np}$.

 <!--- __ --->

A MPI-process is in charge of each subset. Let $n:=\#{\mathcal N}$ be the number of degrees of freedom of the global finite element space. Let $R_i$ denote the restriction operator from $\R^n$ onto $\R^{\#{\mathcal N}_i}$. We have also defined local diagonal matrices $D_i\in \R^{\#{\mathcal N}_i}\times \R^{\#{\mathcal N}_i}$ so that we have a partition of unity at the algebraic level:

\begin{equation}
	\label{eq:hpddm:14}
  {\mathbf U} = \sum_{i=1}^{np} R_i^T\,D_i\,R_i\,{\mathbf U}\ \ \ \ \forall\ {\mathbf U}\in\R^n\,.
\end{equation}

A global vector ${\mathbf U}\in\R^n$ is actually not stored. Rather, it is stored in a distributed way. Each process $i$, $1\le i\le N$, stores the local vector ${\mathbf U}_i:=R_i {\mathbf U}\in \R^{\#{\mathcal N}_i}$.

It is important to ensure that the result of all linear algebra operators applied to this representation are coherent.

As an example, consider the scalar product of two distributed vectors ${\mathbf U}, {\mathbf V} \in \mathbb{R}^{n}$. Using the partition of unity \eqref{eq:hpddm:14}, we have:

\begin{align*}({\mathbf U}, {\mathbf V}) = \left({\mathbf U}, \sum_{i=1}^{np} R_i^T D_i R_i {\mathbf V}\right) &= \sum_{i=1}^{np} (R_i {\mathbf U}, D_i R_i {\mathbf V})\\
&=\sum_{i=1}^{np} \left({\mathbf U}_i, D_i {\mathbf V}_i\right)\,.
\end{align*}

Thus, the formula for the scalar product is:

\begin{equation*}
({\mathbf U}, {\mathbf V}) = \sum_{i = 1}^{np} (R_i {\mathbf U}, D_i R_i {\mathbf V})\,.
\end{equation*}

Local scalar products are performed concurrently. Thus, the implementation is parallel except for the sum which corresponds to a `:::freefem MPI_Reduce` call across the $np$ MPI processes. Note also that the implementation relies on the knowledge of a partition of unity so that the FreeFem++ syntax is `:::freefem dscalprod(D, u, v)`.

A `:::freefem axpy` procedure $y \leftarrow \alpha\,x+y$ for $x,y\in \mathbb{R}^{n}$ and $\alpha\in\R$ is easily implemented concurrently for distributed vectors in the form:

\[
y_i \leftarrow \alpha\,x_i+y_i\,, \forall\ 1\le i \le np\,.
\]

The matrix vector product is more involved and details are given in the SIAM book  [An Introduction to Domain Decomposition Methods: algorithms, theory and parallel implementation](https://www.ljll.math.upmc.fr/nataf/OT144DoleanJolivetNataf_full.pdf) and even more details are given in [P. Jolivet's PhD manuscrit](http://jolivet.perso.enseeiht.fr/thesis.pdf).

## References

<a name="KARYPIS1995">[KARYPIS1995]</a> KARYPIS, George et KUMAR, Vipin. METIS--unstructured graph partitioning and sparse matrix ordering system, version 2.0. 1995.

<a name="CAI1989">[CAI1989]</a> CAI, Xiao-Chuan. Some domain decomposition algorithms for nonselfadjoint elliptic and parabolic partial differential equations. 1989.

<a name="SAAD2003">[SAAD2003]</a> SAAD, Yousef. Iterative methods for sparse linear systems. siam, 2003.

<a name="SMITH1996">[SMITH1996]</a> SMITH, B. P. Bj rstad and W. Gropp, Domain Decomposition. 1996.
