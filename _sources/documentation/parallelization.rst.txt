.. role:: freefem(code)
    :language: freefem

.. role:: cpp(code)
    :language: cpp

.. role:: c(code)
    :language: c

Parallelization
===============

A first attempt of parallelization of **FreeFEM** is made here with **MPI**.
An extended interface with MPI has been added to **FreeFEM** version 3.5, (see the `MPI documentation <https://www.mpi-forum.org/docs/>`__ for the functionality of the language).

MPI
---

MPI Keywords
~~~~~~~~~~~~

The following keywords and concepts are used:

-  :freefem:`mpiComm` to defined a *communication world*
-  :freefem:`mpiGroup` to defined a group of *processors* in the communication world
-  :freefem:`mpiRequest` to defined a request to wait for the end of the communication

MPI Constants
~~~~~~~~~~~~~

-  :freefem:`mpisize` The total number of *processes*,
-  :freefem:`mpirank` the id-number of my current process in ``{0, ..., mpisize-1}``,
-  :freefem:`mpiUndefined` The :cpp:`MPI_Undefined` constant,
-  :freefem:`mpiAnySource` The :cpp:`MPI_ANY_SOURCE` constant,
-  :freefem:`mpiCommWorld` The :cpp:`MPI_COMM_WORLD` constant,
-  [ … ] and all the keywords of :freefem:`MPI_Op` for the *reduce* operator: :freefem:`mpiMAX`, :freefem:`mpiMIN`, :freefem:`mpiSUM`, :freefem:`mpiPROD`, :freefem:`mpiLAND`, :freefem:`mpiLOR`, :freefem:`mpiLXOR`, :freefem:`mpiBAND`, :freefem:`mpiBXOR`.

MPI Constructor
~~~~~~~~~~~~~~~

.. code-block:: freefem
   :linenos:

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

   mpiComm ncomm1(mpiCommWorld, grp); //set the MPI_Comm form grp

   mpiComm ncomm2(comm, color, key); //MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *ncomm)

   mpiRequest rq; //defined an MPI_Request
   mpiRequest[int] arq(10); //defined an array of 10 MPI_Request

MPI Functions
~~~~~~~~~~~~~

.. code-block:: freefem
   :linenos:

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

where a :freefem:`processor` is just a integer rank, pointer to a :cpp:`MPI_comm` and pointer to a :cpp:`MPI_Request`, and :freefem:`processorblock` with a special :cpp:`MPI_Request`.

MPI Communicator operator
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: freefem
   :linenos:

   int status; //to get the MPI status of send / recv
   real a, b;

   mpiComm comm(mpiCommWorld, 0, 0);
   mpiRequest req;

   //send a,b asynchronously to the process 1
   processor(1) << a << b;
   //receive a,b synchronously from the process 10
   processor(10) >> a >> b;

   //broadcast from processor of comm to other comm processor
   // broadcast(processor(10, comm), a);
   //send synchronously to the process 10 the data a
   status = Send(processor(10, comm), a);
   //receive synchronously from the process 10 the data a
   status = Recv(processor(10, comm), a);

   //send asynchronously to the process 10 the data a without request
   status = Isend(processor(10, comm), a);
   //send asynchronously to the process 10 the data a with request
   status = Isend(processor(10, comm, req), a);
   //receive asynchronously from the process 10 the data a
   status = Irecv(processor(10, req), a);
   //Error asynchronously without request.
   // status = Irecv(processor(10), a);

where the data type of :freefem:`a` can be of type of :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`int[int]`, :freefem:`real[int]`, :freefem:`complex[int]`, :freefem:`int[int,int]`, :freefem:`double[int,int]`, :freefem:`complex[int,int]`, :freefem:`mesh`, :freefem:`mesh3`, :freefem:`mesh[int]`, :freefem:`mesh3[int]`, :freefem:`matrix`, :freefem:`matrix<complex>`

.. code-block:: freefem
   :linenos:

   //send asynchronously to the process 10 the data a with request
   processor(10, req) << a ;
   //receive asynchronously from the process 10 the data a with request
   processor(10, req) >> a ;

If :freefem:`a, b` are arrays or full matrices of :freefem:`int`, :freefem:`real`, or :freefem:`complex`, we can use the following MPI functions:

.. code-block:: freefem
   :linenos:

   mpiAlltoall(a, b, [comm]);
   mpiAllgather(a, b, [comm]);
   mpiGather(a, b, processor(..) );
   mpiScatter(a, b, processor(..));
   mpiReduce(a, b, processor(..), mpiMAX);
   mpiAllReduce(a, b, comm, mpiMAX);

Thank you to Guy-Antoine Atenekeng Kahou for his help to code this interface.

Schwarz example in parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example is a rewritting of example :ref:`Schwarz overlapping <domainDecompositionSchwarzOverlapping>`.

.. code-block:: bash
   :linenos:

   ff-mpirun -np 2 SchwarzParallel.edp
   # OR
   mpirun -np 2 FreeFem++-mpi SchwarzParallel.edp

.. code-block:: freefem
   :linenos:

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
       processor(1-mpirank) << u[];
       processor(1-mpirank) >> U[];

       // Error
       real err0, err1;
       err0 = int1d(Th[mpirank],interior)(square(U - u));
       // send err0 to the other proc, receive in err1
       processor(1-mpirank) << err0;
       processor(1-mpirank) >> err1;
       real err = sqrt(err0 + err1);
       cout << " err = " << err << " - err0 = " << err0 << " - err1 = " << err1 << endl;
       if (err < 1e-3) break;
   }
   if (mpirank == 0)
       plot(u, U);

.. todo:: script freeze in the loop

True parallel Schwarz example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Thank you to F. Nataf*

This is a explanation of the two examples :ref:`MPI-GMRES 2D <exampleMPIGMRES2D>` and :ref:`MPI-GMRES 3D <exampleMPIGMRES3D>`, a Schwarz parallel with a complexity almost independent of the number of process (with a coarse grid preconditioner).

To solve the following Poisson problem on domain :math:`\Omega` with boundary :math:`\Gamma` in :math:`L^2(\Omega)` :

.. math::
   \begin{array}{rcll}
       -\Delta u &=& f & \mbox{ in } \Omega\\
       u &=& g & \mbox{ on } \Gamma
   \end{array}

where :math:`f` and :math:`g` are two given functions of :math:`L^2(\Omega)` and of :math:`H^{\frac12}(\Gamma)`,

Lets introduce :math:`(\pi_i)_{i=1,.., N_p}` a regular partition of the unity of :math:`\Omega`, q-e-d:

.. math::
   \pi_i \in \mathcal{C}^0(\Omega) : \quad \pi_i\ge 0 \mbox{ and } \sum_{i=1}^{N_p} \pi_i =1 .

Denote :math:`\Omega_i` the sub domain which is the support of :math:`\pi_i` function and also denote :math:`\Gamma_i` the boundary of :math:`\Omega_i`.

The parallel Schwarz method is:

Let :math:`\ell=0` the iterator and a initial guest :math:`u^0` respecting the boundary condition (i.e. :math:`u^0_{|\Gamma} = g`).

.. math::
   \begin{array}{rcll}
       \forall i = 1 .., N_p:&\nonumber\\
       \displaystyle -\Delta u_i^\ell &=& f &\mbox{ in } \Omega_i\\
       u_i^\ell &=& u^\ell & \mbox{ on } \Gamma_i \setminus \Gamma\\
       u_i^\ell &=& g & \mbox{ on } \Gamma_i \cap \Gamma
   \end{array}
   :label: eq:lapl

.. math::
   u^{\ell+1} = \sum_{i=1}^{N_p} \pi_i u_i^\ell
   :label: eq:pu1

After discretization with the Lagrange finite element method, with a compatible mesh :math:`{\mathcal{T}_h}_i` of :math:`\Omega_i`, i. e., the exist a global mesh :math:`{\mathcal{T}_h}` such that :math:`{\mathcal{T}_h}_i` is include in :math:`{\mathcal{T}_h}`.

Let us denote:

-  :math:`{V_h}_i` the finite element space corresponding to domain :math:`\Omega_i`,
-  :math:`{\mathcal{N}_h}_i` is the set of the degree of freedom :math:`\sigma_i^k`,
-  :math:`{\mathcal{N}^{\Gamma_i}_{hi}}` is the set of the degree of freedom of :math:`{V_h}_i` on the boundary :math:`\Gamma_i` of :math:`\Omega_i`,
-  :math:`\sigma_i^k({v_h})` is the value the degree of freedom :math:`k`,
-  :math:`{V_{0h}}_i= \{ {v_h} \in {V_h}_i :\forall k \in {\mathcal{N}^{\Gamma_i}_{hi}}, \quad \sigma_i^k({v_h})=0 \}`,
-  The conditional expression :math:`a\;?\;b:c` is defined like in :c`C` of :cpp:`C++` language by

   .. math::
        a?b: c \equiv
        \left\{
        \begin{array}{l}
        \mbox{if } a \mbox{ is true then return b}\\
        \mbox{else return } c\\
        \end{array}
        \right..

.. note:: We never use finite element space associated to the full domain :math:`\Omega` because it is too expensive.

We have to defined to operator to build the previous algorithm:

We denote :math:`{u_h^{\ell}}_{|i}` the restriction of :math:`u_h^\ell` on :math:`{V_h}_i`, so the discrete problem on :math:`\Omega_i` of problem :eq:`eq:lapl` is find :math:`{u_h^{\ell}}_{i}\in {V_h}_i` such that:

.. math::
   \forall {v_h}_i\in V_{0i}:
   \int_{\Omega_i} \nabla {v_h}_i \cdot \nabla {u_h}^{\ell}_{i}
   = \int_{\Omega_i} f {v_h}_i ,\quad \forall k \in {\mathcal{N}^{\Gamma_i}_{hi}}\;:\; \sigma_i^k({u_h}^\ell_i) = (k\in \Gamma) \; ? \; g_i^k : \sigma_i^k({u_h}^{\ell}_{|i})

where :math:`g_i^k` is the value of :math:`g` associated to the degree of freedom :math:`k\in {\mathcal{N}^{\Gamma_i}_{hi}}`.

In **FreeFEM**, it can be written has with :freefem:`U` is the vector corresponding to :math:`{u_h^{\ell}}_{|i}` and the vector :freefem:`U1` is the vector corresponding to :math:`{u_h^{\ell}}_{i}` is the solution of:

.. code-block:: freefem
   :linenos:

   real[int] U1(Ui.n);
   real[int] b = onG .* U;
   b = onG ? b : Bi ;
   U1 = Ai^-1*b;

where :math:`\mathtt{onG}[i] =(i \in \Gamma_i\setminus\Gamma) ? 1 : 0`, and :math:`\mathtt{Bi}` the right of side of the problem, are defined by

.. code-block:: freefem
   :linenos:

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

where the **FreeFEM** label of :math:`\Gamma` is 1 and the label of :math:`\Gamma_i\setminus \Gamma` is :math:`10`.

To build the transfer/update part corresponding to :eq:`eq:pu1` equation on process :math:`i`, let us call :freefem:`njpart` the number the neighborhood of domain of :math:`\Omega_i` (i.e: :math:`\pi_j` is none :math:`0` of :math:`\Omega_i`), we store in an array :freefem:`jpart` of size :freefem:`njpart` all this neighborhood.

Let us introduce two array of matrix, :freefem:`Smj[j]` to defined the vector to send from :math:`i` to :math:`j` a neighborhood process, and the matrix :math:`rMj[j]` to after to reduce owith neighborhood :math:`j` domain.

So the tranfert and update part compute :math:`v_i= \pi_i u_i + \sum_{j\in J_i} \pi_j u_j` and can be write the **FreeFEM** function Update:

.. code-block:: freefem
   :linenos:

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

where the buffer are defined by:

.. code-block:: freefem
   :linenos:

   InitU(njpart, Whij, Thij, aThij, Usend) //defined the send buffer
   InitU(njpart, Whij, Thij, aThij, Vrecv) //defined the revc buffer

with the following macro definition:

.. code-block:: freefem
   :linenos:

   macro InitU(n, Vh, Th, aTh, U) Vh[int] U(n); for (int j = 0; j < n; ++j){Th = aTh[j]; U[j] = 0;}

*First GMRES algorithm:* you can easily accelerate the fixed point algorithm by using a parallel GMRES algorithm after the introduction the following affine :math:`\mathcal{A}_i` operator sub domain :math:`\Omega_i`.

.. code-block:: freefem
   :linenos:

   func real[int] DJ0 (real[int]& U){
       real[int] V(U.n), b = onG .* U;
       b = onG ? b : Bi ;
       V = Ai^-1*b;
       Update(V, U);
       V -= U;
       return V;
   }

Where the parallel :freefem:`MPIGMRES` or :freefem:`MPICG` algorithm is just a simple way to solve in parallel the following :math:`A_i x_i = b_i, i = 1, .., N_p` by just changing the dot product by reduce the local dot product of all process with the following MPI code:

.. code-block:: cpp
    :linenos:

    template<class R> R ReduceSum1(R s, MPI_Comm *comm){
        R r = 0;
        MPI_Allreduce(&s, &r, 1, MPI_TYPE<R>::TYPE(), MPI_SUM, *comm );
        return r;
    }

This is done in :freefem:`MPIGC` dynamics library tool.

*Second GMRES algorithm:* Use scharwz algorithm as a preconditioner of basic GMRES method to solving the parallel problem.

.. code-block:: freefem
   :linenos:

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

*Third GMRES algorithm:* Add a coarse solver to the previous algorithm

First build a coarse grid on processor 0, and the

.. code-block:: freefem
   :linenos:

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

The New preconditionner

.. code-block:: freefem
   :linenos:

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

The code of the 4 algorithms:

.. code-block:: freefem
   :linenos:

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

We have all ingredient to solve in parallel if we have et the partitions of the unity.
To build this partition we do:

The initial step on process :math:`1` to build a coarse mesh, :math:`{\mathcal{T}_h}^*` of the full domain, and build the partition :math:`\pi` function constant equal to :math:`i` on each sub domain :math:`\mathcal{O}_i, i =1 ,.., N_p`, of the grid with the :freefem:`metis` graph partitioner [KARYPIS1995]_ and on each process :math:`i` in :math:`1..,N_p` do

1. Broadcast from process :math:`1`, the mesh :math:`{\mathcal{T}_h}^*` (call :freefem:`Thii` in **FreeFEM** script), and :math:`\pi` function,
2. remark that the characteristic function :math:`\mathrm{1\!\!I}_{\mathcal{O}_i}` of domain :math:`\mathcal{O}_i`, is defined by :math:`(\pi=i)?1:0`,
3. Let us call :math:`\Pi^2_P` (resp. :math:`\Pi^2_V`) the :math:`L^2` on :math:`P_h^*` the space of the constant finite element function per element on :math:`{\mathcal{T}_h}^*` (resp. :math:`V_h^*` the space of the affine continuous finite element per element on :math:`{\mathcal{T}_h}^*`) and build in parallel the :math:`\pi_i` and :math:`\Omega_i`, such that :math:`\mathcal{O}_i\ \subset \Omega_i` where :math:`\mathcal{O}_i= supp ((\Pi^2_V \Pi^2_C)^m \mathrm{1\!\!I}_{O_i})`, and :math:`m` is a the overlaps size on the coarse mesh (generally one), (this is done in function :freefem:`AddLayers(Thii,suppii[],nlayer,phii[]);` We choose a function :math:`\pi^*_i = (\Pi^2_1 \Pi^2_0)^m \mathrm{1\!\!I}_{\mathcal{O}_i}` so the partition of the unity is simply defined by

    .. math::
        \pi_i = \frac{\pi_i^*}{\sum_{j=1}^{N_p} \pi_j^*}

    The set :math:`J_i` of neighborhood of the domain :math:`\Omega_i`, and the local version on :math:`V_{hi}` can be defined the array :freefem:`jpart` and :freefem:`njpart` with:

    .. code-block:: freefem
        :linenos:

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

4. We call :math:`{\mathcal{T}_h}^*_{ij}` the sub mesh part of :math:`{\mathcal{T}_h}_i` where :math:`\pi_j` are none zero.
   And thanks to the function :freefem:`trunc` to build this array,

    .. code-block:: freefem
        :linenos:

        for(int jp = 0; jp < njpart; ++jp)
            aThij[jp] = trunc(Thi, pij[jp] > 1e-10, label=10);

5. At this step we have all on the coarse mesh, so we can build the fine final mesh by splitting all meshes: :freefem:`Thi, Thij[j], Thij[j]` with **FreeFEM** :freefem:`trunc` mesh function which do restriction and slipping.

6. The construction of the send/recv matrices :freefem:`sMj` and `freefem:`rMj`: can done with this code:

    .. code-block:: freefem
        :linenos:

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

To buil a not too bad application, all variables come from parameters value with the following code

.. code-block:: freefem
    :linenos:

    include "getARGV.idp"
    verbosity = getARGV("-vv", 0);
    int vdebug = getARGV("-d", 1);
    int ksplit = getARGV("-k", 10);
    int nloc = getARGV("-n", 25);
    string sff = getARGV("-p, ", "");
    int gmres = getARGV("-gmres", 3);
    bool dplot = getARGV("-dp", 0);
    int nC = getARGV("-N", max(nloc/10, 4));

And small include to make graphic in parallel of distribute solution of vector :math:`u` on mesh :math:`T_h` with the following interface:

.. code-block:: freefem
    :linenos:

    include "MPIplot.idp"
    func bool plotMPIall(mesh &Th, real[int] &u, string cm){
        PLOTMPIALL(mesh, Pk, Th, u, {cmm=cm, nbiso=20, fill=1, dim=3, value=1});
        return 1;
    }

.. note:: The :freefem:`cmm=cm, ...` in the macro argument is a way to quote macro argument so the argument is :freefem:`cmm=cm, ...`.

.. _parallelSparseSolvers:

Parallel sparse solvers
-----------------------

Parallel sparse solvers use several processors to solve linear systems of equation. Like sequential, parallel linear solvers can be direct or iterative. In **FreeFEM** both are available.

Using parallel sparse solvers in **FreeFEM**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recall that the :freefem:`solver` parameters are defined in the following commands: :freefem:`solve`, :freefem:`problem`, :freefem:`set` (setting parameter of a matrix) and in the construction of the matrix corresponding to a bilinear form.
In these commands, the parameter :freefem:`solver` must be set to :freefem:`sparsesolver` for parallel sparse solver.
We have added specify parameters to these command lines for parallel sparse solvers.
These are:

-  :freefem:`lparams` : vector of integer parameters (``l`` is for the ``C++`` type :cpp:`long`)
-  :freefem:`dparams` : vector of real parameters
-  :freefem:`sparams` : string parameters
-  :freefem:`datafilename` : name of the file which contains solver parameters

The following four parameters are only for direct solvers and are vectors.
These parameters allow the user to preprocess the matrix (see the section on `sparse direct solver <#sparse-direct-solver>`__ for more information).

-  :freefem:`permr` : row permutation (integer vector)
-  :freefem:`permc` : column permutation or inverse row permutation (integer vector)
-  :freefem:`scaler` : row scaling (real vector)
-  :freefem:`scalec` : column scaling (real vector)

There are two possibilities to control solver parameters.
The first method defines parameters with :freefem:`lparams`, :freefem:`dparams` and :freefem:`sparams` in ``.edp`` file.

The second one reads the solver parameters from a data file. The name of this file is specified by :freefem:`datafilename`.
If :freefem:`lparams`, :freefem:`dparams`, :freefem:`sparams` or :freefem:`datafilename` is not provided by the user, the solver’s default values are used.

To use parallel solver in **FreeFEM**, we need to load the dynamic library corresponding to this solver.
For example to use `MUMPS <http://mumps.enseeiht.fr/>`__ solver as parallel solver in **FreeFEM**, write in the ``.edp`` file :freefem:`load "MUMPS_FreeFem"`.

If the libraries are not loaded, the default sparse solver will be loaded (default sparse solver is :freefem:`UMFPACK`). The :numref:`tabParallelizationSparseSolver` gives this new value for the different libraries.


.. table:: Default sparse solver for real and complex arithmetics when we load a parallel sparse solver library
    :name: tabParallelizationSparseSolver

    +------------------------------+-----------------------------------+
    | Libraries                    | Default sparse solver             |
    +------------------------------+-----------------+-----------------+
    |                              | real            | complex         |
    +==============================+=================+=================+
    | MUMPS_FreeFem                | mumps           | mumps           |
    +------------------------------+-----------------+-----------------+
    | real_SuperLU_DIST_FreeFem    | SuperLU_DIST    | previous solver |
    +------------------------------+-----------------+-----------------+
    | complex_SuperLU_DIST_FreeFem | previous solver | SuperLU_DIST    |
    +------------------------------+-----------------+-----------------+
    | real_pastix_FreeFem          | PaStiX          | previous solver |
    +------------------------------+-----------------+-----------------+
    | complex_pastix_FreeFem       | previous solver | PaStiX          |
    +------------------------------+-----------------+-----------------+
    | hips_FreeFem                 | hips            | previous solver |
    +------------------------------+-----------------+-----------------+
    | hypre_FreeFem                | hypre           | previous solver |
    +------------------------------+-----------------+-----------------+
    | parms_FreeFem                | parms           | previous solver |
    +------------------------------+-----------------+-----------------+

We also add functions (see :numref:`tabParallelizationFunction`) with no parameter to change the default sparse solver in the ``.edp`` file.
To use these functions, we need to load the library corresponding to the solver.
An example of using different parallel sparse solvers for the same problem is given in :ref:`Direct solvers example <exampleDirectSolvers>`.

.. table:: Functions that allow to change the default sparse solver for real and complex arithmetics and the result of these functions
    :name: tabParallelizationFunction

    +-------------------------------+-----------------------------------+
    | Function                      | default sparse solver             |
    +-------------------------------+-----------------+-----------------+
    |                               | real            | complex         |
    +===============================+=================+=================+
    | defaulttoMUMPS()              | mumps           | mumps           |
    +-------------------------------+-----------------+-----------------+
    | realdefaulttoSuperLUdist()    | SuperLU_DIST    | previous solver |
    +-------------------------------+-----------------+-----------------+
    | complexdefaulttoSuperLUdist() | previous solver | SuperLU_DIST    |
    +-------------------------------+-----------------+-----------------+
    | realdefaultopastix()          | pastix          | previous solver |
    +-------------------------------+-----------------+-----------------+
    | complexdefaulttopastix()      | previous solver | pastix          |
    +-------------------------------+-----------------+-----------------+
    | defaulttohips()               | hips            | previous solver |
    +-------------------------------+-----------------+-----------------+
    | defaulttohypre()              | hypre           | previous solver |
    +-------------------------------+-----------------+-----------------+
    | defaulttoparms()              | parms           | previous solver |
    +-------------------------------+-----------------+-----------------+


.. tip:: Test direct solvers

    .. code-block:: freefem
        :linenos:

        load "MUMPS_FreeFem"
        //default solver: real-> MUMPS, complex -> MUMPS
        load "real_SuperLU_DIST_FreeFem"
        //default solver: real-> SuperLU_DIST,
        complex -> MUMPS load "real_pastix_FreeFem"
        //default solver: real-> pastix, complex -> MUMPS

        // Solving with pastix
        {
            matrix A =
                [[1, 2, 2, 1, 1],
                [ 2, 12, 0, 10, 10],
                [ 2, 0, 1, 0, 2],
                [ 1, 10, 0, 22, 0.],
                [ 1, 10, 2, 0., 22]];

            real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
            b = A*xx;
            cout << "b =" << b << endl; cout << "xx =" << xx << endl;

            set(A, solver=sparsesolver, datafilename="ffpastix_iparm_dparm.txt");
            cout << "solve" << endl;
            x = A^-1*b;
            cout << "b =" << b << endl;
            cout << "x =" << endl;
            cout << x << endl;
            di = xx - x;
            if (mpirank == 0){
                cout << "x-xx =" << endl;
                cout << "Linf =" << di.linfty << ", L2 =" << di.l2 << endl;
            }
        }

        // Solving with SuperLU_DIST
        realdefaulttoSuperLUdist();
        //default solver: real-> SuperLU_DIST, complex -> MUMPS
        {
            matrix A =
                [[1, 2, 2, 1, 1],
                [ 2, 12, 0, 10, 10],
                [ 2, 0, 1, 0, 2],
                [ 1, 10, 0, 22, 0.],
                [ 1, 10, 2, 0., 22]];

            real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
            b = A*xx;
            cout << "b =" << b << endl;
            cout << "xx =" << xx << endl;

            set(A, solver=sparsesolver, datafilename="ffsuperlu_dist_fileparam.txt");
            cout << "solve" << endl;
            x = A^-1*b;
            cout << "b =" << b << endl;
            cout << "x =" << endl;
            cout << x << endl;
            di = xx - x;
            if (mpirank == 0){
                cout << "x-xx =" << endl;
                cout << "Linf =" << di.linfty << ", L2 =" << di.l2 << endl;
            }
        }

        // Solving with MUMPS
        defaulttoMUMPS();
        //default solver: real-> MUMPS, complex -> MUMPS
        {
            matrix A =
                [[1, 2, 2, 1, 1],
                [ 2, 12, 0, 10, 10],
                [ 2, 0, 1, 0, 2],
                [ 1, 10, 0, 22, 0.],
                [ 1, 10, 2, 0., 22]];

            real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
            b = A*xx;
            cout << "b =" << b << endl;
            cout << "xx =" << xx << endl;

            set(A, solver=sparsesolver, datafilename="ffmumps_fileparam.txt");
            cout << "solving solution" << endl;
            x = A^-1*b;
            cout << "b =" << b << endl;
            cout << "x =" << endl;
            cout << x << endl;
            di = xx - x;
            if (mpirank == 0){
                cout << "x-xx =" << endl;
                cout << "Linf =" << di.linfty << ", L2" << di.l2 << endl;
            }
        }

Sparse direct solver
~~~~~~~~~~~~~~~~~~~~

In this section, we present the sparse direct solvers interfaced with **FreeFEM**.

MUMPS solver
^^^^^^^^^^^^

MUltifrontal Massively Parallel Solver (`MUMPS <http://mumps.enseeiht.fr/>`__) is an open-source library.

This package solves linear system of the form :math:`A \: x = b` where :math:`A` is a square sparse matrix with a direct method.
The square matrix considered in MUMPS can be either unsymmetric, symmetric positive definite or general symmetric.

The method implemented in MUMPS is a direct method based on a multifrontal approach.
It constructs a direct factorization :math:`A \:= \: L\:U`, :math:`A\: = \: L^t \: D \: L` depending of the symmetry of the matrix :math:`A`.

MUMPS uses the following libraries :
    * `BLAS <http://www.netlib.org/blas/>`__,
    * `BLACS <http://www.netlib.org/blacs/>`__,
    * `ScaLAPACK <http://www.netlib.org/scalapack/>`__.

.. warning:: MUMPS does not solve linear system with a rectangular matrix.

**MUMPS parameters:**

There are four input parameters in `MUMPS <http://mumps.enseeiht.fr/index.php?page=doc>`__.
Two integers :cpp:`SYM` and :cpp:`PAR`, a vector of integer of size 40 :cpp:`INCTL` and a vector of real of size 15 :cpp:`CNTL`.

The first parameter gives the type of the matrix: 0 for unsymmetric matrix, 1 for symmetric positive matrix and 2 for general symmetric.

The second parameter defined if the host processor work during the factorization and solves steps : :cpp:`PAR=1` host processor working and :cpp:`PAR=0` host processor not working.

The parameter :cpp:`INCTL` and :cpp:`CNTL` is the control parameter of MUMPS.
The vectors :cpp:`ICNTL` and :cpp:`CNTL` in MUMPS becomes with index 1 like vector in ``Fortran``.
For more details see the `MUMPS user’s guide <http://mumps.enseeiht.fr/index.php?page=doc>`__.

We describe now some elements of the main parameters of :cpp:`ICNTL` for MUMPS.

-  **Input matrix parameter** The input matrix is controlled by parameters ``ICNTL(5)`` and ``ICNTL(18)``.
    The matrix format (resp. matrix pattern and matrix entries) are controlled by ``INCTL(5)`` (resp. ``INCTL(18)``).

    The different values of ``ICNTL(5)`` are 0 for assembled format and 1 for element format.
    In the current release of **FreeFEM**, we consider that FE matrix or matrix is storage in assembled format.
    Therefore, ``INCTL(5)`` is treated as 0 value.

    The main option for ``ICNTL(18)``: ``INCLTL(18)=0`` centrally on the host processor, ``ICNTL(18)=3`` distributed the input matrix pattern and the entries (recommended option for distributed matrix by developer of MUMPS).
    For other values of ``ICNTL(18)`` see the `MUMPS user’s guide <http://mumps.enseeiht.fr/index.php?page=doc>`__.
    These values can be used also in **FreeFEM**.

    The default option implemented in **FreeFEM** are ``ICNTL(5)=0`` and ``ICNTL(18)=0``.

-  **Preprocessing parameter** The preprocessed matrix :math:`A_{p}` that will be effectively factored is defined by

    .. math::
        A_{p} = P \: D_r \: A \: Q_c \ D_c P^t


    where :math:`P` is the permutation matrix, :math:`Q_c` is the column permutation, :math:`D_r` and :math:`D_c` are diagonal matrix for respectively row and column scaling.

    The ordering strategy to obtain :math:`P` is controlled by parameter ``ICNTL(7)``.
    The permutation of zero free diagonal :math:`Q_c` is controlled by parameter ``ICNTL(6)``.
    The row and column scaling is controlled by parameter ``ICNTL(18)``.
    These option are connected and also strongly related with ``ICNTL(12)`` (see the `MUMPS user’s guide <http://mumps.enseeiht.fr/index.php?page=doc>`__ for more details).

    The parameters :freefem:`permr`, :freefem:`scaler`, and :freefem:`scalec` in **FreeFEM** allow to give permutation matrix(\ :math:`P`), row scaling (:math:`D_r`) and column scaling (:math:`D_c`) of the user respectively.

**Calling MUMPS in FreeFEM**

To call MUMPS in **FreeFEM**, we need to load the dynamic library ``MUMPS_freefem.dylib`` (MacOSX), ``MUMPS_freefem.so`` (Unix) or ``MUMPS_freefem.dll`` (Windows).

This is done in typing :freefem:`load "MUMPS_FreeFem"` in the ``.edp`` file. We give now the two methods to give the option of MUMPS solver in **FreeFEM**.

-  **Solver parameters is defined in .edp file:** In this method, we need to give the parameters :freefem:`lparams` and :freefem:`dparams`.
    These parameters are defined for MUMPS by :

        * :freefem:`lparams[0] = SYM`, :freefem:`lparams[1] = PAR`,
        * :math:`\forall i` = 1,…,40, :freefem:`lparams[i+1] = ICNTL(i)`
        * :math:`\forall i` = 1,…,15, :freefem:`dparams[i-1] = CNTL(i)`

-  **Reading solver parameters on a file:**

    The structure of data file for MUMPS in **FreeFEM** is : first line parameter ``SYM`` and second line parameter ``PAR`` and in the following line the different value of vectors ``ICNTL`` and ``CNTL``.
    An example of this parameter file is given in :freefem:`ffmumpsfileparam.txt`.

    .. code-block:: freefem
        :linenos:

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
        1 /* ICNTL(9) :: 0 solve Ax = b, 1 solve the transposed system A^t x = b : parameter is not considered in the current release of FreeFEM*/
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
        0 /* ICNTL(20) :: right hand side form ( 0 dense form, 1 sparse form) : parameter will be set to 0 for FreeFEM */
        0 /* ICNTL(21) :: 0, 1 kept distributed solution : parameter is not considered in the current release of FreeFEM */
        0 /* ICNTL(22) :: controls the in-core/out-of-core (OOC) facility */
        0 /* ICNTL(23) :: maximum size of the working memory in Megabyte than MUMPS can allocate per working processor */
        0 /* ICNTL(24) :: control the detection of null pivot */
        0 /* ICNTL(25) :: control the computation of a null space basis */
        0 /* ICNTL(26) :: This parameter is only significant with Schur option (ICNTL(19) not zero). : parameter is not considered in the current release of FreeFEM */
        -8 /* ICNTL(27) (Experimental parameter subject to change in next release of MUMPS) :: control the blocking factor for multiple righthand side during the solution phase : parameter is not considered in the current release of FreeFEM */
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

    If no solver parameter is given, we used default option of MUMPS solver.

.. tip:: MUMPS example

    A simple example of calling MUMPS in **FreeFEM** with this two methods is given in the :ref:`Test solver MUMPS example <exampleSolverMUMPS>`.

SuperLU distributed solver
^^^^^^^^^^^^^^^^^^^^^^^^^^

The package `SuperLU_DIST <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__ solves linear systems using LU factorization.
It is a free scientific library

This library provides functions to handle square or rectangular matrix in real and complex arithmetics.
The method implemented in SuperLU_DIST is a supernodal method.
New release of this package includes a parallel symbolic factorization.
This scientific library is written in C and MPI for communications.

**SuperLU_DIST parameters:**

We describe now some parameters of SuperLU_DIST.
The SuperLU_DIST library use a 2D-logical process group.
This process grid is specified by :math:`nprow` (process row) and :math:`npcol` (process column) such that :math:`N_{p} = nprow \: npcol` where :math:`N_{p}` is the number of all process allocated for SuperLU_DIST.

The input matrix parameters is controlled by "matrix=" in :freefem:`sparams` for internal parameter or in the third line of parameters file.
The different value are

-  :freefem:`matrix=assembled` global matrix are available on all process
-  :freefem:`matrix=distributedglobal` The global matrix is distributed among all the process
-  :freefem:`matrix=distributed` The input matrix is distributed (not yet implemented)

The option arguments of SuperLU_DIST are described in the section Users-callable routine of the `SuperLU users’ guide <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ug.pdf>`__.

The parameter ``Fact`` and ``TRANS`` are specified in **FreeFEM** interfaces to SuperLU_DIST during the different steps.
For this reason, the value given by the user for this option is not considered.

The factorization LU is calculated in SuperLU_DIST on the matrix :math:`A_p`.

.. math::
    A_{p} = P_{c} \: P_r \: D_r \: A \: D_{c} \: P_{c}^{t}

where :math:`P_c` and :math:`P_r` is the row and column permutation matrix respectively, :math:`D_r` and :math:`D_c` are diagonal matrix for respectively row and column scaling.

The option argument ``RowPerm`` (resp. ``ColPerm``) control the row (resp. column) permutation matrix.
:math:`D_r` and :math:`D_c` is controlled by the parameter ``DiagScale``.

The parameter :freefem:`permr`, :freefem:`permc`, :freefem:`scaler`, and :freefem:`scalec` in **FreeFEM** is provided to give row permutation, column permutation, row scaling and column scaling of the user respectively.

The other parameters for LU factorization are ``ParSymFact`` and ``ReplaceTinyPivot``.
The parallel symbolic factorization works only on a power of two processes and need the ``ParMetis`` ordering.
The default option argument of SuperLU_DIST are given in the file ``ffsuperlu_dist_fileparam.txt``.

**Calling SuperLU_DIST in FreeFEM**

To call SuperLU_DIST in **FreeFEM**, we need to load the library dynamic correspond to interface.
This done by the following line :freefem:`load "real_superlu _DIST_FreeFem"` (resp. :freefem:`load "complex_superlu_DIST_FreeFem"`) for real (resp. complex) arithmetics in the file ``.edp``.

**Solver parameters is defined in .edp file:**

To call SuperLU_DIST with internal parameter, we used the parameters ``sparams``.
The value of parameters of SuperLU_DIST in ``sparams`` are defined by :

-  ``nprow=1``,
-  ``npcol=1``,
-  ``matrix= distributedgloba``,
-  ``Fact= DOFACT``,
-  ``Equil=NO``,
-  ``ParSymbFact=NO``,
-  ``ColPerm= MMD_AT_PLUS_A``,
-  ``RowPerm= LargeDiag``,
-  ``DiagPivotThresh=1.0``,
-  ``IterRefine=DOUBLE``,
-  ``Trans=NOTRANS``,
-  ``ReplaceTinyPivot=NO``,
-  ``SolveInitialized=NO``,
-  ``PrintStat=NO``,
-  ``DiagScale=NOEQUIL``

This value correspond to the parameter in the file ``ffsuperlu_dist_fileparam.txt``.
If one parameter is not specified by the user, we take the default value of SuperLU_DIST.

**Reading solver parameters on a file:** The structure of data file for SuperLU_DIST in **FreeFEM** is given in the file ``ffsuperlu_dist_fileparam.txt`` (default value of the **FreeFEM** interface).

.. code-block:: freefem
    :linenos:

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

If no solver parameter is given, we used default option of SuperLU_DIST solver.

.. tip:: A simple example of calling SuperLU_DIST in **FreeFEM** with this two methods is given in the :ref:`Solver superLU_DIST example <exampleSolverSuperLUDist>`.

PaStiX solver
^^^^^^^^^^^^^

`PaStiX <http://pastix.gforge.inria.fr/files/README-txt.html>`__ (Parallel Sparse matrix package) is a free scientific library under CECILL-C license.
This package solves sparse linear system with a direct and block ILU(k) iterative methods.
his solver can be applied to a real or complex matrix with a symmetric pattern.

**PaStiX parameters:**

The input :freefem:`matrix` parameter of **FreeFEM** depend on PaStiX interface.
:freefem:`matrix = assembled` for non distributed matrix.
It is the same parameter for SuperLU_DIST.

There are four parameters in PaStiX : ``iparm``, ``dparm``, ``perm`` and ``invp``.
These parameters are respectively the integer parameters (vector of size 64), real parameters (vector of size 64), permutation matrix and inverse permutation matrix respectively.
``iparm`` and ``dparm`` vectors are described in `PaStiX RefCard <https://gforge.inria.fr/docman/?group_id=186&view=listfile&dirid=246>`__.

The parameters :freefem:`permr` and :freefem:`permc` in **FreeFEM** are provided to give permutation matrix and inverse permutation matrix of the user respectively.

**Solver parameters defined in .edp file:**

To call PaStiX in **FreeFEM** in this case, we need to specify the parameters :freefem:`lparams` and :freefem:`dparams`.
These parameters are defined by :

:math:`\forall i` = 0,… ,63, :freefem:`lparams[i] = iparm[i]`.

:math:`\forall i` = 0,… ,63, :freefem:`dparams[i] = dparm[i]`.

**Reading solver parameters on a file:**

The structure of data file for PaStiX parameters in **FreeFEM** is: first line structure parameters of the matrix and in the following line the value of vectors ``iparm`` and ``dparm`` in this order.

.. code-block:: freefem
    :linenos:

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

An example of this file parameter is given in ``ffpastix_iparm_dparm.txt`` with a description of these parameters.
This file is obtained with the example file ``iparm.txt`` and ``dparm.txt`` including in the PaStiX package.

If no solver parameter is given, we use the default option of PaStiX solver.

.. tip:: A simple example of calling PaStiX in **FreeFEM** with this two methods is given in the :ref:`Solver PaStiX example <exampleSolverPastix>`.

In :numref:`tabParallelizationDirectSolver`, we recall the different matrix considering in the different direct solvers.

.. table:: Type of matrix used by the different direct sparse solver
    :name: tabParallelizationDirectSolver

    +---------------+---------------------------+---------------------------+
    | direct solver | square matrix             | rectangular matrix        |
    +---------------+-----+-------------+-------+-----+-------------+-------+
    |               | sym | sym pattern | unsym | sym | sym pattern | unsym |
    +===============+=====+=============+=======+=====+=============+=======+
    | SuperLU_DIST  | yes | yes         | yes   | yes | yes         | yes   |
    +---------------+-----+-------------+-------+-----+-------------+-------+
    | MUMPS         | yes | yes         | yes   | no  | no          | no    |
    +---------------+-----+-------------+-------+-----+-------------+-------+
    | Pastix        | yes | yes         | no    | no  | no          | no    |
    +---------------+-----+-------------+-------+-----+-------------+-------+

Parallel sparse iterative solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Concerning iterative solvers, we have chosen `pARMS <https://www-users.cs.umn.edu/~saad/software/pARMS/>`__, `HIPS <http://hips.gforge.inria.fr/>`__ and `Hypre <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__.

Each software implements a different type of parallel preconditioner.

So, pARMS implements algebraic domain decomposition preconditioner type such as additive Schwartz [CAI1989]_ and interface method; while HIPS implement hierarchical incomplete factorization and finally HYPRE implements multilevel preconditioner are AMG(Algebraic MultiGrid) and parallel approximated inverse.

To use one of these programs in **FreeFEM**, you have to install it independently of **FreeFEM**.
It is also necessary to install the MPI communication library which is essential for communication between the processors and, in some cases, software partitioning graphs like `METIS <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`__ or `Scotch <http://www.labri.fr/perso/pelegrin/scotch/>`__.

All this preconditioners are used with Krylov subspace methods accelerators.

Krylov subspace methods are iterative methods which consist in finding a solution :math:`x` of linear system :math:`Ax=b` inside the affine space :math:`x_0+K_m` by imposing that :math:`b-Ax \bot \mathcal{L}_m`, where :math:`K_m` is Krylov subspace of dimension :math:`m` defined by :math:`K_m=\{r_0, Ar_0, A^2r_0,...,A^{m-1}r_0\}` and :math:`\mathcal{L}_m` is another subspace of dimension :math:`m` which depends on type of Krylov subspace. For example in GMRES, :math:`\mathcal{L}_m=AK_m`.

We realized an interface which is easy to use, so that the call of these different softwares in **FreeFEM** is done in the same way.
You just have to load the solver and then specify the parameters to apply to the specific solvers.
In the rest of this chapter, when we talk about Krylov subspace methods we mean one among GMRES, CG and BICGSTAB.

pARMS solver
^^^^^^^^^^^^

`pARMS <https://www-users.cs.umn.edu/~saad/software/pARMS/>`__ (parallel Algebraic Multilevel Solver) is a software developed by Youssef Saad and al at University of Minnesota.

This software is specialized in the resolution of large sparse non symmetric linear systems of equation.
Solvers developed in pARMS are of type "Krylov’s subspace".

It consists of variants of GMRES like FGMRES (Flexible GMRES), DGMRES (Deflated GMRES) [SAAD2003]_ and BICGSTAB.
pARMS also implements parallel preconditioner like RAS (Restricted Additive Schwarz) [CAI1989]_ and Schur Complement type preconditioner.

All these parallel preconditioners are based on the principle of domain decomposition.
Thus, the matrix :math:`A` is partitioned into sub matrices :math:`A_i`\ (:math:`i=1,...,p`) where p represents the number of partitions one needs.
The union of :math:`A_i` forms the original matrix.
The solution of the overall system is obtained by solving the local systems on :math:`A_i` (see [SMITH1996]_).
Therefore, a distinction is made between iterations on :math:`A` and the local iterations on :math:`A_i`.

To solve the local problem on :math:`A_i` there are several preconditioners as **ilut** (Incomplete LU with threshold), **iluk** (Incomplete LU with level of fill in) and **ARMS** (Algebraic Recursive Multilevel Solver).

.. tip:: Default parameters

    .. code-block:: freefem
        :linenos:

        load "parms_FreeFem" //Tell FreeFem that you will use pARMS

        // Mesh
        border C(t=0, 2*pi){x=cos(t); y=sin(t); label=1;}
        mesh Th = buildmesh (C(50));

        // Fespace
        fespace Vh(Th, P2); Vh u, v;

        // Function
        func f= x*y;

        // Problem
        problem Poisson (u, v, solver=sparsesolver)
            = int2d(Th)(
                  dx(u)*dx(v)
                + dy(u)*dy(v) )
            + int2d(Th)(
                - f*v
            )
            + on(1, u=0) ;

        // Solve
        real cpu = clock();
        Poisson;
        cout << " CPU time = " << clock()-cpu << endl;

        // Plot
        plot(u);

    In line 1, the pARMS dynamic library is loaded with interface **FreeFEM**.
    After this, in line 15 we specify that the bilinear form will be solved by the last sparse linear solver load in memory which, in this case, is pARMS.

    The parameters used in pARMS in this case are the default one since the user does not have to provide any parameter.

    .. note:: In order to see the plot of a parallel script, run the command ``FreeFem++-mpi -glut ffglut script.edp``

Here are some default parameters:

-  ``solver=FGMRES``,
-  ``Krylov dimension=30``,
-  ``Maximum of Krylov=1000``,
-  ``Tolerance for convergence=1e-08`` (see book
   [SAAD2003]_ to understand all this parameters),
-  ``preconditionner=Restricted Additif Schwarz``
   [CAI1989]_,
-  ``Inner Krylov dimension=5``,
-  ``Maximum of inner Krylov dimension=5``,
-  ``Inner preconditionner=ILUK``.

To specify the parameters to apply to the solver, the user can either give an integer vector for **integer parameters** and real vectors for **real parameters** or provide a **file** which contains those parameters.

.. tip:: User specifies parameters inside two vectors

    Lets us consider Navier-Stokes example.
    In this example we solve linear systems coming from discretization of Navier-Stokes equations with pARMS.
    Parameters of solver is specified by user.

    .. code-block:: freefem
        :linenos:

        load "parms_FreeFem"

        // Parameters
        real nu = 1.;
        int[int] iparm(16);
        real[int] dparm(6);
        for (int ii = 0; ii < 16; ii++)
            iparm[ii] = -1;
        for (int ii = 0; ii < 6; ii++)
            dparm[ii] = -1.0; iparm[0]=0;

        // Mesh
        mesh Th = square(10, 10);
        int[int] wall = [1, 3];
        int inlet = 4;

        // Fespace
        fespace Vh(Th, [P2, P2, P1]);

        // Function
        func uc = 1.;

        // Problem
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
            + on(inlet, u=uc, v=0) ;

        matrix AA = Stokes(Vh, Vh);
        set(AA, solver=sparsesolver, lparams=iparm, dparams=dparm); //set pARMS as linear solver
        real[int] bb = Stokes(0, Vh);
        real[int] sol(AA.n);
        sol = AA^-1 * bb;

    We need two vectors to specify the parameters of the linear solver.
    In line 5-6 of the example, we have declared these vectors(:freefem:`int[int] iparm(16); real[int] dparm(6);`).
    In line 7-10 we have initialized these vectors by negative values.

    We do this because all parameters values in pARMS are positive and if you do not change the negative values of one entry of this vector, the default value will be set.

    In :numref:`tabParallelizationLparams` and :numref:`tabParallelizationDparams`, we have the meaning of different entries of these vectors.

    We run this example on a cluster paradent of Grid5000 and report results in :numref:`tabParallelizationConvergenceTime`.

    .. table:: Convergence and time for solving linear system
        :name: tabParallelizationConvergenceTime

        +----------------------------------+
        | :math:`n=471281`                 |
        | :math:`nnz=13\times10^6`         |
        | :math:`Te=571.29`                |
        +----+--------------+--------------+
        | np | add(iluk)    | shur(iluk)   |
        +----+-----+--------+-----+--------+
        |    | nit | time   | nit | time   |
        +====+=====+========+=====+========+
        | 4  | 230 | 637.57 | 21  | 557.8  |
        +----+-----+--------+-----+--------+
        | 8  | 240 | 364.12 | 22  | 302.25 |
        +----+-----+--------+-----+--------+
        | 16 | 247 | 212.07 | 24  | 167.5  |
        +----+-----+--------+-----+--------+
        | 32 | 261 | 111.16 | 25  | 81.5   |
        +----+-----+--------+-----+--------+

    .. table:: Legend of :numref:`tabParallelizationConvergenceTime`
        :name: tabParallelizationLegend

        +----------+---------------------------------------------+
        | ``n``    | matrix size                                 |
        +==========+=============================================+
        | ``nnz``  | number of non null entries inside matrix    |
        +----------+---------------------------------------------+
        | ``nit``  | number of iteration for convergence         |
        +----------+---------------------------------------------+
        | ``time`` | Time for convergence                        |
        +----------+---------------------------------------------+
        | ``Te``   | Time for constructing finite element matrix |
        +----------+---------------------------------------------+
        | ``np``   | number of processor                         |
        +----------+---------------------------------------------+

    In this example, we fix the matrix size (in term of finite element, we fix the mesh) and increase the number of processors used to solve the linear system.
    We saw that, when the number of processors increases, the time for solving the linear equation decreases, even if the number of iteration increases.
    This proves that, using pARMS as solver of linear systems coming from discretization of partial differential equation in **FreeFEM** can decrease drastically the total time of simulation.

.. table:: Meaning of :freefem:`lparams` corresponding variables
    :name: tabParallelizationLparams

    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | Entries of :freefem:`iparm` | Significations of each entries                                                                         |
    +=============================+========================================================================================================+
    | ``iparm[0]``                | Krylov subspace methods                                                                                |
    |                             |                                                                                                        |
    |                             | Different values for this parameters are specify on :numref:`tabParallelizationParmsKrylov`            |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[1]``                | Preconditionner                                                                                        |
    |                             |                                                                                                        |
    |                             | Different preconditionners for this parameters are  specify on :numref:`tabParallelizationParmsPrecon` |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[2]``                | Krylov subspace dimension in outer iteration: default value 30                                         |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[3]``                | Maximum of iterations in outer iteration: default value 1000                                           |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[4]``                | Number of level in arms when used                                                                      |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[5]``                | Krylov subspace dimension in inner iteration: default value 3                                          |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[6]``                | Maximum of iterations in inner iteration: default value 3                                              |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[7]``                | Symmetric(=1 for symmetric) or unsymmetric matrix:                                                     |
    |                             |                                                                                                        |
    |                             | default value 0(unsymmetric matrix)                                                                    |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[8]``                | Overlap size between different subdomain: default value 0(no overlap)                                  |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[9]``                | Scale the input matrix or not: Default value 1 (Matrix should be scaled)                               |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[10]``               | Block size in arms when used: default value 20                                                         |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[11]``               | lfil0 (ilut, iluk, and arms) : default value 20                                                        |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[12]``               | lfil for Schur complement const : default value 20                                                     |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[13]``               | lfil for Schur complement const : default value 20                                                     |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[14]``               | Multicoloring or not in ILU when used : default value 1                                                |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[15]``               | Inner iteration : default value 0                                                                      |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+
    | ``iparm[16]``               | Print message when solving: default 0 (no message print)                                               |
    |                             |                                                                                                        |
    |                             | * 0: no message is print,                                                                              |
    |                             |                                                                                                        |
    |                             | * 1: Convergence informations like number of iteration and residual,                                   |
    |                             |                                                                                                        |
    |                             | * 2: Timing for a different step like preconditioner,                                                  |
    |                             |                                                                                                        |
    |                             | * 3 : Print all informations                                                                           |
    +-----------------------------+--------------------------------------------------------------------------------------------------------+

.. table:: Significations of :freefem:`dparams` corresponding variables
    :name: tabParallelizationDparams

    +-----------------------------+---------------------------------------------------------------------+
    | Entries of :freefem:`dparm` | Significations of each entries                                      |
    +=============================+=====================================================================+
    | ``dparm[0]``                | precision for outer iteration : default value 1e-08                 |
    +-----------------------------+---------------------------------------------------------------------+
    | ``dparm[1]``                | precision for inner iteration: default value 1e-2                   |
    +-----------------------------+---------------------------------------------------------------------+
    | ``dparm[2]``                | tolerance used for diagonal domain: : default value 0.1             |
    +-----------------------------+---------------------------------------------------------------------+
    | ``dparm[3]``                | drop tolerance droptol0 (ilut, iluk, and arms) : default value 1e-2 |
    +-----------------------------+---------------------------------------------------------------------+
    | ``dparm[4]``                | droptol for Schur complement const: default value 1e-2              |
    +-----------------------------+---------------------------------------------------------------------+
    | ``dparm[5]``                | droptol for Schur complement const: default value 1e-2              |
    +-----------------------------+---------------------------------------------------------------------+

.. table:: Krylov Solvers in pARMS
    :name: tabParallelizationParmsKrylov

    +------------------------+-------------------------+
    | Values of ``iparm[0]`` | Krylov subspace methods |
    +========================+=========================+
    | 0                      | FGMRES (Flexible GMRES) |
    +------------------------+-------------------------+
    | 1                      | DGMRES (Deflated GMRES) |
    +------------------------+-------------------------+
    | 2                      | BICGSTAB                |
    +------------------------+-------------------------+

.. table:: Preconditionners in pARMS
    :name: tabParallelizationParmsPrecon

    +------------------------+-------------------------------------------------------------------------+
    | Values of ``iparm[1]`` | Preconditionners type                                                   |
    +========================+=========================================================================+
    | 0                      | additive Schwartz preconditioner with ilu0 as local preconditioner      |
    +------------------------+-------------------------------------------------------------------------+
    | 1                      | additive Schwartz preconditioner with iluk as local preconditioner      |
    +------------------------+-------------------------------------------------------------------------+
    | 2                      | additive Schwartz preconditioner with ilut as local preconditioner      |
    +------------------------+-------------------------------------------------------------------------+
    | 3                      | additive Schwartz preconditioner with arms as local preconditioner      |
    +------------------------+-------------------------------------------------------------------------+
    | 4                      | Left Schur complement preconditioner with ilu0 as local preconditioner  |
    +------------------------+-------------------------------------------------------------------------+
    | 5                      | Left Schur complement preconditioner with ilut as local preconditioner  |
    +------------------------+-------------------------------------------------------------------------+
    | 6                      | Left Schur complement preconditioner with iluk as local preconditioner  |
    +------------------------+-------------------------------------------------------------------------+
    | 7                      | Left Schur complement preconditioner with arms as local preconditioner  |
    +------------------------+-------------------------------------------------------------------------+
    | 8                      | Right Schur complement preconditioner with ilu0 as local preconditioner |
    +------------------------+-------------------------------------------------------------------------+
    | 9                      | Right Schur complement preconditioner with ilut as local preconditioner |
    +------------------------+-------------------------------------------------------------------------+
    | 10                     | Right Schur complement preconditioner with iluk as local preconditioner |
    +------------------------+-------------------------------------------------------------------------+
    | 11                     | Right Schur complement preconditioner with arms as local preconditioner |
    +------------------------+-------------------------------------------------------------------------+
    | 12                     | sch_gilu0, Schur complement preconditioner with global ilu0             |
    +------------------------+-------------------------------------------------------------------------+
    | 13                     | SchurSymmetric GS preconditioner                                        |
    +------------------------+-------------------------------------------------------------------------+

Interfacing with HIPS
^^^^^^^^^^^^^^^^^^^^^

`HIPS <http://hips.gforge.inria.fr/>`__ (*Hierarchical Iterative Parallel Solver*) is a scientific library that provides an efficient parallel iterative solver for very large sparse linear systems.
HIPS is available as free software under the CeCILL-C licence.

HIPS implements two solver classes which are the iteratives class (GMRES, PCG) and the Direct class.
Concerning preconditionners, HIPS implements a type of multilevel ILU.
For further informations on those preconditionners see the `HIPS documentation <http://hips.gforge.inria.fr/doc/hips_user.pdf>`__.

.. tip:: Laplacian 3D solved with HIPS

    Let us consider the 3D Laplacian example inside **FreeFEM** package where after discretization we want to solve the linear equation with HIPS.

    The following example is a Laplacian 3D using Hips as linear solver.
    We first load Hips solver at line 2.
    From line 7 to 18 we specify the parameters for the Hips solver and in line 82 we set these parameters in the linear solver.

    In :numref:`tabParallelizationConvergenceTimeHips` results of running on Cluster Paradent of Grid5000 are reported.
    We can see in this running example the efficiency of parallelism.

    .. code-block:: freefem
        :linenos:

        load "msh3"
        load "hips_FreeFem" //load Hips library

        // Parameters
        int nn = 10;
        real zmin = 0, zmax = 1;
        int[int] iparm(14);
        real[int] dparm(6);
        for (int iii = 0; iii < 14; iii++)
            iparm[iii] = -1;
        for (int iii = 0; iii < 6; iii++)
            dparm[iii] = -1;
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

        mesh3 Th=buildlayers(Th2, nn, zbound=[zmin, zmax], reffacemid=rmid,
            reffaceup = rup, reffacelow = rdown);

        // Fespace
        fespace Vh2(Th2, P2);
        Vh2 ux, uz, p2;

        fespace Vh(Th, P2);
        Vh uhe = ue;
        cout << "uhe min =" << uhe[].min << ", max =" << uhe[].max << endl;
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
                  ue*v + (uex*N.x + uey*N.y + uez*N.z)*v
            )
            + on(1, u=ue);

        varf l (unused, v) = int3d(Th)(f*v);

        real cpu=clock();
        matrix Aa = va(Vh, Vh);

        F[] = va(0, Vh);

        if (mpirank == 0){
            cout << "Size of A =" << Aa.n << endl;
            cout << "Non zero coefficients =" << Aa.nbcoef << endl;
            cout << "CPU TIME FOR FORMING MATRIX =" << clock()-cpu << endl;
        }

        set(Aa, solver=sparsesolver, dparams=dparm, lparams=iparm); //Set hips as linear solver

        // Solve
        u[] = Aa^-1*F[];

        // Plot
        plot(u);

    .. table:: Legend of this table are give in :numref:`tabParallelizationLegend`
        :name: tabParallelizationConvergenceTimeHips

        +-----------------------------+
        | :math:`n=4\times 10^6`      |
        | :math:`nnz=118 \times 10^6` |
        | :math:`Te=221.34`           |
        +--------+---------+----------+
        | ``np`` | ``nit`` | ``time`` |
        +========+=========+==========+
        | 8      | 190     | 120.34   |
        +--------+---------+----------+
        | 16     | 189     | 61.08    |
        +--------+---------+----------+
        | 32     | 186     | 31.70    |
        +--------+---------+----------+
        | 64     | 183     | 23.44    |
        +--------+---------+----------+

.. tip::

    .. table:: Significations of :freefem:`lparams` corresponding to HIPS interface
        :name: tabParallelizationHipsInterface

        +----------------------+--------------------------------------------------------------------+
        | Entries of ``iparm`` | Significations of each entries                                     |
        +======================+====================================================================+
        | ``iparm[0]``         | Strategy use for solving (Iterative=0 or Hybrid=1 or Direct=2).    |
        |                      |                                                                    |
        |                      | Defaults values are : Iterative                                    |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[1]``         | Krylov methods.                                                    |
        |                      |                                                                    |
        |                      | If iparm[0]=0, give type of Krylov methods: 0 for GMRES, 1 for PCG |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[2]``         | Maximum of iterations in outer iteration: default value 1000       |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[3]``         | Krylov subspace dimension in outer iteration: default value 40     |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[4]``         | Symmetric(=0 for symmetric) and 1 for unsymmetricmatrix:           |
        |                      |                                                                    |
        |                      | default value 1 (unsymmetric matrix)                               |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[5]``         | Pattern of matrix are symmetric or not: default value 0            |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[6]``         | Partition type of input matrix: default value 0                    |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[7]``         | Number of level that use the HIPS locally consistentfill-in:       |
        |                      |                                                                    |
        |                      | Default value 2                                                    |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[8]``         | Numbering in indices array will start at 0 or 1: Default value 0   |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[9]``         | Scale matrix. Default value 1                                      |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[10]``        | Reordering use inside subdomains for reducingfill-in:              |
        |                      |                                                                    |
        |                      | Only use for iterative. Default value 1                            |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[11]``        | Number of unknowns per node in the matrix non-zeropattern graph:   |
        |                      |                                                                    |
        |                      | Default value 1                                                    |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[12]``        | This value is used to set the number of time the                   |
        |                      |                                                                    |
        |                      | normalization is applied to the matrix: Default 2.                 |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[13]``        | Level of informations printed during solving: Default 5.           |
        +----------------------+--------------------------------------------------------------------+
        | ``iparm[14]``        | HIPS_DOMSIZE Subdomain size                                        |
        +----------------------+--------------------------------------------------------------------+

    .. table:: Significations of :freefem:`dparams` corresponding to HIPS interface
        :name: tabParallelizationHipsDparms

        +--------------+--------------------------------------------------------------------------+
        | ``dparm[0]`` | HIPS_PREC: Relative residual norm: Default=1e-9                          |
        +--------------+--------------------------------------------------------------------------+
        | ``dparm[1]`` | HIPS_DROPTOL0: Numerical threshold in ILUT for interior domain           |
        |              |                                                                          |
        |              | (important : set 0.0 in HYBRID: Default=0.005)                           |
        +--------------+--------------------------------------------------------------------------+
        | ``dparm[2]`` | HIPS_DROPTOL1 : Numerical threshold in ILUT for Schur preconditioner:    |
        |              |                                                                          |
        |              | Default=0.005                                                            |
        +--------------+--------------------------------------------------------------------------+
        | ``dparm[3]`` | HIPS_DROPTOLE : Numerical threshold for coupling between the interior    |
        |              |                                                                          |
        |              | level and Schur: Default 0.005                                           |
        +--------------+--------------------------------------------------------------------------+
        | ``dparm[4]`` | HIPS_AMALG : Numerical threshold for coupling between the interior level |
        |              |                                                                          |
        |              | and Schur: Default=0.005                                                 |
        +--------------+--------------------------------------------------------------------------+
        | ``dparm[5]`` | HIPS_DROPSCHUR : Numerical threshold for coupling between the interior   |
        |              |                                                                          |
        |              | level and Schur: Default=0.005                                           |
        +--------------+--------------------------------------------------------------------------+


Interfacing with HYPRE
^^^^^^^^^^^^^^^^^^^^^^

`Hypre <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__ (High Level Preconditioner) is a suite of parallel preconditioner developed at Lawrence Livermore National Lab.

There are two main classes of preconditioners developed in HYPRE: AMG (Algebraic MultiGrid) and Parasails (Parallel Sparse Approximate Inverse).

Now, suppose we want to solve :math:`Ax=b`.

At the heart of AMG there is a series of progressively coarser (smaller) representations of the matrix :math:`A`.
Given an approximation :math:`\hat{x}` to the solution :math:`x`, consider solving the residual equation :math:`Ae=r` to find the error :math:`e`, where :math:`r=b-A\hat{x}`.
A fundamental principle of AMG is that it is an algebraically smooth error.
To reduce the algebraically smooth errors further, they need to be represented by a smaller defect equation (coarse grid residual equation) :math:`A_ce_c=r_c`, which is cheaper to solve.
After solving this coarse equation, the solution is then interpolated in fine grid represented here by matrix :math:`A`.
The quality of AMG depends on the choice of coarsening and interpolating operators.

The *sparse approximate inverse* approximates the inverse of a matrix :math:`A` by a sparse matrix :math:`M`.
A technical idea to construct matrix :math:`M` is to minimize the Frobenuis norm of the residual matrix :math:`I-MA`.
For more details on this preconditioner technics see [CHOW1997]_.

HYPRE implement three Krylov subspace solvers: GMRES, PCG and BiCGStab.

.. tip:: Laplacian 3D solved with HYPRE

    Let us consider again the 3D Laplacian example inside **FreeFEM** package where after discretization we want to solve the linear equation with Hypre.
    The following example is a Laplacian 3D using Hypre as linear solver.
    This is the same example as Hips one, so we just show here the lines where we set some Hypre parameters.

    We first load the Hypre solver at line 2.
    From line 6 to 18 we specifies the parameters to set to Hypre solver and in line 22 we set parameters to Hypre solver.

    It should be noted that the meaning of the entries of these vectors is different from those of Hips.
    In the case of HYPRE, the meaning of differents entries of vectors :freefem:`iparm` and :freefem:`dparm` are given in :numref:`tabParallelizationIparmsDparmsHypre` to :numref:`tabParallelizationIparmsDparmsSchwartz`.

    In :numref:`tabParallelizationConvergenceTimeHypre` the results of running on Cluster Paradent of Grid5000 are reported.
    We can see in this running example the efficiency of parallelism, in particular when AMG are use as preconditioner.

    .. code-block:: freefem
        :linenos:

        load "msh3"
        load "hipre_FreeFem" //Load Hipre librairy

        // Parameters
        int nn = 10;
        int[int] iparm(20);
        real[int] dparm(6);
        for (int iii = 0; iii < 20; iii++)
            iparm[iii] = -1;
        for (int iii = 0; iii < 6; iii++)
            dparm[iii] = -1;
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

.. table:: Definitions of common entries of :freefem:`iparms` and :freefem:`dparms` vectors for every preconditioner in HYPRE
    :name: tabParallelizationIparmsDparmsHypre

    +---------------+--------------------------------------------------------------+
    | ``iparms[0]`` | Solver identification:                                       |
    |               |                                                              |
    |               | 0: BiCGStab, 1: GMRES, 2: PCG. Default=1                     |
    +---------------+--------------------------------------------------------------+
    | ``iparms[1]`` | Preconditioner identification:                               |
    |               |                                                              |
    |               | 0: BOOMER AMG, 1: PILUT, 2: Parasails, 3: Schwartz Default=0 |
    +---------------+--------------------------------------------------------------+
    | ``iparms[2]`` | Maximum of iteration: Default=1000                           |
    +---------------+--------------------------------------------------------------+
    | ``iparms[3]`` | Krylov subspace dim: Default= 40                             |
    +---------------+--------------------------------------------------------------+
    | ``iparms[4]`` | Solver print info level: Default=2                           |
    +---------------+--------------------------------------------------------------+
    | ``iparms[5]`` | Solver log: Default=1                                        |
    +---------------+--------------------------------------------------------------+
    | ``iparms[6]`` | Solver stopping criteria only for BiCGStab : Default=1       |
    +---------------+--------------------------------------------------------------+
    | ``dparms[0]`` | Tolerance for convergence: Default=:math:`1.0e-11`           |
    +---------------+--------------------------------------------------------------+

.. table:: Definitions of other entries of :freefem:`iparms` and :freefem:`dparms` if preconditioner is BOOMER AMG
    :name: tabParallelizationIparmsDparmsBoomer

    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[7]``  | AMG interpolation type: Default=6                                                   |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[8]``  | Specifies the use of GSMG - geometrically smooth coarsening and                     |
    |                |                                                                                     |
    |                | interpolation: Default=1                                                            |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[9]``  | AMG coarsen type: Default=6                                                         |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[10]`` | Defines whether local or global measures are used: Default=1                        |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[11]`` | AMG cycle type: Default=1                                                           |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[12]`` | AMG Smoother type: Default=1                                                        |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[13]`` | AMG number of levels for smoothers: Default=3                                       |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[14]`` | AMG number of sweeps for smoothers: Default=2                                       |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[15]`` | Maximum number of multigrid levels: Default=25                                      |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[16]`` | Defines which variant of the Schwartz method isused:                                |
    |                |                                                                                     |
    |                | 0: hybrid multiplicative Schwartz method (no overlap across processor boundaries)   |
    |                |                                                                                     |
    |                | 1: hybrid additive Schwartz method (no overlap across processor boundaries)         |
    |                |                                                                                     |
    |                | 2: additive Schwartz method                                                         |
    |                |                                                                                     |
    |                | 3: hybrid multiplicative Schwartz method (with overlap across processor boundaries) |
    |                |                                                                                     |
    |                | Default=1                                                                           |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[17]`` | Size of the system of PDEs: Default=1                                               |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[18]`` | Overlap for the Schwarz method: Default=1                                           |
    +----------------+-------------------------------------------------------------------------------------+
    | ``iparms[19]`` | Type of domain used for the Schwarz method                                          |
    |                |                                                                                     |
    |                | 0: each point is a domain                                                           |
    |                |                                                                                     |
    |                | 1: each node is a domain (only of interest in "systems" AMG)                        |
    |                |                                                                                     |
    |                | 2: each domain is generated by agglomeration (default)                              |
    +----------------+-------------------------------------------------------------------------------------+
    | ``dparms[1]``  | AMG strength threshold: Default=0.25                                                |
    +----------------+-------------------------------------------------------------------------------------+
    | ``dparms[2]``  | Truncation factor for the interpolation: Default=1e-2                               |
    +----------------+-------------------------------------------------------------------------------------+
    | ``dparms[3]``  | Sets a parameter to modify the definition of strength for                           |
    |                |                                                                                     |
    |                | diagonal dominant portions of the matrix: Default=0.9                               |
    +----------------+-------------------------------------------------------------------------------------+
    | ``dparms[3]``  | Defines a smoothing parameter for the additive Schwartz method. Default=1           |
    +----------------+-------------------------------------------------------------------------------------+

.. table:: Definitions of other entries of :freefem:`iparms` and :freefem:`dparms` if preconditioner is PILUT

    +---------------+-------------------------------------------------+
    | ``iparms[7]`` | Row size in Parallel ILUT: Default=1000         |
    +---------------+-------------------------------------------------+
    | ``iparms[8]`` | Set maximum number of iterations: Default=30    |
    +---------------+-------------------------------------------------+
    | ``dparms[1]`` | Drop tolerance in Parallel ILUT: Default=1e-5   |
    +---------------+-------------------------------------------------+

.. table:: Definitions of other entries of :freefem:`iparms` and :freefem:`dparms` if preconditioner is ParaSails

    +---------------+------------------------------------------------------------------------------------------------+
    | ``iparms[7]`` | Number of levels in Parallel Sparse Approximate inverse: Default=1                             |
    +---------------+------------------------------------------------------------------------------------------------+
    | ``iparms[8]`` | Symmetric parameter for the ParaSails preconditioner:                                          |
    |               |                                                                                                |
    |               | 0: nonsymmetric and/or indefinite problem, and nonsymmetric preconditioner                     |
    |               |                                                                                                |
    |               | 1: SPD problem, and SPD (factored) preconditioner                                              |
    |               |                                                                                                |
    |               | 2: nonsymmetric, definite problem, and SPD (factored) preconditioner                           |
    |               |                                                                                                |
    |               | Default=0                                                                                      |
    +---------------+------------------------------------------------------------------------------------------------+
    | ``dparms[1]`` | Filters parameters. The filter parameter is used to drop small nonzeros in the preconditioner, |
    |               |                                                                                                |
    |               | to reduce the cost of applying the preconditioner: Default=0.1                                 |
    +---------------+------------------------------------------------------------------------------------------------+
    | ``dparms[2]`` | Threshold parameter: Default=0.1                                                               |
    +---------------+------------------------------------------------------------------------------------------------+

.. table:: Definitions of other entries of :freefem:`iparms` and :freefem:`dparms` if preconditionner is Schwartz
    :name: tabParallelizationIparmsDparmsSchwartz

    +---------------+-------------------------------------------------------------------------------------+
    | ``iparms[7]`` | Defines which variant of the Schwartz method isused:                                |
    |               |                                                                                     |
    |               | 0: hybrid multiplicative Schwartz method (no overlap across processor boundaries)   |
    |               |                                                                                     |
    |               | 1: hybrid additive Schwartz method (no overlap across processor boundaries)         |
    |               |                                                                                     |
    |               | 2: additive Schwartz method                                                         |
    |               |                                                                                     |
    |               | 3: hybrid multiplicative Schwartz method (with overlap across processor boundaries) |
    |               |                                                                                     |
    |               | Default=1                                                                           |
    +---------------+-------------------------------------------------------------------------------------+
    | ``iparms[8]`` | Overlap for the Schwartz method: Default=1                                          |
    +---------------+-------------------------------------------------------------------------------------+
    | ``iparms[9]`` | Type of domain used for the Schwartz method                                         |
    |               |                                                                                     |
    |               | 0: each point is a domain                                                           |
    |               |                                                                                     |
    |               | 1: each node is a domain (only of interest in "systems" AMG)                        |
    |               |                                                                                     |
    |               | 2: each domain is generated by agglomeration (default)                              |
    +---------------+-------------------------------------------------------------------------------------+

.. table:: Convergence and time for solving linear system
    :name: tabParallelizationConvergenceTimeHypre

    +-----------------------+--------------------------+---------------------+
    | :math:`n=4\times10^6` | :math:`nnz=13\times10^6` | :math:`Te = 571.29` |
    +-----------------------+--------------------------+---------------------+
    | np                    | AMG                                            |
    +-----------------------+-------+----------------------------------------+
    |                       | `nit` | `time`                                 |
    +=======================+=======+========================================+
    | 8                     | 6     | 1491.83                                |
    +-----------------------+-------+----------------------------------------+
    | 16                    | 5     | 708.49                                 |
    +-----------------------+-------+----------------------------------------+
    | 32                    | 4     | 296.22                                 |
    +-----------------------+-------+----------------------------------------+
    |  64                   | 4     | 145.64                                 |
    +-----------------------+-------+----------------------------------------+

Conclusion
^^^^^^^^^^

With the different runs presented here, we wanted to illustrate the gain in time when we increase the number of processors used for the simulations.
We saw that in every case the time for the construction of the finite element matrix is constant.
This is normal because until now this phase is sequential in **FreeFEM**.
In contrast, phases for solving the linear system are parallel.
We saw on several examples presented here that when we increase the number of processors, in general we decrease the time used for solving the linear systems.
But this is not true in every case.
In several case, when we increase the number of processors the time to convergence also increases.
There are two main reasons for this.
First, the increase of processors can lead to the increase of volume of exchanged data across processors consequently increasing the time for solving the linear systems.

Furthermore, in decomposition domain type preconditioners, the number of processors generally corresponds to the number of sub domains.
In subdomain methods, generally when we increase the number of subdomains we decrease convergence quality of the preconditioner.
This can increase the time used for solving linear equations.

To end this, we should note that good use of the preconditioners interfaced in **FreeFEM** is empiric, because it is difficult to know what is a good preconditioner for some type of problems.
Although, the efficiency of preconditioners sometimes depends on how its parameters are set.
For this reason we advise the user to pay attention to the meaning of the parameters in the user guide of the iterative solvers interfaced in **FreeFEM**.

Domain decomposition
~~~~~~~~~~~~~~~~~~~~

In the previous section, we saw that the phases to construct a matrix are sequential.
One strategy to construct the matrix in parallel is to divide geometrically the domain into subdomains.
In every subdomain we construct a local submatrix and after that we assemble every submatrix to form the global matrix.

We can use this technique to solve PDE directly in domain :math:`\Omega`.
In this case, in every subdomains you have to define artificial boundary conditions to form consistent equations in every subdomains.
After this, you solve equation in every subdomains and define a strategy to obtain the global solution.

In terms of parallel programming for **FreeFEM**, with MPI, this means that the user must be able to divide processors avaible for computation into subgroups of processors and also must be able to realize different type of communications in **FreeFEM** script.
Here is a wrapper of some MPI functions.

Communicators and groups
^^^^^^^^^^^^^^^^^^^^^^^^

**Groups**

:freefem:`mpiGroup grpe(mpiGroup gp, KN_<long>)`: Create MPI_Group from existing group :freefem:`gp` by given vector.

**Communicators**

Communicators is an abstract MPI object which allows MPI user to communicate across group of processors.
Communicators can be Intra-communicators(involves a single group) or Inter-communicators (involves two groups).
When we not specify type of communicator it will be Intra-communicators

**mpiComm cc(mpiComm comm, mpiGroup gp):** Creates a new communicator.

:freefem:`comm` communicator(handle), :freefem:`gp` group which is a subset of the group of :freefem:`comm` (handle).
Return new communicator

**mpiComm cc(mpiGroup gp)**: Same as previous constructor but default :freefem:`comm` here is ``MPI_COMM_WORLD``.

**mpiComm cc(mpiComm comm, int color, int key):** Creates new communicators based on :freefem:`colors` and :freefem:`key`.
This constructor is based on MPI_Comm_split routine of MPI.

**mpiComm cc(MPIrank p, int key):** Same constructor than the last one.

Here :freefem:`colors` and :freefem:`comm` is defined in :freefem:`MPIrank`.
This constructor is based on ``MPI_Comm_split`` routine of MPI.

.. tip:: Split communicator

    .. code-block:: freefem
        :linenos:

        mpiComm comm(mpiCommWorld, 0, 0);
        int color = mpiRank(comm)%2;
        mpiComm ccc(processor(color, comm), 0);
        mpiComm qpp(comm, 0, 0);
        mpiComm cp(ccc, color, 0);

**mpiComm cc(mpiComm comm, int high):** Creates an intracommunicator from an intercommunicator. :freefem:`comm` intercommunicator, :freefem:`high`.

Used to order the groups within :freefem:`comm` (logical) when creating the new communicator.
This constructor is based on ``MPI_Intercomm_merge`` routine of MPI.

**mpiComm cc(MPIrank p1, MPIrank p2, int tag):** This constructor creates an intercommuncator from two intracommunicators.
:freefem:`p1` defined local (intra)communicator and rank in ``local_comm`` of leader (often 0) while :freefem:`p2` defined remote communicator and rank in ``peer_comm`` of remote leader (often 0).
:freefem:`tag` Message tag to use in constructing intercommunicator.
This constructor is based on ``MPI_Intercomm_create``.

.. tip:: Merge

    .. code-block:: freefem
        :linenos:

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

Process
^^^^^^^

In **FreeFEM** we wrap MPI process by function call :freefem:`processor` which create internal **FreeFEM** object call :freefem:`MPIrank`.
This mean that do not use :freefem:`MPIrank` in **FreeFEM** script.

:freefem:`processor(int rk)`: Keep process rank inside object :freefem:`MPIrank`.
Rank is inside ``MPI_COMM_WORLD``.

:freefem:`processor(int rk, mpiComm cc)` and :freefem:`processor(mpiComm cc, int rk)` process rank inside communicator cc.

:freefem:`processor(int rk, mpiComm cc)` and :freefem:`processor(mpiComm cc, int rk)` process rank inside communicator cc.

:freefem:`processorblock(int rk)`: This function is exactlly the same than :freefem:`processor(int rk)` but is use in case of blocking communication.

:freefem:`processorblock(int rk, mpiComm cc)`: This function is exactly the same as :freefem:`processor(int rk, mpiComm cc)` but uses a synchronization point.

Points to Points communicators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In **FreeFEM** you can call MPI points to points communications functions.

:freefem:`Send(processor(int rk, mpiComm cc), Data D)` : Blocking send of :freefem:`Data D` to processor of :freefem:`rank rk` inside communicator :freefem:`cc`.
Note that :freefem:`Data D` can be: :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`int[int]`, :freefem:`real[int]`, :freefem:`complex[int]`, :freefem:`Mesh`, :freefem:`Mesh3`, :freefem:`Matrix`.

:freefem:`Recv(processor(int rk, mpiComm cc), Data D)`: Receive :freefem:`Data D` from process of rank :freefem:`rk` in communicator :freefem:`cc`.

Note that :freefem:`Data D` can be: :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`int[int]`, :freefem:`real[int]`, :freefem:`complex[int]`, :freefem:`Mesh`, :freefem:`Mesh3`, :freefem:`Matrix` and should be the same type than corresponding send.

:freefem:`Isend(processor(int rk, mpiComm cc), Data D)` : Non blocking send of :freefem:`Data D` to processor of :freefem:`rank rk` inside communicator :freefem:`cc`.

Note that :freefem:`Data D` can be: :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`int[int]`, :freefem:`real[int]`, :freefem:`complex[int]`, :freefem:`mesh`, :freefem:`mesh3`, :freefem:`matrix`.

:freefem:`Recv(processor(int rk, mpiComm cc), Data D)`: Receive corresponding to send.

Global operations
^^^^^^^^^^^^^^^^^

In **FreeFEM** you can call MPI global communication functions.

:freefem:`broadcast(processor(int rk, mpiComm cc), Data D)`: Process :freefem:`rk` Broadcast :freefem:`Data D` to all process inside :freefem:`communicator cc`.
Note that :freefem:`Data D` can be: :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`int[int]`, :freefem:`real[int]`, :freefem:`complex[int]`, :freefem:`Mesh`, :freefem:`Mesh3`, :freefem:`Matrix`.

:freefem:`broadcast(processor(int rk), Data D)`: Process :freefem:`rk` Broadcast :freefem:`Data D` to all process inside ``MPI_COMM_WORLD``.
Note that :freefem:`Data D` can be: :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`int[int]`, :freefem:`real[int]`, :freefem:`complex[int]`, :freefem:`Mesh`, :freefem:`Mesh3`, :freefem:`Matrix`.

:freefem:`mpiAlltoall(Data a, Data b)`: Sends :freefem:`data a` from all to all processes.
Receive buffer is :freefem:`Data b`.
This is done inside communicator ``MPI_COMM_WORLD``.

:freefem:`mpiAlltoall(Data a, Data b, mpiComm cc)`: Sends :freefem:`data a` from all to all processes. Receive buffer is :freefem:`Data b`.
This is done inside communicator ``cc``.

:freefem:`mpiGather(Data a, Data b, processor(mpiComm, int rk)`: Gathers together values :freefem:`Data a` from a group of processes.
Process of rank :freefem:`rk` get data on communicator :freefem:`rk`.
This function is like ``MPI_Gather``.

:freefem:`mpiAllgather(Data a, Data b)`: Gathers :freefem:`Data a` from all processes and distribute it to all in :freefem:`Data b`.
This is done inside communicator ``MPI_COMM_WORLD``.
This function is like ``MPI_Allgather``.

:freefem:`mpiAllgather(Data a, Data b, mpiComm cc)`: Gathers :freefem:`Data a` from all processes and distribute it to all in :freefem:`Data b`.
This is done inside :freefem:`communicator cc`.
This function is like ``MPI_Allgather``.

:freefem:`mpiScatter(Data a,Data b,processor(int rk, mpiComm cc))`: Sends :freefem:`Data a` from one process whith rank :freefem:`rk` to all other processes in group represented by communicator :freefem:`mpiComm cc`.

:freefem:`mpiReduce(Data a, Data b, processor(int rk, mpiComm cc), MPI_Op op)` Reduces values :freefem:`Data a` on all processes to a single value :freefem:`Data b` on process of rank :freefem:`rk` and communicator :freefem:`cc`.

Operation use in reduce is: :freefem:`MPI_Op op` which can be: :freefem:`mpiMAX`, :freefem:`mpiMIN`, :freefem:`mpiSUM`, :freefem:`mpiPROD`, :freefem:`mpiLAND`, :freefem:`mpiLOR`, :freefem:`mpiLXOR`, :freefem:`mpiBAND`, :freefem:`mpiBXOR`, :freefem:`mpiMAXLOC`, :freefem:`mpiMINLOC`.

Note that, for all global operations, only :freefem:`int[int]` and :freefem:`real[int]` are data type take in account in **FreeFEM**.
