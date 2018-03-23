



## Eigenvalue Problems

This section depends on your installation of FreeFem++; you need to have compiled (see `:::freefem README\_arpack` $\codered$), ARPACK.
This tool is available in FreeFem++ if the word eigenvalue appears in line Load:, like:

```freefem
-- FreeFem++ v1.28 (date Thu Dec 26 10:56:34 CET 2002)
 file : LapEigenValue.edp
 Load: lg_fem lg_mesh eigenvalue
```

This tool is based on the `:::freefem arpack++` \footnote{\url{http://www.caam.rice.edu/software/ARPACK/}}
the object-oriented version of ARPACK eigenvalue package \cite{arpack} $\codered$.

The function EigenValue computes the generalized eigenvalue
of  $A u = \lambda B u$. The Shift-invert method is used by default, with sigma =$\sigma$ the shift of the method.
The matrix $OP$ is defined with $A - \sigma B$.
The return value is the number of converged eigenvalues (can be greater than the number of requested eigenvalues nev=)

```freefem
int k=EigenValue(OP,B,nev= , sigma= );
```

where the matrix $OP=  A - \sigma B $ with a solver and boundary condition,
and the matrix $B$.

From version 3.31 there is also a functional interface:

```freefem
int k=EigenValue(n,FOP1,FB,nev= , sigma= );
```

where $n$ is the size of the problem,  and the operators are now defined through functions, defining respectively the matrix product of $OP^{-1}$ and $B$, as in

```freefem
int n = OP1.n;
func real[int] FOP1(real[int] & u) { real[int] Au=OP^-1*u;return Au;}
func real[int] FB(real[int] & u) { real[int] Au=B*u;return Au;}
```

If you want finer control over the method employed in ARPACK, you can specify which mode ARPACK will work with (mode= , see ARPACK documentation). The operators necessary for the chosen mode can be passed through the optional parameters A=, A1= , B=, B1= (see below).

* mode=1: Regular mode for solving $A u = \lambda u$

	```freefem
	int k=EigenValue(n,A=FOP,mode=1,nev= );
	```

	where the function FOP defines the matrix product of A

* mode=2: Regular inverse mode for solving $A u = \lambda B u$

	```freefem
	int k=EigenValue(n,A=FOP,B=FB,B1=FB1,mode=2, nev= );
	```

	where the functions FOP, FB and FB1 define respectively the matrix product of $A$, $B$ and $B^{-1}$

* mode=3: Shift-invert mode for solving $A u = \lambda B u$

	```freefem
	int k=EigenValue(n,A1=FOP1,B=FB,mode=3,sigma=sigma, nev= );
	```

	where the functions FOP1 and FB define respectively the matrix product of $OP^{-1} = (A - \sigma B)^{-1}$ and $B$

You can also specify which subset of eigenvalues you want to compute (which= ). The default value is which="LM", for eigenvalues with largest magnitude. "SM" is for smallest magnitude, "LA" for largest algebraic value, "SA" for smallest algebraic value, and "BE" for both ends of the spectrum.

Remark: For complex problems, you need to use the keyword `:::freefem complexEigenValue` instead of `:::freefem EigenValue` when passing operators through functions.

!!! note
	Boundary condition and Eigenvalue Problems

	The locking (Dirichlet) boundary condition is make with exact penalization so we put 1e30=tgv on the diagonal term of the locked degree of freedom (see equation (\ref{eq tgv} $\codered$)). So take Dirichlet boundary condition just on $A$ and not on $B$. because we solve $ w=OP^-1*B*v$.

	If you put locking (Dirichlet ) boundary condition on $B$ matrix (with key work `:::freefem on`) you get small spurious modes $(10^{-30})$, due to boundary condition, but if you forget the locking boundary condition on $B$ matrix (no key work "on") you get huge spurious $(10^{30})$ modes associated to these boundary conditons. We compute only small mode, so we get the good one in this case.

* `:::freefem sym=`
	The problem is symmetric (all the eigen value are real)

* `:::freefem nev=`
	The number desired eigenvalues (nev)  close to the shift.

* `:::freefem value=`
	The array to store the real part of the eigenvalues

* `:::freefem ivalue=`
	The array to store the imag. part of the eigenvalues

* `:::freefem vector=`
	The FE function array to store the eigenvectors

* `:::freefem rawvector=`
	An array of type `:::freefem real[int,int]` to store eigenvectors by column. (up to version 2-17).$\codered$

 	For real non symmetric problems, complex eigenvectors are given as two consecutive vectors, so if eigenvalue $k$ and $k+1$ are complex conjugate eigenvalues, the $k$th vector will contain the real part and the $k+1$th vector the imaginary part of the corresponding complex conjugate eigenvectors.

* `:::freefem tol=`
	The relative accuracy to which eigenvalues are to be determined;

* `:::freefem sigma=`
	The shift value;

* `:::freefem maxit=`
	The maximum number of iterations allowed;

* `:::freefem ncv=`
	The number of Arnoldi vectors generated at each iteration of ARPACK;

* `:::freefem mode=`
	The computational mode used by ARPACK (see above);

* `:::freefem which=`
The requested subset of eigenvalues (see above).

 __Example 9.17__ lapEignenValue.edp

In the first example, we compute the eigenvalues and the eigenvectors of the Dirichlet problem on square $\Omega=]0,\pi[^2$.

The problem is to find: $\lambda$, and $\nabla u_{\lambda}$ in $\mathbb{R}{\times} H^1_0(\Omega)$
$$ \int_\Omega \nabla u_{\lambda} \nabla v = \lambda \int_\Omega u v \quad \forall v \in H^1_0(\Omega)$$

The exact eigenvalues are $\lambda_{n,m} =(n^2+m^2), (n,m)\in {\mathbb{N}_*}^2$ with
the associated eigenvectors are $u_{{m,n}}=sin(nx)*sin(my)$.

We use the generalized inverse shift mode of the `:::freefem arpack++` library, to find
20 eigenvalues and eigenvectors close to the shift value $\sigma=20$.

```freefem
// Computation of the eigen value and eigen vector of the
// Dirichlet problem on square $]0,\pi[^2$
// ----------------------------------------
// we use the inverse shift mode
// the shift is given with the real sigma
// -------------------------------------
// find $\lambda$ and $u_\lambda\in H^1_0(\Omega)$ such that:
// $\displaystyle  \int_{\Omega}  \nabla u_{\lambda} \nabla v = \lambda \int_{\Omega} u_{\lambda}   v , \forall v \in H^1_0(\Omega) $
verbosity=10;
mesh Th=square(20,20,[pi*x,pi*y]);
fespace Vh(Th,P2);
Vh u1,u2;

real sigma = 20; // value of the shift

// OP = A - sigma B ; // the shifted matrix
varf op(u1,u2)= int2d(Th)(  dx(u1)*dx(u2) + dy(u1)*dy(u2) - sigma* u1*u2 )
                    +  on(1,2,3,4,u1=0) ; // Boundary condition

varf b([u1],[u2]) = int2d(Th)(  u1*u2 );//no Boundary condition see note \ref{note BC EV} $\codered$
matrix OP= op(Vh,Vh,solver=Crout,factorize=1); // crout solver because the matrix in not positive
matrix B= b(Vh,Vh,solver=CG,eps=1e-20);

// important remark:
// the boundary condition is make with exact penalization:
// we put 1e30=tgv on the diagonal term of the lock degree of freedom.
// So take Dirichlet boundary condition just on $a$ variational form
// and not on  $b$ variational form.
// because we solve $ w=OP^-1*B*v $

int nev=20; // number of computed eigen value close to sigma

real[int] ev(nev); // to store the nev eigenvalue
Vh[int] eV(nev); // to store the nev eigenvector

int k=EigenValue(OP,B,sym=true,sigma=sigma,value=ev,vector=eV,
                   tol=1e-10,maxit=0,ncv=0);

// tol= the tolerance
// maxit= the maximum iteration see arpack doc.
// ncv see arpack doc. \url{http://www.caam.rice.edu/software/ARPACK/} $\codered$
// the return value is number of converged eigen value.

for (int i=0;i<k;i++)
{
  u1=eV[i];
  real gg = int2d(Th)(dx(u1)*dx(u1) + dy(u1)*dy(u1));
  real mm= int2d(Th)(u1*u1) ;
  cout << " ---- " <<  i<< " " << ev[i]<< " err= "
       <<int2d(Th)(dx(u1)*dx(u1) + dy(u1)*dy(u1) - (ev[i])*u1*u1) << " --- "<<endl;
  plot(eV[i],cmm="Eigen Vector "+i+" valeur =" + ev[i]  ,wait=1,value=1);
}
```

The output of this example is:

```freefem
Nb of edges on Mortars  = 0
Nb of edges on Boundary = 80, neb = 80
Nb Of Nodes = 1681
Nb of DF = 1681
Real symmetric eigenvalue problem: A*x - B*x*lambda

Thanks to ARPACK++ class ARrcSymGenEig
Real symmetric eigenvalue problem: A*x - B*x*lambda
Shift and invert mode sigma=20

Dimension of the system            : 1681
Number of 'requested' eigenvalues  : 20
Number of 'converged' eigenvalues  : 20
Number of Arnoldi vectors generated: 41
Number of iterations taken         : 2

Eigenvalues:
  lambda[1]: 5.0002
  lambda[2]: 8.00074
  lambda[3]: 10.0011
  lambda[4]: 10.0011
  lambda[5]: 13.002
  lambda[6]: 13.0039
  lambda[7]: 17.0046
  lambda[8]: 17.0048
  lambda[9]: 18.0083
  lambda[10]: 20.0096
  lambda[11]: 20.0096
  lambda[12]: 25.014
  lambda[13]: 25.0283
  lambda[14]: 26.0159
  lambda[15]: 26.0159
  lambda[16]: 29.0258
  lambda[17]: 29.0273
  lambda[18]: 32.0449
  lambda[19]: 34.049
  lambda[20]: 34.0492

 ---- 0 5.0002 err= -0.000225891 ---
 ---- 1 8.00074 err= -0.000787446 ---
 ---- 2 10.0011 err= -0.00134596 ---
 ---- 3 10.0011 err= -0.00134619 ---
 ---- 4 13.002 err= -0.00227747 ---
 ---- 5 13.0039 err= -0.004179 ---
 ---- 6 17.0046 err= -0.00623649 ---
 ---- 7 17.0048 err= -0.00639952 ---
 ---- 8 18.0083 err= -0.00862954 ---
 ---- 9 20.0096 err= -0.0110483 ---
 ---- 10 20.0096 err= -0.0110696 ---
 ---- 11 25.014 err= -0.0154412 ---
 ---- 12 25.0283 err= -0.0291014 ---
 ---- 13 26.0159 err= -0.0218532 ---
 ---- 14 26.0159 err= -0.0218544 ---
 ---- 15 29.0258 err= -0.0311961 ---
 ---- 16 29.0273 err= -0.0326472 ---
 ---- 17 32.0449 err= -0.0457328 ---
 ---- 18 34.049 err= -0.0530978 ---
 ---- 19 34.0492 err= -0.0536275 ---
```

|Fig. 9.17: Isovalue of 11th eigenvector $u_{4,3}-u_{3,4}$|
|:----:|
|![eigen11](images/eigen11.png)|

|Fig. 9.18: Isovalue of 12th eigenvector $u_{4,3}+u_{3,4}$|
|:----:|
|![eigen12](images/eigen12.png)|
