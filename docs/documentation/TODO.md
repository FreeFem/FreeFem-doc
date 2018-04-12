<!--- THIS FILE IS AUTOMATICALY GENERATED --->
<!--- DO NOT EDIT --->

# TODO

## Home

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Notations

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Mesh generation

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:99%">
	</div>
	<span class="progress-label">99</span>
</div>

- [ ] line 2597
	``` mmg3d-v4.0```

## Finite element

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Visualization

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Algorithms & Optimization

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:97%">
	</div>
	<span class="progress-label">97</span>
</div>

- [ ] line 56
	```* `:::freefem stop=` `:::freefem stopfunc`  add your test function to stop before (after version 3.31) . The prototype for the function `:::freefem stopfunc` is```
- [ ] line 157
	```then $J(\vec{x})$ is minimized by the solution $\vec{x}$ of $A\vec{x}=\vec{b}$. In this case, we can use the function `:::freefem LinearCG`  LINEAR OR AFFINE ?? (FRANCK)```
- [ ] line 176
	```For detail of these algorithms, refer to \cite{Lucquin} [Chapter IV, 1.3].```
- [ ] line 180
	```Two algorithms of COOOL a package \cite{coool}  are interfaced with the Newton Raphson method (call `:::freefem Newton`) and the `:::freefem BFGS` method. These two ones are directly available in FreeFem (no dynamical link to load). Be careful with these algorithms, because their implementation uses full matrices. We also provide several optimization algorithms from the NLopt library \cite{nlopt}  as well as an interface for Hansen's implementation of CMAES (a MPI version of this one is also available). These last algorithms can be found as dynamical links in the `:::freefem example++-load` folder as the `:::freefem ff-NLopt` and `:::freefem CMA_ES` files (`:::freefem CMA_ES_MPI` from the `:::freefem example++-mpi` folder for the mpi version).```
- [ ] line 211
	```%\end{example} ```
- [ ] line 222
	```This algorithm works with a normal multivariate distribution in the parameters space and try to adapt its covariance matrix using the information provides by the successive function evaluations (see \cite{hansen}  for more details). Thus, some specific parameters can be passed to control the starting distribution, size of the sample generations etc... Named parameters for this are the following :```
- [ ] line 243
	``` * `:::freefem popsize=` Integer value used to change the sample size. The default value is $4+ \lfloor 3\ln (n) \rfloor$, see \cite{hansen}  for more details. Increasing the population size usually improves the global search capabilities at the cost of an at most linear reduction of the convergence speed with respect to `:::freefem popsize`.```
- [ ] line 249
	```The  `:::freefem ff-Ipopt` package is an interface for the IPOPT \cite{ipopt}  optimizer. IPOPT is a software library for large scale, non-linear, constrained optimization. Detailed informations about it are in \cite{ipopt}  and [https://projects.coin-or.org/Ipopt](https://projects.coin-or.org/Ipopt). It implements a primal-dual interior point method along with filter method based line searchs.```
- [ ] line 254
	```In this section, we give a very brief glimpse at the underlying mathematics of IPOPT. For a deeper introduction on interior methods for nonlinear smooth optimization, one may consults \cite{ipintro} , or \cite{ipopt}  for more IPOPT specific elements. IPOPT is designed to perform optimization for both equality and inequality constrained problems. Though, nonlinear inequalities are rearranged before the beginning of the optimization process in order to restrict the panel of nonlinear constraints to those of the equality kind. Each nonlinear inequality ones are transformed into a pair of simple bound inequality and nonlinear equality constraint by the introduction of as many slack variables as is needed : $c_{i}(x)\leq 0$ becomes $c_{i}(x) + s_{i} = 0$ and $s_{i}\leq 0$, where $s_{i}$ is added to the initial variables of the problems $x_{i}$. Thus, for convenience, we will assume that the minimization problem does not contain any nonlinear inequality constraint. It means that, given a function $f:\mathbb{R}^{n}\mapsto\mathbb{R}$, we want to find :```
- [ ] line 266
	```As a barrier method, interior points algorithms try to find a Karush-Kuhn-Tucker point for (\ref{minimproblem})  by solving a sequence of problems, unconstrained with respect to the inequality constraints, of the form :```
- [ ] line 274
	```The remaining equality constraints are handled with the usual Lagrange multipliers method. If the sequence of barrier parameters $\mu$ converge to 0, intuition suggests that the sequence of minimizers of (\ref{barrier})  converge to a local constrained minimizer of (\ref{minimproblem}) . For a given $\mu$, (\ref{barrier})  is solved by finding $(x_{\mu},\lambda_{\mu})\in\R^{n}\times\R^{m}$ such that :```
- [ ] line 294
	```In this equation, the $z_l$ and $z_u$ vectors seems to play the role of Lagrange multipliers for the simple bounds inequalities, and indeed, when $\mu\rightarrow 0$, they converge toward some suitable Lagrange multipliers for the KKT conditions, provided some technical assumptions are fulfilled (see \cite{ipintro} ).```
- [ ] line 296
	```Equation \ref{muproblemlambda} is solved by performing a Newton method in order to find a solution of (\ref{muproblem}) for each of the decreasing values of $\mu$. Some order 2 conditions are also taken into account to avoid convergence to local maximizer, see \cite{ipintro}  for precision about them. In the most classical IP algorithms, the Newton method is directly applied to (\ref{muproblem}). This is in most case inefficient due to frequent computation of infeasible points. These difficulties are avoided in Primal-Dual interior points methods where (\ref{muproblem}) is transformed into an extended system where $z_u$ and $z_l$ are treated as unknowns and the barrier problems are finding $(x,\lambda,z_u,z_l)\in\R^n\times\R^m\times\R^n\times\R^n$ such that :```
- [ ] line 307
	```Where if $a$ is a vector of $\R^n$, $A$ denotes the diagonal matrix $A=(a_i \delta_{ij})_{1\leq i,j\leq n}$ and $e\in\R^{n} = (1,1,\dots,1)$. Solving this nonlinear system by the Newton methods is known as being the _primal-dual_ interior points method. Here again, more details are available in \cite{ipintro} . Most actual implementations introduce features in order to globalize the convergence capability of the method, essentially by adding some line-search steps to the Newton algorithm, or by using trust regions. For the purpose of IPOPT, this is achieved by a _filter line search_ methods, the details of which can be found in \cite{iplinesearch} .```
- [ ] line 309
	```More IPOPT specific features or implementation details can be found in \cite{ipopt} . We will just retain that IPOPT is a smart Newton method for solving constrained optimization problem, with global convergence capabilities due to a robust line search method (in the sense that the algorithm will convergence no matter the initializer). Due to the underlying Newton method, the optimization process requires expressions of all derivatives up to the order 2 of the fitness function as well as those of the constraints. For problems whose hessian matrices are difficult to compute or lead to high dimensional dense matrices, it is possible to use a BFGS approximation of these objects at the cost of a much slower convergence rate.```
- [ ] line 350
	``` PAS SUR DE LA OU SE TROUVE LA FIN DU 2EME WARNING```
- [ ] line 358
	```where $\lambda\in\R^{m}$ and $\sigma\in\R$  FIX EQUATION WHERE=... (FRANCK). Your hessian function should then have the following prototype :```
- [ ] line 776
	```Here are the functions related to the area computation and its shape derivative, according to equations \ref{msarea} and \ref{msdarea}  :```
- [ ] line 812
	```The function returning the hessian of the area for a given shape is a bit blurry, thus we won't show here all of equation \ref{msd2area}  coefficients definition, they can be found in the `:::freefem edp` file.```
- [ ] line 961
	```The `:::freefem ff-NLopt` package provides a FreeFem interface to the free/open-source library for nonlinear optimization, thus easing the use of several different free optimization (constrained or not) routines available online along with the PDE solver. All the algorithms are well documented in \cite{nlopt} , thus no exhaustive informations concerning their mathematical specificities will be found here and we will focus on the way they are called in a FreeFem script. One needing detailed informations about these algorithms should visit the said cite where a description of each of them is given, as well as many bibliographical links.```
- [ ] line 996
	```* `:::freefem IConst`/`:::freefem EConst` : Allows to pass the name of a function implementing some inequality (resp. equality) constraints on the search space. The function type must be `:::freefem real[int]` $\rightarrow$ `:::freefem real[int]` where the size of the returned array is equal to the number of constraints (of the same type - it means that all the constraints are computed in one vectorial function). In order to mix inequality and equality constraints in a same minimization attempt, two vectorial functions have to be defined and passed. See example \ref{varineqex}  for more details about how these constraints have to be implemented.```
- [ ] line 1028
	```The following table sums up the main characteristics of each algorithm, providing the more important information about which features are supported by which algorithm and what are the unavoidable arguments they need. More details can be found in \cite{nlopt} .```
- [ ] line 1034
	``````

## Parallelization

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:99%">
	</div>
	<span class="progress-label">99</span>
</div>

- [ ] line 54
	``` script bug```
- [ ] line 132
	``` script bug```
- [ ] line 154
	``` mpiReduceScatter is commented in parallelempi.cpp```
- [ ] line 156
	```See the `:::freefem examples++-mpi/essai.edp`  to test of all this functionality and thank you to Guy-Antoine Atenekeng Kahou for his help to code this interface.```
- [ ] line 237
	``` script bug```
- [ ] line 241
	``` check this part```

## Plugins

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Developers

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:98%">
	</div>
	<span class="progress-label">98</span>
</div>

- [ ] line 596
	``````
- [ ] line 663
	```It is agood idea to first try the example `load.edp` in directory `example++-load` .```
- [ ] line 671
	```Now, assume that you are in a shell window (a `cygwin` window under Windows) in the directory `example++-load` .```
- [ ] line 854
	```The associed edp file is `examples++-load/convect_dervieux.edp`. ```
- [ ] line 856
	```See `mat_dervieux.cpp`. ```
- [ ] line 860
	```First read the [Adding a new finite element section](#adding-a-new-finite-element), we add two new finite elements examples in the directory `examples++-load`. ```
- [ ] line 878
	```See `BernardiRaugel.cpp`. ```
- [ ] line 927
	```A real example using this finite element, just a small modification of the `NSP2P1.edp`  examples, just the begenning is change to```
- [ ] line 942
	``` PR needed: BernaRdiRaugel.cpp```
- [ ] line 945
	```See the example `bilapMorley.edp` .```
- [ ] line 950
	```Warning the sparse solver interface as been completely rewritten in version 3.2, so the section is obsolete, the example in are correct/ ```

