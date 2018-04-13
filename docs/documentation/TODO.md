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
	<div class="progress-bar" style="width:99%">
	</div>
	<span class="progress-label">99</span>
</div>

- [ ] line 262
	``` * `:::freefem paramFile=` This `:::freefem string` type parameter allows the user to pass all the parameters using an extern file as in Hansen's original code. More parameters related to the CMA-ES algorithm can be changed with this file. A sample of it can be found in the `:::freefem examples++-load/ffCMAES/` folder under the name `:::freefem initials.par` . Note that the parameters passed to the CMAES function in the __`FreeFem++`__ script will be ignored if an input parameters file is given.```
- [ ] line 332
	```Where if $a$ is a vector of $\R^n$, $A$ denotes the diagonal matrix $A=(a_i \delta_{ij})_{1\leq i,j\leq n}$ and $e\in\R^{n} = (1,1,\dots,1)$. Solving this nonlinear system by the Newton methods is known as being the _primal-dual_ interior points method. Here again, more details are available in [FORSGREN2002](#FORSGREN2002). Most actual implementations introduce features in order to globalize the convergence capability of the method, essentially by adding some line-search steps to the Newton algorithm, or by using trust regions. For the purpose of IPOPT, this is achieved by a _filter line search_ methods, the details of which can be found in [](#) .```

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

