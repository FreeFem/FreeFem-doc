# Getting started

## Solving Poisson’s equation  
_For a given function $f(x,y)$, find a function $u(x,y)$ satisfying_

$$\begin{eqnarray}
\label{eqn:Poisson}
-\Delta u(x,y) &=& f(x,y)\quad \mbox{ for all }(x,y)\in\Omega,
 \\ \label{eqn:Dirichlet}
  u(x,y) &=& 0\quad \mbox{ for all }(x,y)\mbox{ on }\p\Omega,.
\end{eqnarray}$$

Here $\p\Omega$ is the boundary of the bounded open set $\Omega\subset \R^2$
and  $\Delta u = \frac{\p^2 u}{\p x^2 } + \frac{\p^2 u}{\p y^2}$.

The following is a Freefem++ program which computes $u$ when
$f(x,y)=xy$  and $\Omega$ is the unit disk. The boundary
$C=\partial\Omega$ is
$$
C=\{(x,y)|\; x=\cos(t),\, y=\sin(t),\, 0\le t\le 2\pi\}
$$

Note that in FreeFem++ the domain $\Omega$ is assumed to described by its boundary
that is on the left side of its boundary oriented by the parameter.

As illustrated in Fig. 2.2,
we can see the isovalue of $u$ by using freefem++'s `plot` command (see line 13
below).

Figure 2.1: mesh Th by `build(C(50))` |  Figure 2.2: isovalue by `plot(u)`
:-------------------------:|:-------------------------:
![mesh TH](images/firstTh.svg)  |  ![isovalue](images/firstU.svg)

## Example 1

```freefem
// Define mesh boundary
border C(t=0, 2*pi){x=cos(t); y=sin(t);}

// The triangulated domain Th is on the left side of its boundary
mesh Th = buildmesh(C(50));

// The finite element space defined over Th is called here Vh
fespace Vh(Th, P1);
Vh u, v;// Define u and v as piecewise-P1 continuous functions

// Define a function f
func f= x*y;

// Get the clock in second
real cpu=clock();

// Define the PDE
solve Poisson(u, v, solver=LU)
	= int2d(Th)(	// The bilinear part
		  dx(u)*dx(v)
		+ dy(u)*dy(v)
	)
	- int2d(Th)(	// The right hand side
		  f*v
	)
	+ on(C, u=0);	// The Dirichlet boundary condition

// Plot the result
plot(u);

// Display the total computational time
cout << "CPU time = " << (clock()-cpu) << endl;
```

Note that the qualifier `solver=LU` is not required and by default a
multi-frontal `LU` would have been used. Note also that the lines
containing `clock` are equally not required. Finally note how
close to the mathematics FreeFem++ input language is. Lines 19 to 24
correspond to the mathematical variational equation
\[
    \int_{T_h}(\frac{\p u}{\p x}\frac{\p v}{\p x}
    +\frac{\p u}{\p y}\frac{\p v}{\p
    y})\d x \d y
    =
   \int_{T_h}f v\d x\d y
\]
for all $v$ which are in the finite element space $V_h$ and zero on
the boundary $C$.


Exercise:
Change `P1` into `P2` and run the program.

