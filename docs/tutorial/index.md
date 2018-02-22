# Getting started

Let's see how FreeFem++ solves **Poisson**’s equation:  
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

As illustrated in Fig. \ref{firstU},
we can see the isovalue of $u$ by using \ttCC{@plot} (see line 13
below).

Figure 2.1: mesh Th by `build(C(50))` |  Figure 2.2: isovalue by `plot(u)`
:-------------------------:|:-------------------------:
![arccos function](images/firstTh.svg)  |  ![arccos function](images/firstU.svg)
