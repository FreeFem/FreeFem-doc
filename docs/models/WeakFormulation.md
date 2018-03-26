In general, to obtain the integral or "weak" statements equivalent to system (\ref{eq:system} $\codered$) and boundary conditions we form a scalar dot product between an arbitrary magnetic field test function $\vec{H}^t=\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}$ and the components of our vectorial equation $A_1,A_2,A_3$, and integrate over the resonator's cross section domain $\Omega$ (and its boundary for the boundary conditions):

\begin{equation}
	\int\limits_{\Omega}(H^t_{\rho}A_1+H^t_{\phi}A_2+H^t_{z}A_3)d\Omega
\end{equation}

We can reduce the order of partial derivatives in this integral by using the Green's formula for integration by parts. For example:

\begin{equation}
\int\limits_{\Omega}H_z^t \frac{\partial^2 H_z}{\partial \rho^2 }d\Omega=
-\int\limits_{\Omega}\frac{\partial H_z^t}{\partial \rho}\frac{\partial H_z}{\partial \rho }d\Omega+\oint H_z^t\frac{\partial H_z}{\partial \rho}n_{\rho}d\Gamma
\end{equation}

Thus converting equations (\ref{eq:system} $\codered$) we obtain a large expression for the weak form (see [1])

### A dielectric sphere example with FreeFem++

We now compute the fundamental mode frequency for a fused silica sphere. The sphere is 36 micrometer in diameter, the refractive index is 1.46, the boundary condition is the magnetic wall (which can actually be omitted as it holds automatically). The example can be found in `:::freefem examples++-eigenvalue/WGM-sphere.edp` $\codered$ in the distribution.
