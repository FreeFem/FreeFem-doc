
## Whispering gallery modes

Written by I. S. Grudinin.

In whispering gallery mode (WGM) resonators, which are typically spheres or disks, electromagnetic field is trapped by total internal reflections from the boundary. Modes of such resonators are distinguished by compact volume and record high quality factors (Q) in a broad range of frequencies. Modern applications of such resonators include microwave and optical cavities for atomic clocks, cavity optomechanics, nonlinear and quantum optics. Analytical solutions for WG modes are only available for a limited number of idealized geometries, such as sphere or ellipsoid. Since resonator dimensions are typically much larger than optical wavelength, direct application of numerical 3D finite difference time domain (FDTD) or finite element methods (FEM) is not practical. It's possible to solve the vectorial wave equation by reducing it to a two dimensional case by taking axial symmetry into account. Such reduction leads to a system of 3 equations to be solved in a 2D "$\rho-z$" section of a resonator. Please refer to \cite{Oxborrow} $\codered$ for a detailed derivation and to \cite{Grudinin} $\codered$ for an example of using FreeFem++ to compute WGMs.

### Wave equation for the WGMs

Since electric field is discontinuous on the surface of a dielectric and magnetic field is typically not, we derive our equations for the magnetic field. The electric field can be easily derived at a later stage from $\vec{E}=\frac{i}{\omega\epsilon_0}\hat{\epsilon}^{-1}\nabla\times\vec{H}$. Following a standard procedure starting with Maxwell equations we derive a wave equation in a single-axis anisotropic medium such as an optical crystal:

\begin{equation}
\nabla\times\left(\hat{\epsilon}^{-1}\nabla\times\vec{H}\right)-k_0^2\vec{H}-\alpha\nabla\left(\nabla\cdot\vec{H}\right)=0
\end{equation}

Here $k_0=\omega/c$ is the wavenumber, $\alpha$ is the penalty term added to fight spurious FEM solutions.  For anisotropic single-axis medium with $\partial\hat{\epsilon}/\partial\phi=0$ in cylindrical system of coordinates we have:

\begin{equation}
\hat{\epsilon}=\begin{pmatrix} \epsilon_{\rho} & 0 & 0 \\ 0 & \epsilon_{\rho} & 0 \\ 0 & 0 & \epsilon_z \end{pmatrix}. \nonumber
\end{equation}

We now assume axial symmetry of our electromagnetic fields and insert an imaginary unity in front of the $H_{\phi}$ to allow all field components to be real numbers and also to account for the phase shift of this component $\vec{H}(\rho,\phi,z)=\left\{H_{\rho}(\rho,z),iH_{\phi}(\rho,z),H_z(\rho,z)\right\}\times e^{im\phi}$.

We write the wave equation (\ref{eq:wave} $\codered$) explicitly in cylindrical coordinates, thus obtaining a set of three differential equations for the domain $\Omega$ given by the resonator's cross section and some space outside:

\begin{eqnarray}
A_1\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}&=&0\\ \nonumber
A_2\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}&=&0\\ \nonumber
A_3\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}&=&0
\end{eqnarray}

 The numerical solutions of these equations and boundary conditions can be found with FreeFem++ if we write the system in the weak, or integral form.
