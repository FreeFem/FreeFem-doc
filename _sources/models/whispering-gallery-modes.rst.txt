Whispering gallery modes
========================

Author: `I. S. Grudinin <http://linkeding.com/in/grudinin>`__

In whispering gallery mode (WGM) resonators, which are typically spheres or disks, electromagnetic field is trapped by total internal reflections from the boundary.
Modes of such resonators are distinguished by compact volume and record high quality factors (Q) in a broad range of frequencies.

Modern applications of such resonators include microwave and optical cavities for atomic clocks, cavity optomechanics, nonlinear and quantum optics.
Analytical solutions for WG modes are only available for a limited number of idealized geometries, such as sphere or ellipsoid.
Since resonator dimensions are typically much larger than optical wavelength, direct application of numerical 3D finite difference time domain (FDTD) or finite element methods (FEM) is not practical.
It’s possible to solve the vectorial wave equation by reducing it to a two dimensional case by taking axial symmetry into account.

Such reduction leads to a system of 3 equations to be solved in a 2D “:math:`\rho-z`" section of a resonator.
Please refer to [OXBORROW2007]_ for a detailed derivation and to [GRUDININ2012]_ for an example of using **FreeFEM** to compute WGMs.

Wave equation for the WGMs
--------------------------

Since electric field is discontinuous on the surface of a dielectric and magnetic field is typically not, we derive our equations for the magnetic field.
The electric field can be easily derived at a later stage from :math:`\vec{E}=\frac{i}{\omega\epsilon_0}\hat{\epsilon}^{-1}\nabla\times\vec{H}`.
Following a standard procedure starting with Maxwell equations we derive a wave equation in a single-axis anisotropic medium such as an optical crystal:

.. math::
   \nabla\times\left(\hat{\epsilon}^{-1}\nabla\times\vec{H}\right)-k_0^2\vec{H}-\alpha\nabla\left(\nabla\cdot\vec{H}\right)=0
   :label: eqn::wave

Here :math:`k_0=\omega/c` is the wavenumber, :math:`\alpha` is the penalty term added to fight spurious FEM solutions.
For anisotropic single-axis medium with :math:`\partial\hat{\epsilon}/\partial\phi=0` in cylindrical system of coordinates we have:

.. math::
   \hat{\epsilon}=\begin{pmatrix} \epsilon_{\rho} & 0 & 0 \\ 0 & \epsilon_{\rho} & 0 \\ 0 & 0 & \epsilon_z \end{pmatrix}. \nonumber

We now assume axial symmetry of our electromagnetic fields and insert an imaginary unity in front of the :math:`H_{\phi}` to allow all field components to be real numbers and also to account for the phase shift of this component :math:`\vec{H}(\rho,\phi,z)=\left\{H_{\rho}(\rho,z),iH_{\phi}(\rho,z),H_z(\rho,z)\right\}\times e^{im\phi}`.

We write the wave equation :eq:`eqn::wave` explicitly in cylindrical coordinates, thus obtaining a set of three differential equations for the domain :math:`\Omega` given by the resonator’s cross section and some space outside:

.. math::
    \begin{array}{rcl}
        A_1\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}&=&0\\ \nonumber
        A_2\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}&=&0\\ \nonumber
        A_3\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}&=&0
    \end{array}
    :label: eqn::system

The numerical solutions of these equations and boundary conditions can be found with **FreeFEM** if we write the system in the weak, or integral form.

Weak formulation
----------------

In general, to obtain the integral or “weak" statements equivalent to system :eq:`eqn::system` and boundary conditions we form a scalar dot product between an arbitrary magnetic field test function :math:`\mathbf{H}^t=\{{H}_{\rho}^t,{H}_{\phi}^t,{H}_{z}^t\}` and the components of our vectorial equation :math:`A_1,A_2,A_3`, and integrate over the resonator’s cross section domain :math:`\Omega` (and its boundary for the boundary conditions):

.. math::
   \int\limits_{\Omega}(H^t_{\rho}A_1+H^t_{\phi}A_2+H^t_{z}A_3)d\Omega

We can reduce the order of partial derivatives in this integral by using the Green’s formula for integration by parts.
For example:

.. math::
   \int\limits_{\Omega}H_z^t \frac{\partial^2 H_z}{\partial \rho^2 }d\Omega=
   -\int\limits_{\Omega}\frac{\partial H_z^t}{\partial \rho}\frac{\partial H_z}{\partial \rho }d\Omega+\oint H_z^t\frac{\partial H_z}{\partial \rho}n_{\rho}d\Gamma

Thus converting equations :eq:`eqn::system` we obtain a large expression for the weak form.

A dielectric sphere example with FreeFEM
----------------------------------------

We now compute the fundamental mode frequency for a fused silica sphere.
The sphere is 36 micrometer in diameter, the refractive index is 1.46, the boundary condition is the magnetic wall (which can actually be omitted as it holds automatically).

.. code-block:: freefem
   :linenos:

   // Parameters
   real radius = 36; //approximate radius of the cavity
   real yb = -10, yt = -yb; //window yb=bottom and yt=top coordinates
   real xl = radius-5, xr = radius+3; //window xl=left and xr=right coordinates
   real angle = asin((yt)/radius); //angle of the sphere segment to model in radians
   int Nm = 60; //number of mesh vertices per border
   real ne = 1.46; //n_e-extraordinary refractive index (root of permittivity parallel to z-axis, epara)
   real no = 1.46; //n_o-ordinary refractive index (root of permittivity orthogonal to z-axis, eorto)
   real nm = 1; //refractive index of surrounding medium (air)

   int nev = 4; // number of eigen values to find

   int M = 213; //azimuthal mode order ~ 2Pi*n*R/lambda
   real alpha = 1; //penalty term

   // Mesh
   border W1l(t=0, 1){x=xl+(radius*cos(angle)-xl)*(1-t); y=yt; label=1;}
   border W1r(t=0, 1){x=xr-(xr-radius*cos(angle))*(t); y=yt; label=1;}
   border W2(t=0, 1){x=xr; y=yb+(yt-yb)*t; label=1;}
   border W3l(t=0, 1){x=xl+(radius*cos(angle)-xl)*(t); y=yb; label=1;}
   border W3r(t=0, 1){x=xr-(xr-radius*cos(angle))*(1-t); y=yb; label=1;}
   border W4(t=0, 1){x=xl; y=yt-(yt-yb)*t; label=1;}
   border S(t=0, 1){x=radius*cos((t-0.5)*2*angle); y=radius*sin((t-0.5)*2*angle); label=2;}
   mesh Th = buildmesh(W1r(Nm/4) + W1l(Nm/4) + W4(Nm) + W3l(Nm/4) + W3r(Nm/4) + W2(Nm) + S(Nm));
   plot(Th, WindowIndex=0);

   // Fespace
   fespace Ph(Th, P0);
   Ph reg = region;

   int ncav = reg(xl+1, 0); // cavity
   int nair = reg(xr-1, 0); //air
   Ph eorto = no^2*(region==ncav) + nm^2*(region==nair); //subdomains for epsilon values inside and outside the resonators
   Ph epara = ne^2*(region==ncav) + nm^2*(region==nair); //subdomains for epsilon values inside and outside the resonators

   //supplementary variables to store eigenvectors, defined on mesh Th with P2 elements - Largange quadratic.
   fespace Supp(Th, P2);
   Supp eHsqr;

   //3d vector FE space
   fespace Vh(Th, [P2, P2, P2]);
   Vh [Hr, Hphi, Hz], [vHr, vHphi, vHz]; //magnetic field components on Vh space and test functions vH

   // Macro
   //boundary condition macros
   macro EWall(Hr, Hphi, Hz) (
         dy(Hr) - dx(Hz) + Hr*N.x + Hz*N.y
       - epara*(Hz*M - dy(Hphi)*x)*N.y
       + eorto*(Hphi - Hr*M+dx(Hphi)*x)*N.x) //
   macro MWall(Hr, Hphi, Hz) (
         Hphi + Hz*N.x - Hr*N.y
       + epara*(Hz*M - dy(Hphi)*x)*N.x
       + eorto*(Hphi - Hr*M+dx(Hphi)*x)*N.y ) //

   // Problem
   real sigma =(M/(ne*radius))^2+2; // value of the shift (k^2), where the modes will be found
   varf b ([Hr, Hphi, Hz], [vHr, vHphi, vHz])
       = int2d(Th)(
             x*(Hr*vHr+Hphi*vHphi+Hz*vHz)
       )
       ;
   // OP = A - sigma B ; // the shifted matrix
   varf op ([Hr, Hphi, Hz], [vHr, vHphi, vHz])=
       int2d(Th)(
             (
                 (eorto*(vHphi*Hphi - M*(vHphi*Hr + Hphi*vHr) + M^2*vHr*Hr) + epara*M^2*vHz*Hz)/x //A/r
               + eorto*(dx(vHphi)*(Hphi - M*Hr) + dx(Hphi)*(vHphi - M*vHr)) - epara*M*(vHz*dy(Hphi) + Hz*dy(vHphi)) //B
               + x*(eorto*dx(vHphi)*dx(Hphi) + epara*((dx(vHz) - dy(vHr))*(dx(Hz) - dy(Hr)) + dy(vHphi)*dy(Hphi))) //C
           )/(eorto*epara)
           + alpha*(
                 (vHr*Hr - M*(vHphi*Hr + Hphi*vHr) + M^2*vHphi*Hphi)/x //D/r
               + (dx(vHr) + dy(vHz))*(Hr - M*Hphi) + (vHr - M*vHphi)*(dx(Hr) + dy(Hz)) //E
               + x*(dx(vHr) + dy(vHz))*(dx(Hr) + dy(Hz)) //F
           )
           -sigma*x*(vHr*Hr + vHphi*Hphi + vHz*Hz)
       )
       //electric wall boundary condition on the boundary of computation domain
       +int1d(Th, 1)(
             EWall(Hr, Hphi, Hz)*EWall(vHr, vHphi, vHz)
       )
       ;
   //setting sparce matrices and assigning the solver UMFPACK to solve eigenvalue problem
   matrix B = b(Vh, Vh, solver=UMFPACK);
   matrix OP = op(Vh, Vh, solver=UMFPACK);

   // Solve
   real[int] ev(nev); //to store the nev eigenvalue
   Vh[int] [eHr, eHphi, eHz](nev); //to store the nev eigenvector
   //calling ARPACK on sparce matrices with the assigned solver UMFPACK:
   int k = EigenValue(OP, B, sym=true, sigma=sigma, value=ev, vector=eHr, tol=1e-10, maxit=0, ncv=0);

   k = min(k, nev); //sometimes the number of converged eigen values
                    //can be greater than nev

   //file to output mode values
   ofstream f("modes.txt");
   //setting number of digits in the file output
   int nold = f.precision(11);

   // Plot & Save
   for (int i = 0; i < k; i++){
       real lambda = 2*pi/sqrt(ev[i]);
       eHsqr = (sqrt(eHr[i]^2 + eHphi[i]^2 + eHz[i]^2)); //intensity from magnetic field components
       plot(eHsqr, WindowIndex=i, value=1, nbiso=20,LabelColors=1, aspectratio=1, cmm="Mode "+i+", lambda="+lambda+", F="+(299792.458/lambda));
       f << "Mode "<< i << ", ka=" << sqrt(ev[i])*radius << endl;
   }
