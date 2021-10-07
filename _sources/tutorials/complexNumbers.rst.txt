An Example with Complex Numbers
===============================

In a microwave oven heat comes from molecular excitation by an electromagnetic field.
For a plane monochromatic wave, amplitude is given by Helmholtz’s equation:

.. math::
    \beta v + \Delta v = 0.

We consider a rectangular oven where the wave is emitted by part of the upper wall.
So the boundary of the domain is made up of a part :math:`\Gamma_1` where :math:`v=0` and of another part :math:`\Gamma_2=[c,d]` where for instance :math:`\displaystyle v=\sin\left(\pi{y-c\over c-d}\right)`.

Within an object to be cooked, denoted by :math:`B`, the heat source is proportional to :math:`v^2`.
At equilibrium, one has :

.. math::
    \begin{array}{rcl}
        -\Delta\theta &=& v^2 I_B\\
        \theta_\Gamma &=& 0
    \end{array}

where :math:`I_B` is :math:`1` in the object and :math:`0` elsewhere.

In the program below :math:`\beta = 1/(1-i/2)` in the air and :math:`2/(1-i/2)` in the object (:math:`i=\sqrt{-1}`):

.. code-block:: freefem
   :linenos:

   // Parameters
   int nn = 2;
   real a = 20.;
   real b = 20.;
   real c = 15.;
   real d = 8.;
   real e = 2.;
   real l = 12.;
   real f = 2.;
   real g = 2.;

   // Mesh
   border a0(t=0, 1){x=a*t; y=0; label=1;}
   border a1(t=1, 2){x=a; y=b*(t-1); label=1;}
   border a2(t=2, 3){ x=a*(3-t); y=b; label=1;}
   border a3(t=3, 4){x=0; y=b-(b-c)*(t-3); label=1;}
   border a4(t=4, 5){x=0; y=c-(c-d)*(t-4); label=2;}
   border a5(t=5, 6){x=0; y=d*(6-t); label=1;}

   border b0(t=0, 1){x=a-f+e*(t-1); y=g; label=3;}
   border b1(t=1, 4){x=a-f; y=g+l*(t-1)/3; label=3;}
   border b2(t=4, 5){x=a-f-e*(t-4); y=l+g; label=3;}
   border b3(t=5, 8){x=a-e-f; y=l+g-l*(t-5)/3; label=3;}

   mesh Th = buildmesh(a0(10*nn) + a1(10*nn) + a2(10*nn) + a3(10*nn) +a4(10*nn) + a5(10*nn)
      + b0(5*nn) + b1(10*nn) + b2(5*nn) + b3(10*nn));
   real meat = Th(a-f-e/2, g+l/2).region;
   real air= Th(0.01,0.01).region;
   plot(Th, wait=1);

   // Fespace
   fespace Vh(Th, P1);
   Vh R=(region-air)/(meat-air);
   Vh<complex> v, w;
   Vh vr, vi;

   fespace Uh(Th, P1);
   Uh u, uu, ff;

   // Problem
   solve muwave(v, w)
      = int2d(Th)(
           v*w*(1+R)
         - (dx(v)*dx(w) + dy(v)*dy(w))*(1 - 0.5i)
      )
      + on(1, v=0)
      + on(2, v=sin(pi*(y-c)/(c-d)))
      ;

   vr = real(v);
   vi = imag(v);

   // Plot
   plot(vr, wait=1, ps="rmuonde.ps", fill=true);
   plot(vi, wait=1, ps="imuonde.ps", fill=true);

   // Problem (temperature)
   ff=1e5*(vr^2 + vi^2)*R;

   solve temperature(u, uu)
      = int2d(Th)(
           dx(u)* dx(uu)+ dy(u)* dy(uu)
      )
      - int2d(Th)(
           ff*uu
      )
      + on(1, 2, u=0)
      ;

   // Plot
   plot(u, wait=1, ps="tempmuonde.ps", fill=true);

Results are shown on :numref:`figComplexReal`, :numref:`figComplexImaginary` and :numref:`figComplexTemperature`.

.. subfigstart::

.. _figComplexReal:

.. figure:: images/real_microwave.png
   :alt: RealMicroWave
   :width: 90%

   Real part

.. _figComplexImaginary:

.. figure:: images/imaginary_microwave.png
   :alt: ImaginaryMicrowave
   :width: 90%

   Imaginary part

.. _figComplexTemperature:

.. figure:: images/temperature_microwave.png
   :alt: TemperatureMicrowave
   :width: 90%

   Temperature

.. subfigend::
   :width: 0.49
   :alt: Microwave
   :label: Microwave

   Microwave
