.. role:: freefem(code)
  :language: freefem

.. _visualization:

Visualization
=============

Results created by the finite element method can be a huge set of data, so it is very important to render them easy to grasp.

There are two ways of visualization in **FreeFEM**:

-  One, the default view, which supports the drawing of meshes, isovalues of real FE-functions, and of vector fields, all by the command :freefem:`plot` (see :ref:`Plot section <plot>` below).
   For publishing purpose, **FreeFEM** can store these plots as postscript files.

-  Another method is to use external tools, for example, gnuplot (see :ref:`Gnuplot section <gnuplot>`, :ref:`medit section <medit>`, :ref:`Paraview section <paraview>`, :ref:`Matlab/Octave section <matlab>`) using the command :freefem:`system` to launch them and/or to save the data in text files.

.. _plot:

Plot
----

With the command :freefem:`plot`, meshes, isovalues of scalar functions, and vector fields can be displayed.

The parameters of the plot command can be meshes, real FE functions, arrays of 2 real FE functions, arrays of two double arrays, to plot respectively a mesh, a function, a vector field, or a curve defined by the two double arrays.

.. note:: The length of an arrow is always bound to be in [5‰, 5%] of the screen size in order to see something.

The :freefem:`plot` command parameters are listed in the :ref:`Reference part <referencePlot>`.

The keyboard shortcuts are:

-  **enter** tries to show plot
-  **p** previous plot (10 plots saved)
-  **?** shows this help
-  **+,-** zooms in/out around the cursor 3/2 times
-  **=** resets the view
-  **r** refreshes plot
-  **up, down, left, right** special keys to tanslate
-  **3** switches 3d/2d plot keys :

   -  **z,Z** focal zoom and zoom out
   -  **H,h** increases or decreases the Z scale of the plot

-  **mouse motion**:

   -  **left button** rotates
   -  **right button** zooms (ctrl+button on mac)
   -  **right button +alt** tanslates (alt+ctrl+button on mac)

-  **a,A** increases or decreases the arrow size
-  **B** switches between showing the border meshes or not
-  **i,I** updates or not: the min/max bound of the functions to the window
-  **n,N** decreases or increases the number of iso value arrays
-  **b** switches between black and white or color plotting
-  **g** switches between grey or color plotting
-  **f** switches between filling iso or iso line
-  **l** switches between lighting or not
-  **v** switches between show or not showing the numerical value of colors
-  **m** switches between show or not showing the meshes
-  **w** window dump in file ffglutXXXX.ppm
-  **\*** keep/drop viewpoint for next plot
-  **k** complex data / change view type
-  **ESC** closes the graphics process before version 3.22, after no way to close
-  **otherwise** does nothing

For example:

.. code-block:: freefem
   :linenos:

   real[int] xx(10), yy(10);

   mesh Th = square(5,5);

   fespace Vh(Th, P1);

   //plot scalar and vectorial FE function
   Vh uh=x*x+y*y, vh=-y^2+x^2;
   plot(Th, uh, [uh, vh], value=true, ps="three.eps", wait=true);

   //zoom on box defined by the two corner points [0.1,0.2] and [0.5,0.6]
   plot(uh, [uh, vh], bb=[[0.1, 0.2], [0.5, 0.6]],
      wait=true, grey=true, fill=true, value=true, ps="threeg.eps");

   //compute a cut
   for (int i = 0; i < 10; i++){
      x = i/10.;
      y = i/10.;
      xx[i] = i;
      yy[i] = uh; //value of uh at point (i/10., i/10.)
   }
   plot([xx, yy], ps="likegnu.eps", wait=true);

.. subfigstart::

.. _figVisuMesh:

.. figure:: images/Visualization_Plot.png
   :alt: Visualization_Plot
   :width: 90%

   Mesh, isovalue and vector

.. _figVisuGrey:

.. figure:: images/Visualization_Plot_Grey.png
   :alt: Visualization_Plot_Grey
   :width: 90%

   Enlargement in grey of isovalue and vector

.. _figVisuCut:

.. figure:: images/Visualization_Plot_Gnuplot.png
   :alt: Visualization_Plot_Gnuplot
   :width: 90%

   Plots a cut of :freefem:`uh`. Note that a refinement of the same can be obtained in combination with gnuplot

.. subfigend::
   :width: 0.49
   :alt: Plot
   :label: Plot

   Plot

To change the color table and to choose the value of iso line you can do:

.. code-block:: freefem
   :linenos:

   // from: \url{http://en.wikipedia.org/wiki/HSV_color_space}
   // The HSV (Hue, Saturation, Value) model defines a color space
   // in terms of three constituent components:
   // HSV color space as a color wheel
   // Hue, the color type (such as red, blue, or yellow):
   // Ranges from 0-360 (but normalized to 0-100% in some applications, like here)
   // Saturation, the "vibrancy" of the color: Ranges from 0-100%
   // The lower the saturation of a color, the more "grayness" is present
   // and the more faded the color will appear.
   // Value, the brightness of the color: Ranges from 0-100%

   mesh Th = square(10, 10, [2*x-1, 2*y-1]);

   fespace Vh(Th, P1);
   Vh uh=2-x*x-y*y;

   real[int] colorhsv=[ // color hsv model
      4./6., 1 , 0.5, // dark blue
      4./6., 1 , 1, // blue
      5./6., 1 , 1, // magenta
      1, 1. , 1, // red
      1, 0.5 , 1 // light red
      ];
    real[int] viso(31);

    for (int i = 0; i < viso.n; i++)
      viso[i] = i*0.1;

    plot(uh, viso=viso(0:viso.n-1), value=true, fill=true, wait=true, hsv=colorhsv);

.. subfigstart::

.. _figVisuHSV:

.. figure:: images/Visualization_HSV_Space.png
   :alt: Visualization_HSV_Space
   :width: 90%

   HSV color cylinder

.. _figVisuIsoColorTable:

.. figure:: images/Visualization_HSV.png
   :alt: Visualization_HSV
   :width: 90%

   Isovalue with an other color table

.. subfigend::
   :width: 0.49
   :alt: HSV
   :label: HSV

   HSV

.. note:: See :ref:`HSV example <exampleHSV>` for the complete script.

.. _gnuplot:

Link with gnuplot
-----------------

Example :ref:`Membrane <tutorialMembrane>` shows how to generate a gnuplot from a **FreeFEM** file.
Here is another technique which has the advantage of being online, i.e. one doesn’t need to quit **FreeFEM** to generate a gnuplot.

However, this works only if `gnuplot <http://www.gnuplot.info>`__ is installed, and only on an Unix-like computer.

Add to the previous example:

.. code-block:: freefem
   :linenos:

   {// file for gnuplot
      ofstream gnu("plot.gp");
      for (int i = 0; i < n; i++)
         gnu << xx[i] << " " << yy[i] << endl;
   }

   // to call gnuplot command and wait 5 second (due to the Unix command)
   // and make postscript plot
   exec("echo 'plot \"plot.gp\" w l \n pause 5 \n set term postscript \n set output \"gnuplot.eps\" \n replot \n quit' | gnuplot");

.. figure:: images/Visualization_Gnuplot.png
   :name: figVisuGnuplot
   :width: 50%

   Plots a cut of uh with gnuplot

.. note:: See :ref:`Plot example <examplePlot>` for the complete script.

.. _medit:

Link with medit
---------------

As said above, ``medit`` is a freeware display package by Pascal Frey using OpenGL. Then you may run the following example.

Now ``medit`` software is included in **FreeFEM** under ``ffmedit`` name.

The :freefem:`medit` command parameters are listed in the :ref:`Reference part <referenceMedit>`.

.. figure:: images/Visualization_Medit.png
   :name: figVisuMedit
   :width: 50%

   :freefem:medit` plot

With version 3.2 or later

.. code-block:: freefem
   :linenos:

   load "medit"

   mesh Th = square(10, 10, [2*x-1, 2*y-1]);

   fespace Vh(Th, P1);
   Vh u=2-x*x-y*y;

   medit("u", Th, u);

Before:

.. code-block:: freefem
   :linenos:

   mesh Th = square(10, 10, [2*x-1, 2*y-1]);

   fespace Vh(Th, P1);
   Vh u=2-x*x-y*y;

   savemesh(Th, "u", [x, y, u*.5]); //save u.points and u.faces file
   // build a u.bb file for medit
   {
      ofstream file("u.bb");
      file << "2 1 1 " << u[].n << " 2 \n";
      for (int j = 0; j < u[].n; j++)
         file << u[][j] << endl;
   }
   //call medit command
   exec("ffmedit u");
   //clean files on unix-like OS
   exec("rm u.bb u.faces u.points");

.. note:: See :ref:`Medit example <exampleMedit>` for the complete script.

.. _paraview:

Link with Paraview
------------------

One can also export mesh or results in the ``.vtk`` format in order to post-process data using `Paraview <https://www.paraview.org/>`__.

.. code-block:: freefem
   :linenos:

   load "iovtk"

   mesh Th = square(10, 10, [2*x-1, 2*y-1]);

   fespace Vh(Th, P1);
   Vh u=2-x*x-y*y;

   int[int] Order = [1];
   string DataName = "u";
   savevtk("u.vtu", Th, u, dataname=DataName, order=Order);

.. figure:: images/Visualization_Paraview.png
   :name: figVisuParaview
   :width: 50%

   Paraview plot

.. warning:: Finite element variables saved using paraview **must be in P0 or P1**

.. note:: See :ref:`Paraview example <exampleParaview>` for the complete script.

.. _matlab:

Link with Matlab© and Octave
----------------------------

In order to create a plot from a **FreeFEM** simulation in `Octave <https://www.gnu.org/software/octave/>`__ and `Matlab <https://www.mathworks.com/>`__ the mesh, the finite element space connectivity and the simulation data must be written to files:

.. code-block:: freefem
   :linenos:

   include "ffmatlib.idp"

   mesh Th = square(10, 10, [2*x-1, 2*y-1]);
   fespace Vh(Th, P1);
   Vh u=2-x*x-y*y;

   savemesh(Th,"export_mesh.msh");
   ffSaveVh(Th,Vh,"export_vh.txt");
   ffSaveData(u,"export_data.txt");

Within Matlab or Octave the files can be plot with the `ffmatlib library <https://github.com/samplemaker/freefem_matlab_octave_plot>`__:

.. code-block:: matlab
   :linenos:

   addpath('path to ffmatlib');
   [p,b,t]=ffreadmesh('export_mesh.msh');
   vh=ffreaddata('export_vh.txt');
   u=ffreaddata('export_data.txt');
   ffpdeplot(p,b,t,'VhSeq',vh,'XYData',u,'ZStyle','continuous','Mesh','on');
   grid;

.. figure:: images/Visualization_Matlab_Octave.png
   :name: figVisuMatlab
   :width: 50%

   Matlab / Octave plot

.. note:: For more Matlab / Octave plot examples have a look at the tutorial section :ref:`Matlab / Octave Examples <tutorialMatlabOctavePlot>` or visit the `ffmatlib library <https://github.com/samplemaker/freefem_matlab_octave_plot>`__ on github.
