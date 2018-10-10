# Visualization

Results created by the finite element method can be a huge set of data, so it is very important to render them easy to grasp.

There are two ways of visualization in FreeFem++:

* One, the default view, which supports the drawing of meshes, isovalues of real FE-functions, and of vector fields, all by the command `:::freefem plot` (see [Plot section](#plot) below). For publishing purpose, FreeFem++ can store these plots as postscript files.

* Another method is to use external tools, for example, gnuplot (see [Gnuplot section](#link-with-gnuplot), [medit section](#link-with-medit), [Paraview section](#link-with-paraview)) using the command `:::freefem system` to launch them and/or to save the data in text files.

## Plot

With the command `:::freefem plot`, meshes, isovalues of scalar functions, and vector fields can be displayed.

The parameters of the plot command can be meshes, real FE functions, arrays of 2 real FE functions, arrays of two double arrays, to plot respectively a mesh, a function, a vector field, or a curve defined by the two double arrays.

!!! note
	The length of an arrow is always bound to be in [5‰, 5%] of the screen size in order to see something.

The `:::freefem plot` command parameters are listed in the [Reference part](../reference/Functions/#plot).

The keyboard shortcuts are :

* __enter__ tries to show plot
* __p__ previous plot (10 plots saved)
* __?__ shows this help
* __+,-__ zooms in/out around the cursor 3/2 times
* __=__ resets the view
* __r__ refreshes plot
* __up, down, left, right__ special keys to tanslate
* __3__ switches 3d/2d plot keys :

	__- z,Z__ focal zoom and zoom out
	__- H,h__ increases or decreases the Z scale of the plot
	__- mouse motion__
	__- left button__ rotates
	__- right button__ zooms (ctrl+button on mac)
	__- right button +alt__ tanslates (alt+ctrl+button on mac)

* __a,A__ increases or decreases the arrow size
* __B__ switches between showing the border meshes or not
* __i,I__ updates or not: the min/max bound of the functions to the window
* __n,N__ decreases or increases the number of iso value arrays
* __b__ switches between black and white or color plotting
* __g__ switches between grey or color plotting
* __f__ switches between filling iso or iso line
* __l__ switches between lighting or not
* __v__ switches between show or not showing the numerical value of colors
* __m__ switches between show or not showing the meshes
* __w__ window dump in file ffglutXXXX.ppm
* __*__ keep/drop viewpoint for next plot
* __k__ complex data / change view type
* __ESC__ closes the graphics process before version 3.22, after no way to close
* __otherwise__ does nothing

For example:

```freefem
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
```

|Fig. 1: mesh, isovalue, and vector|Fig. 2: Enlargement in grey of isovalue, and vector|
|:----:|:----:|
|![Three](images/Visualization_Plot.png)|![Threeg](images/Visualization_Plot_Grey.png)|

<center>

|Fig. 3: Plots a cut of uh. Note that a refinement of the same can be obtained in combination with gnuplot|
|:----:|
|![likegnu](images/Visualization_Plot_Gnuplot.png)|

</center>

To change the color table and to choose the value of iso line you can do :

```freefem
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
```

|Fig. 4: hsv color cylinder|Fig. 5: isovalue with an other color table|
|:----:|:----:|
|![hsv](images/Visualization_HSV_Space.png)|![threehsv](images/Visualization_HSV.png)|

!!!note
	See [HSV.edp](../examples/#visualization-hsv) for the complete script.

## Link with gnuplot

Example [Membrane](../tutorial/Membrane) shows how to generate a gnuplot from a FreeFem++ file. Here is another technique which has the advantage of being online, i.e. one doesn't need to quit FreeFem++ to generate a gnuplot.

However, this works only if [gnuplot](http://www.gnuplot.info) is installed, and only on an Unix-like computer.

Add to the previous example:

```freefem
{// file for gnuplot
	ofstream gnu("plot.gp");
	for (int i = 0; i < n; i++)
		gnu << xx[i] << " " << yy[i] << endl;
}

// to call gnuplot command and wait 5 second (due to the Unix command)
// and make postscript plot
exec("echo 'plot \"plot.gp\" w l \n pause 5 \n set term postscript \n set output \"gnuplot.eps\" \n replot \n quit' | gnuplot");
```

|Fig. 6: Plots a cut of uh with gnuplot|
|:----:|
|![gnuplot](images/Visualization_Gnuplot.png)|

!!!note
	See [Plot.edp](../examples/#visualization-plot) for the complete script.

## Link with medit

As said above, `medit` is a freeware display package by Pascal Frey using OpenGL. Then you may run the following example.

Now `medit` software is included in FreeFem++ under `ffmedit` name.

The `:::freefem medit` command parameters are listed in the [Reference part](../reference/ExternalLibraries/#medit).

<center>

|Fig. 7: medit plot|
|:----:|
|![medit2](images/Visualization_Medit.png)|

</center>

With version 3.2 or later

```freefem
load "medit"

mesh Th = square(10, 10, [2*x-1, 2*y-1]);

fespace Vh(Th, P1);
Vh u=2-x*x-y*y;

medit("u", Th, u);
```

Before:

```freefem
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
```

!!!note
	See [Medit.edp](../examples/#visualization-medit) for the complete script.

## Link with Paraview

One can also export mesh or results in the `.vtk` format in order to post-process data using [Paraview](https://www.paraview.org/).

```freefem
load "iovtk"

mesh Th = square(10, 10, [2*x-1, 2*y-1]);

fespace Vh(Th, P1);
Vh u=2-x*x-y*y;

int[int] Order = [1];
string DataName = "u";
savevtk("u.vtu", Th, u, dataname=DataName, order=Order);
```

|Fig. 8: Paraview plot|
|:----:|
|![Paraview](images/Visualization_Paraview.png)|

!!!note
	See [Paraview.edp](../examples/#visualization-paraview) for the complete script.

## Link with Matlab<sup>&copy;</sup> and Octave

In order to create plots from FreeFem++ simulations in [Octave](https://www.gnu.org/software/octave/) and [Matlab](https://www.mathworks.com/) the FEM mesh and the FE function must be exported to text files:

```freefem
mesh Th = square(10, 10, [2*x-1, 2*y-1]);

fespace Vh(Th, P1);
Vh u=2-x*x-y*y;

savemesh(Th,"export_mesh.msh");

ofstream file("export_data.txt"); 
for (int j=0; j<u[].n; j++)
   file << u[][j] << endl;
```

Within Matlab or Octave the files can be processed with the [ffmatlib library](https://github.com/samplemaker/freefem_matlab_octave_plot):

```Matlab
addpath('path to ffmatlib');
[p,b,t]=ffreadmesh('export_mesh.msh');
u=ffreaddata('export_data.txt');
ffpdeplot(p,b,t,'XYData',u,'ZStyle','continuous','Mesh','on');
grid;
```

<center>

|Fig. 9: Matlab / Octave plot|
|:----:|
|![Matlab / Octave](images/Visualization_Matlab_Octave.png)|

</center>

!!!note
    For more Matlab / Octave plot examples have a look at the tutorial section [Matlab / Octave Examples](../tutorial/MatlabOctavePlot/) or visit the [ffmatlib library](https://github.com/samplemaker/freefem_matlab_octave_plot) at github.

