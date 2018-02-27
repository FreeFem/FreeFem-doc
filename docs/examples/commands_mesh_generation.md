# Commands for Mesh Generation

Let us begin with the two important keywords `:::freefem border` and `:::freefem buildmesh`.

All examples in this section come from the files `:::freefem mesh.edp` and `:::freefem tablefunction.edp`.

## Square

The command `:::freefem square` triangulates the unit square.
The following generates a 4*5 grid in the unit square $[0,1]^2$. The labels of the boundaries
are shown in Fig. 5.1.

```freefem
mesh Th = square(4,5);
```

|Fig 5.1: Boundary labels of the mesh by `:::freefem square(10,10)`|
|:----|
|![Square](images/square.svg)|

To construct a $n\times m$ grid in the rectangle \index{rectangle} $[x_0,x_1]\times [y_0,y_1]$, proceeds as follows:

```freefem
real x0=1.2,x1=1.8;
real y0=0,y1=1;
int n=5,m=20;
mesh Th=square(n,m,[x0+(x1-x0)*x,y0+(y1-y0)*y]);

```

!!! note
	Adding the  named  parameter `flags=icase` with icase:

	0. will produce a mesh where all quads are split with diagonal $ x-y=cte$
	1. will produce Union Jack flag type of mesh.
	2. will produce a mesh where all quads are split with diagonal $ x+y=cte$ (v 3.8)
	3. same as case 1  except in two corners such that no triangle with 3 vertices on boundary (v 3.8)
	4. same as case 3  except in two corners such that no triangle with 3 vertices on boundary (v 3.8)

	````freefem
	mesh Th=square(n,m,[x0+(x1-x0)*x,y0+(y1-y0)*y],flags=icase);
	````

!!! note
	Adding the  named  parameter `:::freefem label=labs` will change the 4 default label numbers to `:::freefem labs[i-1]`, for example `:::freefem int[int] labs=[11,12,13,14]`, and adding the  named parameter `:::freefem region=10` will change  the region number  to $10$, for instance (v 3.8).

	To see all these flags at work,  try the file `examples++/square-mesh.edp` :

	```freefem
	for (int i=0;i<5;++i)
  	{
    	int[int] labs=[11,12,13,14];
    	mesh Th=square(3,3,flags=i,label=labs,region=10);
    	plot(Th,wait=1,cmm=" square flags = "+i );
  	}
	```

## Border

Boundaries are defined piecewise by parametrized curves.
The pieces can only intersect at their endpoints, but it is possible to 
join more than two endpoints. This can be used to structure the mesh
if an area thouches a border and create new regions by dividing larger ones:

```freefem
int upper = 1;
int others = 2;
int inner = 3;

border C01(t=0,1){x = 0;         y = -1+t;        label = upper;}
border C02(t=0,1){x = 1.5-1.5*t; y = -1;          label = upper;}
border C03(t=0,1){x = 1.5;       y = -t;          label = upper;}
border C04(t=0,1){x = 1+0.5*t;   y = 0;           label = others;}
border C05(t=0,1){x = 0.5+0.5*t; y = 0;           label = others;}
border C06(t=0,1){x = 0.5*t;     y = 0;           label = others;}
border C11(t=0,1){x = 0.5;       y = -0.5*t;      label = inner;}
border C12(t=0,1){x = 0.5+0.5*t; y = -0.5;        label = inner;}
border C13(t=0,1){x = 1;         y = -0.5+0.5*t;  label = inner;}

int n = 10;
plot(C01(-n)+C02(-n)+C03(-n)+C04(-n)+C05(-n)+C06(-n)+
      C11(n)+C12(n)+C13(n), wait=true);

mesh Th = buildmesh(C01(-n)+C02(-n)+C03(-n)+C04(-n)+C05(-n)+C06(-n)+
      C11(n)+C12(n)+C13(n));

plot(Th, wait=true); // Fig. 5.3

cout << "Part 1 has region number " << Th(0.75, -0.25).region << endl;
cout << "Part 2 has redion number " << Th(0.25, -0.25).region << endl;
```

|Fig. 5.2: Multiple border ends intersect|Fig. 5.3: Generated mesh|
|:----:|:----:|
|![Multiple border ends intersect](images/multiendborder.svg)|![Generated Mesh](images/multiendmesh.svg)|

Triangulation keywords assume that the domain is defined as being on the _left_ (resp right) of its oriented parameterized boundary

$$
\Gamma_j=\{(x,y)\left|\; x=\varphi_x(t),\, y=\varphi_y(t),\, a_j\le t\le b_j\right.\}
$$

To check the orientation plot $t\mapsto (\varphi_x(t),\varphi_y(t)),\, t_0\le t\le t_1$.
If it is as in Fig. 5.4, then the domain lies on the shaded area, otherwise it lies on the opposite side.

|Fig. 5.4: Orientation of the boundary defined by $(\phi_x(t),\phi_y(t))$|
|:----|
|![Border](images/border.svg)|

The general expression to define a triangulation with `:::freefem buildmesh` is

```freefem
mesh Mesh_Name = buildmesh(Gamma1(m1)+...+GammaJ(mj),OptionalParameter);
```

where $m_j$ are positive or negative numbers to indicate how many vertices should be on $\Gamma_j,\,
\Gamma=\cup_{j=1}^J \Gamma_J$, and the optional parameter (separed with comma) can be

* `:::freefem nbvx=<int value>`,  to set the maximal number of  vertices in the mesh.
* `:::freefem fixedborder=<bool value>`, to say if the mesh generator can change the boundary mesh
or not (by default the boundary mesh can change; beware that with periodic boundary conditions 
(see. \ref{periodic BC}), it can be dangerous.


The orientation of boundaries can be changed by changing the sign of $m_j$.
The following example shows how to change the orientation.
The example generates the unit disk with a small circular hole, and assign "1" to the unit disk ("2" to the circle inside).
The boundary label must be non-zero, but it can also be omitted.

```freefem
border a(t=0,2*pi){ x=cos(t); y=sin(t);label=1;}
border b(t=0,2*pi){ x=0.3+0.3*cos(t); y=0.3*sin(t);label=2;}
plot(a(50)+b(+30)); // To see a plot of the border mesh \index{plot!border}
mesh Thwithouthole= buildmesh(a(50)+b(+30));
mesh Thwithhole   = buildmesh(a(50)+b(-30));
plot(Thwithouthole,wait=1,ps="Thwithouthole.eps"); //fig. 5.5
plot(Thwithhole,wait=1,ps="Thwithhole.eps"); // fig. 5.6
```

!!! note
	Notice that the orientation is changed by `:::freefem b(-30)` in 5th line. In 7th line, `:::freefem ps="fileName"` is used to generate a postscript file with identification shown on the figure.

|Fig. 5.5: Mesh without hole |Fig. 5.6: Mesh with hole |
|:----:|:----:|
|![Mesh without hole](images/Th_without_hole.svg)|![Mesh with hole](images/Th_with_hole.svg)|

!!! note
	Borders are evaluated only at the time `:::freefem plot` or `:::freefem buildmesh` is called so the global variable are defined at this time andhere since $r$ is changed between the two border calls the following code will not work because the first border will be computed with r=0.3:

	```freefem
	   real r=1;  border a(t=0,2*pi){ x=r*cos(t); y=r*sin(t);label=1;}
	   r=0.3;     border b(t=0,2*pi){ x=r*cos(t); y=r*sin(t);label=1;}
	   mesh Thwithhole = buildmesh(a(50)+b(-30)); // bug (a trap) because
	   // the two circle have the same radius = $0.3$
	```

## Multi-Border

Sometime it can be useful to make an array of border, but unfortunately it is incompatible with the FreeFem++ syntax. 
So to bypass this problem, the idea is small, if the number of segment of the discretization $n$ is a array, we make  a implicit loop on all the value of the array, and
the index variable $i$ of the loop  is defined after  parameter definition, like in `:::freefem border a(t=0,2*pi;i)` ...

A first very small example: 

```freefem
1: border a(t=0,2*pi;i){ x=(i+1)*cos(t); y=(i+1)*sin(t);label=1;} 
2: int[int] nn=[10,20,30]; 
3: plot(a(nn)); // Plot 3 circles with 10,20,30 points ..
``` 

And  more complex exemple (taken from `:::freefem mesh.edp` example) to define a square with small circles: 

```freefem
// multi border syntax version 3.30 avril 2014 ... 
real[int] xx=[0,1,1,0],
          yy=[0,0,1,1];
// radius, centre of the 4 circles .. 
real[int] RC=[ 0.1, 0.05, 0.05, 0.1],
          XC= [0.2,0.8,0.2,0.8],
          YC= [0.2,0.8,0.8,0.2];
int[int]  NC=[-10,-11,-12,13]; //list number of $\pm$ segments
// of the 4 circles borders  

border bb(t=0,1;i) 
{
// i is the the index variable of the multi border loop 
  int ii = (i+1)%4; real t1 = 1-t;
  x =  xx[i]*t1 + xx[ii]*t;
  y =  yy[i]*t1 + yy[ii]*t;
  label = 0; ; 
}

border cc(t=0,2*pi;i) 
{
  x = RC[i]*cos(t)+XC[i];
  y = RC[i]*sin(t)+YC[i];
  label = i+1; 
}
int[int] nn=[4,4,5,7]; // 4 border, with 4,4,5,7 segment respectively . 
plot(bb(nn),cc(NC),wait=1);
mesh th= buildmesh(bb(nn)+cc(NC)) ; 
plot(th,wait=1); 
```

## Data Structures and Read/Write Statements for a Mesh

Users who want to read a triangulation made elsewhere should see the structure
of the file generated below:

```freefem
border C(t=0,2*pi) { x=cos(t); y=sin(t); }
mesh Th = buildmesh(C(10));
savemesh("mesh_sample.msh");
```

the mesh is shown on Fig. 5.7.

The informations about `:::freefem Th` are saved in the file `:::freefem mesh_sample.msh`
whose structure is shown on Table 5.1.

There $n_v$ denotes the number of vertices, $n_t$ number of triangles and $n_s$ the number of edges on boundary.

For each vertex $q^i,\, i=1,\cdots,n_v$, denote by $(q^i_x,q^i_y)$ the $x$-coordinate and $y$-coordinate.

Each triangle $T_k, k=1,\cdots,10$ has three vertices $q^{k_1},\, q^{k_2},\,q^{k_3}$
that are oriented counterclockwise.

The boundary consists of 10 lines $L_i,\, i=1,\cdots,10$ whose end points are
$q^{i_1},\, q^{i_2}$.

|Fig. 5.7: Mesh by `:::freefem buildmesh(C(10))`||
|:----|:----|
|![Mesh Sample](images/mesh_sample.svg)|In the left figure, we have the following.<br>$n_v=14, n_t=16, n_s=10$<br>$q^1=(-0.309016994375, 0.951056516295)$<br>$\vdots\qquad \vdots\qquad \vdots$<br>$q^{14}=(-0.309016994375, -0.951056516295)$<br>The vertices of $T_1$ are $q^9, q^{12},\, q^{10}$.<br>$\vdots\qquad \vdots\qquad \vdots$<br>The vertices of $T_{16}$ are $q^9, q^{10}, q^{6}$.<br>The edge of 1st side $L_1$ are $q^6, q^5$.<br>$\vdots\qquad \vdots\qquad \vdots$<br>The edge of 10th side $L_{10}$ are $q^{10}, q^6$.|

|Table. 5.1: The structure of `:::freefem mesh_sample.msh`||
|:----|:----|
|Content of the file|Explanation|
|14 16 10<br>-0.309016994375 0.951056516295 1<br>0.309016994375 0.951056516295 1<br>$\cdots$  $\cdots$ $\vdots$<br>-0.309016994375 -0.951056516295 1|$n_v\qquad n_t\qquad n_e$<br>$q^1_x\qquad q^1_y\qquad$ boundary label=1<br>$q^2_x\qquad q^2_y\qquad$ boundary label=1<br><br>$q^{14}_x\qquad q^{14}_y\quad$ boundary label=1|
|9 12 10 0<br>5 9 6 0<br>$\cdots$<br>9 10 6 0|$1_1\qquad 1_2\qquad 1_3\qquad$ region label=0<br>$2_1\qquad 2_2\qquad 2_3\qquad$ region label=0<br><br>$16_1\quad 16_2\qquad 16_3\qquad$ region label=0|
|6 5 1<br>5 2 1<br>$\cdots$<br>10 6 1|$1_1\qquad 1_2\qquad$ boundary label=1<br>$2_1\qquad 2_2\qquad$ boundary label=1<br><br>$10_1\quad 10_2\qquad$ boundary label=1|


In FreeFem++ there are many mesh file formats available for communication with other tools such as emc2, modulef.. $\codered$ (see \refSec{Mesh Files}), The extension of a file implies its format. More details can be found on the file format .msh in the article by F. Hecht "bamg : a bidimensional anisotropic mesh generator".

A mesh file can be read into FreeFem++ except that the names of the borders are lost and only their reference numbers are kept. So these borders have to be referenced by the number which corresponds to their order of appearance in the program, unless this number is overwritten by the keyword `:::freefem label`.  Here are some examples:

```freefem
border floor(t=0,1){ x=t; y=0; label=1;}; // The unit square
border right(t=0,1){ x=1; y=t; label=5;};
border ceiling(t=1,0){ x=t; y=1; label=5;};
border left(t=1,0){ x=0; y=t; label=5;};
int n=10;
mesh th= buildmesh(floor(n)+right(n)+ceiling(n)+left(n));
savemesh(th,"toto.am_fmt"); // "formatted Marrocco" format \index{file!am\_fmt}
savemesh(th,"toto.Th");     // "bamg"-type mesh   \index{file!bamg}
savemesh(th,"toto.msh");    // freefem format \index{file!mesh}
savemesh(th,"toto.nopo");   // modulef format \index{file!nopo} see \cite{modulef}
mesh th2 = readmesh("toto.msh"); // Read the mesh
```

```freefem
// file readmesh.edp
border floor(t=0,1){ x=t; y=0; label=1;}; // The unit square
border right(t=0,1){ x=1; y=t; label=5;};
border ceiling(t=1,0){ x=t; y=1; label=5;};
border left(t=1,0){ x=0; y=t; label=5;};
int n=10;
mesh th= buildmesh(floor(n)+right(n)+ceiling(n)+left(n));
savemesh(th,"toto.am_fmt"); // format "formated Marrocco"
savemesh(th,"toto.Th");     // format database  db mesh "bamg"
savemesh(th,"toto.msh");    // format freefem
savemesh(th,"toto.nopo");   // modulef format see \cite{modulef}
mesh th2 = readmesh("toto.msh");
fespace femp1(th,P1);
femp1 f = sin(x)*cos(y),g;
{ // Save solution
ofstream file("f.txt");
file << f[] << endl;
}  // Close the file (end block)
{  // Read
ifstream file("f.txt");
file >> g[] ;
} // Close reading file (end block)
fespace Vh2(th2,P1);
Vh2 u,v;
plot(g);
//  Find $u$ such that \hfilll
// $ u + \Delta u = g $ in $\Omega $ , \hfilll
// $ u=0$ on $\Gamma_1$ and $\frac{\p u }{\p n} = g$ on $\Gamma_2$  \hfilll
solve pb(u,v) =
    int2d(th)( u*v - dx(u)*dx(v)-dy(u)*dy(v) )
  + int2d(th)(-g*v)
  + int1d(th,5)( g*v) //  $\frac{\p u }{\p n} = g$ on $\Gamma_2$
  + on(1,u=0) ;
plot (th2,u);
```


## Mesh Connectivity and data

The following example explains methods to obtain mesh information.

```freefem
{ // Get mesh information (version 1.37)
mesh Th=square(2,2);
// Get data of the mesh
int nbtriangles=Th.nt;
real area = Th.measure, borderlen = Th.bordermeasure; // Version 3.56
cout << " nb of Triangles = " << nbtriangles << endl;
for (int i=0;i<nbtriangles;i++)
  for (int j=0; j <3; j++)
    cout << i << " " << j << " Th[i][j] = "
         << Th[i][j] << "  x = "<< Th[i][j].x  << " , y= "<< Th[i][j].y
         << ",  label=" << Th[i][j].label << endl;

// Th(i) return the vextex i of Th
// Th[k] return the triangle k of Th

fespace femp1(Th,P1);
femp1 Thx=x,Thy=y; // Hack of get vertex coordinates
// Get vertices information :
int nbvertices=Th.nv;
cout << " nb of vertices = " << nbvertices << endl;
for (int i=0;i<nbvertices;i++)
      cout << "Th(" <<i  << ") : "   // << endl;
           << Th(i).x << " " << Th(i).y  << " " << Th(i).label //v 2.19
           << "       old method: " << Thx[][i] << " " << Thy[][i] << endl;

// Method to find information of point (0.55,0.6)

int it00 = Th(0.55,0.6).nuTriangle; // Then triangle number
int nr00 = Th(0.55,0.6).region;

// Info of a triangle
real area00 = Th[it00].area;
real nrr00 = Th[it00].region;
real nll00 = Th[it00].label; // Same as region in this case.

// Hack  to get a triangle containing point x,y
// or region number (old method)
// -------------------------------------------------------
fespace femp0(Th,P0);
femp0 nuT; // a P0 function  to get triangle numbering
  for (int i=0;i<Th.nt;i++)
   nuT[][i]=i;
femp0 nuReg=region; // A P0 function to get the region number
//  inquire
int it0=nuT(0.55,0.6); // Number of triangle Th's containing (0.55,0,6);
int nr0=nuReg(0.55,0.6); // Number of region of Th's containing (0.55,0,6);


// dump
// -------------------------------------------------------

cout << "  point (0.55,0,6) :triangle number " << it00 << " " << it00
     << ", region = " << nr0 << " == " << nr00 << ",  area K " << area00 << endl;

// New method to get boundary information and mesh adjacent

int k=0,l=1,e=1;
Th.nbe; // Return the number of boundary element
Th.be(k); // return the boundary element k $\in \{0,...,Th.nbe-1\}$
Th.be(k)[l]; // return the vertices l $\in \{0,1\}$ of  boundary elmt k
Th.be(k).Element; // Return the triangle containing the  boundary elmt k
Th.be(k).whoinElement ;   // Return the edge number of triangle containing
//the  boundary elmt k
Th[k].adj(e); // Return adjacent triangle to k by edge e, and change \index{mesh!adj}\hfill
// The value of e to the corresponding edge in the adjacent triangle
Th[k] == Th[k].adj(e) // Non adjacent triangle return the same
Th[k] != Th[k].adj(e) // True adjacent triangle

cout << " print mesh connectivity " << endl;
int nbelement = Th.nt;
for (int k=0;k<nbelement;++k)
  cout << k << " :  " << int(Th[k][0]) << " " << int(Th[k][1])
       << " " <<  int(Th[k][2])
       << " , label  " << Th[k].label << endl;
//

for (int k=0;k<nbelement;++k)
  for (int e=0,ee;e<3;++e)
    //  remark FH hack:  set ee to e, and ee is change by method adj,
    //  in () to make difference with  named parameters.
	    cout << k <<  " " << e << " <=>  " << int(Th[k].adj((ee=e))) << " " << ee
	     << "  adj: " << ( Th[k].adj((ee=e)) != Th[k]) << endl;
    // note :     if k == int(Th[k].adj(ee=e)) not adjacent element


int nbboundaryelement = Th.nbe;

for (int k=0;k<nbboundaryelement;++k)
    cout << k << " : " <<  Th.be(k)[0] << " " << Th.be(k)[1] << " , label "
         << Th.be(k).label <<  " tria  " << int(Th.be(k).Element)
         << " " << Th.be(k).whoinElement <<  endl;
real[int] bb(4);
boundingbox(Th,bb); //  \index{mesh!boundingbox}  \index{boundingbox}
// bb[0] = xmin, bb[1] = xmax, bb[2] = ymin, bb[3] =ymax 
cout << "\n boundingbox  xmin: " << bb[0] << " xmax: " << bb[1] 
                  << " ymin: " << bb[2] << " ymax: " << bb[3] << endl; 
}
```

The output is:

```freefem
 -- square mesh : nb vertices  =9 ,  nb triangles = 8 ,  nb boundary edges 8
    Nb of Vertices 9 ,  Nb of Triangles 8
    Nb of edge on user boundary  8 ,  Nb of edges on true boundary  8
 number of real boundary edges 8
 nb of Triangles = 8
0 0 Th[i][j] = 0  x = 0 , y= 0,  label=4
0 1 Th[i][j] = 1  x = 0.5 , y= 0,  label=1
0 2 Th[i][j] = 4  x = 0.5 , y= 0.5,  label=0
...
6 0 Th[i][j] = 4  x = 0.5 , y= 0.5,  label=0
6 1 Th[i][j] = 5  x = 1 , y= 0.5,  label=2
6 2 Th[i][j] = 8  x = 1 , y= 1,  label=3
7 0 Th[i][j] = 4  x = 0.5 , y= 0.5,  label=0
7 1 Th[i][j] = 8  x = 1 , y= 1,  label=3
7 2 Th[i][j] = 7  x = 0.5 , y= 1,  label=3
 Nb Of Nodes = 9
 Nb of DF = 9
 -- vector function's bound  0 1
 -- vector function's bound  0 1
 nb of vertices = 9
Th(0) : 0 0 4       old method: 0 0
Th(1) : 0.5 0 1       old method: 0.5 0
...
Th(7) : 0.5 1 3       old method: 0.5 1
Th(8) : 1 1 3       old method: 1 1
 Nb Of Nodes = 8
 Nb of DF = 8

 print mesh connectivity
0 :  0 1 4 , label  0
1 :  0 4 3 , label  0
...
6 :  4 5 8 , label  0
7 :  4 8 7 , label  0
0 0 <=>  3 1  adj: 1
0 1 <=>  1 2  adj: 1
0 2 <=>  0 2  adj: 0
...
6 2 <=>  3 0  adj: 1
7 0 <=>  7 0  adj: 0
7 1 <=>  4 0  adj: 1
7 2 <=>  6 1  adj: 1
0 : 0 1 , label 1 tria  0 2
1 : 1 2 , label 1 tria  2 2
...
6 : 0 3 , label 4 tria  1 1
7 : 3 6 , label 4 tria  5 1

 boundingbox  xmin: 0 xmax: 1 ymin: 0 ymax: 1
```

The real characteristic function of a mesh `:::freefem Th` is  `:::freefem chi(Th)`
in 2d and 3d where


`:::freefem chi(Th)(P)=1` if $P\in Th;\qquad$ `:::freefem chi(Th)(P)=0` if $P\not\in Th;$



## The keyword "triangulate"

FreeFem++ is able to build a triangulation from a set of points. This
triangulation is a Delaunay mesh of the convex hull of the set of points.
It can be useful to build a mesh form a table function.

The coordinates of the points and the value of the table function
are defined separately with rows of the form: `:::freefem x  y  f(x,y)`
in a file such as:

```freefem
0.51387 0.175741 0.636237
0.308652 0.534534 0.746765
0.947628 0.171736 0.899823
0.702231 0.226431 0.800819
0.494773 0.12472 0.580623
0.0838988 0.389647 0.456045
...............
```

|Fig. 5.8: Delaunay mesh of the convex hull of point set in file xy|Fig. 5.9: Isovalue of table function|
|:----:|:----:|
|![Th xy](images/Thxy.svg)|![xyf](images/xyf.svg)

The third column of each line is left untouched by the
`:::freefem triangulate` command. But you can use this third value to
define a table function with rows of the form: `:::freefem x  y  f(x,y)`.

The following example shows how to make a mesh from the file "xyf" with the format stated just above.
The command `:::freefem triangulate` command use only use 1st and 2nd rows.

```freefem
mesh Thxy=triangulate("xyf"); // Build the Delaunay mesh of the convex hull
// Points are defined by the first 2 columns of file `xyf}
plot(Thxy,ps="Thxyf.ps"); // (see figure  \ref{Thxy})

fespace Vhxy(Thxy,P1); // create a P1 interpolation
Vhxy fxy; // the function

// Reading the 3rd row to define the function
{ ifstream file("xyf");
   real xx,yy;
   for(int i=0;i<fxy.n;i++)
   file >> xx >>yy >> fxy[][i]; // To read third row only.
   // xx and yy are just skipped
}
plot(fxy,ps="xyf.eps"); // Plot the function (see figure  \ref{xyf})
```

One  new way to build a mesh is to have two arrays one  the $x$ values and the other for the $y$ values:

```freefem
Vhxy xx=x,yy=y; // To set two arrays for the x's and y's
mesh Th=triangulate(xx[],yy[]);
```

# Boundary FEM Spaces Built as Empty Meshes

To define a Finite Element space on a boundary,
we came up with the idea of a mesh with no internal points (call empty mesh).
It can be useful to handle Lagrange multipliers in mixed and mortar methods.

So the function `:::freefem emptymesh` remove all the internal points of a mesh except
points  on  internal boundaries.

```freefem
{  //  new stuff 2004 emptymesh (version 1.40)
 // -- useful to build Multiplicator space
 //  build a mesh without internal point
 // with the same boundary
 //  -----
assert(version>=1.40);
border a(t=0,2*pi){ x=cos(t); y=sin(t);label=1;}
mesh Th=buildmesh(a(20));
Th=emptymesh(Th);
plot(Th,wait=1,ps="emptymesh-1.eps");//see figure \ref{fig emptymesh-1}
}
```

It is also possible to build an empty mesh of a pseudo subregion
with `:::freefem emptymesh(Th,ssd)` using the set of edges of the mesh `:::freefem Th`;
a edge $e$ is in  this set  if with the two adjacent triangles $e =t1\cap t2$
and  $ ssd[T1] \neq ssd[T2]$ where $ssd$  refers to the pseudo region
numbering of triangles, when they are stored in an `:::freefem int[int]` array of size the number of triangles.

```freefem
{  //  new stuff 2004 emptymesh (version 1.40) \hfilll
// -- useful to build Multiplicator space \hfilll
//  build a mesh without internal point \hfilll
// of peusdo sub domain  \hfilll
//  ----- \hfilll
assert(version>=1.40);
mesh Th=square(10,10);
int[int] ssd(Th.nt);
for(int i=0;i<ssd.n;i++) // build the  pseudo region numbering
 {  int iq=i/2;   // because 2 triangle per quad
    int ix=iq%10; //
    int iy=iq/10; //
  ssd[i]= 1 + (ix>=5) +  (iy>=5)*2;
 }
Th=emptymesh(Th,ssd); // build emtpy  with
//  all edge $e = T1 \cap T2$ and $ ssd[T1] \neq ssd[T2]$
plot(Th,wait=1,ps="emptymesh-2.eps");//see figure \ref{fig emptymesh-2}
savemesh(Th,"emptymesh-2.msh");
}
```

|Fig. 5.10: The empty mesh with boundary|Fig. 5.11: An empty mesh defined from a pseudo region numbering of triangle|
|:----:|:----:|
|![Empty mesh 1](images/emptymesh-1.svg)|![Empty mesh 2](images/emptymesh-2.svg)|

# Remeshing
## Movemesh

Meshes can be translated, rotated and deformed by ':::freefem movemesh`; this is useful for elasticity to watch the deformation due to the displacement
$\vec\Phi(x,y)=(\Phi_1(x,y),\Phi_2(x,y))$ of shape. It is also useful to
handle free boundary  problems or optimal shape problems.

If $\Omega$ is triangulated as $T_h(\Omega)$,
and $\Phi$ is a displacement vector then $\Phi(T_h)$ is obtained by

```freefem
mesh  Th=movemesh(Th,[Phi1,Phi2]);
```

Sometimes the transformed mesh is invalid because some triangle
have flip over (now has negative area). To spot such problems one may check the
minimum triangle area in the transformed mesh with
':::freefem checkmovemesh` before any real transformation.

$\codered$
$\Phi_1(x,y)=x+k*\sin(y*\pi)/10)$, $\Phi_2(x,y)=y+k*\cos(y\pi)/10)$ for a big number $k>1$.

```freefem
verbosity=4;
border a(t=0,1){x=t;y=0;label=1;};
border b(t=0,0.5){x=1;y=t;label=1;};
border c(t=0,0.5){x=1-t;y=0.5;label=1;};
border d(t=0.5,1){x=0.5;y=t;label=1;};
border e(t=0.5,1){x=1-t;y=1;label=1;};
border f(t=0,1){x=0;y=1-t;label=1;};
func uu= sin(y*pi)/10;
func vv= cos(x*pi)/10;

mesh Th = buildmesh ( a(6) + b(4) + c(4) +d(4) + e(4) + f(6));
plot(Th,wait=1,fill=1,ps="Lshape.eps");// see figure \ref{lshape}
real coef=1;
real minT0= checkmovemesh(Th,[x,y]); // the min triangle area
while(1) // find a correct move mesh
{
  real minT=checkmovemesh(Th,[x+coef*uu,y+coef*vv]);//the min triangle area
  if (minT > minT0/5) break ; // if big enough
  coef/=1.5;
}

Th=movemesh(Th,[x+coef*uu,y+coef*vv]);
plot(Th,wait=1,fill=1,ps="movemesh.eps");// see figure \ref{movemesh}
```

|Fig. 5.12: L-shape|Fig. 5.13: moved L-shape|
|:----:|:----:|
|![L-shape](images/L-shape.svg)|![moved L shaped](images/moved-L-shape.svg)|

!!! note
	Consider a function $u$ defined on a mesh `:::freefem Th`. A statement like `:::freefem Th=movemesh(Th...)` does not change $u$ and so the old mesh still exists. It will be destroyed when no function use it. A statement like $u=u$ redefines $u$ on the new mesh `:::freefem Th` with interpolation and therefore destroys the old `:::freefem Th` if $u$ was the only function using it.

Now, we give an example of moving mesh with a lagrangian function $u$ defined on the moving mesh.

```freefem
// Simple movemesh example
mesh Th=square(10,10);
fespace Vh(Th,P1);
real t=0;
// ---
// The problem is how to build data without interpolation
// So the data u is moving with the mesh as you can see in the plot
// ---
Vh u=y;
for (int i=0;i<4;i++)
{
 t=i*0.1;
 Vh f= x*t;
 real minarea=checkmovemesh(Th,[x,y+f]);
 if (minarea >0 ) // movemesh will be ok
   Th=movemesh(Th,[x,y+f]);

 cout << " Min area  " << minarea << endl;

 real[int] tmp(u[].n);
 tmp=u[];  // save the value
 u=0;        // to change the FEspace and mesh associated with u
 u[]=tmp;  // set the value of u without any mesh update
 plot(Th,u,wait=1);
};
// In this program, since u is only defined on the last mesh, all the
// previous meshes are deleted from memory.
```


# Regular Triangulation: `:::freefem hTriangle`

For a set $S$, we define the diameter of $S$ by

\[
\textrm{diam}(S)=\sup\{|\mathbf{x}-\mathbf{y}|; \; \mathbf{x},\, \mathbf{y}\in S\}
\]

The sequence $\{\mathcal{T}_h\}_{h\downarrow 0}$ of $\Omega$ is called
_regular_ if they satisfy the following:


1. \[\lim_{h\downarrow 0}\max\{\textrm{diam}(T_k)|\; T_k\in \mathcal{T}_h\}=0\]

2. There is a number $\sigma>0$ independent of $h$ such that
\[\frac{\rho(T_k)}{\textrm{diam}(T_k)}\ge \sigma\qquad \textrm{for all }T_k\in \mathcal{T}_h\]
	where $\rho(T_k)$ are the diameter of the inscribed circle of $T_k$.

We put $h(\mathcal{T}_h)=\max\{\textrm{diam}(T_k)|\; T_k\in \mathcal{T}_h\}$,
which is obtained by

```freefem
mesh Th = ......;
fespace Ph(Th,P0);
Ph h = hTriangle;
cout << "size of mesh = " << h[].max << endl;
```

# Adaptmesh

The function
\[
f(x,y) = 10.0x^3+y^3+\tan^{-1}[\varepsilon/(\sin(5.0y)-2.0x)]
\qquad \varepsilon =  0.0001
\]
sharply varies in value and the initial mesh given by one of the commands of Section \ref{sec:InitialMesh}
cannot reflect its sharp variations.

```freefem
real eps =  0.0001;
real h=1;
real hmin=0.05;
func f = 10.0*x^3+y^3+h*atan2(eps,sin(5.0*y)-2.0*x);

mesh Th=square(5,5,[-1+2*x,-1+2*y]);
fespace Vh(Th,P1);
Vh fh=f;
plot(fh);
for (int i=0;i<2;i++)
 {
   Th=adaptmesh(Th,fh);
   fh=f;  // old mesh is deleted
   plot(Th,fh,wait=1);
 }
```

|Fig. 5.14: 3D graphs for the initial mesh and 1st and 2nd mesh adaptation|
|:----|
|![Mesh adaptation](images/adaptmesh.svg)|

FreeFem++ uses a variable metric/Delaunay automatic meshing
algorithm.

The command:

```freefem
mesh ATh = adaptmesh(Th, f);
```
create the new mesh `:::freefem ATh` adapted to the Hessian

$$
D^2f=(\p^2 f/\p x^2,\, \p^2 f/\p x\p y,
\p^2 f/\p y^2)
$$

of a function (formula or FE-function).
Mesh adaptation is a very powerful tool when the solution of a problem
varies locally and sharply.

Here we solve the problem $\codered$ (\ref{eqn:Poisson})-(\ref{eqn:Dirichlet}),
when $f=1$ and $\Omega$ is a L-shape domain.

|Fig. 5.15: L-shape domain and its boundary name|Fig. 5.16: Final solution after 4-times adaptation|
|:----|:----|
|![L-shape2](images/L-shape2.svg)|![L Shape solution](images/lshapesol.svg)|

**example** (Adapt.edp) The solution has the singularity $r^{3/2},\, r=|x-\gamma|$
at the point $\gamma$ of the intersection of two lines $bc$ and $bd$ (see Fig. 5.15).

```freefem
border ba(t=0,1.0){x=t;   y=0;  label=1;};
border bb(t=0,0.5){x=1;   y=t;  label=1;};
border bc(t=0,0.5){x=1-t; y=0.5;label=1;};
border bd(t=0.5,1){x=0.5; y=t;  label=1;};
border be(t=0.5,1){x=1-t; y=1;  label=1;};
border bf(t=0.0,1){x=0;   y=1-t;label=1;};
mesh Th = buildmesh ( ba(6)+bb(4)+bc(4)+bd(4)+be(4)+bf(6) );
fespace Vh(Th,P1); // set FE space
Vh u,v;            // set unknown and test function
func f = 1;
real error=0.1;    // Level of error
problem Poisson(u,v,solver=CG,eps=1.0e-6) =
    int2d(Th)(  dx(u)*dx(v) + dy(u)*dy(v))
  - int2d(Th) ( f*v )
  + on(1,u=0)  ;
for (int i=0;i< 4;i++)
{
  Poisson;
  Th=adaptmesh(Th,u,err=error);
  error = error/2;
} ;
plot(u);
```

To speed up the adaptation
the default parameter `:::freefem err` of `:::freefem adaptmesh` is changed by hand; it specifies the required precision, so as to make the new mesh finer or coarser.

The problem is coercive and symmetric,
so the linear system can be solved with the conjugate gradient method (parameter `:::freefem solver=CG`) with the stopping criteria on the residual, here `:::freefem eps=1.0e-6`).
By `:::freefem adaptmesh`, the slope of the final solution is correctly computed near
the point of intersection of $bc$ and $bd$ as in Fig. 5.16.

This method is described in detail in $\codered$ \cite{bamg}. It has a number of
default parameters which can be modified :

Si `:::freefem f1,f2` sont des functions  et `:::freefem thold, Thnew` des maillages.

```freefem
    Thnew = adaptmesh(Thold, f1  ...  );
    Thnew = adaptmesh(Thold, f1,f2  ...  ]);
    Thnew = adaptmesh(Thold, [f1,f2]  ...  );
```

The additional parameters of adaptmesh are not written here, hence the  "..."

* `:::freefem hmin=` Minimum edge size (`:::freefem val` is a real. Its default is related to the size of the domain to be meshed and the precision of the mesh generator).

* `:::freefem hmax=` Maximum edge size (`:::freefem val` is a real. It defaults to the diameter of the domain to be meshed)

* `:::freefem err=` $P_1$ interpolation error level (0.01 is the default).  

* `:::freefem errg=` Relative geometrical error. By default this error is 0.01, and in any case it must be lower than $1/\sqrt{2}$.  Meshes created with this option may have some edges smaller than the `:::freefem -hmin` due to geometrical constraints.  

* `:::freefem nbvx=` Maximum number of vertices generated by the mesh generator (9000 is the default).

* `:::freefem nbsmooth=` number of iterations of the smoothing procedure (5 is the default).

* `nbjacoby=` number of iterations in a smoothing procedure during the metric construction, 0 means no smoothing (6 is the default).

* `:::freefem ratio=` ratio for a prescribed smoothing on the metric. If the value is 0 or less than 1.1 no smoothing is done on the metric (1.8 is the default).

	If `:::freefem ratio > 1.1`, the speed of mesh size variations is bounded by $log(\mathtt{ratio})$.  Note: As `:::freefem ratio` gets closer to 1, the number of generated vertices increases. This may be useful to control the thickness of refined regions near shocks or boundary layers .  

* `:::freefem omega=` relaxation parameter for the smoothing procedure (1.0 is the default).

* `:::freefem iso=` If true, forces the metric to be isotropic (false is the default).  

* `:::freefem abserror=` If false, the metric is evaluated using the criterium of equi-repartion of relative error (false is the default). In this case the metric is defined by
\begin{equation}
  \mathcal{M} = \left({1\over\mathtt{err}\,\, \mathtt{coef}^2} \quad {
  |\mathcal{H}| \over max(\mathtt{CutOff},|\eta|)}\right)^p
  \label{eq err rel}
\end{equation}
	otherwise, the metric is evaluated using the criterium of equi-distribution of errors. In this case the metric is defined by
\begin{equation}
  \mathcal{M} = \left({1\over \mathtt{err}\,\,\mathtt{coef}^2} \quad
  {|{\mathcal{H}|} \over
  {\sup(\eta)-\inf(\eta)}}\right)^p.\label{eq err abs}
\end{equation}

* `:::freefem cutoff=` lower limit for the relative error evaluation (1.0e-6 is the default).

* `:::freefem verbosity=` informational messages level (can be chosen between 0 and $\infty$). Also changes the value of the global variable verbosity (obsolete).  

* `:::freefem inquire=` To inquire graphically about the mesh (false is the default).

* `:::freefem splitpbedge=` If true, splits all internal edges in half with two boundary vertices (true is the
default).

* `:::freefem maxsubdiv=` Changes the metric such that the maximum subdivision of a background edge is bound by `:::freefem val` (always limited by 10, and 10 is also the default).

* `:::freefem rescaling=` if true, the function with respect to which the mesh is adapted is rescaled to be between 0 and 1 (true is the default).

* `:::freefem keepbackvertices=` if true, tries to keep as many vertices from the original mesh as possible (true is the default).

* `:::freefem isMetric=` if true, the metric is defined explicitly (false is the default).  If the 3 functions $m_{11}, m_{12}, m_{22}$ are given, they directly define a symmetric matrix field whose Hessian is computed to define a metric. If only one function is given, then it represents the isotropic mesh size at every point.

	For example, if the partial derivatives `:::freefem fxx` ($=\p^2 f/\p x^2$), `:::freefem fxy` ($=\p^2 f/\p x\p y$), `:::freefem fyy` ($=\p^2 f/\p y^2$) are given, we can set `:::freefem Th=adaptmesh(Th,fxx,fxy,fyy,IsMetric=1,nbvx=10000,hmin=hmin);`

* `:::freefem power=` exponent power of the Hessian used to compute the metric (1 is the default).

* `:::freefem thetamax=` minimum corner angle of in degrees (default is $10^\circ$) where the corner is $ABC$ and the angle is the angle of the two vectors ${AB}, {BC}$, ($0$ imply no corner, $90$ imply perp. corner, ...).

* `:::freefem splitin2=` boolean value. If true, splits all triangles of the final mesh into 4 sub-triangles.

* `:::freefem metric=` an array of 3 real arrays to set or get metric data information. The size of these three arrays must be the number of vertices. So if `:::freefem m11,m12,m22` are three P1 finite elements related to the mesh to adapt, you can write: `:::freefem metric=[m11[],m12[],m22[]]` (see file `:::freefem convect-apt.edp` for a full
example)

* `:::freefem nomeshgeneration=` If true, no adapted mesh is generated (useful to compute
only a metric).

* `:::freefem periodic=` Writing `:::freefem periodic=[[4,y],[2,y],[1,x],[3,x]];` builds an adapted periodic mesh. The sample build a biperiodic mesh of a square. (see periodic finite element spaces $\codered$ \ref{periodic BC}, and see `:::freefem sphere.edp` for a  full example)

We can use the command `:::freefem adaptmesh` to build uniform mesh with a contant mesh size. So to build a mesh with a constant mesh size equal to $\frac{1}{30}$ try:

```freefem
// file uniformmesh.edp
mesh Th=square(2,2); // To have initial mesh
plot(Th,wait=1,ps="square-0.eps");
Th= adaptmesh(Th,1./3As writing
0.,IsMetric=1,nbvx=10000);
plot(Th,wait=1,ps="square-1.eps");
Th= adaptmesh(Th,1./30.,IsMetric=1,nbvx=10000); // More the one time du to
Th= adaptmesh(Th,1./30.,IsMetric=1,nbvx=10000); // Adaptation bound `maxsubdiv=`
plot(Th,wait=1,ps="square-2.eps");
```

|Fig. 5.17: Initial mesh|Fig. 5.18: First iteration|Fig. 5.19: Last iteration|
|:----|:----|:----|
|![Initial mesh](images/square-0.svg)|![First iteration](images/square-1.svg)|![Last iteration](images/square-2.svg)|

# Trunc

Two operators have been introduce to remove triangles from a mesh or to divide them.
Operator `:::freefem trunc` has two parameters :

  * `:::freefem label=` sets the label number of new boundary item (one by default)
  * `:::freefem split=` sets the level $n$ of triangle splitting. Each triangle is splitted in  $n\times n$ (one by default).

To create the mesh `:::freefem Th3`
where alls  triangles of a mesh `:::freefem Th`  are splitted in $3{\times}3$, just write:

```freefem
  mesh Th3 = trunc(Th,1,split=3);
```

The  `:::freefem truncmesh.edp` example construct all "trunc" mesh to the support of the basic function  of the space `:::freefem Vh` (cf. `:::freefem abs(u)>0`), split all the  triangles in $5{\times} 5$, and put a label number to $2$ on new boundary.

```freefem
mesh Th=square(3,3);
fespace Vh(Th,P1);
Vh u;
int i,n=u.n;
u=0;
for (i=0;i<n;i++) // All degree of freedom
 {
  u[][i]=1;       // The basic function i
  plot(u,wait=1);
  mesh Sh1=trunc(Th,abs(u)>1.e-10,split=5,label=2);
  plot(Th,Sh1,wait=1,ps="trunc"+i+".eps");// plot the mesh of
  // the function's support
  u[][i]=0; // reset
 }
```

|Fig. 5.20: mesh of support the function P1  number 0, splitted in $5{\times}5$|Fig. 5.21: Mesh of support the function P1  number 6, splitted in $5{\times}5$|
|:----|:----|
|![Trunc0](images/trunc0.svg)|![Trunc6](images/trunc6.svg)|

# Splitmesh

Another way to split mesh triangles is to use {\tt splitmesh}, for example:
```freefem
{  //  new stuff 2004 splitmesh (version 1.37)
  assert(version>=1.37);
  border a(t=0,2*pi){ x=cos(t); y=sin(t);label=1;}
  mesh Th=buildmesh(a(20));
  plot(Th,wait=1,ps="nosplitmesh.eps"); // see figure \ref{fig nosplitmesh}
  Th=splitmesh(Th,1+5*(square(x-0.5)+y*y));
  plot(Th,wait=1,ps="splitmesh.eps"); // see figure \ref{fig splitmesh}
}
```

\twoplot[height=6cm]{nosplitmesh}{splitmesh}{\label{fig nosplitmesh}initial mesh}{\label{fig splitmesh}all left mesh triangle is split  conformaly in `int(1+5*(square(x-0.5)+y*y)\^2} triangles.}


# Meshing Examples

\begin{example}[Two rectangles touching by a side]~
\index{mesh!beam}
```freefem
border a(t=0,1){x=t;y=0;};
border b(t=0,1){x=1;y=t;};
border c(t=1,0){x=t ;y=1;};
border d(t=1,0){x = 0; y=t;};
border c1(t=0,1){x=t ;y=1;};
border e(t=0,0.2){x=1;y=1+t;};
border f(t=1,0){x=t ;y=1.2;};
border g(t=0.2,0){x=0;y=1+t;};
int n=1;
mesh th = buildmesh(a(10*n)+b(10*n)+c(10*n)+d(10*n));
mesh TH = buildmesh ( c1(10*n) + e(5*n) + f(10*n) + g(5*n) );
plot(th,TH,ps="TouchSide.esp"); // Fig. \ref{TouchSide}
```


\begin{example}[NACA0012 Airfoil]~
\index{mesh!NACA0012}
```freefem
border upper(t=0,1) { x = t;
     y = 0.17735*sqrt(t)-0.075597*t
  - 0.212836*(t^2)+0.17363*(t^3)-0.06254*(t^4); }
border lower(t=1,0) { x = t;
     y= -(0.17735*sqrt(t)-0.075597*t
  -0.212836*(t^2)+0.17363*(t^3)-0.06254*(t^4)); }
border c(t=0,2*pi) { x=0.8*cos(t)+0.5;  y=0.8*sin(t); }
mesh Th = buildmesh(c(30)+upper(35)+lower(35));
plot(Th,ps="NACA0012.eps",bw=1);  // Fig. \ref{NACA0012}
```


\twoplot[height=5cm]{TouchSide}{NACA0012}{Two rectangles touching by a side}
{NACA0012 Airfoil}

\begin{example}[Cardioid]~
\index{mesh!Cardioid}
```freefem
real b = 1, a = b;
border C(t=0,2*pi) { x=(a+b)*cos(t)-b*cos((a+b)*t/b);
                        y=(a+b)*sin(t)-b*sin((a+b)*t/b); }
mesh Th = buildmesh(C(50));
plot(Th,ps="Cardioid.eps",bw=1); // Fig. \ref{Cardioid}
```

\begin{example}[Cassini Egg]~
\index{mesh!Cassini Egg}
```freefem
border C(t=0,2*pi) { x=(2*cos(2*t)+3)*cos(t);
                      y=(2*cos(2*t)+3)*sin(t); }
mesh Th = buildmesh(C(50));
plot(Th,ps="Cassini.eps",bw=1); // Fig. \ref{Cassini}
```

\twoplot[height=5cm]{Cardioid}{Cassini}{Domain with Cardioid curve boundary}
{Domain with Cassini Egg curve boundary}

\begin{example}[By cubic Bezier curve]~
\index{mesh!Bezier curve}
```freefem
// A cubic Bezier curve connecting two points with two control points
func real bzi(real p0,real p1,real q1,real q2,real t)
{
  return p0*(1-t)^3+q1*3*(1-t)^2*t+q2*3*(1-t)*t^2+p1*t^3;
}

real[int] p00=[0,1], p01=[0,-1], q00=[-2,0.1], q01=[-2,-0.5];
real[int] p11=[1,-0.9], q10=[0.1,-0.95], q11=[0.5,-1];
real[int] p21=[2,0.7], q20=[3,-0.4], q21=[4,0.5];
real[int] q30=[0.5,1.1], q31=[1.5,1.2];
border G1(t=0,1) { x=bzi(p00[0],p01[0],q00[0],q01[0],t);
                   y=bzi(p00[1],p01[1],q00[1],q01[1],t); }
border G2(t=0,1) { x=bzi(p01[0],p11[0],q10[0],q11[0],t);
                   y=bzi(p01[1],p11[1],q10[1],q11[1],t); }
border G3(t=0,1) { x=bzi(p11[0],p21[0],q20[0],q21[0],t);
                   y=bzi(p11[1],p21[1],q20[1],q21[1],t); }
border G4(t=0,1) { x=bzi(p21[0],p00[0],q30[0],q31[0],t);
                   y=bzi(p21[1],p00[1],q30[1],q31[1],t); }
int m=5;
mesh Th = buildmesh(G1(2*m)+G2(m)+G3(3*m)+G4(m));
plot(Th,ps="Bezier.eps",bw=1);  // Fig \ref{Bezier}
```


\begin{example}[Section of Engine]~
\index{mesh!Section of Engine}
```freefem
real a= 6., b= 1., c=0.5;
border L1(t=0,1) { x= -a; y= 1+b - 2*(1+b)*t; }
border L2(t=0,1) { x= -a+2*a*t; y= -1-b*(x/a)*(x/a)*(3-2*abs(x)/a );}
border L3(t=0,1) { x= a; y=-1-b + (1+ b )*t; }
border L4(t=0,1) { x= a - a*t;   y=0; }
border L5(t=0,pi) { x= -c*sin(t)/2; y=c/2-c*cos(t)/2; }
border L6(t=0,1) { x= a*t;  y=c; }
border L7(t=0,1) { x= a;  y=c + (1+ b-c )*t; }
border L8(t=0,1) { x= a-2*a*t; y= 1+b*(x/a)*(x/a)*(3-2*abs(x)/a); }
mesh Th = buildmesh(L1(8)+L2(26)+L3(8)+L4(20)+L5(8)+L6(30)+L7(8)+L8(30));
plot(Th,ps="Engine.eps",bw=1); // Fig. \ref{Engine}
```


\begin{figure}[hbt]
\begin{multicols}{2}
\begin{center}
\includegraphics*[height=5cm]{Bezier}
\caption{\label{Bezier} Boundary drawed by Bezier curves}
\end{center}
\begin{center}
\vspace{3cm}~~\\
\includegraphics*[height=2.8cm]{Engine}
\caption{\label{Engine} Section of Engine}
\end{center}
\end{multicols}


\begin{example}[Domain with U-shape channel]~
\index{mesh!U-shape channel}
```freefem
real d = 0.1; // width of U-shape
border L1(t=0,1-d) { x=-1; y=-d-t; }
border L2(t=0,1-d) { x=-1; y=1-t; }
border B(t=0,2) { x=-1+t; y=-1; }
border C1(t=0,1) { x=t-1; y=d; }
border C2(t=0,2*d) { x=0; y=d-t; }
border C3(t=0,1) { x=-t; y=-d; }
border R(t=0,2) { x=1; y=-1+t; }
border T(t=0,2) { x=1-t; y=1; }
int n = 5;
mesh Th = buildmesh (L1(n/2)+L2(n/2)+B(n)+C1(n)+C2(3)+C3(n)+R(n)+T(n));
plot(Th,ps="U-shape.eps",bw=1); // Fig \ref{U-shape}
```

\begin{example}[Domain with V-shape cut]~
\index{mesh!V-shape cut}
```freefem
real dAg = 0.01; // angle of V-shape
border C(t=dAg,2*pi-dAg) { x=cos(t); y=sin(t); };
real[int] pa(2), pb(2), pc(2);
pa[0] = cos(dAg); pa[1] = sin(dAg);
pb[0] = cos(2*pi-dAg); pb[1] = sin(2*pi-dAg);
pc[0] = 0; pc[1] = 0;
border seg1(t=0,1) { x=(1-t)*pb[0]+t*pc[0]; y=(1-t)*pb[1]+t*pc[1]; };
border seg2(t=0,1) { x=(1-t)*pc[0]+t*pa[0]; y=(1-t)*pc[1]+t*pa[1]; };
mesh Th = buildmesh(seg1(20)+C(40)+seg2(20));
plot(Th,ps="V-shape.eps",bw=1);  // Fig. \ref{V-shape}
```

\twoplot[height=5cm]{U-shape}{V-shape}{Domain with U-shape channel changed by `d}}
{Domain with V-shape cut changed by `dAg}}

\begin{example}[Smiling face]~
\index{mesh!Smiling face}
```freefem
real d=0.1;
int m=5;
real a=1.5, b=2, c=0.7, e=0.01;
border F(t=0,2*pi) { x=a*cos(t); y=b*sin(t); }
border E1(t=0,2*pi) { x=0.2*cos(t)-0.5; y=0.2*sin(t)+0.5; }
border E2(t=0,2*pi) { x=0.2*cos(t)+0.5; y=0.2*sin(t)+0.5; }
func real st(real t) {
   return sin(pi*t)-pi/2;
}
border C1(t=-0.5,0.5) { x=(1-d)*c*cos(st(t)); y=(1-d)*c*sin(st(t)); }
border C2(t=0,1){x=((1-d)+d*t)*c*cos(st(0.5));y=((1-d)+d*t)*c*sin(st(0.5));}
border C3(t=0.5,-0.5) { x=c*cos(st(t)); y=c*sin(st(t)); }
border C4(t=0,1) { x=(1-d*t)*c*cos(st(-0.5)); y=(1-d*t)*c*sin(st(-0.5));}

border C0(t=0,2*pi) { x=0.1*cos(t); y=0.1*sin(t); }
mesh Th=buildmesh(F(10*m)+C1(2*m)+C2(3)+C3(2*m)+C4(3)
                  +C0(m)+E1(-2*m)+E2(-2*m));
plot(Th,ps="SmileFace.eps",bw=1);  // see Fig. \ref{SmileFace}
}```


\begin{example}[3point bending]~
\index{mesh!3point bending}
```freefem
// Square for Three-Point Bend Specimens fixed on \ttCC{Fix1, Fix2}
// It will be loaded on \ttCC{Load}.
real a=1, b=5, c=0.1;
int n=5, m=b*n;
border Left(t=0,2*a) { x=-b; y=a-t; }
border Bot1(t=0,b/2-c) { x=-b+t; y=-a; }
border Fix1(t=0,2*c) { x=-b/2-c+t; y=-a; }
border Bot2(t=0,b-2*c) { x=-b/2+c+t; y=-a; }
border Fix2(t=0,2*c) { x=b/2-c+t; y=-a; }
border Bot3(t=0,b/2-c) { x=b/2+c+t; y=-a; }
border Right(t=0,2*a) { x=b; y=-a+t; }
border Top1(t=0,b-c) { x=b-t; y=a; }
border Load(t=0,2*c) { x=c-t; y=a; }
border Top2(t=0,b-c) { x=-c-t; y=a; }
mesh Th = buildmesh(Left(n)+Bot1(m/4)+Fix1(5)+Bot2(m/2)+Fix2(5)+Bot3(m/4)
                    +Right(n)+Top1(m/2)+Load(10)+Top2(m/2));
plot(Th,ps="ThreePoint.eps",bw=1); // Fig. \ref{ThreePoint}
```


\begin{figure}[hbt]
\begin{multicols}{2}
\begin{center}
\includegraphics*[height=5cm]{SmileFace}
\caption{\label{SmileFace} Smiling face (Mouth is changeable)}
\end{center}
\begin{center}
\vspace{2cm}~~\\
\includegraphics*[height=2.8cm]{ThreePoint}
\caption{\label{ThreePoint} Domain for three-point bending test}
\end{center}
\end{multicols}


%%%  3d


# How to change the label of elements and border elements of a mesh

Changing the label of elements and border elements will be done using the keyword {\bf{change}}. The parameters for this
command line are for a two dimensional and dimensional case:
\begin{description}
\item [`label =}] is a vector of integer that contains successive pair of the old label number to  the new label number .
\item [`region =}] is a vector of integer that contains successive pair of the old region number to new region number.
\item [`flabel =}]  is a integer function with given the new value of the label (version 3.21).
\item [`fregion =}] is a integer function with given the new value of the region .
\end{description}
%and for a three dimensional case:
%\begin{description}
%\item [`region =}] is a vector of integer that contains the old region number at index $2i$ and the new region number at index  $2i+1$ of tetrahedra.
%\item [`label =}] is a vector of integer that contains the old labels number at index $2i$ and the new labels number at index $2i+1$ of traingles.
%\end{description}

These vectors are composed of $n_{l}$ successive pair of number $O,N$  where $n_{l}$ is the number (label or region)
that we want to change.
For example, we have
\begin{eqnarray}
\label{eq.org.vector.change.label}
\mathtt{label} &= &[ O_{1}, N_{1},  ..., O_{n_{l}},N_{n_{l}} ] \\
\mathtt{region} & =& [ O_{1}, N_{1},  ..., O_{n_{l}},N_{n_{l}} ] 
\end{eqnarray}

%%%ALH-25/2/10-compilation error
%%%where $O_{i}$ is the $i^\mathrm{nd}$ old number (label or region) to change in new number $N_{i}$.

An example of using this function is given in "glumesh2D.edp":  \index{mesh!change}\index{change}\index{label!change}\index{region!change}\index{flabel!change}\index{fregion!change}
\begin{example}[glumesh2D.edp]
\label{changelabel}~
```freefem

 1:
 2: mesh Th1=square(10,10);
 3: mesh Th2=square(20,10,[x+1,y]);
 4: verbosity=3;
 5: int[int] r1=[2,0],  r2=[4,0];
 6: plot(Th1,wait=1);
 7: Th1=change(Th1,label=r1);//Change the label of Edges  2 in 0.
 8: plot(Th1,wait=1);
 9: Th2=change(Th2,label=r2);//Change the label of Edges  4 in 0.
10: mesh Th=Th1+Th2;         //  ``gluing together'' of meshes Th1 and Th2
11: cout << " nb lab = " << int1d(Th1,1,3,4)(1./lenEdge)+int1d(Th2,1,2,3)(1./lenEdge)
12:          << " == " << int1d(Th,1,2,3,4)(1./lenEdge) <<" == " << ((10+20)+10)*2 << endl;
13: plot(Th,wait=1);
14: fespace Vh(Th,P1);
15: macro Grad(u) [dx(u),dy(u)]; // definition of a macro
16: Vh u,v;
17: solve P(u,v)=int2d(Th)(Grad(u)'*Grad(v))-int2d(Th)(v)+on(1,3,u=0);
18: plot(u,wait=1);

```


\paragraph{``gluing'' different mesh}
In line 10 of previous file, the method to ``gluing'' different mesh of the same dimension in FreeFem++ is using.
This function is the operator "+" between meshes. \index{mesh!+}
The method implemented need that the point in adjacent mesh are the same.
%%The method implemented need that the result's mesh of ``gluing'' meshes is conformal.


# Mesh in three dimensions

## cube
From version (3.38-2), a new function `cube} like the function `square} in 2d is the simple way to build cubic object, in plugin `msh3} (need `load "msh3")}.

The following code
```freefem
mesh3 Th=cube(3,4,5); 
```
generates a $3\times 4 \times 5$ grid in the unit cube $[0, 1]^3$.\index{cube}

By defaults the label (after version 3.56-2 and to correct otherwise add `label=l6} with `int[int] l6=[1,2,3,4,5,6];}) are : 

\begin{enumerate}[topsep=0pt,itemsep=-1ex,partopsep=1ex,parsep=1ex]
 \item face $y=0$; \item  face $x=1$, \item  face $y=1$, \item  face $x=0$, \item   face $z=0$, \item  face $z=1$,
\end{enumerate}
and the region number is $0$. 

A full examples of the this function to build a mesh of cube $]-1,1[^3$ with face label
given by $(ix + 4*(iy+1) + 16*(iz+1) ) $ where $(ix,iy,iz)$ is coordinate of the barycenter of the current face. 
```freefem
load "msh3"
int[int] l6=[37,42,45,40,25,57];
int r11=11;
mesh3 Th=cube(4,5,6,[x*2-1,y*2-1,z*2-1],label=l6,flags =3,region=r11); 
cout << " volume " << Th.measure << ", border area "<< Th.bordermeasure <<endl; // \index{mesh3!measure} \index{mesh3!bordermeasure}

// Check label dans  region numbering 
int err =0; 
for(int i=0; i<100; ++i)
{
    real s =int2d(Th,i)(1.);
    real sx=int2d(Th,i)(x);
    real sy=int2d(Th,i)(y);
    real sz=int2d(Th,i)(z);
    
    if( s )
    {
     int ix = (sx/s+1.5), iy=(sy/s+1.5), iz=(sz/s+1.5), 
              ii=(ix + 4*(iy+1) + 16*(iz+1) ) ;    
      //  value of ix,iy,iz =>  face min 0 ,  face max 2  , no face 1
      cout <<" label="<< i << " s " << s << " " << ix << iy << iz 
           << " : " << ii << endl; 
      if( i != ii ) err++;
    }
}   
real volr11 = int3d(Th,r11)(1.) ;
cout << " vol region " << 11 << ": " << volr11 << endl; 
if( (volr11 - Th.measure )>1e-8) err++;
plot(Th,fill=0); 
cout << " nb err= " << err <<endl;
assert(err==0); 
```
\index{cube!flags=}\index{cube!label=}\index{cube!region=}
the output of this script is:
\begin{verbatim}
  Enter: BuildCube: 3
    kind = 3 n tet Cube = 6 / n slip 6 19
  Cube  nv=210 nt=720 nbe=296
  Out:  BuildCube
 volume 8, border area 24
 label=25 s 4 110 : 25
 label=37 s 4 101 : 37
 label=40 s 4 011 : 40
 label=42 s 4 211 : 42
 label=45 s 4 121 : 45
 label=57 s 4 112 : 57
 vol region 11: 8
 nb err= 0
times: compile 0.005363s, execution 0.00218s,  mpirank:0
 CodeAlloc : nb ptr  2856,  size :352744 mpirank: 0
\end{verbatim}

\begin{figure}[htbp]
\begin{center}
  \includegraphics[height=6cm]{func-cube}
\end{center}
  \caption{The mesh 3d  of  function `cube(4,5,6,flags =3)}
  \label{fig:cube}} \index{cube}


## Read/Write Statements for a Mesh in 3D

In three dimensions, the file mesh format supported for input and output files by FreeFem++ are the extension .msh and .mesh.
These formats are described in the chapter on Mesh Files in two dimensions.

\paragraph{extension file .msh}
The structure of the files with extension .msh in 3D is given in Table \ref{tab:mesh3DSample}.
In this structure, $n_v$ denotes the number of vertices, $n_{tet}$ the number of tetrahedra and $n_{tri}$ the number of triangles
For each vertex $q^i,\, i=1,\cdots,n_v$, we denote by $(q^i_x,q^i_y,q^i_z)$ the $x$-coordinate, the $y$-coordinate and the $z$-coordinate.
Each tetrahedra $T_k, k=1,\cdots,n_{tet}$ has four vertices $q^{k_1},\, q^{k_2},\,q^{k_3}, \,q^{k_4}$.
The boundary consists of an union of triangles. Each triangle $be_j, j=1,\cdots,n_{tri}$ has three vertices $q^{j_1},\, q^{j_2},\,q^{j_3}$.
%that are oriented counterclockwise par rapport \`{a} la normale sortante.

\begin{table}[htbp]
\hspace*{3cm}
\begin{tabular}{|ccccc|}
\hline
$n_v$&  $n_{tet}$& $n_{tri}$ & &\\
$q^1_x$& $q^1_y$& $q^1_z$ & Vertex label &\\
$q^2_x$& $q^2_y$&  $q^2_z$ & Vertex label &\\
$\vdots$  &$\vdots$ &$\vdots$ &$\vdots$ &\\
$q^{n_v}_x$& $q^{n_v}_y$&  $q^{n_v}_z$ & Vertex label&\\
$1_1$& $1_2$& $1_3$& $1_4$ & region label \\
$2_1$& $2_2$& $2_3$& $2_4$ & region label  \\
$\vdots$  &$\vdots$ &$\vdots$ &$\vdots$  &$\vdots$ \\
$(n_{tet})_1$& $(n_{tet})_2$& $(n_{tet})_3  $& $(n_{tet})_4$ & region label \\
$1_1$ & $1_2$& $1_3$& boundary label & \\
$2_1$ & $2_2$& $2_3$& boundary label & \\
$\vdots$&  $\vdots$ &$\vdots$ &$\vdots$ &\\
$(n_tri)_{1}$ & $(n_{tri})_2$& $(n_{tri})_3$ & boundary label &\\
\hline
\end{tabular}
 \caption{The structure of mesh file format ``.msh'' in three dimensions.}
\label{tab:mesh3DSample}
\end{table}


\paragraph{extension file .mesh}
\def\Int#1{ {\tt(I)} #1}
\def\Vertex#1{{{\tt Vertex}#1}}
\def\Loop#1#2{{\bf\Large(}\,#1\,{\bf\Large{,\,\,}}\,#2\,{\bf\Large)}}

The data structure for a three dimensional mesh is composed of the data structure presented in Section \ref{meshformatfile.mesh}
and a data structure for tetrahedra. The tetrahedra of a three dimensional mesh are refereed using the following field:
\small
\begin{itemize}
\item {\tt{Tetrahedra}}\\
  \Int{NbOfTetrahedrons} \\
    \Loop{\Loop{\Vertex{$_i^j$}}{j=1,4}\,,\,\Int{$Ref \phi_i^{tet}$} }{ i=1\,,\,NbOfTetrahedrons}
\end{itemize}
This field is express with the notation of Section \ref{meshformatfile.mesh}.

## TeGen: A tetrahedral mesh generator

\paragraph{TetGen}

TetGen is a software developed by Dr. Hang Si of Weierstrass Institute for Applied Analysis and Stochastics of
Berlin in Germany \cite{tetgen}. TetGen is a free for research and non-commercial uses. For any commercial
licence utilization, a commercial licence is available upon request to Hang Si.

This software is a tetrahedral mesh generator of a three dimensional domain defined by its boundary.
The input domain take into account a polyhedral or a piecewise linear complex.
This tetrahedralization is a constrained Delaunay tetrahedralization.

The method used in TetGen to control the quality of the mesh is a Delaunay refinement due to 
Shewchuk \cite{tetgenshewchuk}. The quality measure of this algorithm is the Radius-Edge
Ratio (see Section 1.3.1 \cite{tetgen} for more details). A theoretical bounds of this ratio of the algorithm
of Shewchuk is obtained for a given complex of vertices, constrained segments and facets of surface mesh,
with no input angle less than 90 degree. This theoretical bounds is 2.0.\\

\index{tetg}
The launch of Tetgen is done with the keyword `tetg}. The parameters of this command line is:

\begin{description}
%%\item [`reftet  =}] set the label of tetrahedra.
\item [`label =}] is a vector of integer that contains the old labels number at index $2i$  and the new labels number at index $2i+1$ of Triangles.
This parameters is initialized as label for the keyword change (\ref{eq.org.vector.change.label}).
\index{tetg!switch=}
\item [`switch  =}] A string expression. This string corresponds to the command line switch of Tetgen see Section 3.2 of \cite{tetgen}.
\index{tetg!nbofholes=}
\item [`nbofholes=}] Number of holes (default value \verb!size of holelist/3! (version 3.11) ).
\index{tetg!holelist=}
\item [`holelist =}] This array correspond to {\bf{holelist}} of tetgenio data structure \cite{tetgen}.
A real vector of size $3\times `nbofholes}$. In TetGen, each hole is associated with a point inside this domain.
This vector is $x_{1}^{h}, y_{1}^{h}, z_{1}^{h}, x_{2}^{h}, y_{2}^{h}, z_{2}^{h}, \cdots,$ where $x_{i}^{h},y_{i}^{h},z_{i}^{h}$
is the associated point with the $i^{\mathrm{th}}$ hole.
\index{tetg!nbofregions=}
\item [`nbofregions =}] Number of regions (\verb!size of regionlist/5! (version 3.11) ). 
\index{tetg!regionlist=}
\item [`regionlist =}] This array corresponds to {\bf{regionlist}} of tetgenio data structure \cite{tetgen}.
The attribute and the volume constraint of region are given in this real vector of size $5\times `nbofregions}$.
The $i^{\mathrm{th}}$ region is described by five elements: $x-$coordinate, $y-$coordinate and $z-$coordinate of
a point inside this domain ($x_{i},y_{i},z_{i}$); the attribute ($at_{i}$) and the maximum volume for tetrahedra ($mvol_{i}$) for this region.
The `regionlist} vector is: $x_{1}, y_{1}, z_{1}, at_{1}, mvol_{1}, x_{2}, y_{2}, z_{2}, at_{2}, mvol_{2}, \cdots  $.
\index{tetg!nboffacetcl=}
\item [`nboffacetcl=}] Number of facets constraints \verb!size of facetcl/2! (version 3.11) ).
\index{tetg!facetcl=}
\item [`facetcl=}] This array corresponds to {\bf{facetconstraintlist}} of tetgenio data structure \cite{tetgen}.
The $i^{th}$ facet constraint is defined by the facet marker $Ref_{i}^{fc}$ and the maximum area for faces $marea_{i}^{fc}$.
The `facetcl} array is: $Ref_{1}^{fc}, marea_{1}^{fc}, Ref_{2}^{fc}, marea_{2}^{fc}, \cdots$.
This parameters has no effect if switch `q} is not selected.
%\item [`nbofsegcl=}] Number of segments constraints.
%\item [`segcl =}] This array correspond to {\bf{segmentconstraintlist}} of tetgenio data structure \cite{tetgen}.
\end{description}


Principal switch parameters in TetGen:
\begin{itemize}
\item [`p}] Tetrahedralization of boundary.
\item [`q}] Quality mesh generation. The bound of Radius-Edge Ratio will be given after the option q. By default, this value is 2.0.
\item [`a}] Construct with the volumes constraints on tetrahedra. These volumes constraints are defined with the bound of the previous
switch `q} or in the parameter `regionlist}.
\item [`A}] Attributes reference to region given in the `regionlist}. The other regions have label 0.
The option `AA} gives a different label at each region. This switch work with the option 'p'. If option 'r' is used, this switch has no effect.
\index{tetgreconstruction}
\item [`r}] Reconstructs and Refines a previously generated mesh. This character is only used with the command line `tetgreconstruction}.
\item [`Y}] This switch allow to preserve the mesh on the exterior boundary.
This switch must be used to ensure conformal mesh between two adjacents mesh.
\item [`YY}] This switch allow to preserve the mesh on the exterior and interior boundary.
\item [`C}] The consistency of the result's mesh is testing by TetGen.
\item [`CC}] The consistency of the result's mesh is testing by TetGen and also checks constrained delaunay mesh
(if 'p' switch is selected) or the consistency of Conformal Delaunay (if 'q' switch is selected).
\item [`V}] Give information of the work of TetGen. More information can be obtained in specified 'VV' or 'VVV'.
\item [`Q}] Quiet: No terminal output except errors
\item [`M}] The coplanar facets are not merging.
\item [`T}] Set a tolerance for coplanar test. The default value is $1e-8$.
\item [`d}] Itersections of facets are detected.
\end{itemize}

To obtain a tetrahedral mesh generator with tetgen, we need the surface mesh of three dimensional domain.
We give now the command line in FreeFem++ to construct these meshes.

\index{movemesh23}
\paragraph{keyword: ``movemesh23''}

\index{mesh3}
A simple method to construct a surface is to place a two dimensional domain in a three dimensional space.
This corresponding to move the domain by a displacement vector of this form $\Phi(x,y) = ( \Phi1(x,y), \Phi2(x,y), \Phi3(x,y) )$.
The result of moving a two dimensional mesh Th2 by this three dimensional displacement is obtained using:

```freefem
mesh3 Th3 = movemesh23(Th2,transfo=[$\Phi$1,$\Phi$2,$\Phi$3]);
```

The parameters of this command line are:
\begin{description}
\item [`transfo   =}] `[$\Phi$1, $\Phi$2, $\Phi$3]} set the displacement vector of transformation $\Phi(x,y) = [\Phi1(x,y), \Phi2(x,y), \Phi3(x,y) ]$.
\index{movemesh23!transfo=}
\item [`label   =}] set integer label of triangles
\index{movemesh23!orientation=}
\item [`orientation=}] set integer orientation of mesh.
\index{movemesh23!ptmerge=}
\item [`ptmerge =}] A real expression. When you transform a mesh, some points can be merged. This parameters is the criteria to define two merging points. By default, we use
$$
ptmerge \: = \: 1e-7 \: \:Vol( B ),
$$
where $B$ is the smallest axis parallel boxes containing the discretized domain of $\Omega$ and $Vol(B)$ is the volume of this box.
\end{description}

We can do a ``gluing'' of surface meshes using the process given in Section \ref{sec.changelab.gluemesh}. An example to obtain a three dimensional
mesh using the command line `tetg} and `movemesh23} is given in the file tetgencube.edp.

\index{movemesh}
\begin{example}[tetgencube.edp]
\label{tetgenboxedp}~
```freefem
// file tetgencube.edp
load "msh3"
load "tetgen"

real x0,x1,y0,y1;
x0=1.; x1=2.; y0=0.; y1=2*pi;
mesh Thsq1 = square(5,35,[x0+(x1-x0)*x,y0+(y1-y0)*y]);

func ZZ1min = 0;
func ZZ1max = 1.5;
func XX1 = x;
func YY1 = y;

mesh3 Th31h = movemesh23(Thsq1,transfo=[XX1,YY1,ZZ1max]);
mesh3 Th31b = movemesh23(Thsq1,transfo=[XX1,YY1,ZZ1min]);

/////////////////////////////////
x0=1.; x1=2.; y0=0.; y1=1.5;
mesh Thsq2 = square(5,8,[x0+(x1-x0)*x,y0+(y1-y0)*y]);

func ZZ2 = y;
func XX2 = x;
func YY2min = 0.;
func YY2max = 2*pi;

mesh3 Th32h = movemesh23(Thsq2,transfo=[XX2,YY2max,ZZ2]);
mesh3 Th32b = movemesh23(Thsq2,transfo=[XX2,YY2min,ZZ2]);

/////////////////////////////////
x0=0.; x1=2*pi; y0=0.; y1=1.5;
mesh Thsq3 = square(35,8,[x0+(x1-x0)*x,y0+(y1-y0)*y]);
func XX3min = 1.;
func XX3max = 2.;
func YY3 = x;
func ZZ3 = y;

mesh3 Th33h = movemesh23(Thsq3,transfo=[XX3max,YY3,ZZ3]);
mesh3 Th33b = movemesh23(Thsq3,transfo=[XX3min,YY3,ZZ3]);

////////////////////////////////
mesh3 Th33 = Th31h+Th31b+Th32h+Th32b+Th33h+Th33b; // "gluing" surface meshs to obtain the surface of cube
savemesh(Th33,"Th33.mesh");

// build a mesh of a axis parallel box with TetGen
real[int] domain =[1.5,pi,0.75,145,0.0025];
mesh3 Thfinal = tetg(Th33,switch="paAAQY",regionlist=domain);    // Tetrahelize the interior of the cube with tetgen
savemesh(Thfinal,"Thfinal.mesh");

// build a mesh of a half cylindrical shell of interior radius 1. and exterior radius 2 and heigh 1.5
func mv2x = x*cos(y);
func mv2y = x*sin(y);
func mv2z = z;
mesh3 Thmv2 = movemesh3(Thfinal, transfo=[mv2x,mv2y,mv2z]);
savemesh(Thmv2,"halfcylindricalshell.mesh")
```

The command `movemesh} is describe in the following section.

%% description des options de Tetgen dans ce cas.

\paragraph{The keyword ``tetgtransfo''} \index{tetgtransfo}

This keyword correspond to a composition of command line `tetg} and `movemesh23}:
```freefem
tetgtransfo( Th2, transfo= [$\Phi$1, $\Phi$2, $\Phi$3] ), $\cdots$ ) = tetg( Th3surf, $\cdots$ ),
```

where Th3surf = `movemesh23}( Th2,tranfo=[$\Phi$1, $\Phi$2, $\Phi$3] ) and Th2 is the input two dimensional mesh of `tetgtransfo}.

\index{tetgtransfo!refface=}
\index{tetgtransfo!switch=}
\index{tetgtransfo!regionlist=}
\index{tetgtransfo!nboffacetcl=}
\index{tetgtransfo!facetcl=}
\index{tetgtransfo!ptmerge=}
The parameters of this command line are on the one hand the parameters: \\
\hspace*{2cm} `label}, `switch}, `regionlist} `nboffacetcl} `facetcl}\\
of keyword `tetg} and on the other hand the parameter `ptmerge} of keyword `movemesh23}.

\paragraph{Remark:} To use `tetgtransfo}, the result's mesh of `movemesh23} must be an closed surface and define one region only.
Therefore, the parameter `regionlist} is defined for one region.

An example of this keyword can be found in line  of file ``buildlayers.edp''

\paragraph{The keyword "tetgconvexhull"}\index{tetgconvexhull}

FreeFem++, using tetgen, is able to build a tetrahedralization from a set of points. This
tetrahedralization is a Delaunay mesh of the convex hull of the set of points.

The coordinates of the points can be initialized in two ways. The first is a file that contains
the coordinate of points $X_{i}=(x_{i}, y_{i}, z_{i})$. This files is organized as follows:
$$
\begin{array}{ccc}
n_{v} & & \\
x_{1} & y_{1} & z_{1}  \\
x_{2} & y_{2} & z_{2} \\
\vdots &\vdots & \vdots \\
x_{n_v} & y_{n_v} & z_{n_v}
\end{array}
$$
The second way is to give three arrays that correspond respectively to the
$x-$coordinates, $y-$coordinates and $z-$coordinates.\\

The parameters of this command line are
\begin{description}
\item [`switch  =}] A string expression. This string corresponds to the command line {\it{switch}} of TetGen see Section 3.2 of \cite{tetgen}.
\item [`reftet  =}] An integer expression. set the label of tetrahedra.
\item [`label =}] An integer expression. set the label of triangles.
\end{description}

In the string switch, we can't used the option 'p' and 'q' of tetgen.

## Reconstruct/Refine a three dimensional mesh with TetGen

Meshes in three dimension can be refined using TetGen with the command line `tetgreconstruction}.

The parameter of this keyword are
\begin{description}
\item [`region=}] an integer array that allow to change the region number  of tetrahedra.
This array is defined as the parameter `reftet} in the keyword `change}.
\item [`label=}] an integer array that allow to change the label of boundary triangles.
This array is defined as the parameter `label} in the keyword `change}.
\item [`sizevolume=}] a reel function. This function allows to constraint volume size of tetrahedra in the domain. (see example \ref{ex:tetg-adap} to build 3d adapt mesh
\index{adaptation})
\end{description}

The parameter `switch} `nbofregions}, `regionlist},
`nboffacetcl} and `facetcl} of the command line which call TetGen (tetg)
is used for `tetgrefine}.

In the parameter `switch=}, the character 'r' should be used without the character 'p'.
For instance, see the manual of TetGen \cite{tetgen} for effect of 'r' to other character.

The parameter `regionlist} allows to define a new volume constraint in the region.
The label in the `regionlist} will be the previous label of region.
This parameter and `nbofregions} can't be used with parameter `sizevolume}.

Example:
\begin{example}[refinesphere.edp]
```freefem
// file refinesphere.edp

load "msh3"
load "tetgen"
load "medit"

mesh Th=square(10,20,[x*pi-pi/2,2*y*pi]);  //  $]\frac{-pi}{2},frac{-pi}{2}[\times]0,2\pi[ $
//  a parametrization of a sphere
func f1 =cos(x)*cos(y);
func f2 =cos(x)*sin(y);
func f3 = sin(x);
//  partiel derivative of the parametrization DF
func f1x=sin(x)*cos(y);
func f1y=-cos(x)*sin(y);
func f2x=-sin(x)*sin(y);
func f2y=cos(x)*cos(y);
func f3x=cos(x);
func f3y=0;
// $  M = DF^t DF $
func m11=f1x^2+f2x^2+f3x^2;
func m21=f1x*f1y+f2x*f2y+f3x*f3y;
func m22=f1y^2+f2y^2+f3y^2;

func perio=[[4,y],[2,y],[1,x],[3,x]];
real hh=0.1;
real vv= 1/square(hh);
verbosity=2;
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
plot(Th,wait=1);

verbosity=2;

// construction of the surface of spheres
real Rmin  = 1.;
func f1min = Rmin*f1;
func f2min = Rmin*f2;
func f3min = Rmin*f3;

mesh3 Th3=movemesh23(Th,transfo=[f1min,f2min,f3min]);

real[int] domain = [0.,0.,0.,145,0.01];
mesh3 Th3sph=tetg(Th3,switch="paAAQYY",nbofregions=1,regionlist=domain);

int[int] newlabel = [145,18];
real[int] domainrefine = [0.,0.,0.,145,0.0001];
mesh3 Th3sphrefine=tetgreconstruction(Th3sph,switch="raAQ",reftet=newlabel,
nbofregions=1,regionlist=domain,refinesizeofvolume=0.0001);

int[int] newlabel2 = [145,53];
func fsize = 0.01/(( 1 + 5*sqrt( (x-0.5)^2+(y-0.5)^2+(z-0.5)^2) )^3);
mesh3 Th3sphrefine2=tetgreconstruction(Th3sph,switch="raAQ",reftet=newlabel2,
sizeofvolume=fsize);

medit(``sphere'',Th3sph);
medit(``isotroperefine'' ,Th3sphrefine);
medit(``anisotroperefine'',Th3sphrefine2);

```




## Moving mesh in three dimensions

Meshes in three dimensions can be translated rotated and deformed using the command line movemesh as in the 2D case
(see section movemesh in chapiter 5). If $\Omega$ is tetrahedrized as $T_{h}(\Omega)$, and $\Phi(x,y)=(\Phi1(x,y,z), \Phi1(x,y,z), \Phi3(x,y,z))$
is a displacement vector then $\Phi(T_{h})$ is obtained by

```freefem
mesh3 Th = movemesh( Th, [$\Phi$1, $\Phi$2, $\Phi$3], ... );
```

The parameters of movemesh in three dimensions are
\begin{description}
%\item [`transfo =}] `[$\Phi$1,$\Phi$2, $\Phi$3]} set the displacement vector of 3D transformation $[\Phi1(x,y,z), \Phi2(x,y,z), \Phi3(x,y,z) ]$.
\item [`region  =}] set integer label of tetrahedra. 0 by default.
\item [`label =}] set the label of faces of border. This parameters is initialized as label for the keyword change (\ref{eq.org.vector.change.label}).
\item [`facemerge =}] An integer expression. When you transform a mesh, some faces can be merged. This parameters equals to one if merge's faces is considered.
Otherwise equals to zero. By default, this parameter is equals to 1.
\item [`ptmerge =}] A real expression. When you transform a mesh, some points can be merged. This parameters is the criteria to define two merging points.
By default, we use
$$
ptmerge \: = \: 1e-7 \: \:Vol( B ),
$$
where $B$ is the smallest axis parallel boxes containing the discretion domain of $\Omega$ and $Vol(B)$ is the volume of this box.
\item[`orientation=}]An integer expression ( 1 by by default) , to reverse or not the orientation of tet if not positive.  
\end{description}


%%\paragraph{Remark:} This command line can be used also to move a surface mesh into an other surface mesh.

An example of this command can be found in the file ''Poisson3d.edp'' located in the directory examples++-3d.


## Layer mesh

In this section, we present the command line to obtain a Layer mesh: `buildlayermesh}.
This mesh is obtained by extending a two dimensional mesh in the z-axis.

The domain $\Omega_{3d}$ defined by the layer mesh is equal to $\Omega_{3d} = \Omega_{2d} \times [zmin, zmax]$
where $\Omega_{2d}$ is the domain define by the two dimensional mesh, $zmin$ and $zmax$ are function
of $\Omega_{2d}$ in $R$ that defines respectively the lower surface and upper surface of $\Omega_{3d}$.

\begin{figure}
\hspace*{4cm} \includegraphics[height=5cm]{buillayermesh}
\caption{Example of Layer mesh in three dimension.}
\label{fig-layermeshextend}


For a vertex of a two dimensional mesh $V_{i}^{2d} = (x_{i},y_{i})$, we introduce the number of associated vertices in the $z-$axis $M_{i}+1$.
We denote by $M$ the maximum of $M_{i}$ over the vertices of the two dimensional mesh. This value are called the number of layers
(if $\forall i, \; M_{i}=M$ then there are $M$ layers in the mesh of $\Omega_{3d}$). $V_{i}^{2d}$ generated $M+1$ vertices which are defined by
$$
\forall j=0, \ldots, M, \qquad  V_{i,j}^{3d} = ( x_{i}, y_{i}, \theta_{i}(z_{i,j})  ),
$$
where $(z_{i,j})_{j=0,\ldots,M}$ are the $M+1$ equidistant points on the interval $[zmin( V_{i}^{2d} ), zmax( V_{i}^{2d})]$:
\begin{eqnarray*}
z_{i,j} =  j \: \delta \alpha + zmin(V_{i}^{2d}), \qquad \delta \alpha= \frac{ zmax( V_{i}^{2d} ) - zmin( V_{i}^{2d}) }{M}.
\end{eqnarray*}
The function $\theta_{i}$, defined on  $[zmin( V_{i}^{2d} ), zmax( V_{i}^{2d} )]$, is given by
$$
\theta_{i}(z) = \left \{
\begin{array}{cl}
\theta_{i,0} & \mbox{if} \: z=zmin(V_{i}^{2d}), \\
\theta_{i,j} & \mbox{if} \: z \in ] \theta_{i,j-1}, \theta_{i,j}],\\
\end{array}
\right.
$$
with $(\theta_{i,j})_{j=0,\ldots,M_{i}}$ are the $M_{i}+1$ equidistant points on the interval $[zmin( V_{i}^{2d} ), zmax( V_{i}^{2d} )]$.\\

Set a triangle $K=(V_{i1}^{2d}$, $V_{i2}^{2d}$, $V_{i3}^{2d})$ of the two dimensional mesh.
$K$ is associated with a triangle on the upper surface (resp. on the lower surface) of layer mesh:
$( V_{i1,M}^{3d}, V_{i2,M}^{3d}, V_{i3,M}^{3d} )$ (resp. $( V_{i1,0}^{3d}, V_{i2,0}^{3d}, V_{i3,0}^{3d})$).

Also $K$ is associated with $M$ volume prismatic elements which are defined by
$$
\forall j=0,\ldots,M, \quad H_{j} = ( V_{i1,j}^{3d}, V_{i2,j}^{3d}, V_{i3,j}^{3d}, V_{i1,j+1}^{3d}, V_{i2,j+1}^{3d}, V_{i3,j+1}^{3d} ).
$$

Theses volume elements can have some merged point:
\begin{itemize}
\item 0 merged point : prism
\item 1 merged points : pyramid
\item 2 merged points : tetrahedra
\item 3 merged points : no elements
\end{itemize}

The elements with merged points are called degenerate elements. To obtain a mesh with tetrahedra, we decompose
the pyramid into two tetrahedra and the prism into three tetrahedra. These tetrahedra are obtained by cutting the quadrilateral
face of pyramid and prism with the diagonal which have the vertex with the maximum index (see \cite{hdrHecht} for the reaspn of this choice).\\

The triangles on the middle surface obtained with the decomposition of the volume prismatic elements are the triangles generated by the edges
on the border of the two dimensional mesh. The label of triangles on the border elements and tetrahedra are defined with the label of these
associated elements.\\


The arguments of `buildlayermesh} is a two dimensional mesh and the number of layers $M$.

The parameters of this command are:
\begin{description}
\item [`zbound  =}] [zmin,zmax] where zmin and zmax are functions expression. Theses functions define the lower surface mesh and upper mesh of surface mesh.
\item [`coef    =}] A function expression between [0,1]. This parameter is used to introduce degenerate element in mesh.
The number of associated points or vertex $V_{i}^{2d}$ is the integer part of $coef(V_{i}^{2d}) M$.
\item [`region  =}] This vector is used to initialized the region of tetrahedra. This vector contain  successive pair of  the  2d region number at index $2i$ and the corresponding    3d region number at index $2i+1$, like (\ref{eq.org.vector.change.label}).
become the
\item [`labelmid =}] This vector is used to initialized the 3d labels number  of the vertical face or mid face form the 2d
label number.   This vector contains  successive pair of the  2d label number
at index $2i$ and the corresponding   3d label number at index $2i+1$, like (\ref{eq.org.vector.change.label}).

\item [`labelup  =}] This vector is used to initialized the 3d label numbers  of the upper/top face form the 2d
region number.   This vector contains  successive pair of the  2d region number
at index $2i$ and the corresponding  3d label number at index $2i+1$, like (\ref{eq.org.vector.change.label}).

\item [`labeldown =}] Same as the previous case but for the lower/down face label .
\end{description}

Moreover, we also add post processing parameters that allow to moving the mesh. These parameters correspond to parameters
`transfo}, `facemerge} and `ptmerge} of the command line `movemesh}.

The vector `region}, `labelmid}, `labelup} and `labeldown} These vectors are composed of $n_{l}$ successive pairs of number $O_i,N_l$  where $n_{l}$ is the number (label or region)
that we want to get.

%%% Example a couper entre les differentes commandes
An example of this command line is given in `buildlayermesh.edp}.

\begin{example}[cube.idp]\index{cube}
\label{cube.idp}~
```freefem
load "medit"
load "msh3"
func mesh3 Cube(int[int] & NN,real[int,int] &BB ,int[int,int] & L)
{
  //  first  build the 6 faces of the hex.
  real x0=BB(0,0),x1=BB(0,1);
  real y0=BB(1,0),y1=BB(1,1);
  real z0=BB(2,0),z1=BB(2,1);

  int nx=NN[0],ny=NN[1],nz=NN[2];
  mesh Thx = square(nx,ny,[x0+(x1-x0)*x,y0+(y1-y0)*y]);

  int[int] rup=[0,L(2,1)],  rdown=[0,L(2,0)],
    rmid=[1,L(1,0),  2,L(0,1),  3, L(1,1),  4, L(0,0) ];
  mesh3 Th=buildlayers(Thx,nz,   zbound=[z0,z1],
                       labelmid=rmid,   labelup = rup,
                       labeldown = rdown);

  return Th;
}
```

The unit cube example:
```freefem
 include "Cube.idp"
 int[int]  NN=[10,10,10]; //  the number of step in each  direction
 real [int,int]  BB=[[0,1],[0,1],[0,1]]; // bounding box
 int [int,int]  L=[[1,2],[3,4],[5,6]]; // the label of the 6 face left,right,
//  front, back, down, right
mesh3 Th=Cube(NN,BB,L);
medit("Th",Th); // see figure \ref{figs-cube}
```



The cone example (an axisymtric mesh on a triangle with degenerateness).
\begin{example}[cone.edp]\index{cone}~
```freefem
load "msh3"
load "medit"
// cone using buildlayers with a triangle
real RR=1,HH=1;
border Taxe(t=0,HH){x=t;y=0;label=0;};
border Hypo(t=1,0){x=HH*t;y=RR*t;label=1;};
border Vert(t=0,RR){x=HH;y=t;label=2;};
int nn=10;   real h= 1./nn;
mesh Th2=buildmesh(  Taxe(HH*nn)+ Hypo(sqrt(HH*HH+RR*RR)*nn) + Vert(RR*nn) ) ;
plot(Th2,wait=1); // the 2d mesh

int MaxLayersT=(int(2*pi*RR/h)/4)*4;// number of layers
real zminT = 0, zmaxT = 2*pi; // height $2*pi$
func fx= y*cos(z); func fy= y*sin(z); func fz= x;
int[int] r1T=[0,0], r2T=[0,0,2,2], r4T=[0,2];
// trick  function:
func deg= max(.01,y/max(x/HH,0.4) /RR); // the function defined the proportion
// of number layer close to axis with reference MaxLayersT
mesh3 Th3T=buildlayers(Th2,coef=  deg, MaxLayersT,
           zbound=[zminT,zmaxT],transfo=[fx,fy,fz],
           facemerge=0, region=r1T, labelmid=r2T);
medit("cone",Th3T); // see figure \ref{figs-cone}
```



\twoplot[height=8cm]{cube}{cone}{the mesh of a  cube made with cube.edp \label{figs-cube}}
{the mesh of a cone made with cone.edp \label{figs-cone}}



\index{buildlayers}
\begin{example}[buildlayermesh.edp]
\label{buildlayermesh}~
```freefem

// file buildlayermesh.edp

load "msh3"
load "tetgen"


// Test 1

int C1=99, C2=98; // could be anything
border C01(t=0,pi){ x=t;  y=0;      label=1;}
border C02(t=0,2*pi){ x=pi; y=t;  label=1;}
border C03(t=0,pi){ x=pi-t;  y=2*pi;    label=1;}
border C04(t=0,2*pi){ x=0;    y=2*pi-t; label=1;}

border C11(t=0,0.7){ x=0.5+t;  y=2.5;      label=C1;}
border C12(t=0,2){ x=1.2;    y=2.5+t;  label=C1;}
border C13(t=0,0.7){ x=1.2-t;  y=4.5;     label=C1;}
border C14(t=0,2){ x=0.5;    y=4.5-t; label=C1;}

border C21(t=0,0.7){ x= 2.3+t;     y=2.5;  label=C2;}
border C22(t=0,2){        x=3;   y=2.5+t;  label=C2;}
border C23(t=0,0.7){   x=3-t;     y=4.5;  label=C2;}
border C24(t=0,2){       x=2.3;   y=4.5-t; label=C2;}

mesh Th=buildmesh(    C01(10)+C02(10)+ C03(10)+C04(10)
                    + C11(5)+C12(5)+C13(5)+C14(5)
                    + C21(-5)+C22(-5)+C23(-5)+C24(-5));

mesh Ths=buildmesh(    C01(10)+C02(10)+ C03(10)+C04(10)
                    + C11(5)+C12(5)+C13(5)+C14(5) );

// construction of a box with one hole and two regions
func zmin=0.;
func zmax=1.;
int MaxLayer=10;

func XX = x*cos(y);
func YY = x*sin(y);
func ZZ = z;

int[int] r1=[0,41], r2=[98,98,  99,99, 1,56];
int[int] r3=[4,12];    //  The triangles of uppper surface mesh
// generated by the triangle in the 2D region of mesh Th of label 4 as label 12.
int[int] r4=[4,45];    //  The triangles of lower surface mesh
// generated by the triangle in the 2D region of mesh Th of label 4 as label 45.

mesh3 Th3=buildlayers( Th, MaxLayer, zbound=[zmin,zmax], region=r1,
                labelmid=r2, labelup = r3, labeldown = r4 );
savemesh(Th3,"box2region1hole.mesh");
// construction of a sphere with TetGen
func XX1 = cos(y)*sin(x);
func YY1 = sin(y)*sin(x);
func ZZ1 = cos(x);
string test="paACQ";
cout << "test=" << test << endl;
mesh3 Th3sph=tetgtransfo(Ths,transfo=[XX1,YY1,ZZ1],switch=test,nbofregions=1,
                           regionlist=domain);
savemesh(Th3sph,"sphere2region.mesh");

```


# Meshing examples


\begin{example}[lac.edp]
// file "lac.edp"
```freefem
load ``msh3''
int nn=5;
border cc(t=0,2*pi){x=cos(t);y=sin(t);label=1;}
mesh Th2 = buildmesh(cc(100));
fespace Vh2(Th2,P2);
Vh2 ux,uy,p2;
int[int] rup=[0,2], rdlow=[0,1], rmid=[1,1,2,1,3,1,4,1];
func zmin = 2-sqrt(4-(x*x+y*y));
func zmax = 2-sqrt(3.);

mesh3 Th = buildlayers(Th2,nn,
  coeff = max((zmax-zmin)/zmax, 1./nn),
  zbound=[zmin,zmax],
  labelmid=rmid;
  labelup=rup;
  labeldown=rlow);
savemesh(Th,''Th.meshb'');
exec(``medit Th; Th.meshb'');
```


\begin{example}[tetgenholeregion.edp]
```freefem
// file ``tetgenholeregion.edp''
load "msh3''
load "tetgen"

mesh Th=square(10,20,[x*pi-pi/2,2*y*pi]);  //  $]\frac{-pi}{2},\frac{-pi}{2}[\times]0,2\pi[ $
//  a parametrization of a sphere
func f1 =cos(x)*cos(y);
func f2 =cos(x)*sin(y);
func f3 = sin(x);
//  partiel derivative of the parametrization DF
func f1x=sin(x)*cos(y);
func f1y=-cos(x)*sin(y);
func f2x=-sin(x)*sin(y);
func f2y=cos(x)*cos(y);
func f3x=cos(x);
func f3y=0;
// $  M = DF^t DF $
func m11=f1x^2+f2x^2+f3x^2;
func m21=f1x*f1y+f2x*f2y+f3x*f3y;
func m22=f1y^2+f2y^2+f3y^2;

func perio=[[4,y],[2,y],[1,x],[3,x]];
real hh=0.1;
real vv= 1/square(hh);
verbosity=2;
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
plot(Th,wait=1);

verbosity=2;

// construction of the surface of spheres
real Rmin  = 1.;
func f1min = Rmin*f1;
func f2min = Rmin*f2;
func f3min = Rmin*f3;

mesh3 Th3sph = movemesh23(Th,transfo=[f1min,f2min,f3min]);

real Rmax  = 2.;
func f1max = Rmax*f1;
func f2max = Rmax*f2;
func f3max = Rmax*f3;

mesh3 Th3sph2 = movemesh23(Th,transfo=[f1max,f2max,f3max]);

cout << "addition" << endl;
mesh3 Th3 = Th3sph+Th3sph2;

real[int] domain2 = [1.5,0.,0.,145,0.001,0.5,0.,0.,18,0.001];
cout << "==============================" << endl;
cout << " tetgen call without hole " << endl;
cout << "==============================" << endl;
mesh3 Th3fin = tetg(Th3,switch="paAAQYY",nbofregions=2,regionlist=domain2);
cout << "=============================" << endl;
cout << "finish tetgen call without hole" << endl;
cout << "=============================" << endl;
savemesh(Th3fin,"spherewithtworegion.mesh");

real[int] hole = [0.,0.,0.];
real[int] domain = [1.5,0.,0.,53,0.001];
cout << "=============================" << endl;
cout << "  tetgen call with hole   " << endl;
cout << "=============================" << endl;
mesh3 Th3finhole=tetg(Th3,switch="paAAQYY",nbofholes=1,holelist=hole,
nbofregions=1,regionlist=domain);
cout << "=============================" << endl;
cout << "finish tetgen call with hole   " << endl;
cout << "=============================" << endl;
savemesh(Th3finhole,"spherewithahole.mesh");
```


## Build a 3d mesh of a cube with a balloon

First the `MeshSurface.idp} file to build boundary mesh of a Hexaedra and of a Sphere.


```freefem
func mesh3 SurfaceHex(int[int] & N,real[int,int] &B ,int[int,int] & L,int orientation)
{
    real x0=B(0,0),x1=B(0,1);
    real y0=B(1,0),y1=B(1,1);
    real z0=B(2,0),z1=B(2,1);

    int nx=N[0],ny=N[1],nz=N[2];

    mesh Thx = square(ny,nz,[y0+(y1-y0)*x,z0+(z1-z0)*y]);
    mesh Thy = square(nx,nz,[x0+(x1-x0)*x,z0+(z1-z0)*y]);
    mesh Thz = square(nx,ny,[x0+(x1-x0)*x,y0+(y1-y0)*y]);

    int[int] refx=[0,L(0,0)],refX=[0,L(0,1)];//  Xmin, Ymax faces labels renumbering
    int[int] refy=[0,L(1,0)],refY=[0,L(1,1)];//  Ymin, Ymax faces labesl renumbering
    int[int] refz=[0,L(2,0)],refZ=[0,L(2,1)];//  Zmin, Zmax faces labels renumbering

    mesh3 Thx0 = movemesh23(Thx,transfo=[x0,x,y],orientation=-orientation,label=refx);
    mesh3 Thx1 = movemesh23(Thx,transfo=[x1,x,y],orientation=+orientation,label=refX);
    mesh3 Thy0 = movemesh23(Thy,transfo=[x,y0,y],orientation=+orientation,label=refy);
    mesh3 Thy1 = movemesh23(Thy,transfo=[x,y1,y],orientation=-orientation,label=refY);
    mesh3 Thz0 = movemesh23(Thz,transfo=[x,y,z0],orientation=-orientation,label=refz);
    mesh3 Thz1 = movemesh23(Thz,transfo=[x,y,z1],orientation=+orientation,label=refZ);
    mesh3 Th= Thx0+Thx1+Thy0+Thy1+Thz0+Thz1;
    return Th;
}


func mesh3 Sphere(real R,real h,int L,int orientation)
{
  mesh  Th=square(10,20,[x*pi-pi/2,2*y*pi]);  //  $]\frac{-pi}{2},frac{-pi}{2}[\times]0,2\pi[ $
  //  a parametrization of a sphere
  func f1 =cos(x)*cos(y);
  func f2 =cos(x)*sin(y);
  func f3 = sin(x);
  //    partiel derivative
  func f1x=sin(x)*cos(y);
  func f1y=-cos(x)*sin(y);
  func f2x=-sin(x)*sin(y);
  func f2y=cos(x)*cos(y);
  func f3x=cos(x);
  func f3y=0;
  // the metric on the sphere  $  M = DF^t DF $
  func m11=f1x^2+f2x^2+f3x^2;
  func m21=f1x*f1y+f2x*f2y+f3x*f3y;
  func m22=f1y^2+f2y^2+f3y^2;

  func perio=[[4,y],[2,y],[1,x],[3,x]];  // to store the periodic condition

  real hh=h/R;// hh  mesh size on unite sphere
  real vv= 1/square(hh);
  Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
  Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
  Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
  Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
  int[int] ref=[0,L];

  mesh3  ThS= movemesh23(Th,transfo=[f1*R,f2*R,f3*R],orientation=orientation,refface=ref);
  return ThS;
}

```

The test of the two functions and the call to `tetgen} mesh generator
```freefem
 load "tetgen"
 include "MeshSurface.idp"
    real hs = 0.1;  // mesh size on sphere
    int[int]  N=[20,20,20];
    real [int,int]  B=[[-1,1],[-1,1],[-1,1]];
    int [int,int]  L=[[1,2],[3,4],[5,6]];
    mesh3 ThH = SurfaceHex(N,B,L,1);
    mesh3 ThS =Sphere(0.5,hs,7,1); // "gluing" surface meshs to tolat boundary meshes

    mesh3 ThHS=ThH+ThS;
    savemesh(ThHS,"Hex-Sphere.mesh");
    exec("ffmedit Hex-Sphere.mesh;rm Hex-Sphere.mesh");// see \ref{figs-Hex-Sphere}

    real voltet=(hs^3)/6.;
    cout << " voltet = " << voltet << endl;
    real[int] domaine = [0,0,0,1,voltet,0,0,0.7,2,voltet];

    mesh3 Th = tetg(ThHS,switch="pqaAAYYQ",nbofregions=2,regionlist=domaine);
    medit("Cube-With-Ball",Th);// see \ref{Cube-With-Ball}

```
\twoplot[height=8cm]{Hex-Sphere}{Cube-With-Ball}{The surface mesh of the Hex with internal Sphere \label{figs-Hex-Sphere}}
{The tet mesh of the cube with internal ball\label{figs-Cube-With-Ball}}

# The output solution formats .sol and .solb

With the keyword savesol, we can store a scalar functions, a scalar FE functions,
a vector fields, a vector FE fields, a symmetric tensor and a symmetric FE tensor..
Such format is used in medit.

\paragraph{extension file .sol}
\def\Int#1{ {\tt(I)} #1}
\def\Loop#1#2{{\bf\Large(}\,#1\,{\bf\Large{,\,\,}}\,#2\,{\bf\Large)}}

The first two lines of the file are
\small
\begin{itemize}
\item {\tt MeshVersionFormatted 0}
\end{itemize}
\normalsize

\small
\begin{itemize}
\item {\tt Dimension}
  \Int{dim}
\end{itemize}

The following fields begin with one of the following keyword:
SolAtVertices, SolAtEdges, SolAtTriangles, SolAtQuadrilaterals,
SolAtTetrahedra, SolAtPentahedra, SolAtHexahedra.

In each field, we give then in the next line the number of elements in the solutions
(SolAtVertices: number of vertices, SolAtTriangles: number of triangles, ...). In other lines, we give
 the number of solutions , the type of solution (1: scalar, 2: vector, 3: symmetric tensor).
 And finally,  we give the values of the solutions on the elements.

The file must be ended with the keyword End.

The real element of symmetric tensor
\begin{eqnarray}
\label{savesol.def.symtensor}
ST^{3d}=\left(
\begin{array}{ccc}
ST_{xx}^{3d} & ST_{xy}^{3d} & ST_{xz}^{3d}\\
ST_{yx}^{3d} & ST_{yy}^{3d} & ST_{yz}^{3d} \\
ST_{zx}^{3d} & ST_{zy}^{3d} & ST_{zz}^{3d}
\end{array}
\right)
\qquad \qquad
ST^{2d}= \left(
\begin{array}{cc}
ST_{xx}^{2d} & ST_{xy}^{2d} \\
ST_{yx}^{2d} & ST_{yy}^{2d}
\end{array}
\right)
\end{eqnarray}
stored in the extension .sol are respectively $ST_{xx}^{3d}, ST_{yx}^{3d}, ST_{yy}^{3d}, ST_{zx}^{3d}, ST_{zy}^{3d}, ST_{zz}^{3d}$
and  $ST_{xx}^{2d}, ST_{yx}^{2d}, ST_{yy}^{2d}$

An example of field with the keyword SolAtTetrahedra:
\small
\begin{itemize}
\item {\tt{SolAtTetrahedra}}\\
  \Int{NbOfTetrahedrons}
  {\tt \obeylines
  $ \mathtt{ \quad nbsol \quad typesol^1 \quad ... \quad typesol^n }  $
  $\left(\left(\left( \mathtt{U}_{ij}^k, \quad \forall i \in \{1,...,\mathtt{nbrealsol}^k\}\right), %
\quad \forall k \in \{1,...\mathtt{nbsol}\}\right) %
 \quad \forall j \in \{1,...,\mathtt{NbOfTetrahedrons}\}\right)$
 }
\end{itemize}
where
\begin{itemize}
\item {\tt  nbsol} is an integer equal to the number of solutions
\item $\mathtt{  typesol^k}$, type of the solution  number $ k$, is
  \begin{itemize}
   \item $\mathtt{typesol^k = 1}$ the solution {\tt k} is scalar.
   \item $\mathtt{typesol^k = 2}$ the solution {\tt k} is vectorial.
   \item $\mathtt{typesol^k = 3}$ the solution {\tt k} is a symmetric tensor or symmetric matrix.
   \end{itemize}
\item  $\mathtt{  nbrealsol^k}$ number of real to discribe solution number $k$ is
  \begin{itemize}
   \item $\mathtt{nbrealsol^k = 1}$ the solution {\tt k} is scalar.
   \item $\mathtt{nbrealsol^k = dim}$ the solution {\tt k} is vectorial ($dim$ is the dimension of the solution).
   \item $\mathtt{nbrealsol^k = dim*(dim+1)/2}$ the solution {\tt k} is a symmetric tensor or symmetric matrix.
   \end{itemize}
\item  {\tt U$_{ij}^k$} is a real equal to the value of the component  $i$ of the solution  $k$ at tetrahedra $j$
on the associated mesh.
\end{itemize}

This field is written with the notation of Section \ref{meshformatfile.mesh}.
The format .solb is the same as format .sol but in binary (read/write is faster, storage is less).

A real scalar functions $f1$, a vector fields $\Phi=[\Phi1,\Phi2,\Phi3]$ and a symmetric tensor $ST^{3d}$
(\ref{savesol.def.symtensor}) at the vertices of the three dimensional mesh Th3 is stored in the file "f1PhiTh3.sol" using
\index{savesol!order=}
```freefem
savesol("f1PhiST3dTh3.sol",Th3, $f1$, [$\Phi$1, $\Phi$2, $\Phi$3], VV3, order=1);
```
where $VV3=[ ST_{xx}^{3d}, ST_{yx}^{3d}, ST_{yy}^{3d}, ST_{zx}^{3d}, ST_{zy}^{3d}, ST_{zz}^{3d}]$.
For a two dimensional mesh Th, A real scalar functions $f2$, a vector fields $\Psi=[\Psi1,\Psi2]$ and a symmetric tensor $ST^{2d}$
(\ref{savesol.def.symtensor}) at triangles is stored in the file "f2PsiST2dTh3.solb" using
```freefem
savesol("f2PsiST2dTh3.solb",Th, $f2$, [$\Psi$1, $\Psi$2], VV2, order=0);
```
where $VV2=[ST_{xx}^{2d}, ST_{yx}^{2d}, ST_{yy}^{2d}]$
The arguments of `savesol} functions are the name of a file, a mesh and solutions. These arguments must be given in this order.

The parmameters of this keyword are
\begin{description}
%\item [`nbscalar =}] number of scalar solutions ($nbs$).
%\item [`scalar   =}]  `[$f_{1}$, $f_{2}$, $f_{3}$, ..., $f_{nbs}$]} set the $nbs$ scalar solutions.
%\item [`nbvector =}] number of vector solutions ($nbv$).
%\item [`nbsymtensor =}] number of symmetric tensor solutions ($nbst$).
%\item [`vector    =}] set the real element of $nbv$ vector solutions.
%\item [`symtensor =}] set the real element of $nbst$ symmetric tensor solutions.
\item [`order =}] 0 is the solution is given at the center of gravity of elements.
1 is the solution is given at the vertices of elements.
\end{description}

%%In two dimensional case for vector solutions $v^{j}=(v^{j}_{x},v^{j}_{y}), \quad j=1,\ldots nbv$, the parameter `vector} is equal to
%%`[$v_{x}^{1}$, $v_{y}^{1}$, ..., $v_{x}^{nbv}$, $v_{y}^{nbv}$]}.

%%In two dimensional case for symtensor solutions $st^{j}= \left (
%%\begin{array}{ccc}
%%st_{xx}^{j} & st_{xy}^{j} & st_{xz}^{j} \\
%%st_{yx}^{j} & st_{yy}^{j} & st_{yz}^{j} \\
%%st_{zx}^{j} & st_{zy}^{j} & st_{zz}^{j}
%%\end{array}
%%\right) \quad j=1,\ldots nbst$, the parameter `vector} is equal to
%%`[ $st_{xx}^{1}$, $st_{yx}^{1}$, $st_{yy}^{1}$, ..., $st_{xx}^{nbst}$, $st_{yx}^{nbst}$, $st_{yy}^{nbst}$]}.
%%\item [`symtensor =}] `[ $st_{xx}^{1}$, $st_{yx}^{1}$, $st_{yy}^{1}$, ..., $st_{xx}^{nbst}$, $st_{yx}^{nbst}$, $st_{yy}^{nbst}$]}
%%set the real element of $nbst$ symmetric tensor solutions  $st^{j}=(st^{j}_{x},st^{j}_{y}), \quad j=1,\ldots nbst$.
%%Three dimensional case:
%%\item [`vector   =}] `[$f_{x}^{1}$, $f_{y}^{1}$, $f_{z}^{1}$, ..., $f_{x}^{nbv}$, $f_{y}^{nbv}$, $f_{z}^{nbv}$  ]}
%%set the $nbv$ vector solutions.
%%\item [`symtensor   =}] `[ $f_{xx}^{1}$, $f_{yx}^{1}$, $f_{yy}^{1}$, $f_{zx}^{1}$, $f_{zy}^{1}$, $f_{zz}^{1}$,
%%..., $f_{xx}^{nbst}$, $f_{yx}^{nbst}$, $f_{yy}^{nbst}$, $f_{zx}^{nbst}$, $f_{zy}^{nbst}$, $f_{zz}^{nbst}$ ]}
%%set the $nbst$ symmetric tensor solutions.

In the file, solutions are stored in this order : scalar solutions, vector solutions and finally symmetric
tensor solutions.

# medit

The keyword medit allows to dipslay a mesh alone or a mesh and one or several functions defined on the mesh using the Pascal Frey's freeware medit.
Medit opens its own window and uses OpenGL extensively.
Naturally to use this command  medit must be installed.
%%This command allow to represent one solution.

A vizualisation with medit of scalar solutions $f1$ and $f2$ continuous, piecewise linear and known at the vertices of the mesh Th is obtained using
```freefem
medit("sol1 sol2",Th, $f1$, $f2$, order=1);
```
The first plot  named ``sol1'' display f1. The second plot names ``sol2'' display f2.

The arguments of function `medit} are the name of the differents scenes (separated by a space) of medit, a mesh and solutions.
Each solution is associated with one scene. The scalar, vector and symmetric tensor solutions are specified in the format described in the section dealing with the keyword `savesol}.

%%These arguments must be given in this order.

The parameters of this command line are
\begin{description}
%%\item [`solution =}] An integer parameter. This parameter is equal to 0 if there is no solution to display and otherwise is equal to 1.
\item [`order =}] 0 is the solution is given at the center of gravity of elements.
1 is the solution is given at the vertices of elements.
\index{medit!order=}
\index{medit!meditff=}
\item [`meditff =}] set the name of execute command of medit. By default, this string is medit.
\index{medit!save=}
\item [`save =}] set the name of a file .sol or .solb to save solutions.
%%\item [`scalar =}] set the scalar solutions to display.
%%two dimensional case:
%%\item [`vector =}] set the vector field solution $v=(v_{x},v_{y})$. This parameter is equal to `[$v_{x}$,$v_{y}$]}.
%%\item [`symtensor =}]  set the symmetric tensor solution
%%$
%%f= \left (
%%\begin{array}{cc}
%%f_{xx} & f_{xy} \\
%%f_{yx} & f_{yy}
%%\end{array}
%%\right)
%%$. This parmameter is equal to `[$f_{xx}$, $f_{yx}$, $f_{yy}$]}
%%three dimensional case:
%%\item [`vector =}] set the vector field solution $v=(v_{x},v_{y},v_{z})$. This parameter is equal to  `[$v_{x}$,$v_{y}$,$v_{z}$]}.
%%\item [`symtensor =}] set the symmetric tensor solution
%%$
%%f= \left (
%%\begin{array}{ccc}
%%f_{xx} & f_{xy} & f_{xz}\\
%%f_{yx} & f_{yy} & f_{yz} \\
%%f_{zx} & f_{zy} & f_{zz}
%%\end{array}
%%\right)
%%$. This parameter is equal to `[$f_{xx}$, $f_{yx}$, $f_{yy}$, $f_{zx}$, $f_{zy}$, $f_{zz}$]}
\end{description}

This command line allows also to represent two differents meshes and solutions on them in the same windows.
The nature of solutions must be the same. Hence, we can vizualize in the same window the different
domains in a domain decomposition method for instance. A vizualisation with medit of scalar solutions $h1$ and $h2$
at vertices of the mesh Th1 and Th2 respectively  are obtained using
```freefem
medit("sol2domain",Th1, $h1$, Th2, $h2$, order=1);
```

\begin{example}[meditddm.edp]
```freefem
// meditddm.edp
load "medit"


// Initial Problem:
//Resolution of the following EDP:
//$- \Delta u_s = f$ on   $\Omega =\{ (x,y) |  1 \leq sqrt(x^2+y^2) \geq 2 \}$\hfilll
//$- \Delta u_1 = f1$ on  $\Omega_{1}=\{ (x,y) |  0.5 \leq sqrt(x^2+y^2) \geq 1. \}$\hfilll
//$u = 1$ on $\Gamma$  +  Null Neumman condition on $\Gamma_{1}$ and on $\Gamma_{2}$\hfilll
//We find the solution $u$ in solving two EDP defined on domain $\Omega$ and $\Omega_{1}$\hfilll
//This solution is visualize with medit

verbosity=3;

border Gamma(t=0,2*pi){x=cos(t); y=sin(t); label=1;};
border Gamma1(t=0,2*pi){x=2*cos(t); y=2*sin(t); label=2;};
border Gamma2(t=0,2*pi){x=0.5*cos(t); y=0.5*sin(t); label=3;};

// construction of mesh of domain $\Omega$
mesh Th=buildmesh(Gamma1(40)+Gamma(-40));

fespace Vh(Th,P2);
func f=sqrt(x*x+y*y);
Vh us,v;
macro Grad2(us) [dx(us),dy(us)]  // EOM

problem Lap2dOmega(us,v,init=false)=int2d(Th)(Grad2(v)' *Grad2(us)) 
   - int2d(Th)(f*v)+on(1,us=1) ;

//   Definition of EDP defined on the domain $\Omega$\hfilll
// $- \Delta u_s = f_1$ on $\Omega_{1}$,     $u_s = 1$ on $\Gamma_1$, $\frac{\partial u_s}{\partial n} =0 $ on $\Gamma_{2}$\hfilll
Lap2dOmega;

// construction of mesh of domain $\Omega_{1}$
mesh Th1=buildmesh(Gamma(40)+Gamma2(-40));

fespace Vh1(Th1,P2);
func f1=10*sqrt(x*x+y*y);
Vh1 u1,v1;
macro Grad21(u1) [dx(u1),dy(u1)]  // EOM

problem Lap2dOmega1(u1,v1,init=false)=int2d(Th1)(Grad21(v1)' *Grad21(u1)) 
            - int2d(Th1)(f1*v1)+on(1,u1=1) ;
//   Resolution of EDP defined on the domain $\Omega_{1}$\hfilll
// $- \Delta u_1 = f_1$ on $\Omega_{1}$,     $u-1 = 1$ on $\Gamma_1$, $\frac{\partial u_1}{\partial n} =0 $ on $\Gamma_{2}$\hfilll
Lap2dOmega1;

// vizualisation of solution of the initial problem
medit("solution",Th,us,Th1,u1,order=1,save="testsavemedit.solb");	
```


\begin{example}[StockesUzawa.edp]
```freefem
//  signe of pressure is correct
assert(version>1.18);
real s0=clock();
mesh Th=square(10,10);
fespace Xh(Th,P2),Mh(Th,P1);
Xh u1,u2,v1,v2;
Mh p,q,ppp;


varf bx(u1,q) = int2d(Th)( (dx(u1)*q));
varf by(u1,q) = int2d(Th)( (dy(u1)*q));
varf a(u1,u2)=  int2d(Th)(  dx(u1)*dx(u2) + dy(u1)*dy(u2) )
                    +  on(1,2,4,u1=0)  +  on(3,u1=1) ;

Xh bc1; bc1[] = a(0,Xh);
Xh b;

matrix A= a(Xh,Xh,solver=CG);
matrix Bx= bx(Xh,Mh);
matrix By= by(Xh,Mh);
Xh bcx=1,bcy=0;

func real[int] divup(real[int] & pp)
{
  int verb=verbosity;
   verbosity=0;
   b[]  = Bx'*pp; b[] += bc1[] .*bcx[];
   u1[] = A^-1*b[];
   b[]  = By'*pp; b[] += bc1[] .*bcy[];
   u2[] = A^-1*b[];
   ppp[] =   Bx*u1[];
   ppp[] +=  By*u2[];
   verbosity=verb;
   return ppp[] ;
};
p=0;q=0;u1=0;v1=0;


LinearCG(divup,p[],q[],eps=1.e-6,nbiter=50);

divup(p[]);

plot([u1,u2],p,wait=1,value=true,coef=0.1);
medit("velocity pressure",Th,[u1,u2],p,order=1);
```



%%% frey et al. software


# Mshmet

Mshmet is a software developped by P. Frey that allows to compute an anisotropic metric based on solutions (i.e. Hessian-based). This sofware can return also an isotropic metric. Moreover,  mshmet can construct also a metric suitable for level sets interface capturing. The solution can be defined on 2D or 3D structured/unstructured meshes. For example, the solution can be an error estimate of a FE solutions.

%%Error estimate for 2d and 3d unstructured meshes.
%%Compute anisotropic metric based on solution variations (i.e. Hessian-based).
%%One option allows to construct a metric suitable for level set interface capturing.

Solutions for mshmet are given as an argument. The solution can be a func, a vector func, a symmetric tensor, a FE func, a FE vector func and a FE symmetric tensor. The symmetric tensor argument is defined as this type of data for datasol argument. This software accepts more than one solution.

For example, the metric $M$ computed with mshmet  for the solution $u$ defined on the mesh $Th$ is obtained by writing.
```freefem
fespace Vh(Th,P1);
Vh u; // a scalar FE func
real[int] M = mshmet(Th,u);
```

The parameters of the keyword mshmet are :
\begin{itemize}\parskip=0cm
%\item	`metric =  <3KN_IdE>} 
\item	`normalization =  <b>} do a normalisation of all solution in $[0,1]$.
\item	`aniso =  <b>} build aniso metric if 1 ( delault 0: iso) 
\item	`levelset =  <b>} {build metric for level set method (default: false)}
\item	`verbosity =  <l>}
\item	`nbregul =  <l>} number of regularization's iteration of solutions given (default 0).
\item	`hmin =  <d>}
\item	`hmax =  <d>}
\item	`err =  <d>} level of error. 
\item	`width =  <d>} the width 
\item `metric}= a vector of double. This vector contains an initial metric given to mshmet. The structure of the metric vector is described in the next paragraph.

\item `loptions}=]a vector of integer of size 7. This vector contains the integer parameters of mshmet(for expert only).
\begin{itemize} \item loptions(0): normalization (default 1).
\item loptions(1): isotropic parameters (default 0). 1 for isotropic metric results otherwise 0.
\item loptions(2): level set parameters (default 0). 1 for building level set metric otherwise 0.
\item loptions(3): debug parameters (default 0). 1 for turning on debug mode otherwise 0.
\item loptions(4): level of verbosity (default 10).
\item loptions(5): number of regularization's iteration of solutions given (default 0). 
\item loptions(6): previously metric parameter (default 0). 1 for using previous metric otherwise 0.
\end{itemize}

\item `doptions}= a vector of double of size 4. This vector contains the real parameters of mshmet (for expert only).
\begin{itemize}
\item doptions(0):  hmin : min size parameters  (default 0.01).
\item doptions(1):  hmax : max size parameters (default 1.0).
\item doptions(2):  eps : tolerance parameters ( default 0.01).
\item doptions(2):  width : relative width for Level Set ($0<w<1$) ( default 0.05).
\end{itemize}
\end{itemize}
The result of the keyword `mshmet} is a `real[int]} which contains the metric computed by  `mshmet}  at the different vertices $V_{i}$ of the mesh.

With $nv$ is the number of vertices, the structure of this vector is
$$ M_{iso}= ( m(V_0), m(V_1), \ldots, m(V_{nv}) )^t$$  for a isotropic metric $m$. For a symmetric tensor metric
$
h=\left(
\begin{array}{ccc}
m_{1 1} & m_{1 2} & m_{1 3}\\
m_{2 1} & m_{2 2} & m_{2 3} \\
m_{3 1} & m_{3 2} & m_{3 3}
\end{array}
\right)$, the parameters `metric}  is $$M_{aniso}= ( H(V_{0}), \ldots, H(V_{nv}) )^t $$
where $H(V_{i})$ is the vector of size 6 defined by \verb![m11,m21,m22,m31,m32,m33]!


\begin{example}[mshmet.edp]
\label{mshmet}~
```freefem
load "mshmet"
load "medit"
load "msh3"

border a(t=0,1.0){x=t;   y=0;  label=1;};
border b(t=0,0.5){x=1;   y=t;  label=2;};
border c(t=0,0.5){x=1-t; y=0.5;label=3;};
border d(t=0.5,1){x=0.5; y=t;  label=4;};
border e(t=0.5,1){x=1-t; y=1;  label=5;};
border f(t=0.0,1){x=0;   y=1-t;label=6;};
mesh Th = buildmesh (a(6) + b(4) + c(4) +d(4) + e(4) + f(6));
savemesh(Th,"th.msh");
fespace Vh(Th,P1);
Vh u,v;
real error=0.01;
problem Problem1(u,v,solver=CG,eps=1.0e-6) =
    int2d(Th,qforder=2)( u*v*1.0e-10+  dx(u)*dx(v) + dy(u)*dy(v))
  +int2d(Th,qforder=2)( (x-y)*v);

func zmin=0;
func zmax=1;
int MaxLayer=10;
mesh3 Th3 = buildlayers(Th,MaxLayer,zbound=[zmin,zmax]);
fespace Vh3(Th3,P2);
fespace Vh3P1(Th3,P1);
Vh3 u3,v3;
Vh3P1 usol;
problem Problem2(u3,v3,solver=sparsesolver) =
   int3d(Th3)( u3*v3*1.0e-10+ dx(u3)*dx(v3) + dy(u3)*dy(v3) + dz(u3)*dz(v3))
  - int3d(Th3)( v3) +on(0,1,2,3,4,5,6,u3=0);
Problem2;
cout << u3[].min << " " << u3[].max << endl;
savemesh(Th3,"metrictest.bis.mesh");
savesol("metrictest.sol",Th3,u3);

real[int] bb=mshmet(Th3,u3);
cout << bb << endl;
for(int ii=0; ii<Th3.nv; ii++)
  usol[][ii]=bb[ii];
savesol("metrictest.bis.sol",Th3,usol);
```




# FreeYams

FreeYams is a surface mesh adaptation software which is developed by P. Frey. This software is a new version of yams. The adapted surface mesh is constructed with a geometric metric tensor field. This  field is based on the intrinsic properties of the discrete surface. Also this software allows to construct a simplification of a mesh. This decimation is  based on the Hausdorff distance between the initial and the current triangulation. Compared to the software yams, FreeYams  can be used also to produce anisotropic triangulations adapted to level set simulations. A technical report on FreeYams is not available yet but a documentation on yams exists at  http://www.ann.jussieu.fr/$\sim$frey/software.html \cite{tech.freeyams}.

To call FreeYams in FreeFem++, we used the keyword freeyams. The arguments of this function are the initial mesh and/or metric. The metric with freeyams are a function, a FE function, a symmetric tensor function, a symmetric tensor FE function or a vector of double. If the metric is vector of double, this data must be given in `metric} parameter. Otherwise, the metric is given in the argument.

For example, the  adapted mesh of $Thinit$ defined by the metric $u$ defined as FE function is obtained in writing.
```freefem
fespace Vh(Thinit,P1);
Vh u;
mesh3 Th=freeyams(Thinit,u);
```
The symmetric tensor argument for freeyams keyword is defined as this type of data for datasol argument.
%\let\oldparskip=\parskip
\begin{itemize}
\parskip=0pt
\item	`aniso =  <b>}  aniso or iso metric  (default 0, iso)
\item	`mem =  <l>}  memory of for freeyams in Mb (delaulf -1, freeyams choose)
\item	`hmin =  <d>}  
\item	`hmax =  <d>}
\item	`gradation =  <d>}  
\item	`option =  <l>}
%\parskip=\oldparskip
\begin{description}

 \item  [0] : mesh optimization (smoothing+swapping)
 \item  [1] :  decimation+enrichment adaptated to a metric map.  (default)
 \item [-1]: decimation adaptated to a metric map. 
 \item  [2] : decimation+enrichment with a Hausdorff-like method
 \item [-2]:  decimation  with a Hausdorff-like method
 \item  [4] : split triangles recursively. 
 \item  [9] : No-Shrinkage Vertex Smoothing
 \end{description}

\item	`ridgeangle =  <d>}
\item	`absolute =  <b>}
\item	`verbosity =  <i>}


\item `metric=} vector expression. This parameters contains the metric at the different vertices on the initial mesh. With $nv$ is the number of vertices, this vector is $$ M_{iso}= ( m(V_0), m(V_1), \ldots, m(V_{nv}) )^t$$  for a scalar metric $m$. For a symmetric tensor metric
$
h=\left(
\begin{array}{ccc}
m_{1 1} & m_{1 2} & m_{1 3}\\
m_{2 1} & m_{2 2} & m_{2 3} \\
m_{3 1} & m_{3 2} & m_{3 3}
\end{array}
\right)$, the parameters `metric}  is $$M_{aniso}= ( H(V_{0}), \ldots, H(V_{nv}) )^t $$
where $H(V_{i})$ is the vector of size 6 defined by \verb![m11,m21,m22,m31,m32,m33]!
\item `loptions=} a vector of integer of size 13. This vectors contains the integer options of FreeYams. (just for the expert )
\begin{itemize}
\item loptions(0):  anisotropic parameter (default 0). If you give an anisotropic metric 1 otherwise 0.
\item loptions(1):  Finite Element correction parameter (default 0). 1 for {\it{no}} Finite Element correction otherwise 0.
\item loptions(2):  Split multiple connected points parameter (default 1). 1 for splitting multiple connected points otherwise 0.
\item loptions(3):  maximum value of memory size in Mbytes (default -1: the size is given by freeyams). 	
\item loptions(4):  set the value of the connected component which we want to obtain. (Remark: freeyams give an automatic value at each connected component).
\item loptions(5):  level of verbosity
\item loptions(6):  Create point on straight edge (no mapping) parameter  (default 0). 1 for creating point on straight edge otherwise 0.
\item loptions(7):  validity check during smoothing parameter. This parameter is only used with No-Shrinkage Vertex Smoothing optimization (optimization option parameter 9). 1 for No validity checking during smoothing otherwise 0.
\item loptions(8):  number of desired's vertices  (default -1).
\item loptions(9):  number of iteration of  optimizations (default 30).
\item loptions(10): no  detection parameter (default 0) . 1 for detecting the ridge on the mesh otherwise 0. The ridge definition is given in the parameter doptions(12).
\item loptions(11): no vertex smoothing parameter (default 0). 1 for smoothing the vertices otherwise 0.
\item loptions(12):  Optimization level parameter (default 0). 
\begin{itemize}

 \item\hspace*{0.3cm}  0 : mesh optimization (smoothing+swapping)
 \item\hspace*{0.3cm}  1 :  decimation+enrichment adaptated to a metric map. 
 \item\hspace*{0.3cm} -1: decimation adaptated to a metric map. 
 \item\hspace*{0.3cm}  2 : decimation+enrichment with a Hausdorff-like method
 \item\hspace*{0.3cm} -2:  decimation  with a Hausdorff-like method
 \item\hspace*{0.3cm}  4 : split triangles recursively. 
 \item\hspace*{0.3cm}  9 : No-Shrinkage Vertex Smoothing
 \end{itemize}

 \end{itemize}
\item [`doptions}=] a vector of double of size 11. This vectors contains the real options of freeyams.
\begin{itemize}
%%doptions(0): 0  !! iso (default 0.0). ????
\item doptions(0):  Set  the geometric approximation (Tangent plane deviation)  (default 0.01).
\item doptions(1):  Set the lamda parameter (default -1. ).
\item doptions(2):  Set the mu parmeter (default  -1. ).
\item doptions(3):  Set the 	 gradation value  (Mesh density control)  (default 1.3).
\item doptions(4):  Set the minimal size(hmin) (default -2.0: the size is automatically computed).  
\item doptions(5):  Set the maximal size(hmax)	(default -2.0: the size is automatically computed). 
\item doptions(6):  Set the tolerance of the control of Chordal deviation (default -2.0). 	
\item doptions(7):  Set the quality of degradation  (default 0.599).
\item doptions(8):  Set the declic parameter (default 2.0).
\item doptions(9): Set the angular walton limitation parameter (default 45 degree).	
\item doptions(10):  Set the angular ridge detection (default 45 degree). 
\end{itemize}

\end{itemize}


\begin{example}[freeyams.edp]
\label{freeyams}~
```freefem
load "msh3"
load "medit"
load "freeyams"
int nn=20;
mesh Th2=square(nn,nn);
fespace Vh2(Th2,P2);
Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1; 
mesh3 Th=buildlayers(Th2,nn, zbound=[zmin,zmax],
                          reffacemid=rmid, reffaceup = rup, reffacelow = rdown);

mesh3 Th3 = freeyams(Th);
medit("maillagesurfacique",Th3,wait=1);

```


# mmg3d

Mmg3d is a 3D remeshing software developed by C. Dobrzynski and P. Frey 

(http://www.math.u-bordeaux1.fr/$\sim$dobj/logiciels/mmg3d.php). To obtain a version of this library send an e-mail at : 

cecile.dobrzynskimath.ubordeaux1.fr or pascal.freyupmc.fr.

This software allows to remesh an initial mesh made of tetrahedra. This initial mesh is adapted to a geometric metric tensor field or to a displacement vector (moving rigid body). The metric can be obtained with mshmet (see section \ref{sec:mshmet}).\\

\begin{remark} : 
\begin{description}
\item[(a)] If no metric is given, an isotropic metric is computed by analyzing the size of the edges in the initial mesh.
\item[(b)] if a displacement is given, the vertices of the surface triangles are moved without verifying the geometrical structure of the new surface mesh.
\end{description}
\end{remark}
The parameters of  mmg3d are :
\begin{itemize}\parskip=0cm
%\item [`nvmax}=] integer expresion. It's correspond to the number of maximum vertices in the solution mesh.
%\item [`ntrimax}=] integer expresion. It's correspond to the number of maximum triangles in the solution mesh.
%\item [`ntetmax}=] integer expresion. It's correspond to the number of maximum triangles in the solution mesh.
\item `options}= vector expression. This vector contains the option parameters of `mmg3d}. It is a vector of $6$ values, with the following meaning:
\begin{description}
\item[(0)]   optimization parameters : (default 1) \\
                     \hspace*{0.3cm}  0 : mesh optimization. \\
                     \hspace*{0.3cm}  1 : adaptation with metric (deletion and insertion vertices) and optimization. \\
                     \hspace*{0.3cm} -1: adaptation with metric (deletion and insertion vertices) without optimization. \\
                     \hspace*{0.3cm}  4 : split tetrahedra (be careful modify the surface). \\
                     \hspace*{0.3cm}  9 : moving mesh with optimization. \\
                     \hspace*{0.3cm} -9: moving mesh without optimization.

\item[(1)]  debug mode :  (default 0)\\
		 \hspace*{0.3cm} 1 : turn on debug mode.\\
		 \hspace*{0.3cm} 0 : otherwise.
		
\item[(2)] Specify the size of bucket per dimension ( default 64)

\item[(3)] swapping mode : (default 0)\\
		\hspace*{0.3cm} 1 : no edge or face flipping. \\
		\hspace*{0.3cm} 0 : otherwise.

\item[(4)] insert points mode : (default 0)\\
		\hspace*{0.3cm} 1 : no edge splitting or collapsing and no insert points. \\
		\hspace*{0.3cm} 0 : otherwise.		
		
\item[(5)] verbosity level (default 3)
\end{description}
\item `memory}= integer expression. Set the maximum memory size of new mesh in Mbytes. By default the number of maximum vertices, tetrahedra and triangles are respectively 500 000,  3000 000, 100000 which represent approximately a memory of 100 Mo.
\item `metric}= vector expression. This vector contains the metric given at mmg3d. It is a vector of size $nv$ or 6 $nv$ respectively for an istropic and anisotropic metric where $nv$ is the number of vertices in the initial mesh. The structure of `metric} vector is described in the mshmet's section(section \ref{sec:mshmet}).
\item `displacement}= `[$\Phi$1, $\Phi$2, $\Phi$3]} set the displacement vector of the initial mesh \\
$\Phi(x,y) = [\Phi1(x,y), \Phi2(x,y), \Phi3(x,y) ]$.
\item `displVect=} sets the vector displacement in a vector expression. This vector contains the displacement at each point of the initial mesh. It is a vector of size 3 $nv$.
\end{itemize}

An example using this function is given in "mmg3d.edp":
\begin{example}[mmg3d.edp]
\label{mmg3dsimple}~
```freefem
// test mmg3d
load "msh3"
load "medit"
load "mmg3d"
include "../examples++-3d/cube.idp"

int n=6;
int[int]  Nxyz=[12,12,12];
real [int,int]  Bxyz=[[0.,1.],[0.,1.],[0.,1.]];
int [int,int]  Lxyz=[[1,1],[2,2],[2,2]];
mesh3 Th=Cube(Nxyz,Bxyz,Lxyz);

real[int] isometric(Th.nv);{
  for( int ii=0; ii<Th.nv; ii++)
    isometric[ii]=0.17;
}

mesh3 Th3=mmg3d(  Th,  memory=100, metric=isometric);
				
medit("init",Th);
medit("isometric",Th3);
```


An example of a moving mesh is given in `fallingspheres.edp}":
\begin{example}[fallingspheres.edp]
```freefem
load "msh3"  load "tetgen"  load "medit"  load "mmg3d"                                                               
include "MeshSurface.idp"

// build mesh of a box (311)  wit 2 holes  (300,310)

real hs = 0.8; 
int[int]  N=[4/hs,8/hs,11.5/hs];
real [int,int]  B=[[-2,2],[-2,6],[-10,1.5]];
int [int,int]  L=[[311,311],[311,311],[311,311]];
mesh3 ThH = SurfaceHex(N,B,L,1);
mesh3 ThSg =Sphere(1,hs,300,-1); 
mesh3 ThSd =Sphere(1,hs,310,-1);   ThSd=movemesh3(ThSd,transfo=[x,4+y,z]);
mesh3 ThHS=ThH+ThSg+ThSd;// gluing surface meshes 
medit("ThHS", ThHS); // see surface mesh

real voltet=(hs^3)/6.;
real[int] domaine = [0,0,-4,1,voltet];
real [int] holes=[0,0,0,0,4,0];
mesh3 Th = tetg(ThHS,switch="pqaAAYYQ",nbofregions=1,regionlist=domaine, nbofholes=2,holelist=holes);    
medit("Box-With-two-Ball",Th);
// End build mesh 

int[int] opt=[9,0,64,0,0,3];   // options  of mmg3d see freeem++ doc 
real[int] vit=[0,0,-0.3];
func zero = 0.;
func dep  = vit[2];

fespace Vh(Th,P1); 
macro Grad(u) [dx(u),dy(u),dz(u)] //

Vh uh,vh; //  to compute the displacemnt field 
problem Lap(uh,vh,solver=CG) = int3d(Th)(Grad(uh)'*Grad(vh))  //') for emacs
				  + on(310,300,uh=dep) +on(311,uh=0.); 

for(int it=0; it<29; it++){ 
  cout<<"  ITERATION       "<<it<<endl;
  Lap;
  plot(Th,uh);
  Th=mmg3d(Th,options=opt,displacement=[zero,zero,uh],memory=1000); 
 }
```



%%%  fin 3d

# A first 3d isotope mesh adaptation process

\begin{example}[Laplace-Adapt-3d.edp]
\label{ex:tetg-adap}~\hfill\break
```freefem
load "msh3" load "tetgen" load "mshmet" load "medit"
//build initial mesh
int nn  = 6;
int[int] l1111=[1,1,1,1],l01=[0,1],l11=[1,1];//label numbering to have all label to 1 
mesh3 Th3=buildlayers(square(nn,nn,region=0,label=l1111),
      nn,  zbound=[0,1],  labelmid=l11,   labelup = l01,  labeldown = l01);
Th3 = trunc(Th3,(x<0.5) | (y < 0.5) | (z < 0.5) ,label=1);// remove the $]0.5,1[^3 cube$
//end of build initial mesh
fespace Vh(Th3,P1);
Vh u,v,usol,h;

macro Grad(u) [dx(u),dy(u),dz(u)] // EOM

problem Poisson(u,v,solver=CG) = int3d(Th3)( Grad(u)'*Grad(v) )  
                                 -int3d(Th3)( 1*v ) + on(1,u=0);

real errm=1e-2;// level of error 
for(int ii=0; ii<5; ii++)
{
  Poisson;// solve Poisson equation. 
  cout <<" u min, max = " <<  u[].min << " "<< u[].max << endl;
  h=0. ;// for resizing h[] because the mesh change 
  h[]=mshmet(Th3,u,normalization=1,aniso=0,nbregul=1,hmin=1e-3,hmax=0.3,err=errm);
  cout <<" h min, max = " <<  h[].min << " "<< h[].max 
       << " " << h[].n << " " << Th3.nv << endl;
  plot(u,wait=1);
  errm*= 0.8;// change the level of error
  cout << " Th3" << Th3.nv < " " << Th3.nt << endl;
  Th3=tetgreconstruction(Th3,switch="raAQ",sizeofvolume=h*h*h/6.);//rebuild mesh
  medit("U-adap-iso-"+ii,Th3,u,wait=1);}
```


# Build a 2d mesh from a isoline

 The idea is get the discretization of a isoline to fluid meshes, this tool  can be useful to
 construct meshes from image. First, we give an example of the isovalue meshes $0.2$  of analytical function $ \sqrt{(x-1/2)^2 +(y-1/2)^2}$,
 on unit square. 
\begin{example}[isoline.edp]
\label{ex:isoline.edp}~\hfill\break
```freefem
load "isoline" // load the plugin "isoline"

real[int,int] xy(3,1); // to store the isoline points 
int[int] be(1);// to store the begin , end couple of lines
{// a block for memory management 
  mesh Th=square(10,10);//,[x*.5,y*0.5]);
  fespace Vh(Th,P1);
  Vh u= sqrt(square(x-0.5)+square(y-0.5));
  real iso= 0.2 ;
  real[int] viso=[iso];
  plot(u,viso=viso,Th);// to see the iso line 

  int nbc= isoline(Th,u,xy,close=1,iso=iso,beginend=be,smoothing=0.1);
```

The isoline parameters are `Th} the mesh, the expression $u$ , the  bidimentionnal array `xy}
to store the list coordinate of the points. The list of named parameter are:
\begin{description}
\item[iso=]  value of the isoline  to compute ($\bm{0}$ is the default value)
\item[close=] close the iso line with the border (def. true), we add the part of the mesh border such the value
 is less than the iso value 
\item[smoothing=] nb of smoothing process  is  the ${l} ^{r} {s} $ where
$l$ is the length of the current  line component, $r$ the ratio, $s$ is smoothing  value. The smoothing default value is $\bm{0}$.
\item[ratio=] the ratio ( $1$ by default). 
\item[eps=] relative $\varepsilon$  (see code ??)  (def 1e-10 )
\item[beginend=]  array to get begin, end couple of each of sub line  (resize automatically)
\item[file=] to save the data  curve in data file for gnu plot
\end{description}

In the array `xy} you get the list of vertices of the isoline,
each connex line go from $i= i_0^c ,\dots, i_1^c-1$ with $i_0^c =be(2*c)$  $i_1^c =be(2*c+1)$, and
  where $x_i= xy(0,i), y_i=yx( 1,i), l_i=xy(2,i) $. 
Here  $l_i$ is the length of the line (the origin of the line is point  $i_0^c$).

The  sense  of the  isoline is such that the upper part is at the left  size of the  isoline.
So here : the minimum is a point $0.5,05$ so  the curve 1 turn in the clockwise  sense, 
the order of each component are sort such the the number of point by component is decreasing .

```freefem   
  cout << " nb of the line  component   = " << nbc << endl; 
  cout << " n = " << xy.m << endl; // number  of points 
  cout << "be = " << be << endl; //  begin end of the each componant

  // show the lines component 
  for( int c=0;c<nbc; ++c) 
  {
    int i0 = be[2*c], i1 = be[2*c+1]-1;//begin,end of the line component
    cout << " Curve " << c << endl; 
    for(int i=i0; i<= i1; ++i)
       cout << " x= " << xy(0,i) <<" y= " << xy(1,i) << " s= " 
             << xy(2,i) << endl; 
    plot([xy(0,i0:i1),xy(1,i0:i1)],wait=1,viso=viso,cmm = " curve "+c);
  }
}// end of block for  memory management 

cout << " len of  last  curve = " << xy(2,xy.m-1) << endl;; 
```

We also have a  new function to parametrize easly  a
discret `Curve} defined by couple $be, xy$. \index{Curve}

```freefem
border Curve0(t=0,1) // the extern boundary 
{ int c =0; // component 0
  int i0 = be[2*c], i1 = be[2*c+1]-1;   
  P=Curve(xy,i0,i1,t); // Curve 0
  label=1; 
} 

border Curve1(t=0,1) 
{ int c =1; // component 1
  int i0 = be[2*c], i1 = be[2*c+1]-1;   
  P=Curve(xy,i0,i1,t);  // Curve 1
  label=1; 
} 

plot(Curve1(100)); // show curve. 
mesh Th= buildmesh(Curve1(-100));// because 
plot(Th,wait=1);// 
```


Secondly, we use this idea to build meshes from image, we use the plugins
`ppm2rnm} to read `pgm} gray scale image, and we extract the 
gray contour at level $0.25$.  

\begin{example}[Leman-mesh.edp]
\label{ex:Leman-mesh.edp}~\hfill\break
```freefem
load "ppm2rnm" load "isoline"
string leman="lg.pgm"; //see figure \ref{fig:lg}
real AreaLac =  580.03; // in $Km^2$
real hsize= 5; // mesh sir in pixel ..
real[int,int] Curves(3,1);
int[int] be(1);
int nc;// nb of curve 
{  
  real[int,int] ff1(leman); // read  image (figure \ref{fig:lg}) 
  // and set to an rect.  array \index{ppm2rnm}
  int nx = ff1.n, ny=ff1.m; // grey value between 0 to 1 (dark)
  // build a Cartesian mesh such that the origin is qt the right place.
  mesh Th=square(nx-1,ny-1,[(nx-1)*(x),(ny-1)*(1-y)]);   
   // warning  the numbering is of the vertices (x,y) is 
   // given by $  i = x/nx + nx* y/ny $
  fespace Vh(Th,P1);
   Vh f1; f1[]=ff1; //  transforme array in finite element function.
  nc=isoline(Th,f1,iso=0.25,close=1,Curves,beginend=be,smoothing=.1,ratio=0.5); 
}
// the longest isoline : the lac .. 
int ic0=be(0), ic1=be(1)-1;		
plot([Curves(0,ic0:ic1),Curves(1,ic0:ic1)], wait=1);
int NC= Curves(2,ic1)/hsize;
border G(t=0,1) {  P=Curve(Curves,ic0,ic1,t);  label= 1 + (x>xl)*2 + (y<yl);} 	

plot(G(-NC),wait=1); 
mesh Th=buildmesh(G(-NC));
plot(Th,wait=1);
real scale = sqrt(AreaLac/Th.area);
Th=movemesh(Th,[x*scale,y*scale]);//resize the mesh  
cout << " Th.area = " << Th.area << " Km^2 " << " == " << AreaLac <<  "   Km^2 " << endl ;
plot(Th,wait=1,ps="leman.eps");//  see figure \ref{fig:leman}
```


\twoplot[width=8cm]{lg}{leman}{\label{fig:lg}The image of the leman lac meshes}{\label{fig:leman} the mesh of lac}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{\setS{Finite Elements}} \index{finite element space}\label{finite elements}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

As stated in Section \ref{sec:example}.
FEM approximates all functions $w$ as
\[
w(x,y)\simeq w_0\phi_0(x,y)+w_1\phi_1(x,y)+\cdots+w_{M-1}\phi_{M-1}(x,y)
\]
with finite element basis functions $\phi_k(x,y)$ and numbers $w_k$ ($k=0,\cdots,M-1$).
The functions $\phi_k(x,y)$ are constructed from the triangle $T_{i_k}$, and called  _shape functions_.
In FreeFem++ the finite element space
$$
V_h=\left\{w\left|\; w_0\phi_0+w_1\phi_1+\cdots+w_{M-1}\phi_{M-1},\,
w_i\in \R\right.\right\}
$$
 is easily created by
```freefem
     fespace IDspace(IDmesh,<IDFE>) ;
```
or with $\ell$ pairs of periodic boundary condition in 2d
```freefem
     fespace IDspace(IDmesh,<IDFE>,
                      periodic=[[la$_1$,sa$_1$],[lb$_1$,sb$_1$],
                                ...
                                [la$_k$,sa$_k$],[lb$_k$,sb$_\ell$]]);
```
and in 3d
```freefem
     fespace IDspace(IDmesh,<IDFE>,
                      periodic=[[la$_1$,sa$_1$,ta$_1$],[lb$_1$,sb$_1$,tb$_1$],
                                ...
                                [la$_k$,sa$_k$,ta$_k$],[lb$_k$,sb$_\ell$,tb$_\ell$]]);
```

where

\index{fespace}\index{periodic}
\ttCC{IDspace} is the name of the space (e.g. \ttCC{Vh}),
\\\\
\ttCC{IDmesh} is the name of the associated mesh and  \ttCC{<IDFE>}
is a identifier of finite element type.
\\\\
In 2D we have a pair of periodic boundary condition, \label{periodic BC}
if \ttCC{[la$_i$,sa$_i$],[lb$_i$,sb$_i$]} is a pair of
`int}, and the 2 labels \ttCC{la$_i$} and \ttCC{lb$_i$}
refer to 2 pieces of boundary to be in equivalence.

If \ttCC{[la$_i$,sa$_i$],[lb$_i$,sb$_i$]} is a pair of `real},
then \ttCC{sa$_i$} and \ttCC{sb$_i$}
give two common abscissa on the two boundary curve, and two points are identified as one
if the two abscissa are equal.
\\\\
In 2D, we have a pair of periodic boundary condition,% \label{periodic BC}
if \ttCC{[la$_i$,sa$_i$,ta$_i$],[lb$_i$,sb$_i$,tb$_i$]} is a pair of
`int}, the 2 labels \ttCC{la$_i$} and \ttCC{lb$_i$}
define the 2 piece of boundary to be in equivalence.

If \ttCC{[la$_i$,sa$_i$,ta$_i$],[lb$_i$,sb$_i$,tb$_i$]} is a pair of `real},
then \ttCC{sa$_i$,ta$_i$} and \ttCC{sb$_i$,tb$_i$}
give two common parameters on the two boundary surface, and two points are identified as one
if the two parameters are equal.


\begin{remark} The 2D  mesh of the two identified borders must be the same, so
to be sure,  use  the parameter  `fixedborder=true} in `buildmesh} command (see \ref{buildmesh fixedborder})
like in example `periodic2bis.edp} (see \ref{exm:periodic4bis}).
\end{remark}

\medskip
 As of today, the known
types of finite element are: \index{type of finite element}
\begin{description}
     \item[P0,P03d]  piecewise constant discontinuous finite element  (2d, 3d), the degrees of freedom are  the barycenter element value.
     \index{P0|textbf}\index{fespace!P0}
    \begin{eqnarray}
    \label{eq:P0}
     P0_{h} = \left\{ v \in L^2(\Omega) \left|\; \textrm{for all }K \in \mathcal{T}_{h}\;\;\textrm{there is }\alpha_{K}\in \R :
        \;\; v_{|K} = \alpha_{K } \right.\right\}
     \end{eqnarray}
     \item[P1,P13d]  piecewise linear  continuous finite element (2d, 3d), the degrees of freedom are the vertices values.
     \index{P1|textbf}\index{fespace!P1}\index{fespace!P13d}
     \begin{eqnarray}
     &&P1_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{1} \right.\right\} \label{eq:P1}
     \end{eqnarray}
     \item[P1dc]  piecewise linear  discontinuous finite element
     \index{P1dc|textbf}\index{fespace!P1dc}
     \begin{equation}
     P1dc_{h} = \left\{ v \in L^{2}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{1} \right.\right\} \label{eq:P1dc}
     \end{equation}
     Warning, due to interpolation problem, the degree of freedom is not the vertices but three vectices  move
     inside with $T(X)= G + .99  (X-G) $ where $G$ is the barycenter, (version 2.24-4).
     \item[P1b,P1b3d]  piecewise linear  continuous finite element plus bubble (2d, 3d) \label{warP1dc}
     \index{P1b|textbf}\index{fespace!P1b}\index{fespace!P1b3d}

     \paragraph{The 2d case:}
     \begin{equation}
     P1b_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{1} \oplus \mathrm{Span}\{  \lambda^{K}_{0} \lambda^{K}_{1} \lambda^{K}_{2} \} \right.\right\} \label{eq:P1b}
     \end{equation}
     \paragraph{The 3d case:}
      \begin{equation}
     P1b_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{1} \oplus \mathrm{Span}\{  \lambda^{K}_{0} \lambda^{K}_{1} \lambda^{K}_{2} \lambda^{K}_{3} \} \right.\right\} \label{eq:P1b-3d}
     \end{equation}
    where $\lambda ^{K}_{i}, i=0,..,d$ are the $d+1$ barycentric  coordinate functions of the element  $K$ (triangle or tetrahedron).
    \item[P1bl,P1bl3d]  piecewise linear  continuous finite element plus linear bubble (2d, 3d) \label{warP1dc}
     \index{P1bl|textbf}\index{fespace!P1bl}\index{fespace!P1bl3d}
     the bulle is build by spliting the $K$ a barycenter in $d+1$ sub element. (need \ttCC{load "Element_P1bl"})
     
     \item[P2,P23d] piecewise $P_{2}$  continuous finite element (2d, 3d),
     \index{P2|textbf}\index{fespace!P2}
     \begin{equation}
     P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{2} \right.\right\}
     \end{equation}
     where
     $P_{2}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 2$.
    
     \item[P2b] piecewise $P_{2} $ continuous finite element  plus bubble,
     \index{P2|textbf}\index{fespace!P2}
     \begin{equation}
     P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{2} \oplus \mathrm{Span}\{  \lambda^{K}_{0} \lambda^{K}_{1} \lambda^{K}_{2} \} \right.\right\}
     \end{equation}

     \item[P2dc] piecewise $P_{2}$  discontinuous finite element,
     \index{P2dc|textbf}\index{fespace!P2dc}
     \begin{equation}
     P2dc_{h} = \left\{ v \in L^{2}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{2} \right.\right\}
     \end{equation}
    Warning, due to interpolation problem, the degree of freedom is not the six  P2 nodes  but six  nodes  move
     inside with $T(X)= G + .99  (X-G) $ where $G$ is the barycenter, (version 2.24-4).

      \item[P3] piecewise $P_{3}$  continuous finite element (2d)  (need \ttCC{load "Element_P3"})
     \index{P3|textbf}\index{fespace!P3}
     \begin{equation}
     P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{3} \right.\right\}
     \end{equation}
     where
     $P_{3}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 3$.
     
      \item[P3dc] piecewise $P_{3}$  discontinuous finite element (2d)  (need \ttCC{load "Element_P3dc"})
     \index{P3dc|textbf}\index{fespace!P3dc}
     \begin{equation}
     P2_{h} = \left\{ v \in L^2(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{3} \right.\right\}
     \end{equation}
     where
     $P_{3}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 3$.

      \item[P4] piecewise $P_{4}$  continuous finite element (2d)  (need \ttCC{load "Element_P4"})
     \index{P4|textbf}\index{fespace!P4}
     \begin{equation}
     P2_{h} = \left\{ v \in H^{1}(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{4} \right.\right\}
     \end{equation}
     where
     $P_{4}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 4$.
     
      \item[P4dc] piecewise $P_{4}$  discontinuous finite element (2d)  (need \ttCC{load "Element_P4dc"})
     \index{P4dc|textbf}\index{fespace!P4dc}
     \begin{equation}
     P2_{h} = \left\{ v \in L^2(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{3} \right.\right\}
     \end{equation}
     where
     $P_{4}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 3$.

     \item[P0Edge] piecewise $P_{0}$  discontinuous finite element (2d)  contant on each edge  of the mesh.
     \index{P0Edge|textbf}\index{fespace!P0Edge}

     \item[P1Edge] piecewise $P_{1}$  discontinuous finite element (2d)  (need \ttCC{load "Element_PkEdge"})  $P_1$ on each edge of the mesh.
     \index{P1Edge|textbf}\index{fespace!P1Edge}
     \item[P2Edge] piecewise $P_{2}$  discontinuous finite element (2d)  (need \ttCC{load "Element_PkEdge"})  $P_2$ on each edge of the mesh.
     \index{P2Edge|textbf}\index{fespace!P2Edge}
     \item[P3Edge] piecewise $P_{3}$  discontinuous finite element (2d)  (need \ttCC{load "Element_PkEdge"})  $P_3$ on each edge of the mesh.
     \index{P1Edge|textbf}\index{fespace!P3Edge}
     \item[P4Edge] piecewise $P_{4}$  discontinuous finite element (2d)  (need \ttCC{load "Element_PkEdge"})  $P_4$ on each edge of the mesh.
     \index{P4Edge|textbf}\index{fespace!P4Edge}
    \item[P5Edge] piecewise $P_{5}$  discontinuous finite element (2d)  (need \ttCC{load "Element_PkEdge"})  $P_5$ on each edge of the mesh.
     \index{P5Edge|textbf}\index{fespace!P4Edge}


     \item[P2Morley] piecewise $P_{2}$  non conform finite element (2d)  (need \ttCC{load "Morley"})
     \index{P2Morley|textbf}\index{fespace!P2Morley}\index{Morley}
     \begin{equation}
     P2_{h} = \left\{ v \in L^2(\Omega) \left|\; \forall K \in \mathcal{T}_{h}
        \quad v_{|K} \in P_{3}, 
        \left\{\begin{array}{c} 
        v \mbox{ continuous  at vertices,}\\
        \p_n{v} \mbox{ continuous  at middle of edge,} 
        \end{array}\right. 
         \right.\right\}
     \end{equation}
     where
     $P_{2}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 2$.

      Warning to build the interplant of a function $u$ (scalar) for this  finite element,
       we need the function and 2 partial derivatives $(u,u_x, u_y)$,
      so  this vectorial finite element with 3 components  $(u,u_x,u_y)$.
      
      See example `bilapMorley.edp} of \verb!examples++-load! for solving BiLaplacien problem : 
      \index{BiLaplacien}
```freefem
         load "Morley" 
         fespace Vh(Th,P2Morley);      // the Morley finite element space
         macro bilaplacien(u,v) ( dxx(u)*dxx(v)+dyy(u)*dyy(v)+2.*dxy(u)*dxy(v)) // fin macro 
         real f=1;
         Vh [u,ux,uy],[v,vx,vy];

         solve bilap([u,ux,uy],[v,vx,vy]) =
             int2d(Th)(  bilaplacien(u,v) )
           - int2d(Th)(f*v)
           + on(1,2,3,4,u=0,ux=0,uy=0)      
```

      \item[HCT]  $P_3$ $C^1$ conform finite element (2d)  (need \ttCC{load "Element_HCT"})  one 3 sub triangles (version 3.40). 
     \index{HCT|textbf}\index{fespace!HCT}
     
     Let call $\mathcal{T}^\triangle_{h}$ the sub mesh of $\mathcal{T}_{h}$ where all triangle are split in 3 at the
     barycenter.
     \begin{equation}
       PHCT_{h} = \left\{ v \in C^1(\Omega) \left|\; \forall K \in \mathcal{T}^\triangle_{h}
        \quad v_{|K} \in P_{3} \right.\right\} 
     \end{equation}
     where
     $P_{3}$ is the set of polynomials of $\R^{2}$ of  degrees $\le 3$. The degree of freedom are 
     the value and derivative at vertices and normal derivative a middle edge point of initial meshes,and thank to \cite{HCT}.
     
         Warning to build the interplant of a function $u$ (scalar) for this  finite element,
       we need the function and 2 partial derivatives $(u,u_x, u_y)$,
      so  this vectorial finite element with 3 components  $(u,u_x,u_y)$ like in previous Finite Element.
        

      \item[P2BR]  \index{P2BR|textbf}\index{fespace!P2BR}(need \ttCC{load "BernadiRaugel"})  the Bernadi Raugel Finite Element is a Vectorial   element (2d)  with 2 components,
 See Bernardi, C., Raugel, G.: Analysis of some finite elements for the Stokes problem. Math. Comp. 44, 71-79 (1985).
 It is  a 2d coupled FE, with 
 the Polynomial space is $ P_1^2$ + 3 normals bubbles edges function $(P_2)$
and  the degre of freedom is 6 values at of the $2$ components at the  $3$ vertices
and the $3$ flux on the $3$ edges  
So the number  degrees of freedom is 9.      
      
      \item[RT0,RT03d]  Raviart-Thomas finite element of degree $0$.
     \index{RT0|textbf}\index{fespace!RT0}

     The 2d case:
     \begin{equation}
         RT0_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{x}{y}  \right.\right\}
         \label{eq:RT0}
     \end{equation}
     The 3d case:
     \begin{equation}
         RT0_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \mathbf{v}_{|K}(x,y,z) =
         \vectthree{\alpha^1_{K}}{\alpha^2_{K}}{\alpha^3_{K}} + \beta_{K}\vectthree{x}{y}{z}  \right.\right\}
        \label{eq:RT03d}
     \end{equation}
      where by writing
      $\textrm{div }\mathbf{w}=\sum_{i=1}^d\p w_i/\p x_i$ with
      $ \mathbf{w}=(w_i)_{i=1}^d$,
      $$
      H(\textrm{div})=\left\{\mathbf{w}\in L^{2}(\Omega)^d\left|
      \textrm{div } \mathbf{w}\in L^{2}(\Omega)
      \right.\right\}
      $$
      and where
      $\alpha^1_{K},\alpha^2_{K},\alpha^3_{K},\beta_{K} $ are real numbers.
      
   \item[RT0Ortho]  Raviart-Thomas Orthogonal, or Nedelec finite element type I of degree $0$ in dimension 2
     \index{RT0Ortho|textbf}\index{fespace!RT0Ortho}
     \begin{equation}
         RT0Ortho{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \mathbf{v}_{|K}(x,y) =
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{-y}{x}  \right.\right\}
         \label{RT0Ortho}
     \end{equation}      
      
     \item[Edge03d]  3d Nedelec finite element or Edge  Element of degree $0$.
     \index{Edge03d|textbf}\index{fespace!Edge03d}

     The 3d case:
     \begin{equation}
         Edge0_{h} = \left\{ \mathbf{v} \in H(\textrm{Curl}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \mathbf{v}_{|K}(x,y,z) =
         \vectthree{\alpha^1_{K}}{\alpha^2_{K}}{\alpha^3_{K}} + \vectthree{\beta^1_{K}}{\beta^2_{K}}{\beta^3_{K}}\times\vectthree{x}{y}{z}  \right.\right\}
         \label{eq:Edge03d}
     \end{equation}
      where by writing
      $\textrm{curl}\mathbf{w}=\vectthree{\p w_2/\p x_3-\p w_3/\p x_2}{\p w_3/\p x_1-\p w_1/\p x_3}{\p w_1/\p x_2-\p w_2/\p x_1}$ with
      $ \mathbf{w}=(w_i)_{i=1}^d$,
      $$
      H(\textrm{curl})=\left\{\mathbf{w}\in L^{2}(\Omega)^d\left|
      \textrm{curl } \mathbf{w}\in L^{2}(\Omega)^d
      \right.\right\}
      $$
      and
      $\alpha^1_{K},\alpha^2_{K},\alpha^3_{K},\beta^1_{K},\beta^2_{K},\beta^3_{K} $ are real numbers.
      \item[Edge13d] (need \ttCC{load "Element_Mixte3d"}, version 3.34) 3d Nedelec finite element or Edge  Element of degree $1$. \index{Edge13d|textbf}\index{fespace!Edge13d}

      \item[Edge23d]  (need \ttCC{load "Element_Mixte3d"}, version 3.34) 3d Nedelec finite element or Edge  Element of degree $2$. \index{Edge23d|textbf}\index{fespace!Edge23d}

      
     \item[P1nc] \index{P1nc|textbf}\index{fespace!P1nc} piecewise linear   element continuous at
     the middle of edge only in 2D (Crouzeix-Raviart Finite Element 2d).
     
         \item[P2pnc] \index{P2pnc|textbf}\index{fespace!P2pnc} piecewise quadratic plus  a bubble  P3    element with the 
     continuity of  the 2 moments on each edge (version 3.59)   (need \ttCC{load "Element_P2pnc"}
 
    \item[RT1] \index{RT1|textbf}\index{fespace!RT1} (need \ttCC{load "Element_Mixte"}, version 3.13)
     \begin{equation}
         RT1_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad  \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_1^2,P_0,  \mathbf{v}_{|K}(x,y) = 
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{x}{y}   \right.\right\}
         \label{eq:RT1}
     \end{equation}
    \item[RT1Ortho] \index{RT1Ortho|textbf}\index{fespace!RT1Ortho} (need \ttCC{load "Element_Mixte"}, version 3.13, dimension 2)
         \begin{equation}
         RT1_{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h},  \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_1^2,P_0,  \mathbf{v}_{|K}(x,y) = 
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{-y}{x}   \right.\right\}
         \label{eq:RT1Ortho}
     \end{equation}
     
       \item[RT2] \index{RT2|textbf}\index{fespace!RT2} (need \ttCC{load "Element_Mixte"}, version 3.59
     \begin{equation}
         RT2_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad   \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_2^2, P_1,  \mathbf{v}_{|K}(x,y) = 
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{x}{y}   \right.\right\}
         \label{eq:RT2}
     \end{equation}
   \item[RT2Ortho] \index{RT2Ortho|textbf}\index{fespace!RT2Ortho} (need \ttCC{load "Element_Mixte"}, version 3.59, dimension 2)
         \begin{equation}
         RT2_{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h} ,  \alpha^1_{K}, \alpha^2_{K}, \beta_{K} \in P_2^2, P_1,  \mathbf{v}_{|K}(x,y) = 
         \vecttwo{\alpha^1_{K}}{\alpha^2_{K}} + \beta_{K}\vecttwo{-y}{x}   \right.\right\}
         \label{eq:RT1Ortho}
     \end{equation}
    


    \item[BDM1] \index{BDM1|textbf}\index{fespace!BDM1} (need \ttCC{load "Element_Mixte"}, version 3.13, dimension 2) the Brezzi-Douglas-Marini finite element 
     \begin{equation}
         BDM1_{h} = \left\{ \mathbf{v} \in H(\textrm{div}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad   \mathbf{v}_{|K} \in P_1^2
         \right.\right\}
         \label{eq:BDM1}
     \end{equation}
        
    \item[BDM1Ortho] \index{BDM1Ortho|textbf}\index{fespace!BDM1Ortho} (need \ttCC{load "Element_Mixte"}, version 3.13, dimension 2) the Brezzi-Douglas-Marini Orthogonal also call
    Nedelec of type II , finite element 
       \begin{equation}
         BDM1Ortho_{h} = \left\{ \mathbf{v} \in H(\textrm{curl}) \left|\; \forall K \in
         \mathcal{T}_{h} \quad   \mathbf{v}_{|K} \in P_1^2
         \right.\right\}
         \label{eq:BDM1Ortho}
     \end{equation}
   \item[FEQF]  \index{FEQF|textbf}\index{fespace!FEQF} (need \ttCC{load "Element_QF"}, version 3.45, dimension 2 or 3) the finite element to store 
   function at default quadrature points (so the  quadrature is `qf5pT} in 2d and is   `qfV5} in 3d).
   
   for over quadrature you have the following correspondance  finite element, quadrature formula. 
    \begin{itemize}
      \item `FEQF1} $\mapsto$ `qf1pT}  , \index{FEQF1|textbf}\index{fespace!FEQF1}
      \item `FEQF2} $\mapsto$ `qf2pT}  , \index{FEQF2|textbf}\index{fespace!FEQF2}
      \item `FEQF5} $\mapsto$ `qf5pT}  , \index{FEQF5|textbf}\index{fespace!FEQF5}
      \item `FEQF7} $\mapsto$ `qf7pT}  , \index{FEQF7|textbf}\index{fespace!FEQF7}
      \item `FEQF9} $\mapsto$ `qf9pT}  , \index{FEQF9|textbf}\index{fespace!FEQF9}
      \item `FEQF13d} $\mapsto$ `qfV1}  , \index{FEQF13d|textbf}\index{fespace!FEQF13d}
      \item `FEQF23d} $\mapsto$ `qfV2}  , \index{FEQF23d|textbf}\index{fespace!FEQF23d}
      \item `FEQF53d} $\mapsto$ `qfV5}   \index{FEQF53d|textbf}\index{fespace!FEQF53d}
    \end{itemize}
   You can use this element element of do optimization to store and reuse  function with long formula in non linear process in integral. 
    
     \end{description}
