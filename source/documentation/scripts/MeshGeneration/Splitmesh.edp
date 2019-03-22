// Mesh
border a(t=0, 2*pi){x=cos(t); y=sin(t); label=1;}
mesh Th = buildmesh(a(20));
plot(Th, wait=true, ps="NotSplittedMesh.eps");

// Splitmesh
Th = splitmesh(Th, 1 + 5*(square(x-0.5) + y*y));
plot(Th, wait=true, ps="SplittedMesh.eps");
