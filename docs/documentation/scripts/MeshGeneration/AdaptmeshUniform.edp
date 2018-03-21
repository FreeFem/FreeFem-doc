mesh Th=square(2, 2); //the initial mesh
plot(Th, wait=true, ps="square-0.eps");

Th = adaptmesh(Th, 1./30., IsMetric=1, nbvx=10000);
plot(Th, wait=true, ps="square-1.eps");

Th = adaptmesh(Th, 1./30., IsMetric=1, nbvx=10000); //More the one time du to
Th = adaptmesh(Th, 1./30., IsMetric=1, nbvx=10000); //Adaptation bound `maxsubdiv=`
plot(Th, wait=true, ps="square-2.eps");
