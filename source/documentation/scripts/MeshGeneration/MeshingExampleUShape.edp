real d = 0.1; //width of U-shape
border L1(t=0, 1-d){x=-1; y=-d-t;}
border L2(t=0, 1-d){x=-1; y=1-t;}
border B(t=0, 2){x=-1+t; y=-1;}
border C1(t=0, 1){x=t-1; y=d;}
border C2(t=0, 2*d){x=0; y=d-t;}
border C3(t=0, 1){x=-t; y=-d;}
border R(t=0, 2){x=1; y=-1+t;}
border T(t=0, 2){x=1-t; y=1;}
int n = 5;
mesh Th = buildmesh(L1(n/2) + L2(n/2) + B(n) + C1(n) + C2(3) + C3(n) + R(n) + T(n));
plot(Th, ps="U-shape.eps", bw=true);
