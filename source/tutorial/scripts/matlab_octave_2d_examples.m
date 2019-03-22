%Add ffmatlib to the search path
addpath('ffmatlib');

%Load the mesh
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('capacitor.msh');
%Load the finite element space connectivity
vh=ffreaddata('capacitor_vh.txt');
%Load scalar data
u=ffreaddata('capacitor_potential.txt');
%Load 2D vector field data
[Ex,Ey]=ffreaddata('capacitor_field.txt');

ffpdeplot(p,b,t,'Mesh','on','Boundary','on');

figure;
ffpdeplot(p,b,t,'VhSeq',vh,'XYData',u,'Mesh','on','Boundary','on', ...
          'XLim',[-2 2],'YLim',[-2 2]);

figure;
ffpdeplot(p,b,t,'VhSeq',vh,'XYData',u,'ZStyle','continuous', ...
          'Mesh','off');
lighting gouraud;
view([-47,24]);
camlight('headlight');

figure;
ffpdeplot(p,b,t,'VhSeq',vh,'XYData',u,'Mesh','off','Boundary','on', ...
          'XLim',[-2 2],'YLim',[-2 2],'Contour','on','CColor','b', ...
          'XYStyle','off', 'CGridParam',[150, 150],'ColorBar','off', ...
          'FlowData',[Ex,Ey],'FGridParam',[24, 24]);

