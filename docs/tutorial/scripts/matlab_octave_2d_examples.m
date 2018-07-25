%Where to find the ffmatlib commands
addpath('ffmatlib');

%Loads the mesh
[p,b,t,nv,nbe,nt,labels]=ffreadmesh('capacitor.msh');
%Loads scalar data
[u]=ffreaddata('capacitor_potential.txt');
%Loads vector field data
[Ex,Ey]=ffreaddata('capacitor_field.txt');

ffpdeplot(p,b,t,'Mesh','on','Boundary','on');

figure;
ffpdeplot(p,b,t,'XYData',u,'Mesh','on','Boundary','on', ...
          'XLim',[-2 2],'YLim',[-2 2]);

figure;
ffpdeplot(p,b,t,'XYData',u,'ZStyle','continuous','Mesh','off');
lighting gouraud;
view([-47,24]);
camlight('headlight');

figure;
ffpdeplot(p,b,t,'XYData',u,'XYStyle','off','Mesh','off','Boundary','on', ...
          'Contour','on','CStyle','monochrome','CColor','b', ...
          'CGridParam',[150, 150],'FlowData',[Ex,Ey],'FGridParam',[60, 60], ...
          'ColorBar','off','XLim',[-2 2],'YLim',[-2 2]);

