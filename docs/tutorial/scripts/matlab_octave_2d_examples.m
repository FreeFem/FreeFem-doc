%Plot of a R2->R function

clear all;

addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,C]=ffread2patch('temp_demo1_getstarted.txt','Delimiter',';', ...
                     'Format','%f %f %f');

%Plots the facets and display the mesh
%Sets the the color C equal to the solution of the PDE (Z-value)
patch(X,Y,C,'EdgeColor',[0 0 0],'LineWidth',1);

%Creates a colorbar
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'C');

ylabel('y');
xlabel('x');
title('2D Density Plot');

%Sets the view point specification to 2d
view(2);
%Sets 1:1 aspect ratio
axis tight equal;

pause(2);
close all;
%%%%%% Hide the mesh

clear all;

addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,C]=ffread2patch('temp_demo1_getstarted.txt','Delimiter',';', ...
                     'Format','%f %f %f');

%Plots the facets and display the mesh
%No Edge Color (hide the mesh)
patch(X,Y,C,'EdgeColor','none');

%Creates a colorbar
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'C');

ylabel('y');
xlabel('x');
title('2D Density Plot wo Mesh');

%Sets the view point specification to 2d
view(2);
%Sets 1:1 aspect ratio
axis tight equal;

pause(2);
close all;
%%%%%% Show only the mesh

clear all;

addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,C]=ffread2patch('temp_demo1_getstarted.txt','Delimiter',';', ...
                     'Format','%f %f %f');

%Plots the facets and display the mesh
%Plots white facets with blue edge color in order to show the mesh
patch(X,Y,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);

ylabel('y');
xlabel('x');
title('2D Mesh');

%Sets the view point specification to 2d
view(2);
%Sets 1:1 aspect ratio
axis tight equal;

pause(2);
close all;
%%%%%% Create a 3d surf

clear all;

addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,C]=ffread2patch('temp_demo1_getstarted.txt','Delimiter',';', ...
                     'Format','%f %f %f');

%Plots the facets and display the mesh
%Sets the the color C equal to the solution of the PDE (Z-value)
patch(X,Y,C,C,'EdgeColor',[0 0 0],'LineWidth',1);

%Creates a colorbar
colormap(jet(192));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'C');

ylabel('y');
xlabel('x');
title('3D Surf Plot');

%Sets the view point specification to 2d
view(3);
%Sets 1:1:1 aspect ratio
daspect([1 1 1*(max(max(C))-min(min(C)))]);
grid;

pause(2);
close all;
