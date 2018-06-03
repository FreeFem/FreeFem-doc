%Plot of a R3->R function

clear all;

addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,Z,C] = ffread2patch('temp_demo4_bddata3d_box.txt', ...
                         'Delimiter',';','Format','%f %f %f %f');

%Plots the facets and display the mesh
%Sets the the color C equal to the solution of the PDE
patch(X,Y,Z,C,'EdgeColor',[0 0 0],'LineWidth',1);

%Creates a colorbar
colormap(jet(250));
caxis([min(min(C)) max(max(C))]);
hcb=colorbar;
title(hcb,'C');

zlabel('z');
ylabel('y');
xlabel('x');
title('3D Plot of a Surface Boundary');

%Sets the view point specification to 3d
view(3);
%Sets 1:1:1 aspect ratio
daspect([1 1 1]);

pause(2);
close all;
%%%%%% Hide the mesh

clear all;

addpath('ffmatlib');

%Reads the file content and converts to patch() - plot data
[X,Y,Z,C] = ffread2patch('temp_demo4_bddata3d_box.txt', ...
                         'Delimiter',';','Format','%f %f %f %f');

%Plots the facets and display the mesh
%Sets the the color C equal to the solution of the PDE
patch(X,Y,Z,[1 1 1],'EdgeColor',[0 0 1],'LineWidth',1);

zlabel('z');
ylabel('y');
xlabel('x');
title('3D Mesh');

%Sets the view point specification to 3d
view(3);
%Sets 1:1:1 aspect ratio
daspect([1 1 1]);

pause(2);
close all;
%%%%%% Hide the mesh

clear all;

addpath('ffmatlib');

%Reads boundary and mesh data from within the domain
[bdata,tdata] = ffreadfile('File1','temp_demo4_bddata3d_box.txt', ...
                           'File2','temp_demo4_tetdata3d_box.txt', ...
                           'Delimiter',';','Format','%f %f %f %f');

%Slicing plane definition
S1=[0 0 0]';
S2=[1 0.5 0]';
S3=[1 1.3 1]';

%Cut boundary as well as the complete domain
[BX,BY,BZ,BC] = slicebd2patch(bdata,S1,S2,S3);
[SX,SY,SZ,SC] = slicetet2patch(tdata,S1,S2,S3);

%Plot
patch([SX BX],[SY BY],[SZ BZ],[SC BC]);
%Adjust to correct colors
colormap(jet(250));
caxis([min(min([SC BC])) max(max([SC BC]))]);
hcb=colorbar;
title(hcb,'C');

zlabel('z');
ylabel('y');
xlabel('x');
title('3D Cross Section');

%Sets the view point specification to 3d
view(3);
%Sets 1:1:1 aspect ratio
daspect([1 1 1]);

pause(2);
close all;

