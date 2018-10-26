% make facets for a cube

%https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html

% make a cube
side1 = [0, 1];
[Xm,Ym,Zm] = meshgrid(side1,side1,side1);
x1 = colvec(Xm);
y1 = colvec(Ym);
z1 = colvec(Zm)
% make PoroTomo Box
[Xm,Ym,Zm] = meshgrid([0,500],[0,1500],[0,400]);
x1 = colvec(Xm);
y1 = colvec(Ym);
z1 = colvec(Zm)

% % make a cube
% side2 = [0,1];
% [Xm,Ym,Zm] = meshgrid(side2,side2,side2);
% x2 = colvec(Xm);
% y2 = colvec(Ym);
% z2 = colvec(Zm)

% % try to make 3 coplanar triangles from a triangle plus its centroid
% x2 = colvec([1 0 0 1]);
% y2 = colvec([0 1 0 1]);
% z2 = colvec([0 0 1 1]);

% % make a sphere
% theta = gallery('uniformdata',[100,1],0)*2*pi;
% phi = gallery('uniformdata',[100,1],1)*pi;
% x2 = cos(theta).*sin(phi);
% y2 = sin(theta).*sin(phi);
% z2 = cos(phi);

% use faults
x2 = FAULTS(10).PR.x;
y2 = FAULTS(10).PR.y;
z2 = FAULTS(10).PR.z;
% iok=1:numel(x2);
% iok=intersect(iok,find(x2>=0));
% iok=intersect(iok,find(x2<500));
% iok=intersect(iok,find(y2>=0));
% iok=intersect(iok,find(y2<1500));
% iok=intersect(iok,find(z2>=0));
% iok=intersect(iok,find(z2<450));
% x2=x2(iok);
% y2=y2(iok);
% z2=z2(iok);




% get the triangulation
DT1 = delaunayTriangulation(x1,y1,z1);

%  delaunayTriangulation with properties:
% 
%               Points: [8×3 double]
%     ConnectivityList: [6×4 double]
%          Constraints: []

%DT2 = delaunayTriangulation(x2,y2,z2);
% DT2.Points = [FAULTS(10).PR.x,FAULTS(10).PR.y,FAULTS(10).PR.z];
% DT2.ConnectivityList = FAULTS(10).TRI.ConnectivityList;


% Find the free boundary facets of the triangulation, and use them to create a 2-D triangulation on the surface.
[T1,Xb1] = freeBoundary(DT1);
TR1 = triangulation(T1,Xb1);
% [T2,Xb2] = freeBoundary(DT2);
% TR2 = triangulation(T2,Xb2);
% T2=FAULTS(10).TRI.ConnectivityList;
% TR2 = triangulation(T2,x2,y2,z2);

% Compute the centers and face normals of each triangular facet in TR.
P1 = incenter(TR1);
F1 = faceNormal(TR1);
% P2 = incenter(TR2);
% F2 = faceNormal(TR2);



%function [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2, varargin)
%SURFACEINTERSECTION intersection of 2 surfaces
% [intMatrix, intSurface] = SurfaceIntersection(surface1, surface2)
% calculates the intersection of surfaces 1 and 2. Code can either return
% just the matrix indicating which face of surface1 intersected with face
% of surface2, which is calculated using Tomas Moller algorithm, or can
% also return the actual line of intersection. In case when parts of the
% surface 1 and 2 lay on the same plane the intersection is a 2D area
% instead of 1D edge. In such a case the intersection area will be
% triangulated and intSurface.edges will hold the edges of the
% triangulation surface and intSurface.faces will hold the faces.

surface1.vertices=[x1,y1,z1];
surface1.faces   =T1;
surface1
surface2.vertices=[x2,y2,z2];
surface2.faces   = FAULTS(10).TRI.ConnectivityList;
surface2

[intMatrix, intSurface] = SurfaceIntersection(surface1,surface2);

% surface1 = DT1;
% surface2 = DT2;
% isa(surface1, 'triangulation')
% [intMatrix, intSurface] = SurfaceIntersectionKF(surface1,surface2);


% Plot the triangulation along with the centers and face normals.

figure;
trisurf(T1,Xb1(:,1),Xb1(:,2),Xb1(:,3),'FaceColor','blue',  'FaceAlpha',0.5);
axis equal;
hold on;
%trisurf(T2,Xb2(:,1),Xb2(:,2),Xb2(:,3),'FaceColor','red','FaceAlpha',0.5);

trisurf(intSurface.faces,intSurface.vertices(:,1),intSurface.vertices(:,2),intSurface.vertices(:,3),'FaceColor','green','FaceAlpha',0.5);

% quiver3(P1(:,1),P1(:,2),P1(:,3),F1(:,1),F1(:,2),F1(:,3),0.5,'color','blue');
% quiver3(P2(:,1),P2(:,2),P2(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'color','red');
plot3(intSurface.vertices(:,1),intSurface.vertices(:,2),intSurface.vertices(:,3),'g*-');


