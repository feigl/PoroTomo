function S = triangulatexyz4(X,Y,Z)
%% Generate kfaces of a fault by projecting onto a plane and then triangulating
%function F = triangulatexyz(XYZ,doplots)
% inputs:
%    XYZ        : npoints-by-3 matrix of [X,Y,Z] coordinates
%    titlestr   : if nonzero, then make plots with this title
% output: F stucture with fields:
%    kfaces   : ntriangles-by-3 matrix of indices to points
%    kboundary: vector of indices to points on edge
%
%    derived from: demo for affine_fit by Author: Adrien Leygue
% 20170420, Kurt Feigl
%% TODO return vertices

%% intialize
% indices for figures
nf = 0;


%% check sizes of arrays
X = colvec(X);nx=numel(X);
Y = colvec(Y);ny=numel(Y);
Z = colvec(Z);nz=numel(Z);

if nx ~= ny || nx ~= nz || ny ~= nz
    nx
    ny
    nz
    error('input vectors are not the same sizes');
else
    npoints = nx;
end

if npoints < 3
    npoints;
    %warning('Too few points');
    S.area = 0;
    return;
end

%% find the best-fitting plane and calculate the modeled values of ZM = f(XM,YM) on a mesh
%compute the normal to the plane and a point that belongs to the plane
% normvec : a unit (column) vector normal to the plane
% basisvecs : a 3 by 2 matrix. Its columns form an orthonormal basis of the plane
% point0 : a point belonging to the plane
XYZ = [X,Y,Z];
% whos XYZ
[normvec,basisvecs,point0] = affine_fit(XYZ);

% calculate projections of points into UU,VV plane
UU=nan(npoints,1);
VV=nan(npoints,1);
kount = 0;
for i=1:npoints
    if (sum(isfinite(colvec(basisvecs))) == numel(basisvecs)) && (sum(isfinite(colvec(XYZ(i,1:3)))) == numel((XYZ(i,1:3))))
        kount = kount+1;
        UU(kount) = dot(colvec(XYZ(i,1:3)),colvec(basisvecs(1:3,1)));
        VV(kount) = dot(colvec(XYZ(i,1:3)),colvec(basisvecs(1:3,2)));
    end
end
UU = colvec(UU);
VV = colvec(VV);


% %% find convex hull
% [Khull,S.area] = convhull(UU,VV);
% 
% %% Calculate a triangular mesh using Delaunay triangulation
% TRI = delaunayTriangulation(UU(Khull),VV(Khull));
% 
% [nvertices,n2] = size(TRI.Points)
% S.vertices = nan(nvertices,3);
% 
% for i=1:nvertices
%     k = Khull(i);
%     S.vertices(i,1) = X(k);
%     S.vertices(i,2) = Y(k);    
%     S.vertices(i,3) = Z(k);
% end
% S.faces = TRI.ConnectivityList;
% [nfaces,ncols] = size(S.faces)


%% Calculate a triangular mesh using Delaunay triangulation
warning('off');
TRI = delaunayTriangulation(UU,VV);

S.area = polyarea(UU,VV);

[nvertices,n2] = size(TRI.Points);
S.vertices = nan(nvertices,3);

for i=1:nvertices
    S.vertices(i,1) = X(i);
    S.vertices(i,2) = Y(i);    
    S.vertices(i,3) = Z(i);
end
S.faces = TRI.ConnectivityList;
warning('on');

return
end






