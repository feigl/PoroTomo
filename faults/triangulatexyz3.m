function [TRI,lgood,X,Y,Z] = triangulatexyz3(X,Y,Z,titlestr,make_plots)
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
TRI=struct([]);
lgood = false;


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
    npoints = nx
end

if npoints < 3
    npoints
    warning('Too few points');
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
clear XYZ
save('UUVV.mat','UU','VV');


% remove empty rows
if kount ~= npoints
    warning('pruning');
    UU = UU(1:kount);
    VV = VV(1:kount);
end

%% remove duplicate rows
%    [C,IA,IC] = unique(A) also returns index vectors IA and IC such that
%     C = A(IA) and A = C(IC) (or A(:) = C(IC), if A is a matrix or array).  
[UUVV,iasort,icsort] = unique([UU,VV],'rows');
UU=UU(iasort);
VV=VV(iasort);
X=X(iasort);
Y=Y(iasort);
Z=Z(iasort);

if numel(UU) < 3 || numel(VV) < 3
    UU
    VV
    warning('miscount')
    return
end
    


% theta = nan(3,1);

%% find the boundary in UU,VV plane using fancy routine named alpha shape
%shp = alphaShape(UU,VV);
% choose search radius based on scatter of coordinates
%alpha_radius = 0.1*max([nanstd(UU),nanstd(VV)]); % in meters
%alpha_radius = 2*max([nanstd(UU),nanstd(VV)]); % in meters
% choose search radius based on spacing coordinates
%alpha_radius = 20*max([mean(colvec(diff(UU))),mean(colvec(diff(VV)))]) % in meters
% choose search radius based on extent
alpha_radius = 0.1*min([max(colvec(UU))-min(colvec(UU)),max(colvec(VV))-min(colvec(VV))]); % in meters
shp = alphaShape(UU,VV,alpha_radius);
Kboundary = boundaryFacets(shp);
% find boundary in U,V plane
Ubound = [shp.Points(Kboundary(:,1),1) ; shp.Points(Kboundary(:,2),1)];
Vbound = [shp.Points(Kboundary(:,1),2) ; shp.Points(Kboundary(:,2),2)];

if make_plots == 1
    %% plot points and alpha shape in U-V plane
    nf=nf+1;figure;hold on;
    % plot shape, including vertices and edges
    %plot(shp);
    %plot vertices only
    plot(shp.Points(:,1),shp.Points(:,2),'r*');
    % plot boundary of alpha shape
    plot(Ubound,Vbound,'b-','LineWidth',3);
    title(sprintf('%s\n%s',titlestr,'nNonconvex Alpha Shape'));
end


%%  make a mesh for plotting plane
if make_plots == 1
%     stepsize = 10 % in meters
%     [XM,YM] = meshgrid(  ...
%          [nanmin(X):stepsize:nanmax(X)]...
%         ,[nanmin(Y):stepsize:nanmax(Y)]);
    nsteps = 10; % 
    [XM,YM] = meshgrid(  ...
         linspace(nanmin(X),nanmax(X),nsteps)...
        ,linspace(nanmin(Y),nanmax(Y),nsteps));
    
    %%  evalute the formula for the plane
    ZM = -1*(normvec(1)/normvec(3)*XM+normvec(2)/normvec(3)*YM-dot(normvec,point0)/normvec(3));
end


%% Calculate a triangular mesh using Delaunay triangulation
% use older routine
%Kfaces = delaunay(UU,VV);
% [ntriangles1,ncols1]= size(Kfaces1)

% use newer routine
TRI = delaunayTriangulation(UU,VV);
% use boundary of alpha shape as edge constraint
% TRI = delaunayTriangulation(UU,VV,Kboundary);

if numel(TRI.Points(:,1)) ~= numel(UU)
    error('miscount')
    return
    % Algorithm has added points, update coordinates
    UU = colvec(TRI.Points(:,1));
    VV = colvec(TRI.Points(:,2));
    npoints = numel(UU)
end

Kfaces = TRI.ConnectivityList;
% gives different answers
% https://www.mathworks.com/matlabcentral/answers/252870-why-does-delaunay-sometimes-give-a-different-result-from-delaunaytriangulation
[ntriangles,ncols]= size(Kfaces)


%% evaluate each triangle individually
lgood = false(ntriangles,1);
for i=1:ntriangles
    nok = 0;
    k3 = Kfaces(i,[1:3]);

%     % test if all three vertices of triangle are inside alpha shape
%     for j=1:3
%         k = Kfaces(i,j); % index of vertices in triangle       
%         if inpolygon(UU(k),VV(k),Ubound,Vbound) == 1
%             nok=nok+1;
%         end
%     end

    % test if centroid of triangle is inside alpha shape
    umean=mean(UU(k3));
    vmean=mean(VV(k3));
    if inpolygon(umean,vmean,Ubound,Vbound) == 1
        nok=nok+1;
    end
    
%     % calculate area of triangle
%     area = polyarea(UU(k3),VV(k3)); % in square meters
%     if area > (25.)^2
%         nok = nok+1;
%     end
   
%     % We can use the aspect ratio (AR) of the triangle as a rejection
%     % criterion; the AR being defined as the ratio of the incircle radius to
%     % circumcircle radius.
%     radius_min = min(sqrt((UU(k3)-umean).^2 + (VV(k3)-vmean).^2)); 
%     radius_max = max(sqrt((UU(k3)-umean).^2 + (VV(k3)-vmean).^2));
%     aspect_ratio = radius_min/radius_max; 
%     if aspect_ratio < 0.9 
%         nok=nok+1;
%     end
    
    % must pass all tests
%    if nok==6 
    if nok==1 
        lgood(i) = true;
    end
end

%% prune
igood = find(lgood==true);
ngood = numel(igood)

if make_plots == 1
    %% make a 3D plot
    nf=nf+1;figure;
    plot3(X,Y,Z,'r.');
    hold on;
    axis tight;
    axis equal;
    % plot the best-fitting plane
    surf(XM,YM,ZM,'facecolor','blue','facealpha',0.5);
    xlabel('X');
    xlabel('Y');
    ylabel('Z');
    title(sprintf('%s\n%s',titlestr,'3 D')); 
       
    %% plot points and alpha shape in U-V plane
    nf=nf+1;figure;hold on;
    % plot shape, including vertices and edges
    %plot(shp);
    %plot vertices only
    %plot(shp.Points(:,1),shp.Points(:,2),'r*');
    plot(UU,VV,'k.')
    % plot edges as triangles
    for i=1:numel(igood)
        k4=colvec(Kfaces(igood(i),[1,2,3,1]));
        plot(UU(k4),VV(k4),'k.-');
    end
    % plot boundary of alpha shape
    plot(Ubound,Vbound,'b-','LineWidth',3);
    title(sprintf('%s\n%s',titlestr,'nNonconvex Alpha Shape'));
    
    
    %% plot mesh in 3D
    nf=nf+1;figure;hold on;
    % plot one triangle at a time
    for i=1:numel(igood)
        k4=colvec(Kfaces(igood(i),[1,2,3,1]));
       %fprintf(1,'%5d %5d %5d %5d %5d\n',i,k4);
        plot3(X(k4),Y(k4),Z(k4),'k.-');
    end
    axis tight;
    axis equal;
    title(sprintf('%s\n%s',titlestr,'Mesh'),'Interpreter','none');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    view([1, -1, 1]);
    %printpng(sprintf('%s_%s.png',mfilename,strrep(titlestr,'/','_')));
    close all    
end



return
end






