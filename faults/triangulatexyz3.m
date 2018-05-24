function [Kfaces,lgood] = triangulatexyz3(X,Y,Z,titlestr,make_plots)
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
    npoints
    warning('Too few points');
    Kfaces = NaN;
    Kboundary = NaN;
    return;
end
XYZ = [X,Y,Z];
% whos XYZ

%% find the best-fitting plane and calculate the modeled values of ZM = f(XM,YM) on a mesh
%compute the normal to the plane and a point that belongs to the plane
% normvec : a unit (column) vector normal to the plane
% basisvecs : a 3 by 2 matrix. Its columns form an orthonormal basis of the plane
% point0 : a point belonging to the plane
%
[normvec,basisvecs,point0] = affine_fit(XYZ);

%% calculate projections of points into UU,VV plane
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

% remove empty rows
UU = UU(1:kount);
VV = VV(1:kount);

if numel(UU) < 3 || numel(VV) < 3
    UU
    VV
    Kfaces = nan;
    Kboundary = nan;
    warning('miscount')
    return
end
    
%% try to make our own kfaces
Kfaces = delaunay(UU,VV);
[ntriangles,ncols]= size(Kfaces)
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

%%  make a mesh
[XM,YM] = meshgrid(  ...
     linspace(nanmin(XYZ(:,1)),nanmax(XYZ(:,1)),10)...
    ,linspace(nanmin(XYZ(:,2)),nanmax(XYZ(:,2)),10));

%%  evalute the formula for the plane
ZM = -1*(normvec(1)/normvec(3)*XM+normvec(2)/normvec(3)*YM-dot(normvec,point0)/normvec(3));

%% evaluate each triangle individually
lgood = false(ntriangles,1);
for i=1:ntriangles
    nok = 0;
    k3 = Kfaces(i,[1:3]);

    % test if all three vertices of triangle are inside alpha shape
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
%    
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
ngood = numel(lgood(igood))

if make_plots == 1
    %% plot points and alpha shape in U-V plane
    nf=nf+1;figure;
    plot(shp);
    hold on;
    plot(UU,VV,'k.')
    plot(UU(igood),VV(igood),'ro','MarkerSize',3);
    %plot(shp.Points(:,1),shp.Points(:,2),'r*');
    plot(Ubound,Vbound,'b-','LineWidth',3);
    title(sprintf('%s\n%s',titlestr,'nNonconvex Alpha Shape'));

    %% make a 3D plot    
    nf=nf+1;figure;
   %plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'r.');
    plot3(XYZ(igood,1),XYZ(igood,2),XYZ(igood,3),'r.');
    hold on;
    axis tight;
    axis equal;
    % plot the mesh
    surf(XM,YM,ZM,'facecolor','blue','facealpha',0.5);
    xlabel('X');
    xlabel('Y');
    ylabel('Z');
    title(sprintf('%s\n%s',titlestr,'3 D'));

%     %% make a 2D plot in UU-VV plane
%     nf=nf+1;figure;hold on;
%     % plot points
%     plot(UU,VV,'b.');
%     axis xy;axis tight;axis equal;
%     % plot boundary
%     %plot(UU(kbound),VV(kbound),'r-','linewidth',5);
%     plot(Ubound,Vbound,'r-','linewidth',5);
%     % plot each triangle individually
%     for j=1:ntriangles
%         k3 = Kfaces(j,[1,2,3]); % indices of points around triangle
%         umean = mean(UU(k3));
%         vmean = mean(VV(k3));
%         %plot(xmean,ymean,'b*');
%         %inside = inpolygon(xmean,ymean,UU(kbound),VV(kbound));
%         inside = inpolygon(umean,vmean,Ubound,Vbound);
%         if inside == 1
%             k4 = Kfaces(j,[1,2,3,1]); % indices of points around triangle, closing
%             plot(colvec(UU(k4)),colvec(VV(k4)),'k-');  % ours
%         end
%     end
%     xlabel('U');
%     ylabel('V');
%     title(sprintf('%s\n%s',titlestr,'Projected coordinates with boundary'));
%     %close all
    
    
%% plot mesh in 3D
nf=nf+1;figure;hold on;
for j=1:numel(igood)
    %        k3 = colvec(Kfaces(j,[1,2,3]));
    %         umean = mean(UU(k3));
    %         vmean = mean(VV(k3));
    %         %inside = inpolygon(xmean,ymean,UU(kbound),VV(kbound));
    %         inside = inpolygon(umean,vmean,Ubound,Vbound);
    %         if inside == 1
    %             k4=colvec(Kfaces(j,[1,2,3,1]));
    %             %plot3(colvec(XYZ(k4,1)),colvec(XYZ(k4,2)),colvec(XYZ(k4,3)),'k.-');
    %             plot3(X(k4),Y(k4),Z(k4),'k.-');
    %         end
    k4=colvec(Kfaces(igood(j),[1,2,3,1]));
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


return
end






