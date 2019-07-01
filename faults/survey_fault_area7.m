function FAULTAREAS = survey_fault_area7(FAULTS,MESH3)
% draw faults on planar slices
% input:
% FAULTS == structure containing fault coordinates
% GRIDX,GRIDY,GRIDZ 3-D arrays of coordinates
% 20181015 Kurt Feigl

%% count number of voxels
% [mx,nx,px] = size(GRIDX)
% [my,ny,py] = size(GRIDY)
% [mz,nz,pz] = size(GRIDZ)
% if mx ~= my || my ~= mz
%     error('miscount in X');
% end
% if nx~= ny || ny ~= nz
%     error('miscount in Y');
% end
% if px ~= py || py ~= pz
%     error('miscount in Z');
% end
% nnodes = mx * my * mz

drawplots = 1;
%
%% count number of faults
nfaults = numel(FAULTS);

% find increments
% dx = nanmean(colvec(diff(MESH3.Xp(:,1,1))))
% dy = nanmean(colvec(diff(MESH3.Yp(1,:,1))))
% dz = nanmean(colvec(diff(MESH3.Zp(1,1,:))))
dx = MESH3.dx;
dy = MESH3.dy
dz = MESH3.dz

nnodes = numel(MESH3.Xp);

% define coordinates of centroids of voxels
% FAULTAREAS.X = nan(nnodes,1);
% FAULTAREAS.Y = nan(nnodes,1);
% FAULTAREAS.Z = nan(nnodes,1);
FAULTAREAS.count  = zeros(nnodes,1); % number of faults intersecting voxel [dimless]
FAULTAREAS.area   = zeros(nnodes,1); % surface area of fault inside voxel [m^2]
FAULTAREAS.density= zeros(nnodes,1); % fault area divided by voxel volume [m^2/m^3]
% make plaid array
% FAULTAREAS.X = nan(nx,ny,nz);
% FAULTAREAS.Y = nan(nx,ny,nz);
% FAULTAREAS.Z = nan(nx,ny,nz);
% FAULTAREAS.count  = zeros(nx,ny,nz); % number of faults intersecting voxel [dimless]
% FAULTAREAS.area   = zeros(nx,ny,nz); % surface area of fault inside voxel [m^2]
% FAULTAREAS.density= zeros(nx,ny,nz); % fault area divided by voxel volume [m^2/m^3]
% 

tstart=tic;
% [MESH3.Xp,MESH3.Yp,MESH3.Zp] = meshgrid(xv,yv,zv);
% nnodes = numel(MESH3.Xp)


for knode = 1:nnodes
    %for knode = 33
    % list of coordinates in intersection
    XI=[];
    YI=[];
    ZI=[];
    r1 = [];
    s1 = [];
    t1 = [];
    if knode == 33
        figure;hold on;
    end
    for knorm = 1:3
        switch knorm
            case 1
                kr = 2;
                ks = 3;
                [R1,S1] = meshgrid([MESH3.Yp(knode)-dy/2,MESH3.Yp(knode)+dy/2],[MESH3.Zp(knode)-dz/2,MESH3.Zp(knode)+dz/2]);
            case 2
                kr = 1;
                ks = 3;
                [R1,S1] = meshgrid([MESH3.Xp(knode)-dx/2,MESH3.Xp(knode)+dx/2],[MESH3.Zp(knode)-dz/2,MESH3.Zp(knode)+dz/2]);
            case 3
                kr = 1;
                ks = 2;
                [R1,S1] = meshgrid([MESH3.Xp(knode)-dx/2,MESH3.Xp(knode)+dx/2],[MESH3.Yp(knode)-dy/2,MESH3.Yp(knode)+dy/2]);
            otherwise
                error('unknown knorm');
        end
        
        % loop over left, right sides
        for kside = [-1,1]
            clear SS1;
            switch knorm
                case 1
                    T1 = (MESH3.Xp(knode) + kside*dx/2)*ones(size(R1));
                case 2
                    T1 = (MESH3.Yp(knode) + kside*dy/2)*ones(size(R1));
                case 3
                    T1 = (MESH3.Zp(knode) + kside*dz/2)*ones(size(R1));
                otherwise
                    error('unknown knorm');
            end
            
            r1 = colvec(R1);
            s1 = colvec(S1);
            t1 = colvec(T1);
            DT = delaunayTriangulation([r1,s1]);
            [npoints,n2] = size(DT.Points);
            for i=1:npoints
                SS1.vertices(i,kr) = DT.Points(i,1);
                SS1.vertices(i,ks) = DT.Points(i,2);
                SS1.vertices(i,knorm) = t1(i);
            end
            
            SS1.faces = DT.ConnectivityList;
            
            area_total = 0;
            
            S4.vertices=[];
            kount = 0;
            for j=1:nfaults
                %for j=21
                %for j=20:22
                SS2.faces = FAULTS(j).faces;
                SS2.vertices = FAULTS(j).vertices;
                [nfaces2,n2] = size(SS2.faces);
                %SS2.facevertexcdata=ones(nfaces2,3);
                
                if numel(SS2.faces) > 0 && numel(SS2.vertices) >= 3
                    kount = kount+1;
                    
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
                    
                    [Mintersect, SS3] = SurfaceIntersection(SS1,SS2);
                    if numel(SS3.vertices) > 0 && numel(SS3.faces) > 1 && numel(SS3.edges) > 0
                        
                        %fprintf(1,'Intersection found at j = %d\n',j);
                        
                        if knode == 33
                            trisurf(SS1.faces,SS1.vertices(:,1),SS1.vertices(:,2),SS1.vertices(:,3)...
                                ,'FaceColor','blue','FaceAlpha',0.5);
                            axis equal;
                            hold on;
                            
                            trisurf(SS2.faces,SS2.vertices(:,1),SS2.vertices(:,2),SS2.vertices(:,3)...
                                ,'FaceColor','red','FaceAlpha',0.5);
                            
                            axis([0, 500, 0, 1500, 0, 400]);
                            hold on;
                            xlabel('Xp [m]');
                            ylabel('Yp [m]');
                            zlabel('Zp [m]');
                            
                        end
                        
                        %                         [nvertices,n3] = size(SS3.vertices);
                        %                         [nfaces,n3] = size(SS3.faces);
                        %                         for iface=1:nfaces
                        %                             ivertex=SS3.faces(iface,:);
                        % %                             plot3(colvec(SS3.vertices(ivertex,1))...
                        % %                                 ,colvec(SS3.vertices(ivertex,2))...
                        % %                                 ,colvec(SS3.vertices(ivertex,3))...
                        % %                                 ,'go-','LineWidth',2);
                        %                             XI=[XI;colvec(SS3.vertices(ivertex,1))];
                        %                             YI=[YI;colvec(SS3.vertices(ivertex,2))];
                        %                             ZI=[ZI;colvec(SS3.vertices(ivertex,3))];
                        %                         end
                        
                        % draw edges
                        [n2,nedges] = size(SS3.edges);
                        for iedge = 1:nedges
                            xedge1 = colvec(SS3.vertices(SS3.edges(1,iedge),1));
                            yedge1 = colvec(SS3.vertices(SS3.edges(1,iedge),2));
                            zedge1 = colvec(SS3.vertices(SS3.edges(1,iedge),3));
                            xedge2 = colvec(SS3.vertices(SS3.edges(2,iedge),1));
                            yedge2 = colvec(SS3.vertices(SS3.edges(2,iedge),2));
                            zedge2 = colvec(SS3.vertices(SS3.edges(2,iedge),3));
                            
                            if knode == 33
                                plot3([xedge1,xedge2],[yedge1,yedge2],[zedge1,zedge2],'g-*');
                            end
                            XI=[XI;xedge1;xedge2];
                            YI=[YI;yedge1;yedge2];
                            ZI=[ZI;zedge1;zedge2];
                        end
                    end
                end
            end
        end
    end
    %% project onto plane and calculate area
    SS4 = triangulatexyz4(XI,YI,ZI);
    
    FAULTAREAS.area(knode) = FAULTAREAS.area(knode) + SS4.area;
    FAULTAREAS.density(knode) = FAULTAREAS.area(knode)/(dx*dy*dz);
    
    % for debugging
    if knode == 33
        %figure; hold on;
        plot3(XI,YI,ZI,'ro');
        trisurf(SS4.faces,SS4.vertices(:,1),SS4.vertices(:,2),SS4.vertices(:,3)...
            ,'FaceColor','blue','FaceAlpha',0.5);
        [nvertices,n3] = size(SS4.vertices);
        [nfaces,n3] = size(SS4.faces);
        for iface=1:nfaces
            ivertex=SS4.faces(iface,:);
            plot3(colvec(SS4.vertices(ivertex,1))...
                ,colvec(SS4.vertices(ivertex,2))...
                ,colvec(SS4.vertices(ivertex,3))...
                ,'go-','LineWidth',2);
        end
        
        axis([0, 500, 0, 1500, 0, 400]);
        xlabel('Xp [m]');
        ylabel('Yp [m]');
        zlabel('Zp [m]');
        title(sprintf('node %d area %10.4E [m^2]',knode,area(knode)));
    end
    
    fprintf(1,'After %10.1f seconds, completed %6d of %6d voxels. Area = %10.4E [m^3] Density = %10.4E [m^2/m^3]\n'...
        ,toc(tstart),knode,nnodes,FAULTAREAS.area(knode),FAULTAREAS.density(knode));
    
end
FAULTAREAS.X = MESH3.Xp;
FAULTAREAS.Y = MESH3.Yp;
FAULTAREAS.Z = MESH3.Zp;
return
end