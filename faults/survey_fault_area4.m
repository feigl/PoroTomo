function FAULTAREAS = survey_fault_area3(FAULTS,x0,y0,z0)
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
% nvoxels = mx * my * mz

drawplots = 0;
%
%% count number of faults
nfaults = numel(FAULTS);

% find increments
dX = mean(colvec(diff(x0)))
dY = mean(colvec(diff(y0)))
dZ = mean(colvec(diff(z0)))

nX=numel(x0)
nY=numel(y0)
nZ=numel(z0)
nvoxels = nX*nY*nZ

% define coordinates of centroids of voxels
FAULTAREAS.X = nan(nvoxels,1);
FAULTAREAS.Y = nan(nvoxels,1);
FAULTAREAS.Z = nan(nvoxels,1);
FAULTAREAS.count  = zeros(nvoxels,1); % number of faults intersecting voxel [dimless]
FAULTAREAS.area   = zeros(nvoxels,1); % surface area of fault inside voxel [m^2]
FAULTAREAS.density= zeros(nvoxels,1); % fault area divided by voxel volume [m^2/m^3]


kount=0;
%nvoxels=10;

% nX=2
% nY=2
% nZ=2;

tstart=tic;
for ix=1:nX
    for iy=1:nY
        for iz=1:nZ
            % for ix=3
            %     for iy=3
            %         for iz=3
            kount = kount+1;
            waitbar(kount/nvoxels)
            
            x1 = colvec([x0(ix)-dX/2., x0(ix)+dX/2., x0(ix)+dX/2., x0(ix)-dX/2.]);
            y1 = colvec([y0(iy)-dY/2., y0(iy)-dY/2., y0(iy)+dY/2., y0(iy)+dY/2.]);
            z1 = (z0(iz)-dZ/2.)*ones(size(x1));
            x1 = [x1;x1];
            y1 = [y1;y1];
            z1 = [z1;z1+dZ];
            
            DT = delaunayTriangulation([x1,y1,z1]);
            %S1.faces = DT.ConnectivityList;
            S1.faces=freeBoundary(DT);
            %S1.vertices = [x1,y1,z1];
            S1.vertices = DT.Points;
            
            area3 = 0;
            
            %S4=struct([]);
            S4.vertices=[];
            for j=1:nfaults
                S2.faces = FAULTS(j).faces;
                S2.vertices = FAULTS(j).vertices;
                [nfaces2,n2] = size(S2.faces);
                S2.facevertexcdata=ones(nfaces2,3);
                
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
                
                if numel(S2.faces) > 0 && numel(S2.vertices) >= 3
                    [Mintersect, S3] = SurfaceIntersection(S1,S2);                  
                    if numel(S3.vertices) > 0 && numel(S3.faces) > 1
                        %S3
                        %fprintf(1,'Intersection found at kount = %d\n',kount);
                        
                        
                        % now make another triangulation of the
                        % intersecting surface
                        DT3=delaunayTriangulation(S3.vertices);
                        if numel(DT3.Points) > 3 && numel(DT3.ConnectivityList) > 1
                            S5.vertices = DT3.Points;
                            S5.faces    = freeBoundary(DT3);
                            FAULTAREAS.count(kount)   = FAULTAREAS.count(kount) + 1;
                            if numel(S5.faces) > 0 && numel(S5.vertices) >= 3
                                area3=areaIsosurface(S5.faces,S5.vertices);
                                FAULTAREAS.area(kount)    = FAULTAREAS.area(kount) + area3;
                                FAULTAREAS.density(kount) = FAULTAREAS.area(kount)/dX/dY/dZ;
                            end
                            
                            if drawplots == 1
                                figure(kount);
                                trisurf(S1.faces,S1.vertices(:,1),S1.vertices(:,2),S1.vertices(:,3)...
                                    ,'FaceColor','blue','FaceAlpha',0.5);
                                axis equal;
                                hold on;
                                [nv,n3] = size(S1.vertices);
                                [nf,n3] = size(S1.faces);
                                for p=1:nf
                                    q=[S1.faces(p,:),S1.faces(p,1)];
                                    plot3(colvec(S1.vertices(q,1))...
                                        ,colvec(S1.vertices(q,2))...
                                        ,colvec(S1.vertices(q,3))...
                                        ,'b*-','LineWidth',2);
                                end
                                
                                
                                trisurf(S2.faces,S2.vertices(:,1),S2.vertices(:,2),S2.vertices(:,3)...
                                    ,'FaceColor','red','FaceAlpha',0.5);
                                trisurf(S3.faces,S3.vertices(:,1),S3.vertices(:,2),S3.vertices(:,3)...
                                    ,'FaceColor','green','FaceAlpha',0.5);
                                
                                axis([0, 500, 0, 1500, 0, 400]);
                                %hold on;
                                xlabel('Xp [m]');
                                ylabel('Yp [m]');
                                zlabel('Zp [m]');
                                
                                [nv,n3] = size(S3.vertices);
                                [nf,n3] = size(S3.faces);
                                for p=1:nf
                                    q=[S3.faces(p,:),S3.faces(p,1)];
                                    plot3(colvec(S3.vertices(q,1))...
                                        ,colvec(S3.vertices(q,2))...
                                        ,colvec(S3.vertices(q,3))...
                                        ,'g*-','LineWidth',2);
                                end
                                
                                trisurf(S5.faces,S5.vertices(:,1),S5.vertices(:,2),S5.vertices(:,3)...
                                    ,'FaceColor','yellow','FaceAlpha',0.5);
                            end                           
                        end
                    end
                end
            end
            fprintf(1,'After %10.1f seconds, completed %6d of %6d voxels. Area = %10.4E [m^3] Density = %10.4E [m^2/m^3]\n'...
                ,toc(tstart),kount,nvoxels,FAULTAREAS.area(kount),FAULTAREAS.density(kount));
        end
    end
end
return
end