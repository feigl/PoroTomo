% make facets for a cube

%https://www.mathworks.com/help/matlab/ref/triangulation.facenormal.html
clear all;
close all;
% make a plane
% side1 = [0, 1];
% [Xm,Ym,Zm] = meshgrid(side1,side1,1);
% x1 = colvec(Xm);
% y1 = colvec(Ym);
% z1 = colvec(Zm)
%% make a plane in PoroTomo Box
% [Xm,Ym] = meshgrid([0,500],[0,1500]);
% Zm = 200*ones(size(Xm));
% P1=surf2patch(Xm,Ym,Zm,Zm)

% [x1,y1] = meshgrid([0, 500],[1500, 0]);
% x1 = colvec(x1);
% y1 = colvec(y1);
% x1 = colvec([0 500 500 0 ]);
% y1 = colvec([0 0   1500 1500]);
% z1 = 200*ones(size(x1));
% x1 = [x1;x1];
% y1 = [y1;y1];
% z1 = [z1;z1+100];



dX = 200;
x0 = 0:dX:500;
dY = 300;
y0 = 0:dY:1500;
dZ = 100;
z0 = 0:dZ:400;

nX=numel(x0)
nY=numel(y0)
nZ=numel(z0)
for ix=1:nX
    for iy=1:nY
        for iz=1:nZ
            
            x1 = colvec([x0(ix)-dX/2., x0(ix)+dX/2.-eps, x0(ix)+dX/2.-eps, x0(ix)-dX/2.    ]);
            y1 = colvec([y0(iy)-dY/2., y0(iy)-dY/2.,     y0(iy)+dY/2.-eps, y0(iy)+dY/2.-eps]);
            z1 = (z0(iz)-dZ/2.)*ones(size(x1));
            x1 = [x1;x1];
            y1 = [y1;y1];
            z1 = [z1;z1-eps+dZ];
            
            
            S1.vertices = [x1,y1,z1];
            S1.faces = [...
                1,2,3,4;...
                5,6,7,8;...
                1,2,6,5;...
                2,3,7,6;...
                3,4,8,7;...
                1,4,8,5];
            [nvertices,n1] = size(S1.vertices);
            S1.facevertexcdata=ones(nvertices,3);
            S1
            area1 = areaIsosurface(S1.faces,S1.vertices)
            
            %% use faults
            %% read the faults.
            % S == structure containing fault coordinates
            % /Users/feigl/Box
            % Sync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff
            %load('read_and_mesh_faults9.mat'); % ALL Faults
            %load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff/read_and_mesh_faults10.mat'); % Nick's Top9
            %% choose a fault model
            %geologic_model = 'Jolie9'
            %geologic_model = 'Jolie'
            geologic_model = 'Siler'
            if strcmp(geologic_model,'Jolie9')
                %    load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff/read_and_mesh_faults10.mat'); % Nick's Top9
                load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie9.mat'); % All faults
            elseif strcmp(geologic_model,'Jolie')
                load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie.mat'); % All faults
            elseif strcmp(geologic_model,'Siler')
                load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Siler.mat'); % All in Siler's model
            else
                error(sprintf('Unknown geologic_model: %s\n',geologic_model));
            end
            
            FAULTS = S
            clear S;
            
            nfaults = numel(FAULTS)
            %nfaults = 5;
            
            % get indices of vertices
            for j=1:nfaults
                %for j=[15,59]
                if sum(FAULTS(j).trigood) ~= 0
                    [ntriangles,n3] = size(FAULTS(j).TRI.ConnectivityList);
                    %ntriangles = 10;
                    fprintf(1,'Working on fault number %d named %s\n',j,char(FAULTS(j).faultname));
                    [ntriangles,n3] = size(FAULTS(j).TRI.ConnectivityList);
                    X2=zeros(3*ntriangles,1);
                    Y2=zeros(3*ntriangles,1);
                    Z2=zeros(3*ntriangles,1);
                    FF=zeros(ntriangles,3);
                    kk=0;
                    k1=1;
                    k3=3;
                    for k=1:ntriangles
                        % get indices of vertices
                        kA=colvec(FAULTS(j).TRI.ConnectivityList(k,1));
                        kB=colvec(FAULTS(j).TRI.ConnectivityList(k,2));
                        kC=colvec(FAULTS(j).TRI.ConnectivityList(k,3));
                        
                        % get coordinates of vertices
                        x2=colvec([FAULTS(j).Xtri(kA),FAULTS(j).Xtri(kB),FAULTS(j).Xtri(kC)]);
                        y2=colvec([FAULTS(j).Ytri(kA),FAULTS(j).Ytri(kB),FAULTS(j).Ytri(kC)]);
                        z2=colvec([FAULTS(j).Ztri(kA),FAULTS(j).Ztri(kB),FAULTS(j).Ztri(kC)]);
                        
                        
                        X2(k1:k3,1)=x2;
                        Y2(k1:k3,1)=y2;
                        Z2(k1:k3,1)=z2;
                        kk = kk+1;
                        FF(kk,1:3)=rowvec([k1:k3]);
                        k1 = k1+3;
                        k3 = k3+3;
                        
                        
                        
                        %             S2.vertices = [x2,y2,z2];
                        %             S2.faces = [1,2,3];
                        %             [nfaces2,n2] = size(S2.faces);
                        %             S2.facevertexcdata=ones(nfaces2,3);
                        
                        
                        
                        
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
                        
                        
                        
                        %             [Mintersect, S3] = SurfaceIntersection(S1,S2);
                        %
                        %             if numel(S3.vertices) > 0
                        %
                        %                 trisurf(S2.faces,S2.vertices(:,1),S2.vertices(:,2),S2.vertices(:,3),'FaceColor','blue','FaceAlpha',0.5);
                        %
                        %
                        %                 trisurf(S3.faces,S3.vertices(:,1),S3.vertices(:,2),S3.vertices(:,3),'FaceColor','green','FaceAlpha',0.5);
                        %                 axis equal;
                        %                 hold on;
                        %
                        %                 % quiver3(P1(:,1),P1(:,2),P1(:,3),F1(:,1),F1(:,2),F1(:,3),0.5,'color','blue');
                        %                 % quiver3(S2(:,1),S2(:,2),S2(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'color','red');
                        %                 plot3(S3.vertices(:,1),S3.vertices(:,2),S3.vertices(:,3),'g*-');
                        
                    end
                    S2.vertices = [X2,Y2,Z2];
                    S2.faces = FF;
                    [nfaces2,n2] = size(S2.faces);
                    S2.facevertexcdata=ones(nfaces2,3);
                    %S2
                    
                    %S2 = S1;
                    [Mintersect, S3] = SurfaceIntersection(S1,S2);
                    
                    if numel(S3.vertices) > 0
                        fprintf(1,'Intersection found\n');
                        %close(gcf);
                        figure(j);
                        trisurf(S1.faces,S1.vertices(:,1),S1.vertices(:,2),S1.vertices(:,3),'FaceColor','blue','FaceAlpha',0.5);
                        axis equal;
                        hold on;
                        [nv,n3] = size(S1.vertices);
                        [nf,n3] = size(S1.faces);
                        for p=1:nf
                            pp=[S1.faces(p,:),S1.faces(p,1)];
                            plot3(colvec(S1.vertices(pp,1))...
                                ,colvec(S1.vertices(pp,2))...
                                ,colvec(S1.vertices(pp,3))...
                                ,'b*-','LineWidth',2);
                        end
                        
                        
                        trisurf(S2.faces,S2.vertices(:,1),S2.vertices(:,2),S2.vertices(:,3),'FaceColor','red','FaceAlpha',0.5);
                        
                        
                        trisurf(S3.faces,S3.vertices(:,1),S3.vertices(:,2),S3.vertices(:,3),'FaceColor','green','FaceAlpha',0.5);
                        
                        %axis([0, 500, 0, 1500, 0, 400]);
                        %hold on;
                        xlabel('Xp [m]');
                        ylabel('Yp [m]');
                        zlabel('Zp [m]');
                        
                        %             quiver3(P1(:,1),P1(:,2),P1(:,3),F1(:,1),F1(:,2),F1(:,3),0.5,'color','blue');
                        %             quiver3(S2(:,1),S2(:,2),S2(:,3),F2(:,1),F2(:,2),F2(:,3),0.5,'color','red');
                        [nv,n3] = size(S3.vertices);
                        [nf,n3] = size(S3.faces);
                        % add one more face to close
                        S3.faces=[S3.faces;S3.faces(1,:)]
                        for p=1:nf
                            for q=1:3
                            qq=[S3.faces(p,:),S3.faces(p,1)];
                            line(colvec(S3.vertices(qq,1))...
                                ,colvec(S3.vertices(qq,2))...
                                ,colvec(S3.vertices(qq,3))...
                                ,'g-','LineWidth',2);
                            end
                        end
                        area3 = areaIsosurface(S3.faces,S3.vertices)
                    end
                end
            end
        end
    end
end




