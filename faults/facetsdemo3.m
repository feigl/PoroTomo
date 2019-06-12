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



dX = 25;
x0 = 0:dX:500;
dY = 25;
y0 = 0:dY:1500;
dZ = 25;
z0 = 0:dZ:400;

% %% use faults
% %% read the faults.
% % S == structure containing fault coordinates
% % /Users/feigl/Box
% % Sync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff
% %load('read_and_mesh_faults9.mat'); % ALL Faults
% %load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff/read_and_mesh_faults10.mat'); % Nick's Top9
% %% choose a fault model
% %geologic_model = 'Jolie9'
% %geologic_model = 'Jolie'
% geologic_model = 'Siler'
% if strcmp(geologic_model,'Jolie9')
%     %    load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff/read_and_mesh_faults10.mat'); % Nick's Top9
%     load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie9.mat'); % All faults
% elseif strcmp(geologic_model,'Jolie')
%     load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie.mat'); % All faults
% elseif strcmp(geologic_model,'Siler')
%     load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Siler.mat'); % All in Siler's model
% else
%     error(sprintf('Unknown geologic_model: %s\n',geologic_model));
% end
% 
% FAULTS = S
% clear S;
% 
% nfaults = numel(FAULTS)
% FAULTS = facet_faults(FAULTS)

%FAULTAREAS = survey_fault_area2(FAULTS,x0,y0,z0);
%FAULTAREAS = survey_fault_area4(FAULTS,x0,y0,z0);
%save('FAULTAREAS25m.mat','-struct','FAULTAREAS');

%% recover coordinates neglacted above
load('SAVE20181016.mat');
FAULTAREAS = load('FAULTAREAS25m.mat');
dX = mean(colvec(diff(x0)))
dY = mean(colvec(diff(y0)))
dZ = mean(colvec(diff(z0)))

nX=numel(x0)
nY=numel(y0)
nZ=numel(z0)
nvoxels = nX*nY*nZ



tstart=tic;
kount = 0;
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
            
            % record centroids
            FAULTAREAS.X(kount) = mean(x1);
            FAULTAREAS.Y(kount) = mean(y1);
            FAULTAREAS.Z(kount) = mean(z1);
        end
    end
end
save('FAULTAREAS25mB.mat','-struct','FAULTAREAS');

            


