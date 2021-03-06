%% interpolate fields for PoroTomo
% Build a configuration for PoroTomo project.
% 20171007 Kurt Feigl
% 20171026 with Dante - use constant color scale
% 20171120 with Lesley - add Matzel's latest
% 20171205 Kurt for AGU
% 20180207 Kurt for Stanford
% 20180211 Kurt for GRL 
% 20180507 Kurt for SRL and SSA
% 20180513 Kurt work on lithologic units
% 20180809 Kurt add voxel estimates from Reinisch et al.
% 20180829 Kurt make depth anomalies
% 20181005 Kurt add Avinash model
% 20181008 Kurt try to count faults

%% Examples of Imagemagick
% example how to make an animation
% convert -delay 1 Vp_Thurber20171123*.png VpThurber20171123.gif

% example how to make a montage
% RM MontageZ00400m_Siler.jpg
% montage *ZnormZ00400m_Siler.jpg -geometry +2+2 MontageZ00400m_Siler.jpg



%% initialize
clear all;
close all;
nf=0;
tstart = tic;

save_wells = 1;   % save Z profile
save_meshes = 0;  % save tomograms on REGULAR 3 D mesh?
make_slices = 1;  % plot cross sections?
%show_anomalies = 1; % subtract mean value for depth from all values
show_anomalies = 0; % DO NOT subtract mean value for depth from all values
standardslices = 1


%% set up path for Matlab
%addpath(genpath('/globus/PoroTomo/SOFTWARE'),'-begin');
%addpath('/Users/feigl/gipht/utils','-begin'); % needed for colortables
addpath(genpath('/Users/feigl/PoroTomo'),'-begin');
rmpath(genpath('/Users/feigl/PoroTomo/.git'));

%% before releasing, find dependencies using:
% [fList,pList] = matlab.codetools.requiredFilesAndProducts(mfilename);
% fList{:}
% pList.Name


%% choose the function for printing
%funprint = str2func('printpng'); % images as bitmaps
%funprint = str2func('printpdf'); % PDF with labels
%funprint = str2func('printeps'); % EPS, but rendered ugly by Mac Preview
%funprint = str2func('printfig'); % save as "compact" Matlab figure format
funprint = str2func('printjpg'); % JPG
%funprint = str2func('save2pdf'); % PDF cropped without labels

%% Digital Elevation model in Latitude and Longitude
% TODO: should use webget to ftp://roftp.ssec.wisc.edu/porotomo/PoroTomo/METADATA/
demgrdfilename = '/Users/feigl/BoxSync/PoroTomo/METADATA/brady_dem_srtm_20m.grd'

%% get the coordinates of the study area
utmzone = '11 S' % UTM zone for Brady Hot Springs
% Start with Latitude and Longitude.
% BRADYBOX = read_bradybox(1)
% [BRADYBOX.Xp,BRADYBOX.Yp] = utm2xy_porotomo(BRADYBOX.E,BRADYBOX.N);
% [BRADYBOX.Lat,BRADYBOX.Lon] = utm2deg(colvec(BRADYBOX.E),colvec(BRADYBOX.N),repmat(utmzone,size(colvec(BRADYBOX.E))));
% [BRADYBOX.Xp,BRADYBOX.Yp]
%          -0.10       1533.39
%              0             0
%         460.84          0.00
%         469.25       1531.42
% 20180527 - start with PoroTomo coordinates
xy = [  0., 1500.; ...
        0.,    0.; ...
      500.,    0.; ...
      500., 1500.];
BRADYBOX.Xp = xy(:,1);
BRADYBOX.Yp = xy(:,2);
[BRADYBOX.E,BRADYBOX.N] = xy_porotomo2utm(BRADYBOX.Xp,BRADYBOX.Yp);
[BRADYBOX.Lat,BRADYBOX.Lon] = utm2deg(colvec(BRADYBOX.E),colvec(BRADYBOX.N),repmat(utmzone,size(colvec(BRADYBOX.E))));

%% Set bounding box
% These are Zporotomo coordinates in meters
% %BOUNDS.Zp = [0, 450]; % These are Zporotomo coordinates.
% % Elevation above WGS84 goes from 800 m to 1250 m
% BOUNDS.Zp = [0, 490]; % These are Zporotomo coordinates.
% % Elevation above WGS84 goes from 800 m to 1290 m
ibound = 1
switch ibound
    case 1
        % Official PoroTomo Natural Laboratory
        BOUNDS.Xp = [0,500];
        BOUNDS.Yp = [0,1500];
        BOUNDS.Zp = [0, 450]; 
    case 2
        % Big Boxes for AGU 2017 and Geothermics paper
        BOUNDS.Xp = [-400, 900];
        BOUNDS.Yp = [-300,1800];
        BOUNDS.Zp = [  0, 490];        
    case 3
        % Stay within small area of best resolution - use for MASW
        BOUNDS.Xp = [100,400];
        BOUNDS.Yp = [200,1400];
        BOUNDS.Zp = [390, 450];       
    case 4
        % Stay within small area of best resolution -
        % use for final review and AGU 2018
        BOUNDS.Xp = [0 500];
        BOUNDS.Yp = [0,1500];
        BOUNDS.Zp = [200, 450];       
    otherwise
        ibound
        error(sprintf('Unknown ibound: %d\n',ibound));
end


            
% make a 3-D bounding box
[BOX3.Xp,BOX3.Yp,BOX3.Zp] = meshgrid(BOUNDS.Xp,BOUNDS.Yp,BOUNDS.Zp);          
[BOX3.E, BOX3.N, BOX3.H] = xyz_porotomo2utm(colvec(BOX3.Xp),colvec(BOX3.Yp),colvec(BOX3.Zp));
% % make a 2-D bounding box
% [BOX2.Xp,BOX2.Yp] = meshgrid(BOUNDS.Xp,BOUNDS.Yp,BOUNDS.Zp);          
% [BOX.E,  BOX.N]   = xy_porotomo2utm(colvec(BOX2.Xp),colvec(BOX2.Yp));


% make 2-D grid with 25 m spacing using Brady Box
% [BRADY2DGRID.Xp,BRADY2DGRID.Yp] = meshgrid([min(BRADYBOX.Xp):25:max(BRADYBOX.Xp)]...
%     ,[min(BRADYBOX.Yp):25:max(BRADYBOX.Yp)]);
% make 2-D grid with 25 m spacing using BOUNDS
[BRADY2DGRID.Xp,BRADY2DGRID.Yp] = meshgrid([min(BOUNDS.Xp):25:max(BOUNDS.Xp)]...
    ,[min(BOUNDS.Yp):25:max(BOUNDS.Yp)]);

[BRADY2DGRID.E,BRADY2DGRID.N] = xy_porotomo2utm(colvec(BRADY2DGRID.Xp),colvec(BRADY2DGRID.Yp));
[BRADY2DGRID.Lat,BRADY2DGRID.Lon] = utm2deg(colvec(BRADY2DGRID.E),colvec(BRADY2DGRID.N),repmat(utmzone,size(colvec(BRADY2DGRID.E))));
BRADY2DGRID.ElevWGS84=get_elevation(colvec(BRADY2DGRID.Lon),colvec(BRADY2DGRID.Lat),demgrdfilename);
format bank
elev_mean = nanmean(colvec(BRADY2DGRID.ElevWGS84))
elev_std  = nanstd(colvec(BRADY2DGRID.ElevWGS84))
elev_min  = nanmin(colvec(BRADY2DGRID.ElevWGS84))
elev_max  = nanmax(colvec(BRADY2DGRID.ElevWGS84))
BRADY2DGRID.Zp = reshape(BRADY2DGRID.ElevWGS84 - 800,size(BRADY2DGRID.Xp));

% get vertical coordinates of corners
BRADYBOX.ElevWGS84=get_elevation(colvec(BRADYBOX.Lon),colvec(BRADYBOX.Lat),demgrdfilename);
BRADYBOX.ElevWGS84 = reshape(BRADYBOX.ElevWGS84,size(BRADYBOX.Xp));
BRADYBOX.Zp = BRADYBOX.ElevWGS84 - 800;

save('BRADYBOX.mat','-struct','BRADYBOX');
save('BRADY2DGRID.mat','-struct','BRADY2DGRID');

% make a 3-D mesh
% dx = 10; % meters
% dy = 10; % meters
% dz = 10; % meters
% dx = 100; % meters
% dy = 100; % meters
% dz = 100; % meters
dx = 25; % meters
dy = 25; % meters
dz = 25; % meters


% [MESH3.Xp,MESH3.Yp,MESH3.Zp] = meshgrid(...
%      [min(BOUNDS.Xp):dX:max(BOUNDS.Xp)]...
%     ,[min(BOUNDS.Yp):dY:max(BOUNDS.Yp)]...
%     ,[min(BOUNDS.Zp):dZ:max(BOUNDS.Zp)]);   
% center voxels on nodes
[MESH3.Xp,MESH3.Yp,MESH3.Zp] = meshgrid(...
     [min(BOUNDS.Xp)-dx/2.:dx:max(BOUNDS.Xp)+dx/2.]...
    ,[min(BOUNDS.Yp)-dy/2.:dy:max(BOUNDS.Yp)+dy/2.]...
    ,[min(BOUNDS.Zp)-dz/2.:dz:max(BOUNDS.Zp)+dz/2.]);          
[MESH3.E, MESH3.N, MESH3.H] = xyz_porotomo2utm(colvec(MESH3.Xp),colvec(MESH3.Yp),colvec(MESH3.Zp));
MESH3.E = reshape(MESH3.E,size(MESH3.Xp));
MESH3.N = reshape(MESH3.N,size(MESH3.Xp));
MESH3.H = reshape(MESH3.H,size(MESH3.Xp));
[MESH3.nx,MESH3.ny,MESH3.nz] = size(MESH3.Xp)
MESH3.dx = dx;
MESH3.dy = dy;
MESH3.dz = dz;

% find UTM coordinates
utmzone = '11 S';
[MESH3.Lat,MESH3.Lon] = utm2deg(colvec(MESH3.E),colvec(MESH3.N),utmzone);
MESH3.Lat=reshape(MESH3.Lat,size(MESH3.Xp));
MESH3.Lon=reshape(MESH3.Lon,size(MESH3.Xp));




%% read the file containing coordinates of nodes
dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_4_Adjoint_seismic_tomography/MESH_SPECFEM3D'
fname_nodes = strcat(dirname,filesep,'MESH_topo_xpypzellipsoid_800/nodes_coords_file');
NODES = read_nodes(fname_nodes);

%% read the file containing the list of elements
dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_4_Adjoint_seismic_tomography/MESH_SPECFEM3D'
fname_elements = strcat(dirname,filesep,'MESH_topo_xpypzellipsoid_800/mesh_file')
ELEMENTS = read_elements(fname_elements);

%% Read the coordinates of the centroids from a file
fname_centroids = strcat(dirname,filesep,'MESH_topo_xpypzellipsoid_800/centroids_coords_file')
CENTROIDS = read_centroids(fname_centroids);
%  get UTM coordinates of centroids
[CENTROIDS.E,CENTROIDS.N] = xy_porotomo2utm(CENTROIDS.Xp,CENTROIDS.Yp);
% find latitude and longitude
[CENTROIDS.Lat,CENTROIDS.Lon] = utm2deg(colvec(CENTROIDS.E),colvec(CENTROIDS.N),repmat(utmzone,size(colvec(CENTROIDS.E))));
% get the topographic elevation of the CENTROIDS in the model
CENTROIDS.ElevWGS84 = get_elevation(CENTROIDS.Lon,CENTROIDS.Lat,demgrdfilename) - (CENTROIDS.Zp + 400); % Elevation_in_meters_aboveWGS84ellipsoid
save('CENTROIDS.mat','-struct','CENTROIDS');


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
if fexist('FAULTSwithFacets.mat') == 1
    load('FAULTSwithFacets.mat');
else
    FAULTS = facet_faults(FAULTS);
    save('FAULTSwithFacets.mat','FAULTS');
end







%% get the coordinates for the wells from a file
WELLS = read_wells('brady');
nwells = numel(WELLS.WellName)
for i=1:nwells
    [WELLS.Xp(i),WELLS.Yp(i),WELLS.Zp(i)] =  utm2xyz_porotomo(...
        WELLS.XCoordinate_UTMME_(i)...
        ,WELLS.YCoordinate_UTMMN_(i)...
        ,WELLS.LandSurfElev_mAboveWGS84Ellipsoid_(i));
end
WELLS.Xp=colvec(WELLS.Xp);
WELLS.Yp=colvec(WELLS.Yp);
WELLS.Zp=colvec(WELLS.Zp);

% make an array of depths every 1 m
zprof=colvec([min(BOUNDS.Zp):1:max(BOUNDS.Zp)]);

%% get the coordinates of the fumaroles
T=readtable('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Surface_features/brady_fumerole_ll.txt');
lonlat = table2struct(T,'ToScalar',true);
FUMAROLES.Lon = lonlat.Var1;
FUMAROLES.Lat = lonlat.Var2;
iok = intersect(find(FUMAROLES.Lat > 0.),find(FUMAROLES.Lon < 0.));
FUMAROLES.Lat = FUMAROLES.Lat(iok);
FUMAROLES.Lon = FUMAROLES.Lon(iok);
clear lonlat;
FUMAROLES.H = get_elevation(FUMAROLES.Lon,FUMAROLES.Lat,demgrdfilename); % elevation in meters above WGS84 ellipsoid
[FUMAROLES.E,FUMAROLES.N] = deg2utm(colvec(FUMAROLES.Lat),colvec(FUMAROLES.Lon));
for i=1:numel(FUMAROLES.Lon)
    [FUMAROLES.Xp(i),FUMAROLES.Yp(i),FUMAROLES.Zp(i)] =  utm2xyz_porotomo(...
        FUMAROLES.E(i)...
        ,FUMAROLES.N(i)...
        ,FUMAROLES.H(i));
end
save('FUMAROLES.mat','-struct','FUMAROLES');

%% get the coordinates of the mudpots
T=readtable('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Surface_features/brady_mudpot_ll.txt');
lonlat = table2struct(T,'ToScalar',true);
MUDPOTS.Lon = lonlat.Var1;
MUDPOTS.Lat = lonlat.Var2;
iok = intersect(find(MUDPOTS.Lat > 0.),find(MUDPOTS.Lon < 0.));
MUDPOTS.Lat = MUDPOTS.Lat(iok);
MUDPOTS.Lon = MUDPOTS.Lon(iok);
clear lonlat;
MUDPOTS.H = get_elevation(MUDPOTS.Lon,MUDPOTS.Lat,demgrdfilename); % elevation in meters above WGS84 ellipsoid
[MUDPOTS.E,MUDPOTS.N] = deg2utm(colvec(MUDPOTS.Lat),colvec(MUDPOTS.Lon));
for i=1:numel(MUDPOTS.Lon)
    [MUDPOTS.Xp(i),MUDPOTS.Yp(i),MUDPOTS.Zp(i)] =  utm2xyz_porotomo(...
        MUDPOTS.E(i)...
        ,MUDPOTS.N(i)...
        ,MUDPOTS.H(i));
end
save('MUDPOTS.mat','-struct','MUDPOTS');


%% make figure in PoroTomo coordinates
nf=nf+1;figure(nf);hold on;
imagesc([min(BRADYBOX.Xp):max(BRADYBOX.Xp)],[min(BRADYBOX.Yp):max(BRADYBOX.Yp)],reshape(BRADY2DGRID.ElevWGS84,size(BRADY2DGRID.Xp)));
axis tight
axis equal
xlim(BOUNDS.Xp);
ylim(BOUNDS.Yp);
set(gca,'ClippingStyle','rectangle','Clipping','on');

plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
plot(MUDPOTS.Xp,    MUDPOTS.Yp,'kd','MarkerFaceColor','w','MarkerSize',7);
plot(FUMAROLES.Xp,FUMAROLES.Yp,'kd','MarkerFaceColor','r','MarkerSize',7);

colormap('summer');
% OPTIONS.cmap=colormap('autumn');
% OPTIONS.colorbarlabelstr = 'ElevWGS84 [m]';
colorbar;

xlabel('Xporotomo [m]');
ylabel('Yporotomo [m]');
title('Topographic Elevation (meters above WGS84 ellipsoid');
feval(funprint,'WGS84elevation');

%% Load traces
% On  2018-Nov-30, at 11:27 , Jeremy Patterson
% <jpatterson7@wisc.edu> wrote:
%
% Hi Kurt,
%
% I?ve uploaded the coordinate data for particle traces X / Y
% coordinates are in the rotated PoroTomo coordinate system,
% and the Z / Elevation coordinates are in absolute elevation
% of the WGS84 datum.
%
% Simulated Particle Trace Coordinates
% https://uwmadison.app.box.com/file/359685430200
%

dirname = '/Users/feigl/BoxSync/PoroTomo/Task9_Inverse_Modeling/Subtask9_3_Inverse_Modeling_of_Pressure_Data';
fname25b = 'DmgZoneParticleTraces.csv'; %
TRACES=table2struct(readtable(strcat(dirname,filesep,fname25b)),'ToScalar',true);
% rename the fields in the stucture
TRACES.Xp = TRACES.X;TRACES=rmfield(TRACES,'X');
TRACES.Yp = TRACES.Y;TRACES=rmfield(TRACES,'Y');
TRACES.Zp = TRACES.Z-800; % convert elevation to PoroTomo Z
TRACES=rmfield(TRACES,'Z');
nf=nf+1;figure;hold on;
for i=1:numel(unique(TRACES.ParticleNumber))
    ii=find(TRACES.ParticleNumber == i);
    plot(TRACES.Time(ii));
    xlabel('index');
    ylabel('Time [???]');
end
printpdf(sprintf('TraceTimeTags.pdf'));




%% set observable quantity
% obsq = 0 % Density from Witter et al.
% obsq = 1 % P-wave velocity from Matzel's sweep interferometry
% obsq = 2 % S-wave velocity from Matzel's sweep interferometry
% obsq = 3 % Poisson's ratio from Matzel's sweep interferometry
% obsq = 4 % Young's modulus from Matzel's sweep interferometry
% obsq = 5 % P-wave velocity from body-wave tomography
% obsq = 6 % shear-wave velocity from Zeng MASW
%            Multiscale Analysis of Surface Waves
% obsq = 7 % Quality factor Qp from Matzel's sweep interferometry
% obsq = 8 % Quality factor Qs from Matzel's sweep interferometry
% obsq = 9 % SKIP! Quality factor ratio Qp/Qs from Matzel's sweep interferometry
% obsq = 10 % P-wave velocity from body-wave tomography Thurber20171123
% obsq = 11 % Quality factor ratio Qs/Qp from Matzel's sweep interferometry
% obsq = 12 % Temperature
% obsq = 13 % Thurber20180504 Vp with model resolution
% obsq = 14 % Nayak20180513 Vp with model resolution
% obsq = 15 % Siler lithology
% obsq = 16 % Strain rates in voxels from Reinisch et al.
% obsq = 17 % Nayak2010926 Vp with model resolution from case 13
% obsq = 18 % Nayak2010607 Vp with model resolution from same inversion
% obsq = 19 % Zeng MASW on DAS TODO
% obsq = 20 % Number of faults meshlines that intersect a slice
% obsq = 21 % Modeled Temperature from regression (Training Data set)
% obsq = 22 % Modeled Temperature from regression (Testing Data set)
% obsq = 23 % Pressure from H-T modeling
% obsq = 24 % areal fault density (area per unit volume) [m^2/m^3]
% obsq = 25 % Hydraulic Conductivity

obsqmin = 0;
obsqmax = 25;

% dimension
%obsqs = 10; % for SRL paper
%obsqs = [0:15]
%obsqs = [0,10,17]; % for debugging
%obsqs = 23
%obsqs = [0:8,10:16]; % skip 9
obsqs = [0:8,10:16,18,23,24,25]; % skip 9, 17, 19, and 20, 21 and 22 are testing
%obsqs = 24; % only one
%obsqs = [0,12,16,21,22,23];
%obsqs = [25];

% this fails because of changes to ismember
if ismember(9,obsqs) == true  
    obsqs=setxor(obsqs,9); % case 9 is not defined
end

% number of observable quantities
nobsq = numel(obsqs)
% dimension array for storing vertical profiles at wells
vertprofvals=nan(nwells,nobsq,numel(zprof));
%% dimension array for storing 3-D tomograms, interpolated to regular grid
% named MESH3
%MESHEDTOMO.tomo = nan(nobsq,MESH3.nx,MESH3.ny,MESH3.nz);
% for i=1:numel(obsqs)
%     tomoid = sprintf('tomoid%02d',obsqs(i));
%     MESHEDTOMO.(tomoid) = nan(MESH3.nx,MESH3.ny,MESH3.nz);
% end



%% loop over observable quantity
for obsq = obsqs
    % initialize a few things
    OPTIONS.resolution_threshold = nan;
    %% read the files   
    switch obsq
         case 0
            %% Get Density from Witter et al.
            % Witter, J. B., D. L. Siler, J. E. Faulds, and N. H. Hinz (2016), 3D
            % geophysical inversion modeling of gravity data to test the 3D geologic
            % model of the Bradys geothermal area, Nevada, USA, Geotherm Energy, 4, 14.
            % http://dx.doi.org/10.1186/s40517-016-0056-6
            % Three-dimensional geophysical
            % inversion modeling of gravity data has been performed to test the
            % validity of a 3D geologic model constructed for the Bradys geothermal
            % area. Geophysical modeling was implemented in three different ways: (1)
            % fully unconstrained (i.e., no geologic data included); (2) constrained by
            % the 3D geologic model using homogeneous rock unit densities, and (3)
            % constrained by the 3D geologic model using heterogeneous rock unit
            % densities. We show that the existing 3D geologic model of the Bradys area
            % is broadly consistent with the gravity data. At a more detailed level,
            % however, our analysis suggests that some adjustments to the Bradys 3D
            % geologic model would improve agreement between the observed gravity and
            % the calculated gravity response. The results of the geophysical inversion
            % modeling are important as they serve as a guide to show where and how the
            % boundaries of the 3D geologic model may need to be adjusted to address
            % density excesses and deficiencies. A 3D geologic model that has been
            % independently tested prior to drilling (using a method such as that
            % described in this paper) will be more robust and have less uncertainty
            % than those which have not been tested. Such an approach will facilitate a
            % reduction in drilling risk, lead to more successful drilling programs,
            % and provide valuable geologic input to improve the accuracy of reservoir
            % models.
            %/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_4_Adjoint_seismic_tomography/MESH_SPECFEM3D/
            
            % Read the file containing density values
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/WitterDensity'
            %fname_witter_density = 'Bradys_Output_Density_gcm3_UTMNAD83_zone11_with_NDV.csv'; % with NoDataValues set to 9999
            fname_witter_density = 'Bradys_Output_Density_gcm3_UTMNAD83_zone11.csv'; % omits missing data values
            T=readtable(strcat(dirname,filesep,fname_witter_density));
            DENSITY = table2struct(T,'ToScalar',true);
            
            % Transform coordinates into rotated PoroTomo frame (Xp,Yp,Zp)
            [DENSITY.Xp,DENSITY.Yp,DENSITY.Zp] = utm2xyz_porotomo(DENSITY.X,DENSITY.Y,DENSITY.Z);
            save('DENSITY.mat','-struct','DENSITY');

        case {1,2,3,4,7,8,11}
            if exist('MATZEL.mat','file') == 2
                MATZEL = load('MATZEL.mat');
            else
                
                %% Read the file for Matzel model
                % Vp, Vs from Matzel et al.
                % /Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_3_Seismic_ANT
                
                % PoroTomo/Task8_Analyze_Data_Collected/Subtask8_3_Seismic_ANT/PoroTomo.SweepInterferometry.3D_Vp_Vs.zip
                %
                % This contains the result of 3D tomography. I used the Green's functions
                % estimated using the "sweep interferometry" method and inverted for Vp and
                % Vs from the surface down to several kilometers. Be aware that resolution
                % decreases with depth. I haven't tested the rigorously, yet, but expect
                % the model to be good down to ~ 500 m.
                %
                % When you unzip the file you'll see 40 ascii files with names like:
                %
                % absoluteseis_VP_0.1km.xyz
                %
                % "absoluteseis" simply means this is the actual value of Vp, as opposed to
                % a differential.
                %
                % The file above is for a layer a 100 meters below the surface.
                %
                % Each line of the file is an xyz point (longitude, latitude, velocity
                % (km/s)) at that depth.
                %
                % Note I haven't introduced topography. My default is to use the model as
                % is, with depth relative to the surface topography ( so the layers are not
                % flat ), but that may not be the optimal solution.
                %
                % -Eric
                %
                %dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_3_Seismic_ANT/'
                %fname_matzel_sweep_interf = 'Bradys_Output_Density_gcm3_UTMNAD83_zone11_with_NDV.csv'; % with NoDataValues set to 9999
                %fname_matzel_sweep_interf = 'Porotomo_Sweep_Interferometry_Model.txt'; % omits missing data values
                %         T=readtable(strcat(dirname,fname_matzel_sweep_interf));
                %         MATZEL = table2struct(T,'ToScalar',true);
                
                % Transform coordinates into rotated PoroTomo frame (Xp,Yp,Zp)
                %[MATZEL.E,MATZEL.N,MATZEL.utmzone] = deg2utm(MATZEL.latitude,MATZEL.longitude);
                % should read DEM instead....
                %MATZEL.elevation = 1250-1000*MATZEL.depth_km_;
                %[MATZEL.Xp,MATZEL.Yp,MATZEL.Zp] = utm2xyz_porotomo(MATZEL.E,MATZEL.N,MATZEL.elevation);
                
                %% new version of November 2017, including Qs
                dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_3_Seismic_ANT/'
                fname_matzel_sweep_interf = 'Porotomo_Sweep_Interferometry_Model.version_Nov2017.XYcoordinates.txt'; % omits missing data values
                T=readtable(strcat(dirname,fname_matzel_sweep_interf),'ReadVariableNames',true,'HeaderLines',0);
                MATZEL = table2struct(T,'ToScalar',true);
                
                %         xporo(m)		yporo(m)	z(depth)	vp(m/s)		vs(m/s)		  qp			qs
                %         0.000 	     0.000 	       0.0 	    1245.3 	     275.1 	      21.2 	      33.0
                %         10.000 	     0.000 	       0.0 	    1276.6 	     283.2 	      22.0 	      30.6
                %         20.000 	     0.000 	       0.0 	    1274.2 	     283.0 	      22.0 	      30.7
                %         30.000 	     0.000 	       0.0 	    1257.5 	     278.4 	      21.6 	      31.5
                %         40.000 	     0.000 	       0.0 	    1258.7 	     279.7 	      22.0 	      31.4
                
                % parse columns into fields
                MATZEL=rmfield(MATZEL,'Var8') % remove 8th column full of NaNs
                MATZEL.Z = 1250-MATZEL.z_depth_;
                MATZEL.Xp = MATZEL.xporo_m_;MATZEL=rmfield(MATZEL,'xporo_m_');
                MATZEL.Yp = MATZEL.yporo_m_;MATZEL=rmfield(MATZEL,'yporo_m_');
                MATZEL.Zp = MATZEL.Z - 800.;
                [MATZEL.E,MATZEL.N] = xy_porotomo2utm(MATZEL.Xp,MATZEL.Yp);
                
                MATZEL.Vp = MATZEL.Vp_m_s_;MATZEL=rmfield(MATZEL,'Vp_m_s_');
                MATZEL.Vs = MATZEL.Vs_m_s_;MATZEL=rmfield(MATZEL,'Vs_m_s_');
                
                %% find latitude and longitude
                utmzone = '11 S' % UTM zone for Brady Hot Springs
                [MATZEL.Lat,MATZEL.Lon] = utm2deg(colvec(MATZEL.E),colvec(MATZEL.N),repmat(utmzone,size(colvec(MATZEL.E))));
                
                
                %% make figure in UTM coordinates
                nf=nf+1;figure(nf);hold on;
                plot(MATZEL.E,MATZEL.N,'r.');
                %plot(bradybox_UTM([1,2,3,4,1],1),bradybox_UTM([1,2,3,4,1],2),'k-');
                plot(BRADYBOX.E([1,2,3,4,1]),BRADYBOX.N([1,2,3,4,1]),'k-');
                xlabel('UTM Easting [m]');
                ylabel('UTM Northing [m]');
                axis equal;axis tight
                title('UTM coordinates MATZEL');
                legend('model nodes','PoroTomo Box');
                printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
                
                %% make figure in PoroTomo coordinates
                nf=nf+1;figure(nf);hold on;
                plot(MATZEL.Xp,MATZEL.Yp,'r.');
                %plot(BRADYBOX.Xporotomo([1,2,3,4,1]),BRADYBOX.Yporotomo([1,2,3,4,1]),'k-');
                plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
                
                xlabel('Xporotomo [m]');
                ylabel('Yporotomo [m]');
                axis equal; axis tight;
                title('PoroTomo coordinates MATZEL');
                legend('model nodes','PoroTomo Box');
                printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
                
                               
                %% make figure in UTM coordinates
                nf=nf+1;figure(nf);hold on;
                plot(DENSITY.X,DENSITY.Y,'r.');
                %plot(bradybox_UTM([1,2,3,4,1],1),bradybox_UTM([1,2,3,4,1],2),'k-');
                xlabel('UTM Easting [m]');
                ylabel('UTM Northing [m]');
                axis equal;axis tight
                title('UTM coordinates: DENSITY');
                legend('model nodes','PoroTomo Box');
                printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
                
                % make figure in PoroTomo coordinates
                nf=nf+1;figure(nf);hold on;
                plot(DENSITY.Xp,DENSITY.Yp,'r.');
                plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
                xlabel('Xporotomo [m]');
                ylabel('Yporotomo [m]');
                axis equal; axis tight;
                title('PoroTomo coordinates: DENSITY');
                legend('model nodes','PoroTomo Box');
                printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
                
                % interpolate values of density on Matzel's grid
                fprintf(1,'Starting 3-D scatteredInterpolant at t = %.001f seconds\n',toc(tstart));
                interpolation_method = 'linear';
                extrapolation_method = 'none';
                % Finterp2 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V ...
                %     ,interpolation_method,extrapolation_method);
                % image2d = Finterp2(Q.x,Q.y,Q.z);
                Finterp3 = scatteredInterpolant(DENSITY.Xp,DENSITY.Yp,DENSITY.Zp,1000*DENSITY.Density_gcm3 ...
                    ,interpolation_method,extrapolation_method);
                MATZEL.Density_kg_per_m3 = Finterp3(MATZEL.Xp,MATZEL.Yp,MATZEL.Zp);
                
                %% get the topographic elevation of the nodes in the model
                MATZEL.Zpt = get_elevation(MATZEL.Lon,MATZEL.Lat,demgrdfilename) - MATZEL.z_depth_ - 800;
                
                %% Calculate Elastic moduli
                %https://en.wikipedia.org/wiki/Elastic_modulus
                MATZEL.VpOverVs = MATZEL.Vp./MATZEL.Vs;
                MATZEL.Poisson = 0.5*((MATZEL.VpOverVs.^2) - 2)./((MATZEL.VpOverVs.^2)-1);
                
                % https://en.wikipedia.org/wiki/P-wave_modulus
                % In linear elasticity, the P-wave modulus M {\displaystyle M} M, also
                % known as the longitudinal modulus or the constrained modulus, is one of
                % the elastic moduli available to describe isotropic homogeneous materials.
                %
                % It is defined as the ratio of axial stress to axial strain in a uniaxial
                % strain state
                %
                %     ? z z = M ? z z {\displaystyle \sigma _{zz}=M\epsilon _{zz}}
                %     {\displaystyle \sigma _{zz}=M\epsilon _{zz}}
                %
                % where all the other strains ? ? ? {\displaystyle \epsilon _{**}}
                % {\displaystyle \epsilon _{**}} are zero.
                %
                % This is equivalent to stating that
                %
                %     M = ? V P 2 {\displaystyle M=\rho V_{\mathrm {P} }^{2}}
                %     {\displaystyle M=\rho V_{\mathrm {P} }^{2}}
                %
                % where VP is the velocity of a P-wave.
                %density = MATZEL.Density_kg_per_m3;
                Mp = MATZEL.Density_kg_per_m3 .* (MATZEL.Vp).^2;
                
                % Young's modulus
                MATZEL.EYoung  = Mp .* (1 + MATZEL.Poisson) .* (1 - 2*MATZEL.Poisson) ./ (1 - MATZEL.Poisson);
                clear Mp;
                
                % calculate Qp/Qs ratio
                MATZEL.QpOverQs = MATZEL.Qp ./ MATZEL.Qs;
                % calculate Qs/Qp ratio
                MATZEL.QsOverQp = MATZEL.Qs ./ MATZEL.Qp;
                
                save('MATZEL.mat','-struct','MATZEL');
            end
        case 5
            
            %% P-wave velocity from Lesley Parker and Cliff Thurber
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff'
            tomo_filename = 'velfile20170727.txt'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            iok = find(XYZVthurber(:,4) > 0); % find indices of non-zero velocities
            TOMO.Xthurber = XYZVthurber(iok,1)*1000; % meters long  axis positive toward SW
            TOMO.Ythurber = XYZVthurber(iok,2)*1000; % meters short axis positive toward NW
            TOMO.Zthurber = XYZVthurber(iok,3)*1000; % meters above WGS84 ellipsoid
            TOMO.V        = XYZVthurber(iok,4)*1000; % velocity in meters per second
            
            %Rotate coordinates from Thurber to UTM
            % From: Clifford Thurber <clifft@geology.wisc.edu>
            % Date: 2017- Jul-21 Friday at 11:40
            % To: Kurt Feigl <feigl@wisc.edu>
            % Subject: Origin
            %
            %   origin :  latitude   longitude   rotation
            %              39 48.060   119  0.648    53.580
            %[x,y,utmzone] = deg2utm(Lat,Lon)
            [ORIGIN.Eutm,ORIGIN.Nutm,ORIGIN.utmzone] = deg2utm(39+48.060/60,-119. - 0.648/60);
            ORIGIN.Zutm = 0; % elevation in meters above WGS84 ellipsoid
            rotang = 53.58*pi/180; % rotation angle in degrees
            
            % Rotation matrix
            Rz = [cos(rotang) -sin(rotang); sin(rotang) cos(rotang)];
            
            U=nan(2,1); % vector in Thurber coordinates
            TOMO.Eutm = nan(size(TOMO.Xthurber));
            TOMO.Nutm = nan(size(TOMO.Ythurber));
            for i=1:numel(TOMO.Xthurber)
                U(1) =  -1*TOMO.Xthurber(i);
                U(2) =     TOMO.Ythurber(i);
                V = Rz * U;
                TOMO.Eutm(i) = V(1) + ORIGIN.Eutm;
                TOMO.Nutm(i) = V(2) + ORIGIN.Nutm;
            end
            
            %% rotate coordinates from UTM to PoroTomo
            [TOMO.Xp,TOMO.Yp,TOMO.Zp] = utm2xyz_porotomo(TOMO.Eutm,TOMO.Nutm,-1*TOMO.Zthurber);
            [ORIGIN.Xp,ORIGIN.Yp,ORIGIN.Zp] = utm2xyz_porotomo(ORIGIN.Eutm,ORIGIN.Nutm,ORIGIN.Zutm);
            
            %% make figure in Thurber coordinates
            nf=nf+1;figure(nf);hold on;axis xy;
            %plot(TOMO.Xthurber,TOMO.Ythurber,'r*');
            plot(-1*TOMO.Ythurber,TOMO.Xthurber,'r*');
            axis ij;
            plot(0.,0.,'ob');
            ylabel('Xthurber [m]');
            xlabel('-1*Ythurber [m]');
            axis equal;axis tight;
            title('Thurber coordinates');
            legend('model nodes','Thurber origin');
            printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
            
            %% make figure in UTM coordinates
            nf=nf+1;figure(nf);hold on;
            plot(TOMO.Eutm,TOMO.Nutm,'r*');
            %plot(bradybox_UTM([1,2,3,4,1],1),bradybox_UTM([1,2,3,4,1],2),'k-');
            plot(BRADYBOX.E([1,2,3,4,1]),BRADYBOX.N([1,2,3,4,1]),'k-');
            plot(ORIGIN.Eutm,ORIGIN.Nutm,'ob');
            xlabel('UTM Easting [m]');
            ylabel('UTM Northing [m]');
            axis equal;axis tight
            title('UTM coordinates');
            legend('model nodes','PoroTomo Box','Thurber origin');
            printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
            
            %% make figure in PoroTomo coordinates
            nf=nf+1;figure(nf);hold on;
            plot(TOMO.Xp,TOMO.Yp,'r*');
            %plot(BRADYBOX.Xporotomo([1,2,3,4,1]),BRADYBOX.Yporotomo([1,2,3,4,1]),'k-');
            plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
            plot(ORIGIN.Xp,ORIGIN.Yp,'ob');
            plot(WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,'k^','MarkerSize',2,'MarkerFaceColor','k');
            xlabel('Xporotomo [m]');
            ylabel('Yporotomo [m]');
            axis equal; axis tight;
            title('PoroTomo coordinates');
            legend('model nodes','PoroTomo Box','Thurber origin');
            printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
                        printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
            save('VpThurber20170727.mat','-struct','TOMO');          
         case 10
            
            %% P-wave velocity from Cliff Thurber
            %             From: Clifford Thurber
            % Sent: Thursday, November 23, 2017 7:50 AM To: Kurt Feigl
            % Subject: When you have time
            %
            % Hi Kurt.  I decided to take a step sort of backwards and
            % determine a smoother model using 100 m horizontal node
            % spacing, and adjusting faster the deeper part of the starting
            % model a bit.  The result is attached.  Please plot with your
            % nice faults!  The annoyance is it is in my coordinate system.
            % So Z is w.r.t. WGS84 zero, negative means above WGS84 zero
            % (positive elevation, negative depth), and the signs on X and
            % Y are the opposite of PoroTomo.
            
            
            % From: Clifford Thurber Sent: Thursday, November 23, 2017
            % 10:07 AM To: Kurt Feigl Subject: Fwd: PoroTomo Z
            %
            % Well, if it is correct that PoroTomo Z is positive upwards
            % and PoroTomo Z=0 is at WGS84 elevation of 800 m, then the
            % attached Vp model file is in correct PoroTomo coordinates.
            %
            % Cliff
          
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time'
            tomo_filename = 'VpThurber20171123.vel'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            iok = find(XYZVthurber(:,4) > 0); % find indices of non-zero velocities
            TOMO.Xp = XYZVthurber(iok,1)*1000; % meters short  axis positive toward SE
            TOMO.Yp = XYZVthurber(iok,2)*1000; % meters short axis positive toward  NE
            TOMO.Zp = XYZVthurber(iok,3)*1000; % meters above 800 m
            TOMO.V  = XYZVthurber(iok,4)*1000; % velocity in meters per second
            
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
                                  
%             % make figure in PoroTomo coordinates
%             nf=nf+1;figure(nf);hold on;
%             plot(TOMO.Xp,TOMO.Yp,'r*');
%             %plot(BRADYBOX.Xporotomo([1,2,3,4,1]),BRADYBOX.Yporotomo([1,2,3,4,1]),'k-');
%             plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
% %            plot(ORIGIN.Xp,ORIGIN.Yp,'ob');
% %            plot(WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,'k^','MarkerSize',2,'MarkerFaceColor','k');
%             xlabel('Xporotomo [m]');
%             ylabel('Yporotomo [m]');
%             axis equal; axis tight;
%             title('PoroTomo coordinates');
%             legend('model nodes','PoroTomo Box','Thurber origin');
%             printpdf(sprintf('%s_%03d.pdf',mfilename,nf));
            save('VpThurber20171123.mat','-struct','TOMO');
         case 13           
             %% P-wave velocity with model resolution from Cliff Thurber
             %
             % From: Clifford Thurber Sent: Friday, May 4, 2018 4:00 PM To: Kurt Feigl
             % Subject: Fwd: PoroTomo Z
             %
             % Here is the resolution diagonals, modified to be in the exact same form
             % as the Vp model.  Resolution is not so great so I recommend a threshold
             % of 0.1
             %
             % An inversion with reduced damping result in much better resolution and a
             % somewhat rougher model, but that's not what we submitted.
             %
             % Cliff
             %
             % Begin forwarded message:
             %
             % From: Clifford Thurber <clifft@geology.wisc.edu> Subject: Fwd: PoroTomo Z
             % Date: November 23, 2017 at 10:07:35 AM CST To: Kurt Feigl
             % <feigl@wisc.edu>
             %
             % Well, if it is correct that PoroTomo Z is positive upwards and PoroTomo
             % Z=0 is at WGS84 elevation of 800 m, then the attached Vp model file is in
             % correct PoroTomo coordinates.
             %
             % Cliff
             %

          
            %dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time'
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Thurber20180504PoroTomo_Z'
            % read velocity
            tomo_filename = 'velfile.porotomo'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            % read model resolution
            mres_filename = 'resfile.porotomo'
            XYZRthurber = load(strcat(dirname,filesep,mres_filename));
            
            % check for consistency
            DX = XYZVthurber(:,1) - XYZRthurber(:,1);
            DY = XYZVthurber(:,2) - XYZRthurber(:,2);
            DZ = XYZVthurber(:,3) - XYZRthurber(:,3);
            disdiff = sqrt(DX.^2 + DY.^2 + DZ.^2);
            idiff = find(disdiff > 1.0);
            if numel(idiff) == 0
                clear DX DY DZ;
            else
                figure
                plot(XYZVthurber(:,3),XYZRthurber(:,3),'ko');
                xlabel('Z from velfile');
                ylabel('Z from resfile');
                printpdf('Zdiff.pdf');
                error('coordinate mismatch');
            end

            iok = find(XYZVthurber(:,4) > 0); % find indices of non-zero velocities
            TOMO.Xp = XYZVthurber(iok,1)*1000; % meters short  axis positive toward SE
            TOMO.Yp = XYZVthurber(iok,2)*1000; % meters short axis positive toward  NE
            TOMO.Zp = XYZVthurber(iok,3)*1000; % meters above 800 m
            TOMO.V  = XYZVthurber(iok,4)*1000; % velocity in meters per second
            TOMO.R  = XYZRthurber(iok,4)*1000; % model resolution (dimensionless)
            
         
            
            %            % clip out velocity values where resolution is poor
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %             OPTIONS.resolution_threshold = 0.05;  % threshhold value
            
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
            
            save('VpThurber20171123.mat','-struct','TOMO');
        case 14
            %% P-wave velocity with from Avinash Nayak % From: Avinash
            %              Nayak <avinash07guddu@gmail.com>
            % Sent: Sunday, May 13, 2018 10:56 AM To: Kurt Feigl Cc: Lesley Parker;
            % Clifford Thurber; Elena Reinisch Subject: Re: [porotomo] Routines for
            % plotting cross sections
            %
            % Sorry I didn't find the opportunity to plot with the new code. Here is
            % the how the velocity model slices look now. I was able to plot the
            % injection well locations and the brady box. The injection wells appear to
            % go through the fault zone and tap into the high velocity block between
            % the faults. Thanks to Jeremy for bringing this to attention. Here are the
            % old and new models. But I think I would prefer to not show anything with
            % resolution till we can get the checkerboard tests and the final model
            % (after adding revised DASv picks. I haven't reprocessed the DASv picks
            % yet- do be done after SSA). Thanks. Avinash
            
            
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Nayak20180513'
            % read velocity
            tomo_filename = 'vp_new.txt'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            
            %  Here is the resolution file. SIMUL code at this point doesn't output
            % resolution of the last depth slice. It is a bug. The actual elevation is
            % the Z value + 800 m in these files. Thanks. Avinash
            [nvel,nfour] = size(XYZVthurber)
            
            % read model resolution
            mres_filename = 'resolution_new.txt'
            XYZRthurber = load(strcat(dirname,filesep,mres_filename));
            [nres,nfour] = size(XYZRthurber)

            
            % check for consistency
            if nres < nvel
                warning('padding missing values in resolution');
                for i=nres+1:nvel
                    XYZRthurber(i,1) = XYZVthurber(i,1);
                    XYZRthurber(i,2) = XYZVthurber(i,2);
                    XYZRthurber(i,3) = XYZVthurber(i,3);
                    XYZRthurber(i,4) = 0.0;  % set resolution to zero in "last" slice. Is it top or bottom? 
                end
            end
         
            DX = XYZVthurber(:,1) - XYZRthurber(:,1);
            DY = XYZVthurber(:,2) - XYZRthurber(:,2);
            DZ = XYZVthurber(:,3) - XYZRthurber(:,3);
            disdiff = sqrt(DX.^2 + DY.^2 + DZ.^2);
            idiff = find(disdiff > 1.0);
            if numel(idiff) == 0
                clear DX DY DZ;
            else
                figure
                plot(XYZVthurber(:,3),XYZRthurber(:,3),'ko');
                xlabel('Z from velfile');
                ylabel('Z from resfile');
                printpdf('Zdiff.pdf');
                error('coordinate mismatch');
            end

            iok = find(XYZVthurber(:,4) > 0); % find indices of non-zero velocities
            TOMO.Xp = XYZVthurber(iok,1)*1000; % meters short  axis positive toward SE
            TOMO.Yp = XYZVthurber(iok,2)*1000; % meters short axis positive toward  NE
            TOMO.Zp = XYZVthurber(iok,3)*1000; % meters above 800 m
            TOMO.V  = XYZVthurber(iok,4)*1000; % velocity in meters per second
            TOMO.R  = XYZRthurber(iok,4)*1000; % model resolution (dimensionless)
            
            %clip out velocity values where resolution is poor
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %OPTIONS.resolution_threshold = 0.05;  % threshhold value
                        
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
                                  
            save('VpNayak20180513.mat','-struct','TOMO');
            
        case 17
            %% Avinash 20180926
            %             On  2018-Oct-05, at 10:46 , Clifford             Thurber
            %             <cthurber@wisc.edu> wrote:
            %
            %             Case 14 is what needs to be re-plotted extending not as deep
            %             to match the plots of the older (11/23/17) model, case 13.
            %
            %                         The uncertainty and resolution files were from
            %                         case 13.  They would be approximately correct for
            %                         case 14.
            %
            %**read velocity from case 14 
            % From: Avinash
            %                          Nayak <avinash07guddu@gmail.com>
            %             Sent: Sunday, May 13, 2018 10:56 AM To: Kurt Feigl Cc: Lesley
            %             Parker; Clifford Thurber; Elena Reinisch Subject: Re:
            %             [porotomo] Routines for plotting cross sections
            %
            %             Sorry I didn't find the opportunity to plot with the new
            %             code. Here is the how the velocity model slices look now. I
            %             was able to plot the injection well locations and the brady
            %             box. The injection wells appear to go through the fault zone
            %             and tap into the high velocity block between the faults.
            %             Thanks to Jeremy for bringing this to attention. Here are the
            %             old and new models. But I think I would prefer to not show
            %             anything with resolution till we can get the checkerboard
            %             tests and the final model (after adding revised DASv picks. I
            %             haven't reprocessed the DASv picks yet- do be done after
            %             SSA). Thanks. Avinash
                       
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Nayak20180513'            
            tomo_filename = 'vp_new.txt'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            [nvel,nfour] = size(XYZVthurber)
            
            %** Read Resolution from Case 13
            %  Here is the resolution file. SIMUL code at this point doesn't output
            % resolution of the last depth slice. It is a bug. The actual elevation is
            % the Z value + 800 m in these files. Thanks. Avinash
            
            
            % read model resolution
            %             %% On  2018-Sep-26, at 15:40 , Avinash Nayak
            %             <avinash07guddu@gmail.com> wrote:
            %
            % Prof. Feigl, I have a request. Is it possible for you to generate figures
            % of the some slices of the velocity model (attached) for the final review
            % ? I tried to use interpolate_fields5.m, I was unable to mask the velocity
            % values with low resolution. Of if you have a few tips, on exactly what to
            % change in the code, it will be helpful. Here is the model attached, but I
            % think you have this file. Sorry for this late request. Regards, Thanks.
            % Avinash
            %
            % -- Avinash Nayak, Research Associate, Department of Geoscience,
            % University of Wisconsin, Madison <resolution_new3.txt>
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Nayak20180926'
            mres_filename = 'resolution_new3.txt'
            XYZRthurber = load(strcat(dirname,filesep,mres_filename));
            [nres,nfour] = size(XYZRthurber)

            
            % check for consistency
            if nres < nvel
                warning('padding missing values in resolution');
                for i=nres+1:nvel
                    XYZRthurber(i,1) = XYZVthurber(i,1);
                    XYZRthurber(i,2) = XYZVthurber(i,2);
                    XYZRthurber(i,3) = XYZVthurber(i,3);
                    XYZRthurber(i,4) = 0.0;  % set resolution to zero in "last" slice. Is it top or bottom? 
                end
            end
         
            DX = XYZVthurber(:,1) - XYZRthurber(:,1);
            DY = XYZVthurber(:,2) - XYZRthurber(:,2);
            DZ = XYZVthurber(:,3) - XYZRthurber(:,3);
            disdiff = sqrt(DX.^2 + DY.^2 + DZ.^2);
            idiff = find(disdiff > 1.0);
            if numel(idiff) == 0
                clear DX DY DZ;
            else
                figure
                plot(XYZVthurber(:,3),XYZRthurber(:,3),'ko');
                xlabel('Z from velfile');
                ylabel('Z from resfile');
                printpdf('Zdiff.pdf');
                error('coordinate mismatch');
            end

            iok = find(XYZVthurber(:,4) > 0); % find indices of non-zero velocities
            iok = intersect(iok,find(isfinite(XYZRthurber(:,4))==1)); % finite resolution
            TOMO.Xp = XYZVthurber(iok,1)*1000; % meters short  axis positive toward SE
            TOMO.Yp = XYZVthurber(iok,2)*1000; % meters short axis positive toward  NE
            TOMO.Zp = XYZVthurber(iok,3)*1000; % meters above 800 m
            TOMO.V  = XYZVthurber(iok,4)*1000; % velocity in meters per second
            TOMO.R  = XYZRthurber(iok,4)     ; % model resolution (dimensionless)
            
            %            % clip out velocity values where resolution is poor
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %             OPTIONS.resolution_threshold = 0.05;  % threshhold value
                        
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
            
            % clip out velocity values where resolution is poor
            figure;
            histogram(TOMO.R);
            xlabel('Model Resolution [dimensionless]');
            ylabel('Count');
            % OPTIONS.resolution_threshold = 0.1;  % threshhold value
            % OPTIONS.resolution_threshold = 0.05;  % threshhold value
            OPTIONS.resolution_threshold = quantile(colvec(TOMO.R),0.1);                                            
            save('VpNayak20180926.mat','-struct','TOMO');
        case 18
            %% Avinash 20180607
            % On  2018-Oct-05, at 16:43 , Avinash Nayak <avinash07guddu@gmail.com> wrote:
            %
            % Here the model and resolution files for the model in June 2018 that has
            % revised DASV picks. The output file will be used to extract the
            % uncertainty.
                      
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Nayak20180607'            
            tomo_filename = 'vp_20180607.txt'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            [nvel,nfour] = size(XYZVthurber)
            
            mres_filename = 'resolution_vp_20180607.txt'
            XYZRthurber = load(strcat(dirname,filesep,mres_filename));
            [nres,nfour] = size(XYZRthurber)
          
            % check for consistency
            if nres < nvel
                warning('padding missing values in resolution');
                for i=nres+1:nvel
                    XYZRthurber(i,1) = XYZVthurber(i,1);
                    XYZRthurber(i,2) = XYZVthurber(i,2);
                    XYZRthurber(i,3) = XYZVthurber(i,3);
                    XYZRthurber(i,4) = 0.0;  % set resolution to zero in "last" slice. Is it top or bottom? 
                end
            end
         
            DX = XYZVthurber(:,1) - XYZRthurber(:,1);
            DY = XYZVthurber(:,2) - XYZRthurber(:,2);
            DZ = XYZVthurber(:,3) - XYZRthurber(:,3);
            disdiff = sqrt(DX.^2 + DY.^2 + DZ.^2);
            idiff = find(disdiff > 1.0);
            if numel(idiff) == 0
                clear DX DY DZ;
            else
                figure
                plot(XYZVthurber(:,3),XYZRthurber(:,3),'ko');
                xlabel('Z from velfile');
                ylabel('Z from resfile');
                printpdf('Zdiff.pdf');
                error('coordinate mismatch');
            end

            iok = find(XYZVthurber(:,4) > 0); % find indices of non-zero velocities
            iok = intersect(iok,find(isfinite(XYZRthurber(:,4))==1)); % finite resolution
            TOMO.Xp = XYZVthurber(iok,1)*1000; % meters short  axis positive toward SE
            TOMO.Yp = XYZVthurber(iok,2)*1000; % meters short axis positive toward  NE
            TOMO.Zp = XYZVthurber(iok,3)*1000; % meters above 800 m
            TOMO.V  = XYZVthurber(iok,4)*1000; % velocity in meters per second
            TOMO.R  = XYZRthurber(iok,4)     ; % model resolution (dimensionless)
            
            %            % clip out velocity values where resolution is poor
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %             OPTIONS.resolution_threshold = 0.05;  % threshhold value
                        
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
            
            % clip out velocity values where resolution is poor
            figure;
            histogram(TOMO.R);
            xlabel('Model Resolution [dimensionless]');
            ylabel('Count');
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %OPTIONS.resolution_threshold = 0.05;  % threshhold value
            %OPTIONS.resolution_threshold = quantile(colvec(TOMO.R),0.1);                                            
            save('VpNayak20180926.mat','-struct','TOMO');
        case 6
            %% ZengMASWonDAS
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_1_DAS_Quality_Control/SurfaceWave'
            % read the coordinates of the DAS segments in rotated PoroTomo coordinates
            locf='seg_center_loc.dat';
            fid=fopen(strcat(dirname,filesep,locf),'r');
            all=textscan(fid,'%s %f %f');
            nseg=length(all{1});
            segx=all{2};
            segy=all{3};
            
            % read the locations of the points in the model
            load(strcat(dirname,filesep,'das_loc.mat'));
            % prune
            iok = intersect(find(das_x>0),find(das_y>0));
            % plot as is
            nf=nf+1;figure(nf);
            plot(das_x(iok),das_y(iok),'r+');
            xlabel('UTM Easting from das\_x  [m]');
            ylabel('UTM Northing from das\_y [m]');
            axis equal; axis tight;
            title('das\_loc.mat');
            printpng(sprintf('%s_%03d.png',mfilename,nf));
            
            
            % rotate into PoroTomo coordinate system
            [xx,yy]=utm2xy_porotomo(das_x(51:8671),das_y(51:8671));
            nf=nf+1;figure(nf);
            plot(xx,yy,'b+');
            xlabel('PoroTomo Xp from xx [m]');
            ylabel('PoroTomo Yp from yy [m]');
            axis equal; axis tight;
            title('after rotation');
            printpng(sprintf('%s_%03d.png',mfilename,nf));
            
            
            % find latitude and longitude
            utmzone = '11 S' % UTM zone for Brady Hot Springs
            [das_lat,das_lon] = utm2deg(colvec(das_x),colvec(das_y),repmat(utmzone,size(colvec(das_x))));
            iok = find(abs(das_lat)>0);
            iok = intersect(iok,find(abs(das_lon)> 0));
            iok = intersect(iok,find(abs(das_lat)<80));
            nf=nf+1;figure(nf);
            plot(das_lon(iok),das_lat(iok),'k+');
            xlabel('Longitude');
            ylabel('Latitude');
            axis equal; axis tight;
            title('Geographic coordinates');
            printpng(sprintf('%s_%03d.png',mfilename,nf));
            
            
            % get the topographic elevation of the nodes in the model
            % TODO: should use webget to ftp://roftp.ssec.wisc.edu/porotomo/PoroTomo/METADATA/
            demgrdfilename = '/Users/feigl/BoxSync/PoroTomo/METADATA/brady_dem_srtm_20m.grd'
            das_z = get_elevation(das_lon,das_lat,demgrdfilename);
            elevation_mean = nanmean(das_z)
            nf=nf+1;figure(nf);
            histogram(das_z);
            xlabel('elevation [m]');
            ylabel('number of occurrences');
            
            % get 3-D coordinates of the nodes in the model in PoroTomo coordinates
            [das_xp,das_yp,das_zp] = utm2xyz_porotomo(das_x,das_y,das_z);
            
            % set up interpolation function to return Zp coordinate at a
            % given location (Xp, Yp)
            Fget_zp = scatteredInterpolant(das_xp,das_yp,das_zp);
            
            
            % read the velocity values from the many files
            %dirname = strcat(dirname,filesep,'seg.model');
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_1_DAS_Quality_Control/SurfaceWave2/newmodel';
            postfix='mod.dat';
            %nlay=27;
            nlay=11; % updated 20171211
            
            % here are the first few lines
            % seg_0104_0191.mod.dat
            %
            % 0.0 10.0 309.7
            % 10.0 10.0 420.2
            % 20.0 10.0 499.1
            % 30.0 10.0 527.9
            % 40.0 10.0 536.0
            % 50.0 10.0 537.7
%           here are the first few lines of 20171208
%          0         17.20        241.20
%          17.20         15.10        115.10
%          32.30          9.80        240.40
%          42.10         16.20        303.50
%          58.40         17.30        587.30
%          75.60         12.70        499.10
%          88.30         17.20        333.20
%         105.60         14.00        515.10
%         119.50         13.20        655.00
%         132.70         16.30        741.90
%         149.00         51.00        741.90

% assume   topdepth[m]      bottomdepth[m]  Vs[m]

            
            % loop over segments
            
            tseg1=zeros(nlay,1);
            bseg1=zeros(nlay,1);
            vseg1=zeros(nlay,1);
            dsegs=zeros(nlay,nseg);
            vsegs=zeros(nlay,nseg);
            zsegs=zeros(nlay,nseg);
            
            seg_flag=zeros(1,nseg);
            
            
            for k=1:nseg
                %modf=sprintf('%s.%s',all{1}{k},postfix);
                modf=sprintf('%s%s%s.%s',dirname,filesep,all{1}{k},postfix);
                if(exist(modf,'file') == 2)
                    seg_flag(k)=1;
                    vtmp=load(modf);
%                     %zseg(:,k)=vtmp(:,1); % first column contains depth [ in meters] ** below what elevation?
%                     zsegs(:,k)=elevation_mean-vtmp(:,1)-800.; % PoroTomo coordinate
%                     %tseg(:,k)=vtmp(:,2); % what is in second column ?? **
%                     vsegs(:,k)=vtmp(:,3); % third column contains velocity [in m/s]
%                     update 20171211
                    tseg1=vtmp(:,1); % first  column contains depth [ in meters] to top     of layer
                    bseg1=vtmp(:,2); % second column contains depth [ in meters] to botttom of layer
                    %zsegs(:,k)=elevation_mean - 800.0 - (tseg1+bseg1)/2.; % PoroTomo vertical coordinate
                    dsegs(:,k)=(tseg1+bseg1)/2.; % set depth to be mean of top and bottom
                    vsegs(:,k)=vtmp(:,3); % third column contains velocity [in m/s]
                end
            end
            
            %% find good segments
            isegx=find(seg_flag>0);
            
            %% plot in 3 dimensions
            nf=nf+1;figure(nf);
            hold on;
            %h1=subplot(1,3,1);
            disksize = 200;            
            k1 = 1;
            nn = numel(isegx);
            for iz=1:10
                TOMO.Xp(k1:k1+nn-1) = colvec(segx(isegx));
                TOMO.Yp(k1:k1+nn-1) = colvec(segy(isegx));
                % 20171211 - get Zp coordinate of center of segment              
                zp_cable1 = Fget_zp(colvec(segx(isegx)),colvec(segy(isegx)));
                % subtract depth to obtain PoroTomo Zp coordinate
                TOMO.Zp(k1:k1+nn-1) = zp_cable1-colvec(dsegs(iz,isegx));
                TOMO.V(k1:k1+nn-1)  = colvec(vsegs(iz,isegx));
                k1 = k1+nn;
                scatter3(segx(isegx),segy(isegx),-1*dsegs(iz,isegx),disksize,vsegs(iz,isegx),'filled');
            end
            TOMO.Xp = colvec(TOMO.Xp);
            TOMO.Yp = colvec(TOMO.Yp);
            TOMO.Zp = colvec(TOMO.Zp);
            TOMO.V = colvec(TOMO.V);            
            %axis equal;
            axis tight; grid on;
            colormap(flipud(jet));
            %caxis([200,600]);
            xlabel('Xp [m]');
            ylabel('Yp [m]');
            zlabel('negative depth [m]');
            h1=gca;
            h1.FontSize=20;
            h1=colorbar;
            h1.Label.String='Vs (m/s)';
            
            %         view([1,0,0]);printpng(sprintf('%s_view100.png',mfilename));
            %         view([0,1,0]);printpng(sprintf('%s_view010.png',mfilename));
            %         view([0,0,1]);printpng(sprintf('%s_view001.png',mfilename));
            view([1,1,1]);printpng(sprintf('%s_view111.png',mfilename));
            
            %% make a histogram
            figure;
            histogram(colvec(TOMO.V));
            xlabel('Vs [m/s]');
            ylabel('number of occurences');
            title('Zeng Vs MASW (iobsq = 6) before clipping','Interpreter','none');
            feval(funprint,sprintf('%s_histogram','ZengMASW_before_clipping'));
            
%             %% keep only reasonable values
%             iok = find(TOMO.V > 100); % m/s
%             iok = intersect(iok,find(TOMO.V < 800)); % m/s
%             TOMO.Xp = colvec(TOMO.Xp(iok));
%             TOMO.Yp = colvec(TOMO.Yp(iok));
%             TOMO.Zp = colvec(TOMO.Zp(iok));
%             TOMO.V  = colvec(TOMO.V(iok)); 
            
            %% clip low and high values
            ilo = find(TOMO.V < 100); % m/s         
            TOMO.V(ilo)  = 100; 
            ihi = find(TOMO.V > 800); % m/s         
            TOMO.V(ihi)  = 800; 
            
            %% make a histogram
            figure;
            histogram(colvec(TOMO.V));
            xlabel('Vs [m/s]');
            ylabel('number of occurences');
            title('Zeng Vs MASW (iobsq = 6) after clipping','Interpreter','none');
            feval(funprint,sprintf('%s_histogram','ZengMASW_after_clipping'));

            %% save the thing
            save('ZengMASWonDAS','-struct','TOMO');           
        case {12,23}
            %% Read the file containing Pressure and Temperature Values
            url = 'https://gdr.openei.org/files/939'
            data_filename='PTfield_brady_for_GDR_PT_20150324.csv'
            websave(data_filename,strcat(url,filesep,data_filename));
            table=readtable(data_filename);
            PT=table2struct(table,'ToScalar',true);
            
            % rename the fields in the stucture
            PT.Xp = PT.XPoroTomoInMeters;PT=rmfield(PT,'XPoroTomoInMeters');
            PT.Yp = PT.YPoroTomoInMeters;PT=rmfield(PT,'YPoroTomoInMeters');
            PT.Zp = PT.ZPoroTomoInMeters;PT=rmfield(PT,'ZPoroTomoInMeters');
            PT.E  = PT.UTMeastingInMeters; PT=rmfield(PT,'UTMeastingInMeters');
            PT.N  = PT.UTMnorthingInMeters;PT=rmfield(PT,'UTMnorthingInMeters');
            %
            % make figure in UTM coordinates
            nf=nf+1;figure(nf);hold on;
            plot(PT.E,PT.N,'r.');
            plot(BRADYBOX.E([1,2,3,4,1]),BRADYBOX.N([1,2,3,4,1]),'k-');
            xlabel('UTM Easting [m]');
            ylabel('UTM Northing [m]');
            axis equal;axis tight
            title('UTM coordinates');
            legend('model nodes','PoroTomo Box');
            printpdf(sprintf('PT_UTM.pdf'));
           
            % make figure in PoroTomo coordinates
            nf=nf+1;figure(nf);hold on;
            plot(PT.Xp,PT.Yp,'r.');
            plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
            plot(WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,'k^','MarkerSize',2,'MarkerFaceColor','k');
            xlabel('Xporotomo [m]');
            ylabel('Yporotomo [m]');
            axis equal; axis tight;
            title('PoroTomo coordinates');
            legend('model nodes','PoroTomo Box');
            printpdf(sprintf('PT_XpYp.pdf'));         
            save('PT.mat','-struct','PT');
        case 15
            %             %% Get Lithology from Siler et al.
            %             Folder:
            %
            % https://uwmadison.box.com/s/zjxp79g2aeqssem9jbjszmmwy9yzvxt1
            %
            % /Users/feigl/Box Sync/PoroTomo/Task4_Design_Deployment/Subtask4_9_Incorporate_Geologic_Information/SilerBradysGeoModelData
            %
            % File: SilerLithologyKey.csv
            % contains: 24 lithologies as defined by Drew Siler:
            %
            % ID, Symbol, Lithology
            % 1,  Q, Quaternary sediments, undifferentiated
            % 2,  Tsy, Tertiary sediments
            % 3,  Tsl6, Tertiary lacustrine sediments
            % 4,  Tsl5, Tertiary lacustrine sediments
            % 5,  Tls3, Tertiary limestone
            % 6,  Tbo7, Tertiary basalt flows
            % 7,  Tsl4, Tertiary lacustrine sediments
            % 8,  Tls2, Tertiary limestone
            % 9,  Tbo6, Tertiary basalt flows
            % 10,  Tsl3, Tertiary lacustrine sediments
            % 11,  Tls1, Tertiary limestone
            % 12,  Tbo5, Tertiary basalt flows
            % 13,  Tsl2, Tertiary lacustrine sediments
            % 14,  Tbo4, Tertiary basalt flows
            % 15,  Tsl1, Tertiary lacustrine sediments
            % 16,  Tbo3, Tertiary basalt flows
            % 17,  Tpd, Tertiary Porphrytic (hornblende-biotite) dacite to rhyodacite flows and domes
            % 18,  Tbo2, Tertiary basalt flows
            % 19,  Tlr, Tertiary rhyolite lavas and lesser tuffs
            % 20,  Tbo1, Tertiary basalt flows
            % 21,  Tda, Tertiary andesite to dacite lavas
            % 22,  Tslo, Tertiary lacustrine sediments
            % 23,  Trt, Tertiary (Oligocene) ash-flow tuffs
            % 24,  Mzu, Mesozoic undifferentiated
            %
            % ?????????????
            %
            %
            % File: Lithology.pdat
            %
            % # Type: property scattered data
            % # Version: 80
            % # Description: Query points for Earth Vision (Tabrez Ali, 26 June 2015)
            % # Format: free
            % # Field: 1 x meters
            % # Field: 2 y meters
            % # Field: 3 z meters
            % # Field: 4 Lithology non-numeric
            % # Coordinate_System_Id: CRD589B1
            % # Coordinate_System_Name: NAD83 / UTM Zone 11 North
            % # Projection: Universal Transverse Mercator
            % # Zone: 11
            % # Units: meters
            % # Ellipsoid: GRS 1980/NAD83
            % # End:
            %  some examples: Is ?? air or no information, or both?
            % 329148.327	4408248.871	1287.5	""
            % 329163.17	4408268.987	1287.5	"Q"
            % 329178.013	4408289.104	1287.5	"Tsl2"
            % 329192.856	4408309.221	1287.5	"Tsl2"
            % 329207.699	4408329.338	1287.5	"Tsl2"
            % 329222.542	4408349.454	1287.5	"Q"
            % 329237.386	4408369.571	1287.5	"Tsl2"
            % 329252.229	4408389.688	1287.5	"Tsl2"
            % 329267.072	4408409.804	1287.5	"Tsl2"
            %
             
            % Read the file containing density values
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task4_Design_Deployment/Subtask4_9_Incorporate_Geologic_Information/SilerBradysGeoModelData'
            fname_siler_lithology = 'Lithology.pdat'; % 
            T=readtable(strcat(dirname,filesep,fname_siler_lithology),'HeaderLines',16,'FileType','text');
            LITHOLOGY = table2struct(T,'ToScalar',true)
            
            % rename fields
            LITHOLOGY.Eutm_in_meters      = LITHOLOGY.Var1;LITHOLOGY=rmfield(LITHOLOGY,'Var1');
            LITHOLOGY.Nutm_in_meters      = LITHOLOGY.Var2;LITHOLOGY=rmfield(LITHOLOGY,'Var2');
            LITHOLOGY.ElevNAD83_in_meters = LITHOLOGY.Var3;LITHOLOGY=rmfield(LITHOLOGY,'Var3');
            LITHOLOGY.Code                = LITHOLOGY.Var4;LITHOLOGY=rmfield(LITHOLOGY,'Var4');
            
            % Transform coordinates into rotated PoroTomo frame (Xp,Yp,Zp)
            [LITHOLOGY.Xp,LITHOLOGY.Yp,LITHOLOGY.Zp] = utm2xyz_porotomo(...
                 colvec(LITHOLOGY.Eutm_in_meters)...
                ,colvec(LITHOLOGY.Nutm_in_meters)...
                ,colvec(LITHOLOGY.ElevNAD83_in_meters));
            
            
            %% Convert lithologic units' text codes into numerical ID numbers
            %LITHCODES.Symbols_List = char(unique(LITHOLOGY.Code))
            T=readtable(strcat(dirname,filesep,'SilerLithologyKey.xlsx'));
            %/Volumes/KForange/Box Sync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler
           %LITHCODES.Lithology = strrep(LITHCODES.Lithology,'\t',' ');           
           %[LITHCODES.ncodes,nchar] = size(LITHCODES.Symbol);
%             
%             % add a code for Unknown
%             LITHCODES.ncodes = LITHCODES.ncodes+1;
%             LITHCODES.Symbol{end+1} = sprintf('----');
%             LITHCODES.Lithology{end+1}    = sprintf('unspecified');
%             LITHCODES.ID(LITHCODES.ncodes)    = LITHCODES.ncodes;
%             
            %% make a color table
            % for more values, see Color Palette City http://soliton.vm.bytemark.co.uk/pub/cpt-city/
            % use canned values from http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
%             T=readtable(strcat(dirname,filesep,'Paired_12.gpl'),'HeaderLines',6,'FileType','text');
%             C=table2struct(T,'ToScalar',true);
%             r1      = C.Var1/256.;C=rmfield(C,'Var1');
%             g1      = C.Var2/256.;C=rmfield(C,'Var2');
%             b1      = C.Var3/256.;C=rmfield(C,'Var3');
%                     
%             % add the values from the matlab color table named "vga"
%             LITHCODES.colormap =[r1, g1, b1;vga]; 
%             
%             % place gray at the top
%             LITHCODES.colormap(11,:) = [0, 1, 0]; % green
%             LITHCODES.colormap(LITHCODES.ncodes,:)  = [0.5, 0.5, 0.5];
%             % truncate unused values
%             LITHCODES.colormap=LITHCODES.colormap(1:LITHCODES.ncodes,:);
%             

% 
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler'
            T=readtable(strcat(dirname,filesep,'SilerCorrelations.xlsx'),'Sheet',1);
            LITHCODES=table2struct(T,'ToScalar',true)
            LITHCODES.Lithology = strrep(LITHCODES.Lithology,'\t',' ');
            LITHCODES.Symbols_List = char(unique(LITHOLOGY.Code))
           [LITHCODES.ncodes,nchar] = size(LITHCODES.Symbol);
            

            LITHCODES.Symbols_List = char(LITHCODES.Symbol);
            LITHCODES.ncodes = numel(LITHCODES.ID)
            LITHCODES.colormap = zeros(LITHCODES.ncodes,3);
            for i=1:LITHCODES.ncodes
                LITHCODES.colormap(i,1) = LITHCODES.R(i)/256;
                LITHCODES.colormap(i,2) = LITHCODES.G(i)/256;
                LITHCODES.colormap(i,3) = LITHCODES.B(i)/256;
            end
            
                

            % make a figure with a color key
            figure;hold on;
            colcol = [1:LITHCODES.ncodes]';
            imagesc([0,1],[1,LITHCODES.ncodes],colcol);
            
            colormap(LITHCODES.colormap);
            for i=1:LITHCODES.ncodes
                str1 = char(LITHCODES.Lithology{i});
                % truncate
                if numel(str1) > 35
                    str1 = str1(1:35); 
                end
                text(0.,i+0.5,sprintf('%2d %-4s %-35s'...
                    ,LITHCODES.ID(i)...
                    ,char(LITHCODES.Symbol{i})...
                    ,str1) ...
                    ,'BackgroundColor','white'...
                    ,'FontName','fixedwidth'...
                    ,'VerticalAlignment','top');
            end  
            axis xy;axis tight;
            %axis([-Inf,+Inf,0.5,LITHCODES.ncodes+0.5]);
            set(gca,'XTickLabel','');
            ylabel('Siler ID');
            colorbar;
            %colorbar('YTickLabel',LITHCODES.Lithology);
            feval(funprint,sprintf('lithologic_units'));

            LITHOLOGY.ID = zeros(size(LITHOLOGY.Xp));
            for i=1:numel(LITHOLOGY.Code)
                code1 = char(LITHOLOGY.Code{i});
                if numel(code1) > 0
                    for j=1:LITHCODES.ncodes
                        if strcmp(LITHCODES.Symbol(j,1:nchar),code1)==1
                            LITHOLOGY.ID(i)=j;
                            break;
                        end
                    end
                end
            end
            

            save('LITHOLOGY.mat','-struct','LITHOLOGY');
            save('LITHCODES.mat','-struct','LITHCODES');
        case 16
            %% strain rates from Reinisch et al.
            % On  2018-Jul-03, at 14:49 , Elena Reinisch <ebaluyut@wisc.edu> wrote: 
            % Hi Kurt, I have attached a
            % VoxelMap.mat file which contains a structure with information
            % for making the voxel plot as well as some code and necessary
            % .mat files used to make the plot itself.  Some of the fields
            % in VoxelMap (dV, dT, dE) are structures themselves containing
            % both values and unit information.  Please let me know if I
            % forgot to include anything.
            %
            % Thanks!
            %
            % Elena
            %
            % VoxelMap
            %
            % VoxelMap =
            %
            %                dV: [1x1 struct]
            %                dT: [1x1 struct]
            %                dE: [1x1 struct]
            %        strainrate: [1x1 struct]
            %              mest: [1656x1 double]
            %         centx_utm: [1656x1 double]
            %         centy_utm: [1656x1 double]
            %    centx_porotomo: [1656x1 double]
            %    centy_porotomo: [1656x1 double]
            %
            %
            %
            %
            % <VoxelMap.mat><cmappolargraynan.mat><SYMS.mat><plot_voxel_map.m>
             
            %% Read the file 
%             dirname = '/Users/feigl/BoxSync/PoroTomo/Task9_Inverse_Modeling/Subtask9_2_Inverse_Modeling_of_Geodetic_Data/'
%             fname_reinsich_voxels = 'VoxelMap.mat'; % 
%             load(strcat(dirname,filesep,fname_reinsich_voxels));
%             VoxelMap
            %       On  2018-Aug-13, at 09:57 , Elena Reinisch <ebaluyut@wisc.edu>
            %       wrote:
            %
            % Hi Kurt,
            %
            % I have attached a new version of the VoxelMap.mat file with uncertainties
            % for the strain rate estimates included (under the field name msig).
            %
            % Elena
            %
            % <VoxelMap.mat>

            dirname = '/Users/feigl/BoxSync/PoroTomo/Task9_Inverse_Modeling/Subtask9_2_Inverse_Modeling_of_Geodetic_Data/'
            fname_reinsich_voxels = 'VoxelMapWithUncertainties.mat'; % 
            load(strcat(dirname,filesep,fname_reinsich_voxels));
            VoxelMap

            
            
            % make figure in UTM coordinates
            nf=nf+1;figure;hold on;
            plot(VoxelMap.centx_utm,VoxelMap.centy_utm,'r.');
            plot(BRADYBOX.E([1,2,3,4,1]),BRADYBOX.N([1,2,3,4,1]),'k-');
            xlabel('UTM Easting [m]');
            ylabel('UTM Northing [m]');
            axis equal;axis tight
            title('UTM coordinates');
            legend('model nodes','PoroTomo Box');
            printpdf(sprintf('VoxelMap_UTM.pdf'));
           
            % make figure in PoroTomo coordinates
            nf=nf+1;figure;hold on;
            plot(VoxelMap.centx_porotomo,VoxelMap.centy_porotomo,'r.');
            plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
            plot(WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,'k^','MarkerSize',2,'MarkerFaceColor','k');
            xlabel('Xporotomo [m]');
            ylabel('Yporotomo [m]');
            axis equal; axis tight;
            title('PoroTomo coordinates');
            legend('model nodes','PoroTomo Box');
            printpdf(sprintf('VoxelMap_XpYp.pdf'));         
          
            
            %% load a colormap named cmappolargraynan
            load(strcat(dirname,'cmappolargraynan.mat'))  
            
                       %             TOMO.Xp = VoxelMap.centx_porotomo;
            %             TOMO.Yp = VoxelMap.centy_porotomo;
            %             % what is depth of centroid?
            %             TOMO.Zp = (elev_mean- 100 -800) * ones(size(VoxelMap.centy_porotomo));
            %             TOMO.V  = VoxelMap.dV.value;   % volume change rate [m^3/year?]
            %             TOMO.R  = nan(size(TOMO.V));
            npoints = numel(VoxelMap.centx_utm)
            dside = 100; % length of side of cube [m]
            DELVOL.E=nan(8*npoints,1);
            DELVOL.N=nan(8*npoints,1);
            DELVOL.Z=nan(8*npoints,1);
            DELVOL.V=nan(8*npoints,1);
            % get strain rate in picostrain per second
            epsdot = unit(VoxelMap.strainrate.value,'1/year');
            epsdot=epsdot.value*1e12;
            
            % prune 
            insignificant = find(abs(VoxelMap.mest./VoxelMap.msig) < 1);
            epsdot(insignificant) = NaN;

            k1 = 1;
            k2 = npoints;
            for i=1:4
                for j=1:2
                    % create vertices of cubes
                    switch i
                        case 1
                            de = -1*dside/2; dn = -1*dside/2;
                        case 2
                            de = +1*dside/2; dn = -1*dside/2;
                        case 3
                            de = +1*dside/2; dn = -1*dside/2;
                        case 4
                            de = +1*dside/2; dn = +1*dside/2;
                        otherwise
                            error('bad i');                            
                    end
                    switch j
                        case 1
                            dz = -1*dside/2;
                        case 2
                            dz = +1*dside/2;
                        otherwise
                            error('bad j');    
                    end
                    % UTM coordinates     
                    DELVOL.E(k1:k2,1) = VoxelMap.centx_utm + de;
                    DELVOL.N(k1:k2,1) = VoxelMap.centy_utm + dn;
                    DELVOL.Z(k1:k2,1) = (elev_mean - 100 + dz) * ones(size(VoxelMap.centy_utm));
                    % volume change rate [m^3/year]
                    %DELVOL.V(k1:k2) = VoxelMap.dV.value;   
                    % strain rate [1/year]
                    DELVOL.V(k1:k2) = epsdot;   
                    % update indices
                    k1 = k1 + npoints;
                    k2 = k2 + npoints;
                end
            end
            %add values of zero above and below
            DELVOL.E = repmat(DELVOL.E,3,1);
            DELVOL.N = repmat(DELVOL.N,3,1);
            DELVOL.Z = repmat(DELVOL.Z,3,1);
            DELVOL.V = repmat(DELVOL.V,3,1);
            k1 = 1;
            k2 = 8*npoints;
            k3 = k1 + 8*npoints;
            k4 = k2 + 8*npoints;
            k5 = k3 + 8*npoints;
            k6 = k4 + 8*npoints;

            DELVOL.Z(k1:k2)=DELVOL.Z(k3:k4) + dside + 1.0;
            DELVOL.V(k1:k2)=0;
            DELVOL.Z(k5:k6)=DELVOL.Z(k3:k4) - dside - 1.0;
            DELVOL.V(k5:k6)=0;

            [DELVOL.Xp,DELVOL.Yp, DELVOL.Zp] = utm2xyz_porotomo(DELVOL.E,DELVOL.N,DELVOL.Z);
            
            DELVOL.R  = nan(size(DELVOL.V)); 
            DELVOL
            save('DELVOL.mat','-struct','DELVOL');
        case 20
            % count faults
            dthresh = 50; % threshold distance to nearest fault [m]
%             FAULTCOUNTS = count_faults3(FAULTS,MESH3.Xp,MESH3.Yp,MESH3.Zp,dthresh);
%             
%             save('FAULTCOUNTS3.mat','-struct','FAULTCOUNTS');
            FAULTCOUNTS=load('FAULTCOUNTS.mat');
            iok = find(FAULTCOUNTS.count>0);           
            figure;hold on;
            imagesc(squeeze(FAULTCOUNTS.count(1,:,:)));
            axis equal
            colorbar
        case 21 
            % Modeled temperature from regression (training data set)
            % estimated using routine named plot_correlations.m
            Ttrain=load('Ttraining.mat');
        case 22
            % Modeled temperature from regression (testing data set)
            % estimated using routine named plot_correlations.m
            Ttestr=load('Ttesting.mat');
        case 24
            %% measure fault density
            % survey with 25 m box and ibox = 4
            %FAULTAREAS = survey_fault_area7(FAULTS,MESH3);
            %save('FAULTAREAS25mV7.mat','-struct','FAULTAREAS');
            %save('FAULTAREAS25mV7ibound4.mat','-struct','FAULTAREAS');
            %FAULTAREAS = load('FAULTAREAS25mV7.mat')
            FAULTAREAS = load('FAULTAREAS25mV7ibound4.mat')
            %             FAULTAREAS =
            %
            %             struct with fields:
            %
            %             count: [16368�1 double]
            %             area: [16368�1 double]
            %             density: [16368�1 double]
            %             X: [62�22�12 double]
            %             Y: [62�22�12 double]
            %             Z: [62�22�12 double]
%             FAULTAREAS.count  = zeros(nnodes,1); % number of faults intersecting voxel [dimless]
%             FAULTAREAS.area   = zeros(nnodes,1); % surface area of fault inside voxel [m^2]
%             FAULTAREAS.density= zeros(nnodes,1); % fault area divided by voxel volume [m^2/m^3]
            
            histogram(FAULTAREAS.density);
            xlabel('Fault density [m^2/m^3]');
            ylabel('Number of voxels');
            
        case 25
            % Hydraulic conductivity
            % On  2018-Nov-30, at 11:27 , Jeremy Patterson <jpatterson7@wisc.edu>
            % wrote:
            % I?ve uploaded the coordinate data for hydraulic conductivity estimates to Box.
            % X / Y coordinates are in the rotated PoroTomo coordinate
            % system, and the Z / Elevation coordinates are in absolute elevation
            % of the WGS84 datum.
            %
            % Respectfully, Jeremy
            %
            % Reservoir Hydraulic Conductivity Coordinates
            % https://uwmadison.app.box.com/file/359672784994
            
            dirname = '/Users/feigl/BoxSync/PoroTomo/Task9_Inverse_Modeling/Subtask9_3_Inverse_Modeling_of_Pressure_Data';
            fname25a = 'HydCondDmgZone.txt'; % 
            S25a=table2struct(readtable(strcat(dirname,filesep,fname25a)),'ToScalar',true);
            %             S25 =
            %
            %   struct with fields:
            %
            %             X: [2534700�1 double]
            %             Y: [2534700�1 double]
            %         layer: [2534700�1 double]
            %      sublayer: [2534700�1 double]
            %          node: [2534700�1 double]
            %     Elevation: [2534700�1 double]
            %            Kx: [2534700�1 double]
            %            Ky: [2534700�1 double]
            %            Kz: [2534700�1 double]
            %         Var10: [2534700�1 double]
               
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
            %break;
    end 

   
      

    
    %% build structure for TOMO
    % TOMO.Xp = colvec(RGRID.Xp);
    % TOMO.Yp = colvec(RGRID.Yp);
    % TOMO.Zp = colvec(RGRID.Zp);
    % TOMO.V  = colvec(RGRID.density_kg_per_cubic_meter);
    switch obsq
        case 0
            TOMO.Xp = DENSITY.Xp;
            TOMO.Yp = DENSITY.Yp;
            TOMO.Zp = DENSITY.Zp;
        case {1,2,3,4,7,8,9,11}
            TOMO.Xp = MATZEL.Xp;
            TOMO.Yp = MATZEL.Yp;
            %TOMO.Zp = MATZEL.Zp; % assumes ground is horizontal flat surface
            TOMO.Zp = MATZEL.Zpt; % assumes ground has topographic relief
        case {5,10,13,14,17,18}
            % Thurber Vp
            % 20180526 Cliff wrote:
            %     Going from 0 to surface (where does 490 come from?) does not make
            %     sense because we know resolution is poor except for the top 200 to
            %     250 m. I would only plot top 250 m max.
            %
            % To satisfy the reviewer, for your cross-section plots, please shift the Z
            % coordinate of the top layer of the model upwards by 25 m and leave it at
            % the same velocity, that?s basically what the inversion thinks is the
            % velocity above the top layer anyway.           
            itop = find(TOMO.Zp >= 450.);
            TOMO.Zp(itop) = TOMO.Zp(itop) + 25;
            TOMO.R  = nan(size(TOMO.V));
            fprintf(1,'%s\n','Thurber');
        case 6
            % Zeng MASW on DAS
            %         TOMO.Xp
            %         TOMO.Yp
            %         TOMO.Zp
            %         TOMO.V % velocity in meters per second
            fprintf(1,'%s\n','MASW');
        case 12
            TOMO.Xp = PT.Xp;
            TOMO.Yp = PT.Yp;
            TOMO.Zp = PT.Zp;
            TOMO.V  = PT.TemperatureInDegCelsius;
            TOMO.R  = nan(size(TOMO.V));
        case 15
            TOMO.Xp = LITHOLOGY.Xp;
            TOMO.Yp = LITHOLOGY.Yp;
            TOMO.Zp = LITHOLOGY.Zp;
            TOMO.V  = LITHOLOGY.ID;
            TOMO.R  = nan(size(TOMO.V));
        case 16
            TOMO.Xp = DELVOL.Xp;
            TOMO.Yp = DELVOL.Yp;
            TOMO.Zp = DELVOL.Zp;
%             % perturb coordinates slightly to avoid problems with coplanar
%             % locations
%             TOMO.Xp = DELVOL.Xp + 1*randn(size(DELVOL.Xp));
%             TOMO.Yp = DELVOL.Yp + 1*randn(size(DELVOL.Yp));
%             TOMO.Zp = DELVOL.Zp + 1*randn(size(DELVOL.Zp));
            TOMO.V  = DELVOL.V;
            TOMO.R  = DELVOL.R;
        case 20
            % fault counts
            TOMO.Xp = colvec(FAULTCOUNTS.X);
            TOMO.Yp = colvec(FAULTCOUNTS.Y);
            TOMO.Zp = colvec(FAULTCOUNTS.Z);
            TOMO.V  = colvec(FAULTCOUNTS.count); 
        case 21
            TOMO.Xp = Ttrain.Xp;
            TOMO.Yp = Ttrain.Yp;
            TOMO.Zp = Ttrain.Zp;
            TOMO.V  = Ttrain.Tmod; % temperature in degC
        case 22
            TOMO.Xp = Ttestr.Xp;
            TOMO.Yp = Ttestr.Yp;
            TOMO.Zp = Ttestr.Zp;
            TOMO.V  = Ttestr.Tmod; % temperature in degC
        case 23
            TOMO.Xp = PT.Xp;
            TOMO.Yp = PT.Yp;
            TOMO.Zp = PT.Zp;
            TOMO.V  = PT.PressureInPascal;
            % convert zero values to NaNs
            TOMO.V(abs(TOMO.V) < eps) = nan;
            %TOMO.V=log10(TOMO.V); % take logarithm
            %TOMO.V = TOMO.V - 9.8e3*(elev_mean - TOMO.Zp - 800.0); % subtract rho g
            TOMO.V = TOMO.V/1.e6; % convert from Pa to MPa
            TOMO.R  = nan(size(TOMO.V));
        case 24
            %% fault density
            TOMO.Xp = colvec(FAULTAREAS.X);
            TOMO.Yp = colvec(FAULTAREAS.Y);
            TOMO.Zp = colvec(FAULTAREAS.Z);
            TOMO.V  = colvec(FAULTAREAS.density);
            %figure
            %histogram(TOMO.V);
        case 25
            %% Hydraulic Conductivity
            TOMO.Xp = colvec(S25a.X);
            TOMO.Yp = colvec(S25a.Y);
            TOMO.Zp = colvec(S25a.Elevation)-800; % convert elevation to PoroTomo Z [m]
            TOMO.V  = colvec(log10((S25a.Kx).^2 + (S25a.Ky).^2 + (S25a.Kz).^2)); % [m/s]
            %figure
            %histogram(TOMO.V);
         otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    
    %% where to make slices
    switch obsq
%         case {0} % Witter Density
% %             SLICES.Xp = [-400:100:900];
% %             SLICES.Yp = [-300:100:1700];
% %             SLICES.Zp = [0 100 200 250 300 350 390 400 410 420 430 440 450];
%             % Try "flattened cube" at 185 m depth at location of Well 56-1
%             SLICES.Xp = 24;
%             SLICES.Yp = 119;
%             SLICES.Zp = 1277-195-800;
%         case {2,3,4,7,8,9,11}  % Matzel's area is smaller
%             SLICES.Xp = [0 50 100 200 270 300 400 500];
%             SLICES.Yp = [0 300 600 900 1000 1200 1500];
%             SLICES.Zp = [25   125   225   325   425];
%         case {1,5,12} % Cliff's models go wider
%             SLICES.Xp = [-400:100:900];
%             SLICES.Yp = [-300:100:1700];
%             SLICES.Zp = [0 100 200 250 300 350 390 400 410 420 430 440 450];
        case {10,14,15,17,18} % P-wave tomography for SRL paper by Parker et al.
            SLICES.Xp = [0 250 500];
            SLICES.Yp = [0:300:1500];
            SLICES.Zp = [1200 1150 1100 1050]-800; 
        case 6 % for surface waves, do not go too deep
            SLICES.Xp = [200 400];
            SLICES.Yp = [ 300 600 900 1200];
            SLICES.Zp = [390 400 410 420 430 440 450];
%         case 14 % Debug faults
%             SLICES.Xp = [300];
%             SLICES.Yp = [600];
%             SLICES.Zp = [400];
        case 16 % only one map view is useful
            SLICES.Xp = [0 250 500];
            SLICES.Yp = [0:300:1500];
            SLICES.Zp = [350];
        otherwise
            %error(sprintf('unknown obsq %d\n',obsq));
%             SLICES.Xp = [-400 0 500 900];
%             SLICES.Yp = [1500 1200 900 600 300 0];
%             SLICES.Zp = [1300 1200 1150 1100 1050]-800; 
            SLICES.Xp = [0 250 500];
            SLICES.Yp = [0:300:1500];
            SLICES.Zp = [1200 1150 1100 1050]-800; 
    end

    


    
    % set up options
    OPTIONS.nmin = 20; % minimum number of intersecting points to qualify a fault for plot
    OPTIONS.norder = 2;  % order of polynomial fit (norder = 1 is linear)
    %OPTIONS.plot_points = 1; % plot intersections as black dots
    OPTIONS.plot_points = 0; % do NOT plot intersections as black dots
    OPTIONS.plot_curves = 1; % plot polynomial as black lines
    OPTIONS.draw_topo = 1; % draw topography as thin line and clip
    OPTIONS.draw_mudpots   = 1; % draw mudpots
    OPTIONS.draw_fumaroles = 1; % draw mudpots
    
    %% set up observed values and titles
    switch obsq
        case 0
            title_str = 'Density_from_Witter_et_al';
            TOMO.V = 1000*DENSITY.Density_gcm3;
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('parula'));
            OPTIONS.colorbarlabelstr = 'Density [kg/m^3]';
        case 1
            title_str = 'Vp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Vp;
            TOMO.R  = nan(size(TOMO.V));           
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
        case 2
            title_str = 'Vs_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Vs;
            TOMO.R  = nan(size(TOMO.V));           
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
        case 3
            title_str = 'Poisson_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Poisson;
            TOMO.R  = nan(size(TOMO.V));            
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = '\nu [dimless]';
        case 4
            title_str = 'EYoung_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.EYoung/1.e9;% convert from Pa to GPa
            TOMO.R  = nan(size(TOMO.V));           
            OPTIONS.cmap=colormap('jet');
            OPTIONS.colorbarlabelstr = 'E [GPa]';
        case 5
            title_str = 'Vp_Thurber_velfile20170727';
            % TOMO.V is defined above
            TOMO.R  = nan(size(TOMO.V));           
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.short_labels = 1;
        case 6
            title_str = 'Vs_MASW_on_DAS'
            % TOMO.V is defined above
            TOMO.R  = nan(size(TOMO.V));                      
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
        case 7
            title_str = 'Qp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Qp;%
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('summer'));
            OPTIONS.colorbarlabelstr = 'Qp [dimless]';
        case 8
            title_str = 'Qs_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Qs;           
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('summer'));
            OPTIONS.colorbarlabelstr = 'Qs [dimless]';
        case 9          
            %             title_str = 'QpOverQs_MatzelSweepInterfNov2017';
            %             TOMO.V = MATZEL.QpOverQs;%
            %             OPTIONS.cmap=flipud(colormap('winter'));
            %             OPTIONS.colorbarlabelstr = 'Qp/Qs [dimless]';
            % better to plot reciprocal
            warning(sprintf('skipping case 9\n'));
        case 10
            title_str = 'Vp_Thurber20171123';
            % TOMO.V is defined above
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 0;
            OPTIONS.short_labels = 1;
            OPTIONS.draw_topo = 1; % not fully implemented yet
            OPTIONS.contours = 1;
        case 11
            title_str = 'QsOverQp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.QsOverQp;%
            TOMO.R  = nan(size(TOMO.V));
            %OPTIONS.cmap=colormap('winter');
            OPTIONS.cmap=colormap('jet'); % 2080213
            OPTIONS.colorbarlabelstr = 'Qs/Qp [dimless]';
        case 12
            title_str = 'Temperature20150324';
            TOMO.V = PT.TemperatureInDegCelsius;
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('hot'));
            OPTIONS.colorbarlabelstr = 'T [degC]';
            OPTIONS.contours = 0;
        case 13
            title_str = sprintf('Vp_Thurber20180504');
            % TOMO.V is defined above
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 0; % not fully implemented yet
            OPTIONS.draw_box  = 1;
            %OPTIONS.short_labels = 0;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
        case 14
            title_str = sprintf('Vp_Nayak20180513');
            % TOMO.V is defined above
            % TOMO.R  is defined above
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            %OPTIONS.short_labels = 0;
            OPTIONS.draw_points = 1;
            OPTIONS.draw_curves = 1;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
        case 15
            title_str = sprintf('SilerLithology');
            % TOMO.V is defined above
            TOMO.R  = nan(size(TOMO.V));
            %OPTIONS.cmap=colormap(LITHCODES.colormap);
            OPTIONS.cmap=LITHCODES.colormap;
            OPTIONS.colorbarlabelstr = '[ID]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            %OPTIONS.short_labels = 0;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
        case 16
            %title_str = sprintf('Rate of Volume Change');OPTIONS.colorbarlabelstr = VoxelMap.dV.unit;
            title_str = sprintf('StrainRate');
            OPTIONS.colorbarlabelstr = '[picostrain/second]';
            % TOMO.V is defined above
            OPTIONS.cmap=cmappolargraynan;          
            OPTIONS.draw_wells = 0;
            OPTIONS.draw_topo = 0;
            OPTIONS.draw_box  = 1;
            %OPTIONS.short_labels = 0;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
      case 17
            title_str = sprintf('Vp_Nayak20180926');
            % TOMO.V is defined above
            % TOMO.R  is defined above
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 0;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            OPTIONS.draw_points = 0;
            OPTIONS.draw_curves = 1;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 0;
      case 18
            title_str = sprintf('Vp_Nayak20180607');
            % TOMO.V is defined above
            % TOMO.R  is defined above
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 0;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            OPTIONS.draw_points = 0;
            OPTIONS.draw_curves = 1;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 0;
      case 20
            title_str = sprintf('FAULTCOUNTS');
            % TOMO.V is defined above
            % TOMO.R  is defined above
            %OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.cmap=colormap(flipud(hot));
            OPTIONS.colorbarlabelstr = 'count';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            OPTIONS.draw_points = 0;
            OPTIONS.draw_curves = 1;
            OPTIONS.short_labels = 0; % 
            OPTIONS.contours = 0;
        case 21
            title_str = 'TmodTraining';
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('hot'));
            OPTIONS.colorbarlabelstr = 'T [degC]';  
            OPTIONS.contours = 0;
        case 22
            title_str = 'TmodTesting';
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('hot'));
            OPTIONS.colorbarlabelstr = 'T [degC]';
            OPTIONS.contours = 0;
        case 23
            title_str = 'Pressure';
            TOMO.R  = nan(size(TOMO.V));
            OPTIONS.cmap=colormap('cool');
            OPTIONS.colorbarlabelstr = 'P [MPa]';
            %OPTIONS.colorbarlabelstr = 'log10(P/(1 Pa))';
            OPTIONS.contours = 0;
        case 24
            title_str = sprintf('FaultAreaPerUnitVolume');
            %OPTIONS.cmap=flipud(colormap('jet'));
            %OPTIONS.cmap=colormap(flipud(hot));
            %OPTIONS.cmap=cmap(10:10:50,1:3); % take only 5 colors
            cmap=flipud(hot);
            OPTIONS.cmap=cmap(10:50,1:3);
            OPTIONS.colorbarlabelstr = 'A/V [1/m]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            OPTIONS.draw_points = 0;
            OPTIONS.draw_curves = 1;
            OPTIONS.short_labels = 0; %
            OPTIONS.contours = 0;
        case 25
            title_str = sprintf('HydraulicConductivity');
            OPTIONS.cmap=flipud(colormap('cool'));
            OPTIONS.colorbarlabelstr = 'log_{10}(K) [m/s]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            OPTIONS.draw_points = 0;
            OPTIONS.draw_curves = 1;
            OPTIONS.draw_traces = 1;
            OPTIONS.short_labels = 0; %
            OPTIONS.contours = 0;
         otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    %% standardize for comparison
    if standardslices == 1
        SLICES.Xp = [0 250 500];
        SLICES.Yp = [0:300:1500];
        SLICES.Zp = [1250 1200 1150 1100 1050]-800;
        OPTIONS.draw_wells = 1;
        OPTIONS.draw_topo = 1;
        OPTIONS.draw_box  = 1;
        OPTIONS.draw_points = 0;
        OPTIONS.draw_curves = 1;
        OPTIONS.contours = 0;
        % how close does a well have to be show up in plot?
        %OPTIONS.rsearch = 2*sqrt(dx^2 + dy^2 + dz^2);  
        OPTIONS.rsearch = 100;  
        OPTIONS.short_labels = 0; % show details
        %OPTIONS.short_labels = 1; % no details
        %title_str = strcat(title_str,'_stdslices');
    end

    
    %% prune
    fprintf(1,'Before pruning\n');
%     TOMO
%     figure;histogram(TOMO.Xp); title('TOMO.Xp');
%     figure;histogram(TOMO.Yp); title('TOMO.Yp');
%     figure;histogram(TOMO.Zp); title('TOMO.Zp');
    %dx = 50; % tolerance in meters
    npoints = numel(TOMO.Xp)
    iok=1:npoints;
    iok=find(isfinite(TOMO.V));
    iok=intersect(iok,find(TOMO.Xp <= nanmax(BOUNDS.Xp)+dx));
    iok=intersect(iok,find(TOMO.Xp >= nanmin(BOUNDS.Xp)-dx));
    iok=intersect(iok,find(TOMO.Yp <= nanmax(BOUNDS.Yp)+dy));
    iok=intersect(iok,find(TOMO.Yp >= nanmin(BOUNDS.Yp)-dy));
    iok=intersect(iok,find(TOMO.Zp <= nanmax(BOUNDS.Zp)+dz));
    iok=intersect(iok,find(TOMO.Zp >= nanmin(BOUNDS.Zp)-dz));
    TOMO.Xp = TOMO.Xp(iok);
    TOMO.Yp = TOMO.Yp(iok);
    TOMO.Zp = TOMO.Zp(iok);
    TOMO.V = TOMO.V(iok);
%     if isfield(TOMO,'R') 
%         TOMO.R = TOMO.R(iok);
%     end
    fprintf(1,'After pruning\n');
    TOMO
    npoints = numel(TOMO.Xp)
    if npoints < 10
        error(sprintf('npoints (%d) less than 10\n',npoints));
    end
    
    %% set up plots
    switch obsq
        case 0
            OPTIONS.vmin=nanmin(TOMO.V); %
            OPTIONS.vmax=nanmax(TOMO.V); %
            OPTIONS.interpolation_method = 'nearest';
        case {2,3,4,7,8,9,11}
            OPTIONS.vmin=quantile(TOMO.V,0.05); % scale to include 90 percent of values
            OPTIONS.vmax=quantile(TOMO.V,0.95); %
            OPTIONS.interpolation_method = 'linear';
        case {1,5,10,13,14,17,18}
            OPTIONS.vmin=1000; % m/s % one scale for all plots of Vp
            OPTIONS.vmax=3000; % m/s
            OPTIONS.interpolation_method = 'linear';
        case 6
            OPTIONS.vmin= 200; % Dante's choice
            OPTIONS.vmax= 800; %
%               OPTIONS.vmin = nan; % scale each slice individually and automatically
%               OPTIONS.vmax = nan;
            OPTIONS.interpolation_method = 'linear';
        case {12,21,22} % Temperature 
            OPTIONS.vmin=quantile(TOMO.V,0.05); % scale to include 90 percent of values
            OPTIONS.vmax=quantile(TOMO.V,0.95); %
            %OPTIONS.interpolation_method = 'linear';
            OPTIONS.interpolation_method = 'natural';
            %OPTIONS.interpolation_method = 'nearest';
       case 15 % Lithology
            OPTIONS.vmin=0; % 
            OPTIONS.vmax=LITHCODES.ncodes;
            OPTIONS.interpolation_method = 'nearest';
       case 16 % Rate of Volume change
%             OPTIONS.vmin=quantile(TOMO.V,0.05); % scale to include 90 percent of values
%             OPTIONS.vmax=quantile(TOMO.V,0.95); %
            vmed = nanmedian(TOMO.V);
            vabs = nanmax([abs(nanmin(TOMO.V)-vmed) ...
                          ,abs(nanmax(TOMO.V)-vmed)]);
            OPTIONS.vmin=vmed-vabs; % symmetric scale 
            OPTIONS.vmax=vmed+vabs; %
            %OPTIONS.interpolation_method = 'linear';
            OPTIONS.interpolation_method = 'natural';
            %OPTIONS.interpolation_method = 'nearest';
       case 20 % Fault counts
             OPTIONS.vmin=0; % 
            OPTIONS.vmax=nanmax(TOMO.V)+1;
            OPTIONS.interpolation_method = 'nearest';
        case 23 % Pressure
            % OPTIONS.vmin=nanmin(TOMO.V); %
            % OPTIONS.vmax=nanmax(TOMO.V); %
            OPTIONS.vmin = nan; % scale each slice individually and automatically
            OPTIONS.vmax = nan;
            OPTIONS.interpolation_method = 'nearest';
        case 24 % Fault density
%             Try stretching histogram
%             TOMO.V = histogram_stretch(TOMO.V,5);
%             OPTIONS.vmin=0; %
%             OPTIONS.vmax=nanmax(TOMO.V);
%             OPTIONS.interpolation_method = 'nearest'
%             scale to include 80 percent of values
            OPTIONS.vmin=quantile(TOMO.V,0.10);
            OPTIONS.vmax=quantile(TOMO.V,0.90); %
%             scale to include 68 percent of values
%             cdf('normal',1,0,1) - cdf('normal',-1,0,1)
%             ans =   0.682689492137086
%             OPTIONS.vmin=quantile(TOMO.V,cdf('normal',-1,0,1)); % 1 sigma below mean
%             OPTIONS.vmax=quantile(TOMO.V,cdf('normal',+1,0,1)); % 1 sigma above mean
            %OPTIONS.interpolation_method = 'natural';
            OPTIONS.interpolation_method = 'nearest';
            %OPTIONS.interpolation_method = 'linear';
        case 25 % Hydraulic Conductivity
             OPTIONS.vmin=0; %
             OPTIONS.vmax=nanmax(TOMO.V);
%             scale to include 80 percent of values
            OPTIONS.vmin=quantile(TOMO.V,0.10);
            OPTIONS.vmax=quantile(TOMO.V,0.90); %
            %OPTIONS.interpolation_method = 'natural';
            OPTIONS.interpolation_method = 'nearest';
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    
    %% choose the search radius
%     switch obsq
%         case {0,1,2,3,4,7,8,9,10,11,13,14,15,16,17,18,20,21,22,23}
%             rsearch = 50; % search radius in meters
%         case 12
%             rsearch = 500;
%         case 24
%             rsearch = 200;
%         case 5
%             rsearch = 30; % search radius in meters
%         case 6
%             rsearch = 5; % search radius in meters
%             %rsearch = 500; % search radius in meters
%         otherwise
%             error(sprintf('unknown obsq %d\n',obsq));
%     end
    
    %% choose the Extrapolation method
    switch obsq
        case {0}
            OPTIONS.extrapolation_method = 'nearest';
        case {1,2,3,4,5,6,7,8,9,10,11,13,14,16,17,18,20,21,22,23,24,25}
            OPTIONS.extrapolation_method = 'none';
        case {12,15}
            OPTIONS.extrapolation_method = 'none';
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    %% choose the padding [m] to avoid edge effects when extrapolating
    switch obsq
        case 12
            OPTIONS.padding = 1000;
        otherwise
            OPTIONS.padding = 50;
    end

   %% subtract mean value for depth
   if show_anomalies == 1
       title_str = strcat(title_str,'_depth_anomaly');
       OPTIONS.colorbarlabelstr = strrep(OPTIONS.colorbarlabelstr,'[','anomaly [');
%        zvalues_list = unique(TOMO.Zp);
%        for i=1:numel(zvalues_list)
%            iz = find(abs(TOMO.Zp - zvalues_list(i)) <= eps);
%            TOMO.V(iz) = TOMO.V(iz) - nanmean(TOMO.V(iz));
%        end
       zvalues_list = [nanmin(BOUNDS.Zp):2*dZ:nanmax(BOUNDS.Zp)];
       for i=1:numel(zvalues_list)
           iz = intersect(find(TOMO.Zp >= zvalues_list(i)-dZ)...
                         ,find(TOMO.Zp <  zvalues_list(i)+dZ));
           TOMO.V(iz) = TOMO.V(iz) - nanmean(TOMO.V(iz));
       end
       
%        % scale color table to [min, max]
%        OPTIONS.vmin=nanmin(TOMO.V); %
%        OPTIONS.vmax=nanmax(TOMO.V); %
       
       % scale color table symmetrically
       vavg = nanmean(TOMO.V);
       vabs = nanmax([abs(nanmin(TOMO.V)-vavg) ...
                     ,abs(nanmax(TOMO.V)-vavg)]);
       OPTIONS.vmin=vavg-vabs; 
       OPTIONS.vmax=vavg+vabs; 
   end
   
    
    %% make a histogram
    figure;
    histogram(colvec(TOMO.V));
    xlabel(OPTIONS.colorbarlabelstr);
    ylabel('number of occurences');
    title(title_str,'Interpreter','none');
    feval(funprint,sprintf('%s_histogram',title_str));
    
    %% save profiles in boreholes
    if save_wells == 1
        j = find(obsqs == obsq);
        WELLS.tomonames{j} = title_str;
        WELLS.tomoids(j) = obsqs(j);
        WELLS.tomoZp     = zprof; %  Zcoordinate
        WELLS.units{j} = OPTIONS.colorbarlabelstr; % units
        
        Finterp3 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,OPTIONS.interpolation_method,OPTIONS.extrapolation_method);
        figure; hold on;
        
        for i=1:nwells
            xprof=repmat(WELLS.Xp(i),size(zprof));
            yprof=repmat(WELLS.Yp(i),size(zprof));
            welltomo = Finterp3(xprof,yprof,zprof);
            plot(welltomo,zprof,'LineWidth',2);
            WELLS.vertprofvals(i,j,:) = colvec(welltomo);
        end
        
        
        ylabel('Zporotomo [m]');
        xlabel(OPTIONS.colorbarlabelstr);
        legend(WELLS.WellName,'Location','Southwest');
        title(strrep(title_str,'_',' '));
        feval(funprint,sprintf('%s_vprofile',title_str));
    end
    
    %% interpolate on to regular mesh
    if save_meshes == 1
        if obsq == obsqs(1)
            MESHEDTOMO.Xp = MESH3.Xp;
            MESHEDTOMO.Yp = MESH3.Yp;
            MESHEDTOMO.Zp = MESH3.Zp;
            MESHEDTOMO.E =  MESH3.E;
            MESHEDTOMO.N =  MESH3.N;
            MESHEDTOMO.H =  MESH3.H;
            MESHEDTOMO.Lat =MESH3.Lat;
            MESHEDTOMO.Lon =MESH3.Lon;
        end

        ii = find(obsqs == obsq)
%         MESHEDTOMO.tomonames{ii} = title_str;
%         MESHEDTOMO.tomoids(ii) = obsqs(ii);
%         MESHEDTOMO.units{ii} = OPTIONS.colorbarlabelstr; % units
        Finterp3 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,OPTIONS.interpolation_method,OPTIONS.extrapolation_method);
        %MESHEDTOMO.tomo(ii,:,:,:) = Finterp3(MESH3.Xp,MESH3.Yp,MESH3.Zp);
        tomoid = sprintf('ID%02d_%s',obsq,title_str);
        %tomoid = sprintf('tomoid%02d',obsq);
        MESHEDTOMO.(tomoid) = Finterp3(MESH3.Xp,MESH3.Yp,MESH3.Zp);
    end
    
    
    % add ibound to file name
    if standardslices == 1
       FileNamePrefix = strcat(sprintf('ibound%1d_stdslices',ibound),filesep)
    else
       FileNamePrefix = strcat(sprintf('ibound%1d',ibound),filesep)
    end

    %% make slices
    if make_slices == 1 
        %% make the plots with slices
        % try "flattened cube"
        OPTIONS.flattened_cube = 0;
        %OPTIONS.contours = 0;
        OPTIONS.geologic_model = geologic_model;
        OPTIONS.draw_box = 1;
        %title_str = strcat(title_str,'_with_',geologic_model,'_faults')
        title_str = sprintf('M%02d_%s_%s_faults',obsq,title_str,geologic_model)

        figfilenames = plot_tomo_and_faults9(TOMO,FAULTS,title_str,BOUNDS,OPTIONS,funprint,WELLS,elev_mean,BRADY2DGRID,SLICES,BRADYBOX...
            ,FUMAROLES,MUDPOTS,FileNamePrefix,TRACES);
        %        [figfilenames,FAULTCOUNTS] = plot_tomo_and_faults8(TOMO,FAULTS,title_str,BOUNDS,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID,SLICES,BRADYBOX...
        %             ,FUMAROLES,MUDPOTS);
    end
    
    %% do not make too many windows
    if numel(obsqs) > 1
       close all;
    end
end % loop over obsq
if save_wells == 1
    save('WELLS2.mat','-struct','WELLS');
end
 
%% 
if save_meshes == 1
    MESHEDTOMO
    clear T1;
    whos T1
    field_names=fieldnames(MESHEDTOMO)
    kount = 0;  
    for i=1:numel(field_names)
        field_name1 = field_names{i}
        [nnx,nny,nnz] = size(MESHEDTOMO.(field_name1));
        %if (numel(field_names{i}) == numel('tomoid00')) && (contains(field_names{i},'tomoid') == 1)
        if nnx == MESH3.nx && nny == MESH3.ny && nnz == MESH3.nz 
            kount = kount+1;
            if kount == 1
%                 T1 = table;
%                 S1 = struct([]);
%                 T1.Description = title_str;
%                 T1.height = nnx*nny*nnz
                index = colvec(1:nnx*nny*nnz);
                S1 = struct('index',index);                 
            end
            tomovec = reshape(MESHEDTOMO.(field_name1),nnx*nny*nnz,1);
            S1.(field_name1) = tomovec;
            
            % write the name of the variable and its units
            switch field_name1
                case {'E','N','H','Xp','Yp','Zp'}
                    unitstr = 'm';
                case 'Lat'
                    unitstr = 'degN';
                case 'Lon'
                    unitstr = 'degE';
                otherwise
                    unitstr = OPTIONS.colorbarlabelstr;
                    unitstr = strrep(unitstr,'/','_per_');
                    unitstr = strrep(unitstr,'[','_');
                    unitstr = strrep(unitstr,']','_');
            end
            % remove special characters, spaces, and leading numbers
%             T1.Properties.VariableNames{i} = sprintf('%s_in_%s' ...
%                 ,T1.Properties.VariableNames{i} ...
%                 ,matlab.lang.makeValidName(unitstr,'ReplacementStyle','delete'));
            field_names{i} = sprintf('%s_in_%s' ...
                ,field_name1 ...
                ,matlab.lang.makeValidName(unitstr,'ReplacementStyle','delete'));
        end
    end
    
    if show_anomalies == 1      
       SaveFileNamePrefix = 'MESHEDANOM'; 
    else
       SaveFileNamePrefix = 'MESHEDTOMO'; 
    end
    
    save(sprintf('%s_%s.mat',SaveFileNamePrefix,datestr(now,'yyyymmdd')),'-struct','MESHEDTOMO');
    T1 = struct2table(S1);
    writetable(T1,sprintf('%s_%s.xlsx',SaveFileNamePrefix,datestr(now,'yyyymmdd')));   
    writetable(T1,sprintf('%s_%s.csv', SaveFileNamePrefix,datestr(now,'yyyymmdd')));   
end

        





