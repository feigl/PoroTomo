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
% example how to make an animation
% convert -delay 1 Vp_Thurber20171123*.png VpThurber20171123.gif
% 20180601 Kurt and Avinash adapt for github



%% initialize
clear all;
close all;
nf=0;
tstart = tic;

save_wells  = 0;   % save Z profile
save_meshes = 0;  % save tomograms on REGULAR 3 D mesh
make_slices = 0;  % plot cross sections



%% set up path for Matlab
%addpath(genpath('/globus/PoroTomo/SOFTWARE'),'-begin');
%addpath('/Users/feigl/gipht/utils','-begin'); % needed for colortables
% addpath(genpath('/Users/feigl/PoroTomo'),'-begin');
% rmpath(genpath('/Users/feigl/PoroTomo/.git'));
homedir = getenv('HOME')
addpath(genpath(strcat(homedir,filesep,'Porotomo')),'-begin');
rmpath(genpath(strcat(homedir,filesep,'Porotomo/.git')));

myos = computer
if strcmp(myos, 'GLNXA64')==1
    fprintf(1,'Linux computer\n');
end

[status,hostname] = system('hostname')
if contains(hostname,'porotomo') == 1
    boxdir = '/mnt/u1/feigl/BoxSync'
else
    boxdir = strcat(homedir,filesep,'BoxSync')
end


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
demgrdfilename = strcat(boxdir,'/PoroTomo/METADATA/brady_dem_srtm_20m.grd')

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

%% Limits of cross sections
% 20180224 - set limits separately from slices
% % Elevation above WGS84 goes from 800 m to 1250 m
% BOUNDS.Zp = [0, 490]; % These are Zporotomo coordinates.
% % Elevation above WGS84 goes from 800 m to 1290 m
% Stay within area of best resolution
% BOUNDS.Xp = [0,500];
% BOUNDS.Yp = [0,1500];
% BOUNDS.Zp = [0, 490]; % These are Zporotomo coordinates.
%BOUNDS.Zp = [200, 490]; % These are Zporotomo coordinates.
BOUNDS.Xp = [-400,900];
BOUNDS.Yp = [-300,1700];
BOUNDS.Zp = [0, 490]; % These are Zporotomo coordinates.


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
dX = 25; % meters
dY = 25; % meters
dZ = 25; % meters

[MESH3.Xp,MESH3.Yp,MESH3.Zp] = meshgrid(...
    [min(BOUNDS.Xp):dX:max(BOUNDS.Xp)]...
    ,[min(BOUNDS.Xp):dY:max(BOUNDS.Xp)]...
    ,[min(BOUNDS.Zp):dZ:max(BOUNDS.Zp)]);
[MESH3.E, MESH3.N, MESH3.H] = xyz_porotomo2utm(colvec(MESH3.Xp),colvec(MESH3.Yp),colvec(MESH3.Zp));
MESH3.E = reshape(MESH3.E,size(MESH3.Xp));
MESH3.N = reshape(MESH3.N,size(MESH3.Xp));
MESH3.H = reshape(MESH3.H,size(MESH3.Xp));
[MESH3.nx,MESH3.ny,MESH3.nz] = size(MESH3.Xp)



%% read the file containing coordinates of nodes
%dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_4_Adjoint_seismic_tomography/MESH_SPECFEM3D'
dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_4_Adjoint_seismic_tomography/MESH_SPECFEM3D')
fname_nodes = strcat(dirname,filesep,'MESH_topo_xpypzellipsoid_800/nodes_coords_file');
NODES = read_nodes(fname_nodes);

%% read the file containing the list of elements
dirname =strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_4_Adjoint_seismic_tomography/MESH_SPECFEM3D')
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
    load(strcat(boxdir,filsep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie9.mat')); % All faults
elseif strcmp(geologic_model,'Jolie')
    load(strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie.mat')); % All faults
elseif strcmp(geologic_model,'Siler')
    load(strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Siler.mat')); % All in Siler's model
else
    error(sprintf('Unknown geologic_model: %s\n',geologic_model));
end

FAULTS = S
clear S;



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
T=readtable(strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Surface_features/brady_fumerole_ll.txt'));
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
T=readtable(strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Surface_features/brady_mudpot_ll.txt'));
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
feval(funprint,sprintf('%s_WGS84elevation.pdf',mfilename));



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
% obsq = 16 % Nayak20180607 Vp3 with model resolution
% obsq = 17 % Thurber20180504 empirical Vs

obsqmin = 0;
obsqmax = 17;

% dimension
%obsqs = 10; % for SRL paper
%obsqs = [0:16] % All
obsqs = [0,17]; % for debugging

if ismember(9,obsqs) == true
    obsqs=setxor(obsqs,9); % case 9 is not defined
end

% number of observable quantities
nobsq = numel(obsqs)
% dimension array for storing vertical profiles at wells
vertprofvals=nan(nwells,nobsq,numel(zprof));
% dimension array for storing 3-D tomograms, interpolated to regular grid
% named MESH3
MESHEDTOMO.tomo = nan(nobsq,MESH3.nx,MESH3.ny,MESH3.nz);



%%loop over observable quantity
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
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/WitterDensity')
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
                dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_3_Seismic_ANT')
                fname_matzel_sweep_interf = 'Porotomo_Sweep_Interferometry_Model.version_Nov2017.XYcoordinates.txt'; % omits missing data values
                T=readtable(strcat(dirname,filesep,fname_matzel_sweep_interf),'ReadVariableNames',true,'HeaderLines',0);
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
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff')
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
            
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time')
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
        case {13,17}
            %% obsq = 13 Thurber20180504 P-wave velocity with model resolution from Cliff Thurber
            
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
            
            
            %% obsq = 17 Thurber20180504 empirical Vs
            %              From: Clifford Thurber
            % Sent: Friday, June 8, 2018 3:51 PM To: Kurt Feigl Subject: Empirical Vs
            % model
            %
            % Hi Kurt.  Here is a Vs model that can be given to Christina along with
            % the velfile.porotomo Vp model, interpolated onto her grid, to use as a
            % starting model.  I hope it is reasonably close to reality such that her
            % inversion can converge.
            %
            % FYI, all I did was divide Vp by 2.75 for Vp < 1.8 km/s and divide by 2.0
            % for Vp >= 1.8 km/s.
            %
            % Arbitrary but reasonable.
            %
            % Cliff
            
            
            
            %dirname = '/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time'
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Thurber20180504PoroTomo_Z')
            
            % read P-wave velocity
            tomo_filename = 'velfile.porotomo'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            % read model resolution
            mres_filename = 'resfile.porotomo'
            XYZRthurber = load(strcat(dirname,filesep,mres_filename));
            
            % read empirical S-wave velocity
            vsem_filename = 'velfileS.porotomo'
            XYZSthurber = load(strcat(dirname,filesep,vsem_filename));           
            
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
            TOMO.Vp = XYZVthurber(iok,4)*1000; % P-wave velocity in meters per second
            TOMO.Rm = XYZRthurber(iok,4)*1000; % model resolution (dimensionless)
            TOMO.Vs = XYZSthurber(iok,4)*1000; % S-wave velocity in meters per second
            
            
            
            %% get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
            
            %% interpolate onto centroids of Morency mesh
            TOMOFEMMESH = CENTROIDS;
            FVp = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.Vp);TOMOFEMMESH.Vp = FVp(CENTROIDS.Xp,CENTROIDS.Yp,CENTROIDS.Zp);
            FVs = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.Vs);TOMOFEMMESH.Vs = FVs(CENTROIDS.Xp,CENTROIDS.Yp,CENTROIDS.Zp);
            FRm = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.Rm);TOMOFEMMESH.Rm = FRm(CENTROIDS.Xp,CENTROIDS.Yp,CENTROIDS.Zp);
            % also need density
            FDkg = scatteredInterpolant(DENSITY.Xp,DENSITY.Yp,DENSITY.Zp,1000*DENSITY.Density_gcm3);
            TOMOFEMMESH.Density_kg_per_m3 = FDkg(CENTROIDS.Xp,CENTROIDS.Yp,CENTROIDS.Zp);
                      

            %% Poisson's ratio from Vp/Vs ratio
            %  https://en.wikipedia.org/wiki/Elastic_modulus
            TOMOFEMMESH.VpOverVs = TOMOFEMMESH.Vp./TOMOFEMMESH.Vs;
            TOMOFEMMESH.Poisson = 0.5*((TOMOFEMMESH.VpOverVs.^2) - 2)./((TOMOFEMMESH.VpOverVs.^2)-1);
            
            %%  Young's modulus from M
            % https://en.wikipedia.org/wiki/P-wave_modulus
            % In linear elasticity, the P-wave modulus M {\displaystyle M} M, also
            % known as the longitudinal modulus or the constrained modulus, is one of
            % the elastic moduli available to describe isotropic homogeneous materials.
            % It is defined as the ratio of axial stress to axial strain in a uniaxial
            % strain state
            %     ? z z = M ? z z {\displaystyle \sigma _{zz}=M\epsilon _{zz}}
            %     {\displaystyle \sigma _{zz}=M\epsilon _{zz}}
            % where all the other strains ? ? ? {\displaystyle \epsilon _{**}}
            % {\displaystyle \epsilon _{**}} are zero.
            % This is equivalent to stating that
            %     M = ? V P 2 {\displaystyle M=\rho V_{\mathrm {P} }^{2}}
            %     {\displaystyle M=\rho V_{\mathrm {P} }^{2}}
            % where VP is the velocity of a P-wave.
            Mp = TOMOFEMMESH.Density_kg_per_m3 .* (TOMOFEMMESH.Vp).^2;           
            % Young's modulus
            TOMOFEMMESH.EYoung  = Mp .* (1 + TOMOFEMMESH.Poisson) .* (1 - 2*TOMOFEMMESH.Poisson) ./ (1 - TOMOFEMMESH.Poisson);
            
            
            %% make table to save as ASCII 
            T= struct2table(TOMOFEMMESH);
            %id      Xp         Yp         Zp          E            N          Lat       Lon      ElevWGS84      Vp         Vs       Rm      Density_kg_per_m3    VpOverVs    Poisson       EYoung    

            T.Properties.Description = sprintf('Material properties for PoroTomo on centroids of Morency FEM Mesh');
            for i=1:numel(T.Properties.VariableNames)
                vname = char(T.Properties.VariableNames{i});
                if contains(vname,{'ID'}) == 1
                    T.Properties.VariableUnits{i} = 'dimless';
                    T.Properties.VariableDescriptions{i} = 'Element ID number';
                elseif contains(vname,'EYoung') == 1
                    T.Properties.VariableUnits{i} = 'Pa';
                    T.Properties.VariableDescriptions{i} = 'Young\''s modulus';
                elseif contains(vname,{'Xp','Yp','Zp'}) == 1
                    T.Properties.VariableUnits{i} = 'm';
                    T.Properties.VariableDescriptions{i} = 'Rotated PoroTomo coordinates';
                elseif contains(vname,{'E','N','Elev'}) == 1
                    T.Properties.VariableUnits{i} = 'm';
                    T.Properties.VariableDescriptions{i} = 'coordinates in UTM Zone 11';
                elseif contains(vname,{'Lat','Lon'}) == 1
                    T.Properties.VariableUnits{i} = 'deg';
                elseif contains(vname,'VpOverVs') == 1
                    T.Properties.VariableUnits{i} = 'dimless';
                    T.Properties.VariableDescriptions{i} = 'empircal';
                elseif contains(vname,'Vp') == 1
                    T.Properties.VariableUnits{i} = 'm/s';
                    T.Properties.VariableDescriptions{i} = 'VpThurber20180504';
                elseif contains(vname,'Vs') == 1
                    T.Properties.VariableUnits{i} = 'm/s';
                    T.Properties.VariableDescriptions{i} = 'EmpiricalVs_from_Thurber20180504';
                else
                    T.Properties.VariableUnits{i} = 'dimless';
                end
            end
            T.Properties.Description
            T.Properties.VariableNames
            T.Properties.VariableUnits
            T.Properties.VariableDescriptions
            
            %% save as text file
            %save('VpThurber20171123.mat','-struct','TOMO');           
            writetable(T,strcat(dirname,filesep,'VpThurber20180504onMorencyMesh.csv'));
            
            %% save file as Matlab binary
            %save('VpThurber20180504onMorencyMesh.mat','-struct','TOMOFEMMESH');
            save(strcat(dirname,filesep,'VpThurber20180504onMorencyMesh.mat'),'T');
            
            %% select quantity for subsequent plotting
            switch obsq
                case 13
                    TOMO.V = TOMO.Vp;
                    % set threshold to mask velocity values where resolution is poor
                    OPTIONS.resolution_threshold = 0.1;  % threshhold value
                case 17
                    TOMO.V = TOMO.Vs;
                    OPTIONS.resolution_threshold = NaN;
                otherwise
                    error(sprintf('Unknown obsq %d\n',obsq));
            end
            
            
        case 14
            %              %% P-wave velocity with from Avinash Nayak % From: Avinash
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
            
            
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Nayak20180513')
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
            TOMO.Rm  = XYZRthurber(iok,4)*1000; % model resolution (dimensionless)
            
            %            % clip out velocity values where resolution is poor
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %             OPTIONS.resolution_threshold = 0.05;  % threshhold value
            
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
            
            save('VpNayak20180513.mat','-struct','TOMO');
        case 16
            %% P-wave velocity with from Avinash Nayak
            %  On  2018-Jun-07, at 23:24 , Avinash Nayak <avinash07guddu@gmail.com>
            %  wrote:
            %
            % Sorry for the delay in replying to your email. Unfortunately, I am
            % traveling and I will be back in Madison on Saturday night. I was able to
            % run interpolate_fields5.m with the latest Vp model. The .fig files
            % generated by the code are in the tar file (dropbox link
            % https://www.dropbox.com/s/umggaso225gef5p/nayak_vp_model_20180605.tar?dl=0
            % ). I have also attached the model and resolution data and 4 depth slices
            % that show Vp over the entire extent of the model. Thanks. Avinash
            
            
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_2_Seismic_travel_time/Nayak20180607')
            % read velocity
            tomo_filename = 'vp_new3.txt'
            XYZVthurber = load(strcat(dirname,filesep,tomo_filename));
            
            %  Here is the resolution file. SIMUL code at this point doesn't output
            % resolution of the last depth slice. It is a bug. The actual elevation is
            % the Z value + 800 m in these files. Thanks. Avinash
            [nvel,nfour] = size(XYZVthurber)
            
            % read model resolution
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
            TOMO.Xp = XYZVthurber(iok,1)*1000; % meters short  axis positive toward SE
            TOMO.Yp = XYZVthurber(iok,2)*1000; % meters short axis positive toward  NE
            TOMO.Zp = XYZVthurber(iok,3)*1000; % meters above 800 m
            TOMO.V  = XYZVthurber(iok,4)*1000; % velocity in meters per second
            TOMO.Rm  = XYZRthurber(iok,4)*1000; % model resolution (dimensionless)
            
            %            % clip out velocity values where resolution is poor
            OPTIONS.resolution_threshold = 0.1;  % threshhold value
            %             OPTIONS.resolution_threshold = 0.05;  % threshhold value
            
            % get UTM coordinates
            [TOMO.Eutm,TOMO.Nutm,TOMO.Dutm] = xyz_porotomo2utm(TOMO.Xp,TOMO.Yp,TOMO.Zp);
            
            save('VpNayak20180607.mat','-struct','TOMO');
        case 6
            %% ZengMASWonDAS
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_1_DAS_Quality_Control/SurfaceWave')
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
            demgrdfilename = strcat(boxdir,filesep,'PoroTomo/METADATA/brady_dem_srtm_20m.grd');
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
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_1_DAS_Quality_Control/SurfaceWave2/newmodel');
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
            
            % find good segments
            isegx=find(seg_flag>0);
            
            % plot in 3 dimensions
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
            caxis([200,600]);
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
            save('ZengMASWonDAS','-struct','TOMO');
            
        case 12
            %% Read the file containing Temperature Values
            url = 'https://gdr.openei.org/files/939'
            data_filename='PTfield_brady_for_GDR_PT_20150324.csv'
            websave(data_filename,strcat(url,filesep,data_filename));
            table=readtable(data_filename);
            TEMPERATURE=table2struct(table,'ToScalar',true);
            
            % rename the fields in the stucture
            TEMPERATURE.Xp = TEMPERATURE.XPoroTomoInMeters;TEMPERATURE=rmfield(TEMPERATURE,'XPoroTomoInMeters');
            TEMPERATURE.Yp = TEMPERATURE.YPoroTomoInMeters;TEMPERATURE=rmfield(TEMPERATURE,'YPoroTomoInMeters');
            TEMPERATURE.Zp = TEMPERATURE.ZPoroTomoInMeters;TEMPERATURE=rmfield(TEMPERATURE,'ZPoroTomoInMeters');
            TEMPERATURE.E  = TEMPERATURE.UTMeastingInMeters; TEMPERATURE=rmfield(TEMPERATURE,'UTMeastingInMeters');
            TEMPERATURE.N  = TEMPERATURE.UTMnorthingInMeters;TEMPERATURE=rmfield(TEMPERATURE,'UTMnorthingInMeters');
            %
            % make figure in UTM coordinates
            nf=nf+1;figure(nf);hold on;
            plot(TEMPERATURE.E,TEMPERATURE.N,'r.');
            plot(BRADYBOX.E([1,2,3,4,1]),BRADYBOX.N([1,2,3,4,1]),'k-');
            xlabel('UTM Easting [m]');
            ylabel('UTM Northing [m]');
            axis equal;axis tight
            title('UTM coordinates');
            legend('model nodes','PoroTomo Box');
            printpdf(sprintf('TEMPERATURE_UTM.pdf'));
            
            % make figure in PoroTomo coordinates
            nf=nf+1;figure(nf);hold on;
            plot(TEMPERATURE.Xp,TEMPERATURE.Yp,'r.');
            plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
            plot(WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,'k^','MarkerSize',2,'MarkerFaceColor','k');
            xlabel('Xporotomo [m]');
            ylabel('Yporotomo [m]');
            axis equal; axis tight;
            title('PoroTomo coordinates');
            legend('model nodes','PoroTomo Box');
            printpdf(sprintf('TEMPERATURE_XpYp.pdf'));
            save('TEMPERATURE.mat','-struct','TEMPERATURE');
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
            dirname = strcat(boxdir,filesep,'PoroTomo/Task4_Design_Deployment/Subtask4_9_Incorporate_Geologic_Information/SilerBradysGeoModelData')
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
            dirname = strcat(boxdir,filesep,'PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler')
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
        case {5,10,13,14,16,17}
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
            % 20180607 Model resolution is not defined for early models
            if ismember(obsq,[5,10]) == 1
                TOMO.Rm  = nan(size(TOMO.V));
            end
            fprintf(1,'%s\n','Thurber');
        case 6
            % Zeng MASW on DAS
            %         TOMO.Xp
            %         TOMO.Yp
            %         TOMO.Zp
            %         TOMO.V % velocity in meters per second
            fprintf(1,'%s\n','MASW');
        case 12
            TOMO.Xp = TEMPERATURE.Xp;
            TOMO.Yp = TEMPERATURE.Yp;
            TOMO.Zp = TEMPERATURE.Zp;
            TOMO.V  = TEMPERATURE.TemperatureInDegCelsius;
            TOMO.Rm  = nan(size(TOMO.V));
        case 15
            TOMO.Xp = LITHOLOGY.Xp;
            TOMO.Yp = LITHOLOGY.Yp;
            TOMO.Zp = LITHOLOGY.Zp;
            TOMO.V  = LITHOLOGY.ID;
            TOMO.Rm  = nan(size(TOMO.V));
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
        case {10,14,16,17} % P-wave tomography for SRL paper by Parker et al.
            SLICES.Xp = [0 250 500];
            SLICES.Yp = [0:300:1500];
            SLICES.Zp = [1200 1150 1100 1050]-800;
        case 6 % for surface waves, do not go too deep
            SLICES.Xp = [0 100 200 300 400];
            SLICES.Yp = [0 200 400 600 800 1000 1200 1400 1500];
            SLICES.Zp = [390 400 410 420 430 440 450];
            %         case 14 % Debug faults
            %             SLICES.Xp = [300];
            %             SLICES.Yp = [600];
            %             SLICES.Zp = [400];
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
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('parula'));
            OPTIONS.colorbarlabelstr = 'Density [kg/m^3]';
        case 1
            title_str = 'Vp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Vp;
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
        case 2
            title_str = 'Vs_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Vs;
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
        case 3
            title_str = 'Poisson_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Poisson;
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = '\nu [dimless]';
        case 4
            title_str = 'EYoung_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.EYoung/1.e9;% convert from Pa to GPa
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=colormap('jet');
            OPTIONS.colorbarlabelstr = 'E [GPa]';
        case 5
            title_str = 'Vp_Thurber_velfile20170727';
            % TOMO.V is defined above
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.short_labels = 1;
        case 6
            title_str = 'Vs_MASW_on_DAS'
            % TOMO.V is defined above
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
        case 7
            title_str = 'Qp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Qp;%
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('summer'));
            OPTIONS.colorbarlabelstr = 'Qp [dimless]';
        case 8
            title_str = 'Qs_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Qs;
            TOMO.Rm  = nan(size(TOMO.V));
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
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 0;
            OPTIONS.short_labels = 1;
            OPTIONS.draw_topo = 1; % not fully implemented yet
            OPTIONS.contours = 1;
        case 11
            title_str = 'QsOverQp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.QsOverQp;%
            TOMO.Rm  = nan(size(TOMO.V));
            %OPTIONS.cmap=colormap('winter');
            OPTIONS.cmap=colormap('jet'); % 2080213
            OPTIONS.colorbarlabelstr = 'Qs/Qp [dimless]';
        case 12
            title_str = 'Temperature20150324';
            TOMO.V = TEMPERATURE.TemperatureInDegCelsius;
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=flipud(colormap('hot'));
            OPTIONS.colorbarlabelstr = 'Temperature [degC]';
        case 13
            title_str = sprintf('Vp_Thurber20180504');
            % TOMO.V is defined above
            % TOMO.Rm is defined above
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
            % TOMO.Rm  is defined above
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
            TOMO.Rm  = nan(size(TOMO.V));
            OPTIONS.cmap=colormap(LITHCODES.colormap);
            OPTIONS.colorbarlabelstr = 'Siler Lithologic ID';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 1;
            OPTIONS.draw_box  = 1;
            %OPTIONS.short_labels = 0;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
        case 16
            title_str = sprintf('Vp_Nayak20180607');
            % TOMO.V is defined above
            % TOMO.Rm  is defined above
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
        case 17
            title_str = sprintf('Vs_Thurber20180504');
            % TOMO.V is defined above
            % TOMO.Rm is defined above
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 0; % not fully implemented yet
            OPTIONS.draw_box  = 1;
            %OPTIONS.short_labels = 0;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    %% prune
    dx = 50; % tolerance in meters
    npoints = numel(TOMO.Xp);
    iok=1:npoints;
    iok=intersect(iok,find(TOMO.Xp <= nanmax(BOUNDS.Xp)+dx));
    iok=intersect(iok,find(TOMO.Xp >= nanmin(BOUNDS.Xp)-dx));
    iok=intersect(iok,find(TOMO.Yp <= nanmax(BOUNDS.Yp)+dx));
    iok=intersect(iok,find(TOMO.Yp >= nanmin(BOUNDS.Yp)-dx));
    iok=intersect(iok,find(TOMO.Zp <= nanmax(BOUNDS.Zp)+dx));
    iok=intersect(iok,find(TOMO.Zp >= nanmin(BOUNDS.Zp)-dx));
    TOMO.Xp = TOMO.Xp(iok);
    TOMO.Yp = TOMO.Yp(iok);
    TOMO.Zp = TOMO.Zp(iok);
    TOMO.V = TOMO.V(iok);
    if isfield(TOMO,'R')
        TOMO.Rm = TOMO.Rm(iok);
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
        case {1,5,10,13,14,16}
            OPTIONS.vmin=1000; % m/s % one scale for all plots of Vp
            OPTIONS.vmax=3000; % m/s
            OPTIONS.interpolation_method = 'linear';
        case 17
            OPTIONS.vmin= 200; % m/s % one scale for all plots of Vs
            OPTIONS.vmax=1800; % m/s
            OPTIONS.interpolation_method = 'linear';
        case 6
            OPTIONS.vmin= 200; % Dante's choice
            OPTIONS.vmax= 800; %
            OPTIONS.interpolation_method = 'linear';
        case 12 % Temperature
            OPTIONS.vmin=quantile(TOMO.V,0.05); % scale to include 90 percent of values
            OPTIONS.vmax=quantile(TOMO.V,0.95); %
            %OPTIONS.interpolation_method = 'linear';
            OPTIONS.interpolation_method = 'natural';
            %OPTIONS.interpolation_method = 'nearest';
        case 15 % Lithology
            OPTIONS.vmin=0; %
            OPTIONS.vmax=LITHCODES.ncodes;
            OPTIONS.interpolation_method = 'nearest';
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    
    %% choose the search radius
    switch obsq
        case {0,1,2,3,4,7,8,9,10,11,13,14,15,16,17}
            rsearch = 50; % search radius in meters
        case 12
            rsearch = 500;
        case 5
            rsearch = 30; % search radius in meters
        case 6
            rsearch = 5; % search radius in meters
            %rsearch = 500; % search radius in meters
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    %% choose the Extrapolation method
    switch obsq
        case 0
            OPTIONS.extrapolation_method = 'nearest';
        case {1,2,3,4,5,6,7,8,9,10,11,13,14,16,17}
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
    
    
    
    %% make a histogram
    figure;
    histogram(colvec(TOMO.V));
    xlabel(OPTIONS.colorbarlabelstr);
    ylabel('number of occurences');
    title(title_str,'Interpreter','none');
    feval(funprint,sprintf('%s_histogram.pdf',title_str));
    
    %% Extract the value of the field in a vertical profile at the location the wells
    if save_wells == 1
        j = find(obsqs == obsq)
        WELLS.tomonames{j} = title_str;
        WELLS.tomoids(j) = obsqs(j);
        WELLS.tomoZp     = zprof; %  Zcoordinate
        %WELLS.colorbarlabelstr{j} = OPTIONS.colorbarlabelstr; % units
        WELLS.units{j} = OPTIONS.colorbarlabelstr; % units
        
        %Finterp3 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,'nearest','nearest');
        %Finterp3 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,'linear','none');
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
        ii = find(obsqs == obsq)
        % MESHEDTOMO.tomo = nan(nobsqn,MESH3.nx,MESH3.ny,MESH3.nz);
        MESHEDTOMO.tomonames{ii} = title_str;
        MESHEDTOMO.tomoids(ii) = obsqs(ii);
        MESHEDTOMO.units{ii} = OPTIONS.colorbarlabelstr; % units
        Finterp3 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,OPTIONS.interpolation_method,OPTIONS.extrapolation_method);
        MESHEDTOMO.tomo(ii,:,:,:) = Finterp3(MESH3.Xp,MESH3.Yp,MESH3.Zp);
    end
    
    
    
    
    
    
    
    
    
    %% make the plots with slices
    %    figfilenames = plot_tomo_and_faults3(TOMO,FAULTS,title_str,SLICES,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID);
    % 20180211 for SRL paper only
    %figfilenames = plot_tomo_and_faults4(TOMO,FAULTS,title_str,SLICES,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID);
    
    % try "flattened cube"
    OPTIONS.flattened_cube = 0;
    %OPTIONS.contours = 0;
    OPTIONS.geologic_model = geologic_model;
    OPTIONS.draw_box = 1;
    title_str = strcat(title_str,'_with_',geologic_model,'_faults')
    
    if make_slices == 1
        OPTIONS.short_labels = 0; % show details
        %OPTIONS.short_labels = 1; % no details
        figfilenames = plot_tomo_and_faults7(TOMO,FAULTS,title_str,BOUNDS,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID,SLICES,BRADYBOX...
            ,FUMAROLES,MUDPOTS);
    end
    
    %close all;
end % loop over obsq
if save_wells == 1
    save('WELLS2.mat','-struct','WELLS');
end

if save_meshes == 1
    MESHEDTOMO.Xp = MESH3.Xp;
    MESHEDTOMO.Yp = MESH3.Yp;
    MESHEDTOMO.Zp = MESH3.Zp;
    save('MESHEDTOMO.mat','-struct','MESHEDTOMO');
end







