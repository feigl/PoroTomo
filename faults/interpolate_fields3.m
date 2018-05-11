%% interpolate fields for PoroTomo
% Build a configuration for PoroTomo project.
% 20171007 Kurt Feigl
% 20171026 with Dante - use constant color scale
% 20171120 with Lesley - add Matzel's latest
% 20171205 Kurt for AGU
% 20180207 Kurt for Stanford
% 20180211 Kurt for GRL 
% 20180507 Kurt for SRL and SSA
% example how to make an animation
% convert -delay 1 Vp_Thurber20171123*.png VpThurber20171123.gif



%% initialize
clear all;
close all;
nf=0;
tstart = tic;


%% set up path for Matlab
%addpath(genpath('/globus/PoroTomo/SOFTWARE'),'-begin');
%addpath('/Users/feigl/gipht/utils','-begin'); % needed for colortables
addpath(genpath('/Users/feigl/PoroTomo'),'-begin');

%% before releasing, find dependencies using:
[fList,pList] = matlab.codetools.requiredFilesAndProducts(mfilename);
fList
pList.Name


%% choose the function for printing
%funprint = str2func('printpng'); % images as bitmaps
%funprint = str2func('printpdf'); % PDF, but makes white lines
%funprint = str2func('printeps'); % EPS, but rendered ugly by Mac Preview
funprint = str2func('printfig'); % save as "compact" Matlab figure format
%funprint = str2func('printjpg'); % JPG

%% Digital Elevation model in Latitude and Longitude
% TODO: should use webget to ftp://roftp.ssec.wisc.edu/porotomo/PoroTomo/METADATA/
demgrdfilename = '/Users/feigl/BoxSync/PoroTomo/METADATA/brady_dem_srtm_20m.grd'

%% get the coordinates of the study area
utmzone = '11 S' % UTM zone for Brady Hot Springs
BRADYBOX = read_bradybox(1)
[BRADYBOX.Xp,BRADYBOX.Yp] = utm2xy_porotomo(BRADYBOX.E,BRADYBOX.N);
[BRADYBOX.Lat,BRADYBOX.Lon] = utm2deg(colvec(BRADYBOX.E),colvec(BRADYBOX.N),repmat(utmzone,size(colvec(BRADYBOX.E))));

% make 2-D grid with 25 m spacing
[BRADY2DGRID.Xp,BRADY2DGRID.Yp] = meshgrid([min(BRADYBOX.Xp):25:max(BRADYBOX.Xp)]...
    ,[min(BRADYBOX.Yp):25:max(BRADYBOX.Yp)])
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

% make figure in PoroTomo coordinates
nf=nf+1;figure(nf);hold on;
imagesc([min(BRADYBOX.Xp):max(BRADYBOX.Xp)],[min(BRADYBOX.Yp):max(BRADYBOX.Yp)],reshape(BRADY2DGRID.ElevWGS84,size(BRADY2DGRID.Xp)));
axis xy
axis image
%imagesc(reshape(BRADY2DGRID.ElevWGS84,size(BRADY2DGRID.Xp)));

%plot(BRADYBOX.Xporotomo([1,2,3,4,1]),BRADYBOX.Yporotomo([1,2,3,4,1]),'k-');
%plot(BRADYBOX.Xp([1,2,3,4,1]),BRADYBOX.Yp([1,2,3,4,1]),'k-');
colormap('summer');
% OPTIONS.cmap=colormap('autumn');
% OPTIONS.colorbarlabelstr = 'ElevWGS84 [m]';
colorbar;

xlabel('Xporotomo [m]');
ylabel('Yporotomo [m]');
axis equal; axis tight;
title('Topographic Elevation (meters above WGS84 ellipsoid');
printpdf(sprintf('%s_WGS84elevation.pdf',mfilename));
save('BRADYBOX.mat','-struct','BRADYBOX');
save('BRADY2DGRID.mat','-struct','BRADY2DGRID');


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
geologic_model = 'Jolie9'
%geologic_model = 'Jolie'
%geologic_model = 'Siler'
if strcmp(geologic_model,'Jolie9')
    load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/faults_for_Cliff/read_and_mesh_faults10.mat'); % Nick's Top9
elseif strcmp(geologic_model,'Jolie')
    load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Jolie.mat'); % All faults
elseif strcmp(geologic_model,'Siler')
    load('/Users/feigl/BoxSync/PoroTomo/Task8_Analyze_Data_Collected/Subtask8_9_Geology/Siler/Siler.mat'); % All in Siler's model
else
    error(sprintf('Unknown geologic_model: %s\n',geologic_model));
end

FAULTS = S
clear S;




%% get the coordinates for the wells from a file
WELLS = read_wells('brady');


% get the coordinates of the field in Well 56-1
i561=find(contains(WELLS.WellName,'56-1'));
zp561=colvec(0:5:450); % Vertical coordinate every 5 meters
xp561=repmat(WELLS.XCoordinate_UTMME_(i561),size(zp561));
yp561=repmat(WELLS.YCoordinate_UTMMN_(i561),size(zp561));



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
% obsq = 9 % Quality factor ratio Qp/Qs from Matzel's sweep interferometry
% obsq = 10 % P-wave velocity from body-wave tomography Thurber20171123
% obsq = 11 % Quality factor ratio Qs/Qp from Matzel's sweep interferometry
% obsq = 12 % Temperature
% obsq = 13 % Thurber20180504 Vp with model resolution

%% loop over observable quantity
for obsq = 13
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
        case {1,2,3,4,7,8,9,11}
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
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));          
    end 
    
    
    %% Limits of cross sections
    % 20180224 - set limits separately from slices
    BOUNDS.Xp = [-400,900];
    BOUNDS.Yp = [-300,1700];
    BOUNDS.Zp = [0, 450]; % These are Zporotomo coordinates.
    % Elevation above WGS84 goes from 800 m to 1250 m
    
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
%         case {10} % P-wave tomography for SRL paper by Parker et al.
%             SLICES.Xp = [-400 900];
%             SLICES.Yp = [1500 1200 900 600 300 0];
%             SLICES.Zp = [1200 1150 1100 1050]-800; 
        case 6 % for surface waves, do not go too deep
            SLICES.Xp = [0 100 200 300 400];
            SLICES.Yp = [0 200 400 600 800 1000 1200 1400 1500];
            SLICES.Zp = [390 400 410 420 430 440 450];
        otherwise
            %error(sprintf('unknown obsq %d\n',obsq));
            SLICES.Xp = [-400 0 500 900];
            SLICES.Yp = [1500 1200 900 600 300 0];
            SLICES.Zp = [1300 1200 1150 1100 1050]-800; 

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
        case {5,10,13}
            % Thurber Vp
            %         TOMO.Xp
            %         TOMO.Yp
            %         TOMO.Zp
            %         TOMO.V % velocity in meters per second
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
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    
    % set up options
    OPTIONS.nmin = 20; % minimum number of intersecting points to qualify a fault for plot
    OPTIONS.norder = 2;  % order of polynomial fit (norder = 1 is linear)
    %OPTIONS.plot_points = 1; % plot intersections as black dots
    OPTIONS.plot_points = 0; % do NOT plot intersections as black dots
    OPTIONS.plot_curves = 1; % plot polynomial as black lines
    OPTIONS.draw_topo = 1; % draw topography as thin line and clip
    
    %% set up titles
    switch obsq
        case 0
            title_str = 'Density_from_Witter_et_al';
            TOMO.V = 1000*DENSITY.Density_gcm3;
            OPTIONS.cmap=flipud(colormap('parula'));
            OPTIONS.colorbarlabelstr = 'Density [kg/m^3]';
        case 1
            title_str = 'Vp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Vp;
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
        case 2
            title_str = 'Vs_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Vs;
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
        case 3
            title_str = 'Poisson_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Poisson;
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = '\nu [dimless]';
        case 4
            title_str = 'EYoung_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.EYoung/1.e9;% convert from Pa to GPa
            OPTIONS.cmap=colormap('jet');
            OPTIONS.colorbarlabelstr = 'E [GPa]';
        case 5
            title_str = 'Vp_Thurber_velfile20170727';
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.short_labels = 1;
        case 6
            title_str = 'Vs_MASW_on_DAS'
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vs [m/s]';
        case 7
            title_str = 'Qp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Qp;%
            OPTIONS.cmap=flipud(colormap('summer'));
            OPTIONS.colorbarlabelstr = 'Qp [dimless]';
        case 8
            title_str = 'Qs_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.Qs;%
            OPTIONS.cmap=flipud(colormap('summer'));
            OPTIONS.colorbarlabelstr = 'Qs [dimless]';
%         case 9
%             title_str = 'QpOverQs_MatzelSweepInterfNov2017';
%             TOMO.V = MATZEL.QpOverQs;%
%             OPTIONS.cmap=flipud(colormap('winter'));
%             OPTIONS.colorbarlabelstr = 'Qp/Qs [dimless]';
        case 10
            title_str = 'Vp_Thurber20171123';
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]'; 
            OPTIONS.draw_wells = 0;
            OPTIONS.short_labels = 1;
            OPTIONS.contours = 1;
        case 11
            title_str = 'QsOverQp_MatzelSweepInterfNov2017';
            TOMO.V = MATZEL.QsOverQp;%
            %OPTIONS.cmap=colormap('winter');          
            OPTIONS.cmap=colormap('jet'); % 2080213
            OPTIONS.colorbarlabelstr = 'Qs/Qp [dimless]';
        case 12
            title_str = 'Temperature20150324';
            TOMO.V = TEMPERATURE.TemperatureInDegCelsius;
            OPTIONS.cmap=flipud(colormap('hot'));
            OPTIONS.colorbarlabelstr = 'Temperature [degC]';
        case 13
            title_str = sprintf('Vp_Thurber20180504');
            OPTIONS.cmap=flipud(colormap('jet'));
            OPTIONS.colorbarlabelstr = 'Vp [m/s]';
            OPTIONS.draw_wells = 1;
            OPTIONS.draw_topo = 0; % not fully implemented yet
            OPTIONS.draw_box  = 1;
           %OPTIONS.short_labels = 0;
            OPTIONS.short_labels = 1; % no details
            OPTIONS.contours = 1;
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    
    %% set up plots
    switch obsq
        case 0
            OPTIONS.vmin=nanmin(TOMO.V); %
            OPTIONS.vmax=nanmax(TOMO.V); %
            %         OPTIONS.vmin= nan; % autoscale each plot individually
            %         OPTIONS.vmax= nan; %
            OPTIONS.interpolation_method = 'nearest';
       case {2,3,4,7,8,9,11}
            OPTIONS.vmin=quantile(TOMO.V,0.05); % scale to include 90 percent of values
            OPTIONS.vmax=quantile(TOMO.V,0.95); %
            OPTIONS.interpolation_method = 'linear';
        case {1,5,10,13}
            OPTIONS.vmin=1000; % m/s % one scale for all plots of Vp
            OPTIONS.vmax=3000; % m/s
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
        otherwise
            error(sprintf('unknown obsq %d\n',obsq));
    end
    
    
    %% choose the search radius
    switch obsq
        case {0,1,2,3,4,7,8,9,10,11,13}
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
        case {1,2,3,4,5,6,7,8,9,10,11,13}
            OPTIONS.extrapolation_method = 'none';
        case 12
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
    
    %% Extract the value of the field in a vertical profile at the location of Well 56-1
    Finterp3 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,'nearest','nearest');
    VPROFILE.Xp=xp561;
    VPROFILE.Yp=yp561;
    VPROFILE.Zp=zp561;
    profval = Finterp3(VPROFILE.Xp,VPROFILE.Yp,VPROFILE.Zp);
    VPROFILE.(title_str) = profval;
    figure
    plot(VPROFILE.(title_str),VPROFILE.Zp,'r-','LineWidth',5);
    ylabel('Zporotomo [m]]');
    xlabel(sprintf('%s',OPTIONS.colorbarlabelstr));
    title(strrep(title_str,'_',' '));
    printpdf(sprintf('%s_vprofile561',title_str));
    save(sprintf('%s_vprofile561.mat',title_str),'-struct','TOMO');
    
    
    
    %% make the plots with slices
%    figfilenames = plot_tomo_and_faults3(TOMO,FAULTS,title_str,SLICES,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID);    
    % 20180211 for SRL paper only
    %figfilenames = plot_tomo_and_faults4(TOMO,FAULTS,title_str,SLICES,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID);
    
    % try "flattened cube"
     OPTIONS.flattened_cube = 0;
     OPTIONS.contours = 0;
     OPTIONS.geologic_model = geologic_model;
     OPTIONS.draw_box = 1;
     %title_str = strcat(title_str,'_with_',geologic_model,'_faults');
     figfilenames = plot_tomo_and_faults6(TOMO,FAULTS,title_str,BOUNDS,OPTIONS,funprint,rsearch,WELLS,elev_mean,BRADY2DGRID,SLICES,BRADYBOX); 

    %close all;
end % loop over obsq

