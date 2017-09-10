%% Save mesh metadata text files to structures and a .MAT file
% Based off of configuration_example1.m from Kurt Feigl (20170901)
% 20170907 Elena C. Reinsich 
% Downloads metadata files from the PoroTomo GitHub repo (master branch),
% loads data to structures, and saves all structures to an array
% currently pulls: nodes_coords_file, mesh_file, and centroids_coords_file

%% Initialize
clear all;
close all;

%% Set up path for Matlab
addpath(genpath('~/PoroTomo'),'-begin');

%% Define urls for mesh files
nodes_url = 'https://github.com/feigl/PoroTomo/raw/master/metadata_txt_files/nodes_coords_file';
elements_url = 'https://github.com/feigl/PoroTomo/raw/master/metadata_txt_files/mesh_file';
centroids_url = 'https://github.com/feigl/PoroTomo/raw/master/metadata_txt_files/centroids_coords_file';

%% Set up local directory for files (if doesn't already exist) and define text file short names 
[status,msg,msgID] = mkdir('MESH_topo_xpypzellipsoid_800');
% node short name
url_splt = strsplit(nodes_url, '/');
fname_nodes = sprintf('%s/%s', 'MESH_topo_xpypzellipsoid_800', char(url_splt(end)));
% elements short name
url_splt = strsplit(elements_url, '/');
fname_elements = sprintf('%s/%s', 'MESH_topo_xpypzellipsoid_800', char(url_splt(end)));
% centroids short name
url_splt = strsplit(centroids_url, '/');
fname_centroids = sprintf('%s/%s', 'MESH_topo_xpypzellipsoid_800', char(url_splt(end)));

%% read the file containing coordinates of nodes
% pull most recent, stable version from GitHub and save to text file
websave(fname_nodes, nodes_url);

% read text file
fid_nodes = fopen(fname_nodes,'rt');
n_nodes1 = sscanf(fgetl(fid_nodes),'%d'); % get the number of nodes
C = textscan(fid_nodes,'%d %f %f %f');
fclose(fid_nodes);
NODES.id = C{1}; % ID number of NODE
NODES.xp = C{2}; % X coordinate [m] of node in (rotated) PoroTomo coordinate system
NODES.yp = C{3}; % Y coordinate [m] of node in (rotated) PoroTomo coordinate system
NODES.zp = C{4}; % Z coordinate [m] of node in (rotated) PoroTomo coordinate system (positive upward above 800 m elevation)
n_nodes2 = numel(NODES.id);
if n_nodes1 == n_nodes2
    n_nodes = n_nodes1
else
    n_nodes1
    n_nodes2
    error('Miscount of nodes');
end

%% read the file containing the list of elements
% pull most recent, stable version from GitHub and save to text file
websave(fname_elements, elements_url);

% read text file
fid_elements = fopen(fname_elements,'rt');
n_elements1 = sscanf(fgetl(fid_elements),'%d'); % get the number of elements
C = textscan(fid_elements,'%d %d %d %d %d %d %d %d %d');
fclose(fid_elements);
ELEMENTS.id = C{1}; % ID number of element
n_elements2 = numel(ELEMENTS.id);
if n_elements1 == n_elements2
    n_elements = n_elements1
else
    n_elements1
    n_elements2
    error('Miscount of elements');
end
ELEMENTS.nodes = nan(n_elements,8);
for j=1:8
    ELEMENTS.nodes(:,j) = C{1+j}; % ID numbers of eight nodes (vertices) comprising a (hexahedral) element (i.e. brick with 6 faces)
end

%% read the file containing coordinates of centroids
% pull most recent, stable version from GitHub and save to text file
websave(fname_centroids, centroids_url);

% read text file
fid_centroids = fopen(fname_centroids,'rt');
n_centroids1 = sscanf(fgetl(fid_centroids),'%d'); % get the number of nodes
C = textscan(fid_centroids,'%d %f %f %f');
fclose(fid_centroids);
CENTROIDS.id = C{1}; % ID number of NODE
CENTROIDS.xp = C{2}; % X coordinate [m] of node in (rotated) PoroTomo coordinate system
CENTROIDS.yp = C{3}; % Y coordinate [m] of node in (rotated) PoroTomo coordinate system
CENTROIDS.zp = C{4}; % Z coordinate [m] of node in (rotated) PoroTomo coordinate system (positive upward above 800 m elevation)
n_centroids2 = numel(CENTROIDS.id);
if n_centroids1 == n_centroids2
    n_centroids = n_centroids1
else
    centroids1
    centroids2
    error('Miscount of nodes');
end

save('MESH.mat', 'NODES', 'ELEMENTS', 'CENTROIDS')

% %% draw the mesh in 3-D
% figure; hold on;
% % plot3(ELEMENTS.centroid_xp,ELEMENTS.centroid_yp,ELEMENTS.centroid_zp,'*r');
% plot3(CENTROIDS.xp,CENTROIDS.yp,CENTROIDS.zp,'*r');
% 
% % very slow
% % for i=1:n_elements
% %     plot3(ELEMENTS.centroid_xp(i),ELEMENTS.centroid_yp(i),ELEMENTS.centroid_zp(i),'*r');
% % end