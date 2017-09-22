%% Save mesh metadata text files to structures and a .MAT file
% Based off of configuration_example1.m from Kurt Feigl (20170901)
% 20170907 Elena C. Reinsich 
% Downloads metadata files from the PoroTomo GitHub repo (master branch),
% loads data to structures, and saves all structures to an array
% currently pulls: nodes_coords_file, mesh_file, and centroids_coords_file
% Update ECR 20170920 adding section for downloading xlsx files from GitHub

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

%% Get xlsx files
% Initialize structure
METADATA = struct;

% pull down list of metadata xlsx files
xlsxlist_url = 'https://github.com/feigl/PoroTomo/raw/master/metadata_txt_files/xlsx_list.txt';
websave('xlsx_list.txt', xlsxlist_url)

% load names
xlsx_table = readtable('xlsx_list.txt', 'ReadVariableNames', 0);

% initialize text file for writing short name mapping keys 
fid = fopen('metadata_name_mapping.txt', 'w+');

% figure out path to wget 
% try first to find GNU or FreeBSD version
[test1, wgetpath] = unix('find /usr -name wget -print -quit');
% if neither of those, try NetBSD
if isempty(strfind(wgetpath, '/wget'))
    [test1, wgetpath] = unix('find /usr -name wget -print -exit');
    % if still can't find version with quick method, print out all matches
    % and then take path
    if isempty(strfind(wgetpath, '/wget'))
        [test1, wgetpath] = unix('find /usr -name wget | head -1');
    end
end

% cycle through names, pulling each file down and saving to structure as we go
for k = 1:numel(xlsx_table)
   xlsx_name =  char(table2array(xlsx_table(k,1)))
   
   % get data file from GITHub
   command = strcat(wgetpath, ' --no-check-certificate https://raw.githubusercontent.com/feigl/PoroTomo/master/metadata_txt_files/', xlsx_name);
   system(command);
   
   % read file that has been downloaded
   data = readtable(xlsx_name, 'ReadVariableNames', 1); 
   [ndata nfields] = size(data);
   
   % assign shortname
   if find(regexpi(xlsx_name, 'nodal')) == 1
       data_short_name = 'NODE';
   elseif find(regexpi(xlsx_name, 'reftek')) == 1
       data_short_name = 'REFT';
   elseif find(regexpi(xlsx_name, 'vibroseis')) == 1
       data_short_name = 'VIBRO';
%    elseif find(regexpi(xlsx_name, 'temp')) == 1
%        if find(regexpi(xlsx_name, 'dtsv')) == 1
%            data_short_name = 'DTSV';
%        elseif find(regexpi(xlsx_name, 'dtsh')) == 1
%            data_short_name = 'DTSH';
%        end
%    elseif find(regexpi(xlsx_name, 'surface')) == 1 & find(regexpi(xlsx_name, 'DAS')) == 1 & find(regexpi(xlsx_name, 'DTS')) == 1
%        data_short_name = 'DASsurf_DTSsurf';
   elseif find(regexpi(xlsx_name, 'borehole')) == 1 & find(regexpi(xlsx_name, 'DAS')) == 1 & find(regexpi(xlsx_name, 'DTS')) == 1
       data_short_name = 'DASV_DTSV';
   elseif find(regexpi(xlsx_name, 'horizontal')) == 1 & find(regexpi(xlsx_name, 'DAS')) == 1 & find(regexpi(xlsx_name, 'DTS')) == 1
       data_short_name = 'DASH_DTSH';
   elseif find(regexpi(xlsx_name, 'dtsv')) == 1
       data_short_name = 'DTSV';
   elseif find(regexpi(xlsx_name, 'dtsh')) == 1
       data_short_name = 'DTSH';
   elseif find(regexpi(xlsx_name, 'permaseis')) == 1
       data_short_name = 'PERMS';
   elseif find(regexpi(xlsx_name, 'insar')) == 1
       data_short_name = 'INSAR';
   elseif find(regexpi(xlsx_name, 'pressure')) == 1
       data_short_name = 'PDAT';
   elseif find(regexpi(xlsx_name, 'geol')) == 1
       data_short_name = 'GEOL';
   elseif find(regexpi(xlsx_name, 'pump')) == 1
       data_short_name = 'PUMP';
   elseif find(regexpi(xlsx_name, 'uav')) == 1
       data_short_name = 'UAV';
   else
       data_short_name = input('Manual entry of data short name needed:', 's')
   end

   % save name of short name to file
   fprintf(fid, '%s %s\n', data_short_name, xlsx_name);
   
   % check if values were recorded correctly (indication of  a second
    % header line in file)
    check = 0;
    for k = 1:nfields
        % count number of fields with class double
       check =  check+double(strcmp(class(data.(k)), 'double'));
    end
    
    % if no fields are class double, then 2 header lines are present
    if check == 0
        fprintf('\n')
        disp('NOTE: Possibility of 2 header lines detected.  Assigning information in second line as separate field of units. Please double check result.')
        fprintf('\n')
        
        % pull unit label line from table
        uvals = table2array(data(1,:));
        data(1,:) = [];
        
        % assign units to variable names
        data.Properties.VariableUnits = uvals;
    end
    
    % convert data to structure 
    METADATA.(data_short_name) = table2struct(data,'ToScalar',true);
    
    % add extra field for units if not already included in field names
    if check == 0
        METADATA.(data_short_name).Units = data.Properties.VariableUnits;
    end
   
end

% close name mapping file
fclose(fid);

% save metadata to .mat file with date in name
 save(strcat('METADATA_', date, '.mat'), 'METADATA')
 
 fprintf('Downloads have completed. Data is saved to %s \n', strcat('METADATA_', date, '.mat'))
