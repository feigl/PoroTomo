function NODES = read_nodes(fname_nodes)
% read the file containing coordinates of nodes
% Coordinates are Xp,Yp,Zp in meters in Rotated PoroTomo coordinate system
% 20170910 Kurt Feigl

%fname_nodes = 'MESH_topo_xpypzellipsoid_800/nodes_coords_file';
fid_nodes = fopen(fname_nodes,'rt');
n_nodes1 = sscanf(fgetl(fid_nodes),'%d'); % get the number of nodes
C = textscan(fid_nodes,'%d %f %f %f');
fclose(fid_nodes);
NODES.id = C{1}; % ID number of NODE
NODES.Xp = C{2}; % X coordinate [m] of node in (rotated) PoroTomo coordinate system
NODES.Yp = C{3}; % Y coordinate [m] of node in (rotated) PoroTomo coordinate system
NODES.Zp = C{4}; % Z coordinate [m] of node in (rotated) PoroTomo coordinate system (positive upward above 800 m elevation)
n_nodes2 = numel(NODES.id);
if n_nodes1 == n_nodes2
    n_nodes = n_nodes1
else
    n_nodes1
    n_nodes2
    error('Miscount of nodes');
end
return
end

