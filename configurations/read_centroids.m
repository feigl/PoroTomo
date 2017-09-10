function CENTROIDS = read_centroids(fname_centroids)
%CENTROIDS = read_centroids(fname_centroids) - read centroid coordinates
% Coordinates are Xp,Yp,Zp in meters in Rotated PoroTomo coordinate system

%fname_centroids = 'MESH_topo_xpypzellipsoid_800/centroids_coords_file';
fid_centroids = fopen(fname_centroids,'rt');
n_centroids1 = sscanf(fgetl(fid_centroids),'%d'); % get the number of centroids
C = textscan(fid_centroids,'%d %f %f %f');
CENTROIDS.id = C{1};
n_centroids2 = numel(CENTROIDS.id);
fclose(fid_centroids);
if  n_centroids1 == n_centroids2
    n_centroids = n_centroids1
    CENTROIDS.xp = C{2};  % Centroid coordinate [m] Xp in rotated PoroTomo system
    CENTROIDS.yp = C{3};  % Centroid coordinate [m] Yp in rotated PoroTomo system
    CENTROIDS.zp = C{4};  % Centroid coordinate [m] Zp in rotated PoroTomo system
else
    n_centroids1
    n_centroids2
    error('Miscount of centroids');
end

return
end



