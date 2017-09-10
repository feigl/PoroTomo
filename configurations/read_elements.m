function ELEMENTS = read_elements(fname_elements )
% read the file containing the list of elements for a hexahedral mesh with
% 8 nodes per element
% 2017 Kurt Feigl
fname_elements = 'MESH_topo_xpypzellipsoid_800/mesh_file';
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

return
end
