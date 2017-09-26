function BRADYBOX = read_bradybox(isite)
% return 3-coordinates of study area
% 20170925 Kurt Feigl

switch isite
    case 1 % Brady PoroTomo
        % get UTM coordinates
        i=0;
        i=i+1;bradybox_UTM(i,:) = [ 328761.1411         4408840.1323];
        i=i+1;bradybox_UTM(i,:) = [ 327850.8122         4407606.2053];
        i=i+1;bradybox_UTM(i,:) = [ 328221.6320         4407332.5948];
        i=i+1;bradybox_UTM(i,:) = [ 329137.6472         4408559.8842];
        BRADYBOX.E = colvec(bradybox_UTM(:,1)); % UTM easting in meters
        BRADYBOX.N = colvec(bradybox_UTM(:,2)); % UTM easting in meters
        
    case 2 % Brady GIPhT
        % add bounding box for gipht plots (edit ECR) >
        i=0;
        i=i+1;giphtintbox_UTM(i,:) = [ 326500         4405500];
        i=i+1;giphtintbox_UTM(i,:) = [ 326500         4409500];
        i=i+1;giphtintbox_UTM(i,:) = [ 329500         4409500];
        i=i+1;giphtintbox_UTM(i,:) = [ 329500         4405500];
        BRADYBOX.E = colvec(giphtintbox_UTM(:,1)); % UTM easting in meters
        BRADYBOX.N = colvec(giphtintbox_UTM(:,2)); % UTM easting in meters        
    otherwise
        error(sprintf('Unknown isite = %d\n',isite));
end


return
end



