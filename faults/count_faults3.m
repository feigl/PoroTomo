function FAULTCOUNTS = count_faults3(FAULTS,GRIDX,GRIDY,GRIDZ,dthresh)
% draw faults on planar slices
% input:
% FAULTS == structure containing fault coordinates
% 20181010 Kurt Feigl

warning('off'); % Interpolation function throws warning with duplicate points


%
%% count number of faults
nfaults = numel(FAULTS);

FAULTCOUNTS.X = GRIDX;
FAULTCOUNTS.Y = GRIDY;
FAULTCOUNTS.Z = GRIDZ;
FAULTCOUNTS.count = zeros(size(GRIDX));

for i=1:nfaults
%for i=1:2
    %for i=12 % for debugging
    %figure;
    if sum(FAULTS(i).trigood) ~= 0
        [ntriangles,n3] = size(FAULTS(i).TRI.ConnectivityList);
        fprintf(1,'Working on fault number %d named %s\n',i,char(FAULTS(i).faultname));
        %kount = 0;
        if ntriangles > 3 && n3 == 3
            if isfinitearray(FAULTS(i).TRI.ConnectivityList) == 1              
                % get indices of vertices
                kA=colvec(FAULTS(i).TRI.ConnectivityList(:,1));
                kB=colvec(FAULTS(i).TRI.ConnectivityList(:,2));
                kC=colvec(FAULTS(i).TRI.ConnectivityList(:,3));
                
                % get coordinates of vertices              
                Xs=colvec([FAULTS(i).Xtri(kA),FAULTS(i).Xtri(kB),FAULTS(i).Xtri(kC)]);
                Ys=colvec([FAULTS(i).Ytri(kA),FAULTS(i).Ytri(kB),FAULTS(i).Ytri(kC)]);
                Zs=colvec([FAULTS(i).Ztri(kA),FAULTS(i).Ztri(kB),FAULTS(i).Ztri(kC)]);
              
                
%                 % make interpolation function for Y as function of (X,Z)
%                 Fy=scatteredInterpolant(Xs,Zs,Ys,'nearest','none');
%                 
%                 % evaluate interpolation function at X,Y of grid points
%                 Gy = Fy(GRIDX,GRIDZ);
                
                % make interpolation function for X as function of (Y,Z)
                Fx=scatteredInterpolant(Ys,Zs,Xs,'nearest','none');
                
                % evaluate interpolation function at Y,Z of grid points
                Gx = Fx(GRIDY,GRIDZ);
                
                % is fault less than threshold distance?
                isfault = find(abs(Gx-GRIDX)< dthresh);
                ncount = numel(isfault);
                if ncount > 0
                    fprintf(1,'Counted %d grid points on fault number %d named %s\n',ncount,i,char(FAULTS(i).faultname));
                    FAULTCOUNTS.count(isfault) = FAULTCOUNTS.count(isfault) + 1;
                end
            end
            
        end
    else
        warning(sprintf('bad fault: %s\n',FAULTS(i).faultname));
    end
end

warning('on'); % turn warnings back on again
return
end

