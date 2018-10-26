function FAULTAREAS = survey_fault_area(FAULTS,GRIDX,GRIDY,GRIDZ,dthresh)
% draw faults on planar slices
% input:
% FAULTS == structure containing fault coordinates
% GRIDX,GRIDY,GRIDZ 3-D arrays of coordinates
% 20181015 Kurt Feigl

warning('off'); %
%options = optimoptions(@fminunc,'Display','final');
options = optimoptions(@fminunc,'Display','none');

%% count number of voxels
[mx,nx,px] = size(GRIDX)
[my,ny,py] = size(GRIDY)
[mz,nz,pz] = size(GRIDZ)
if mx ~= my || my ~= mz
    error('miscount in X');
end
if nx~= ny || ny ~= nz
    error('miscount in Y');
end
if px ~= py || py ~= pz
    error('miscount in Z');
end
nvoxels = mx * my * mz


%
%% count number of faults
nfaults = numel(FAULTS);

% find increments
dX = mean(colvec(diff(GRIDX(:,1,1))));
dY = mean(colvec(diff(GRIDX(1,:,1))));
dZ = mean(colvec(diff(GRIDX(1,1,:))));

% define coordinates of centroids of voxels
FAULTAREAS.X = GRIDX;
FAULTAREAS.Y = GRIDY;
FAULTAREAS.Z = GRIDZ;
FAULTAREAS.count  = zeros(size(GRIDX)); % number of faults intersecting voxel [dimless]
FAULTAREAS.area   = zeros(size(GRIDX)); % surface area of fault inside voxel [m^2]
FAULTAREAS.density= zeros(size(GRIDX)); % fault area divided by voxel volume [m^2/m^3]

% linprog Linear programming.
%     X = linprog(f,A,b) attempts to solve the linear programming problem:
%
%              min f'*x    subject to:   A*x <= b
%               x
%
%     X = linprog(f,A,b,Aeq,beq) solves the problem above while additionally
%     satisfying the equality constraints Aeq*x = beq. (Set A=[] and B=[] if
%     no inequalities exist.)
%
%     X = linprog(f,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%     bounds on the design variables, X, so that the solution is in
%     the range LB <= X <= UB. Use empty matrices for LB and UB
%     if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below;
%     set UB(i) = Inf if X(i) is unbounded above.

% % set up objective function, one for each face of voxel
% % define faces by their normal vector
ff=zeros(6,3);
ff(1,:) = [ 1, 0,  0]; % near side of Y-Z face
ff(2,:) = [-1, 0,  0]; % far  side of Y-Z face
ff(3,:) = [ 0,-1,  0]; % near side of X-Z face
ff(4,:) = [ 0,+1,  0]; % far  side of X-Z face
ff(5,:) = [ 0, 0,  1]; % near side of X-Y face
ff(6,:) = [ 0, 0, -1]; % far  side of X-Y face

%
AA = [eye(3);eye(3)];

kount=0;
nvoxels=10;
for i=1:nvoxels
    tstart=tic;
    % set up inequality
    bb = nan(6,1);
    bb(1) = dX/2. - GRIDX(i);     % flip sign to be less than
    bb(2) = dY/2. - GRIDY(i);
    bb(3) = dZ/2. - GRIDZ(i);
    bb(4) = GRIDX(i)+dX/2. + eps; % to mimic strictly less than
    bb(5) = GRIDY(i)+dY/2. + eps;
    bb(6) = GRIDZ(i)+dZ/2. + eps;
    
    nfaults = 2;
    for j=1:nfaults
        %for i=1:2
        %for i=12 % for debugging
        %figure;
        if sum(FAULTS(j).trigood) ~= 0
            [ntriangles,n3] = size(FAULTS(j).TRI.ConnectivityList);
            fprintf(1,'Working on fault number %d named %s\n',j,char(FAULTS(j).faultname));
            %kount = 0;
            %             if ntriangles > 3 && n3 == 3
            %                 if isfinitearray(FAULTS(j).TRI.ConnectivityList) == 1
            for k=1:ntriangles
                % get indices of vertices
                kA=colvec(FAULTS(j).TRI.ConnectivityList(k,1));
                kB=colvec(FAULTS(j).TRI.ConnectivityList(k,2));
                kC=colvec(FAULTS(j).TRI.ConnectivityList(k,3));
                
                % get coordinates of vertices
                Xs=colvec([FAULTS(j).Xtri(kA),FAULTS(j).Xtri(kB),FAULTS(j).Xtri(kC)]);
                Ys=colvec([FAULTS(j).Ytri(kA),FAULTS(j).Ytri(kB),FAULTS(j).Ytri(kC)]);
                Zs=colvec([FAULTS(j).Ztri(kA),FAULTS(j).Ztri(kB),FAULTS(j).Ztri(kC)]);
                
                A3=colvec([FAULTS(j).Xtri(kA),FAULTS(j).Ytri(kA),FAULTS(j).Ztri(kA)]);
                B3=colvec([FAULTS(j).Xtri(kB),FAULTS(j).Ytri(kB),FAULTS(j).Ztri(kB)]);
                C3=colvec([FAULTS(j).Xtri(kC),FAULTS(j).Ytri(kC),FAULTS(j).Ztri(kC)]);
                
                % normal vector orthogonal to plane containing triangle
                N3 = cross(B3-A3,C3-A3);
                Aeq = rowvec(N3);
                
                % constant for equality for plane
                beq = N3'*N3;
                
                for kk=1:6
                    [xyz,fval,exitflag] = linprog(ff(kk,:),AA,bb,Aeq,beq,[],[],options);
                    if exitflag == 1
                        FAULTAREAS.count(i) = FAULTAREAS.count(i) + 1;
                    end
                end
                
                %                     end
                %                 end
                %
                kount = kount+1;
                %fprintf(1,'%d\r',kount);
                waitbar(kount/nvoxels/nfaults/ntriangles/6);
            end
        else
            warning(sprintf('bad fault: %s\n',FAULTS(j).faultname));
        end
        fprintf(1,'After %.1f seconds, counted %d grid points on fault number %d named %s\n'...
            ,toc(tstart),FAULTAREAS.count(i)...
            ,j,char(FAULTS(j).faultname));
    end
end
warning('on'); % turn warnings back on again
return
end

