function draw_faults4(S,knorm,constant_coordinate1,nmin,norder,plot_points,plot_curves,ksor)
% draw faults on planar slices
% input:
% S == structure containing fault coordinates
% knorm  == index to axis of normal to planar slice:
%     knorm == 1 is vertical slice normal to PoroTomo X axis
%     knorm == 2 is vertical slice normal to PoroTomo Y axis
%     knorm == 3 is horizontal slice normal to PoroTomo Z axis
% constant_coordinate == values in meters of constant coordinate
% nmin == minimum number of intersecting points to qualify a fault for plot
% norder == order of polynomial fit (norder = 1 is linear)
% plot_points == 1 to plot intersections as red dots
% plot_curves == 1 to plot polynomial as black lines
% 20170621 Kurt Feigl
% 20170727 include range of values for coordinates
% 20180525 use "line" function

%% arrange coordinates
switch knorm
    case 1
        % slicing plane is normal to X axis
        ksor = 2;
        kfir = 2; % Y is first coordinate
        ksec = 3; % Z is second coordinate
    case 2
        % slicing plane is normal to Y axis
        ksor = 3;
        kfir = 1; % X is first coordinate
        ksec = 3; % Z is second coordinate
    case 3
        % slicing plane is normal to Z axis
        ksor = 2;
        kfir = 1; % X is first coordinate
        ksec = 2; % Y is second coordinate
    otherwise
        error(sprintf('unknown knorm %d\n',knorm));
end

% count number of faults
nfaults = numel(S)

dc = 10 ; % spacing in meters
constant_coordinate = constant_coordinate1

for i=1:nfaults
%for i=12 % for debugging
    %figure;
    if sum(S(i).trigood) ~= 0
        [ntriangles,n3] = size(S(i).TRI.ConnectivityList);
        fprintf(1,'Working on fault number %d named %s\n',i,char(S(i).faultname));
        kount = 0;
        if ntriangles > 3 && n3 == 3
            %ntriangles = 10;
            %         P0 = nan(3*ntriangles,1); % coordinate for sorting
            %         P1 = nan(3*ntriangles,1); % coordinate on horizontal axis of plot
            %         P2 = nan(3*ntriangles,1); % coordinate on vertical axis of plot
            P0 = nan(2,3*ntriangles); % coordinate for sorting
            P1 = nan(2,3*ntriangles); % coordinate on horizontal axis of plot
            P2 = nan(2,3*ntriangles); % coordinate on vertical axis of plot
            
            %for constant_coordinate = constant_coordinate1-dc:dc:constant_coordinate1+dc
            for j=1:ntriangles
                %if isfinitearray(S(i).TRI.ConnectivityList(j,:)) == 1
                % find indices of vertices
                kA=colvec(S(i).TRI.ConnectivityList(j,1));
                kB=colvec(S(i).TRI.ConnectivityList(j,2));
                kC=colvec(S(i).TRI.ConnectivityList(j,3));
                
                % get coordinates of vertices
                %                     A=colvec([S(i).PR.x(kA),S(i).PR.y(kA),S(i).PR.z(kA)]);
                %                     B=colvec([S(i).PR.x(kB),S(i).PR.y(kB),S(i).PR.z(kB)]);
                %                     C=colvec([S(i).PR.x(kC),S(i).PR.y(kC),S(i).PR.z(kC)]);
                A=colvec([S(i).Xtri(kA),S(i).Ytri(kA),S(i).Ztri(kA)]);
                B=colvec([S(i).Xtri(kB),S(i).Ytri(kB),S(i).Ztri(kB)]);
                C=colvec([S(i).Xtri(kC),S(i).Ytri(kC),S(i).Ztri(kC)]);
                
                % define a plane by two vectors along sides of triangle
                U = B - A;
                V = C - A;               
                P = nan(3,1);
                
                % count number of points inside the triangle
                intri = 0;
                
                % first solution: set s = 0 and solve for t
                tt = (constant_coordinate - A(knorm))/V(knorm);
                if 0 <= tt && tt <= 1
                    t = tt;
                    s = 0;
                    P = A + s*U + t*V;
                    kount = kount + 1;
                    intri = intri + 1;
                    P0(intri,kount) = P(ksor);
                    P1(intri,kount) = P(kfir);
                    P2(intri,kount) = P(ksec);
                end
                
                % second solution: set t = 0 and solve for s
                ss = (constant_coordinate - A(knorm))/U(knorm);
                if 0 <= ss && ss <= 1
                    t = 0;
                    s = ss;
                    P = A + s*U + t*V;
                    if intri==0
                        kount = kount + 1;
                    end
                    intri = intri + 1;
                    P0(intri,kount) = P(ksor);
                    P1(intri,kount) = P(kfir);
                    P2(intri,kount) = P(ksec);
                end
                
                if intri == 1
                    intri = intri + 1;
                    % third solution: set s = 1 - t and solve for t
                    tt = (constant_coordinate - U(knorm) - A(knorm)) / (V(knorm) - U(knorm));
                    if 0 <= tt && tt <= 1
                        s = 1 - tt;
                        t = tt;
                        P = A + s*U + t*V;
                        P0(intri,kount) = P(ksor);
                        P1(intri,kount) = P(kfir);
                        P2(intri,kount) = P(ksec);
                    end
                end
            end
            if plot_curves == 1 || plot_points == 1
                kok = intersect(find(isfinite(P1)==1),find(isfinite(P2)==1));
                [iok,jok] = ind2sub(size(P1),kok);
                if numel(iok) > 2
%                     whos P1 P2
%                     P1(iok,jok)
%                     P2(iok,jok)
                    fprintf(1,'plotting %d points in fault: %d %s\n',numel(iok),i,S(i).faultname);
                    if plot_curves == 1
                        line(P1(iok,jok),P2(iok,jok),'color','k','LineWidth',3);
                    end
                    if plot_points == 1
                        plot(P1,P2,'w.','MarkerSize',1);
                    end
                end
            end
        else
            warning(sprintf('bad fault: %s\n',S(i).faultname));
        end
    else
        warning(sprintf('bad triangle %d in fault: %d %s\n',j,i,S(i).faultname));
    end
end
    

% % return handle for current figure
% h = gcf;


return
end

