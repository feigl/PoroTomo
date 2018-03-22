function draw_faults(S,knorm,constant_coordinate1,nmin,norder,plot_points,plot_curves,ksor)
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
nfaults = numel(S);

dc = 10 ; % spacing in meters
for i=1:nfaults
    %figure;
    [ntriangles,n3] = size(S(i).kfaces);
    if ntriangles > 3 && n3 == 3
        %ntriangles = 10;
        P0 = nan(3*ntriangles,1); % coordinate for sorting
        P1 = nan(3*ntriangles,1); % coordinate on horizontal axis of plot
        P2 = nan(3*ntriangles,1); % coordinate on vertical axis of plot
        kount = 0;
        %for constant_coordinate = constant_coordinate1-dc:dc:constant_coordinate1+dc
        for constant_coordinate = constant_coordinate1
            for j=1:ntriangles
                if isfinitearray(S(i).kfaces(j,:)) == 1
                    % find indices of vertices
                    kA=colvec(S(i).kfaces(j,1));
                    kB=colvec(S(i).kfaces(j,2));
                    kC=colvec(S(i).kfaces(j,3));
                    
                    % get coordinates of vertices
                    A=colvec([S(i).PR.x(kA),S(i).PR.y(kA),S(i).PR.z(kA)]);
                    B=colvec([S(i).PR.x(kB),S(i).PR.y(kB),S(i).PR.z(kB)]);
                    C=colvec([S(i).PR.x(kC),S(i).PR.y(kC),S(i).PR.z(kC)]);
                    
                    % define a plane by two vectors along sides of triangle
                    U = B - A;
                    V = C - A;
                    
                    P = nan(3,1);
                    intri = 0;
                    % first solution: set s = 0 and solve for t
                    tt = (constant_coordinate - A(knorm))/V(knorm);
                    if 0 <= tt && tt <= 1
                        t = tt;
                        s = 0;
                        kount = kount + 1;
                        P = A + s*U + t*V;
                        P0(kount) = P(ksor);
                        P1(kount) = P(kfir);
                        P2(kount) = P(ksec);
                        intri = intri + 1;
                    end
                    
                    % second solution: set t = 0 and solve for s
                    ss = (constant_coordinate - A(knorm))/U(knorm);
                    if 0 <= ss && ss <= 1
                        t = 0;
                        s = ss;
                        kount = kount + 1;
                        P = A + s*U + t*V;
                        P0(kount) = P(ksor);
                        P1(kount) = P(kfir);
                        P2(kount) = P(ksec);
                        intri = intri + 1;
                    end
                    
                    % third solution: set s = 1 - t and solve for t
                    tt = (constant_coordinate - U(knorm) - A(knorm)) / (V(knorm) - U(knorm));
                    if 0 <= tt && tt <= 1
                        s = 1 - tt;
                        t = tt;
                        kount = kount + 1;
                        P = A + s*U + t*V;
                        P0(kount) = P(ksor);
                        P1(kount) = P(kfir);
                        P2(kount) = P(ksec);
                        intri = intri + 1;
                    end
                    if intri ~= 2
                        kount = kount - intri;
                    end
                else
                    warning(sprintf('bad triangle %d in fault: %d %s\n',j,i,S(i).faultname));
                end
            end
        end
        
        % prune
        iok = intersect(find(isfinite(P1)),find(isfinite(P2)));
        ndata = numel(iok);
        
        
        % plot only faults with sufficient number of points
        % nmin = 5
        if ndata > nmin
            P0 = P0(iok);
            P1 = P1(iok);
            P2 = P2(iok);
            
            
            %         %         %   fit a polynomial curve
            %         %         %norder = 3;
            %         [Pfit,Sfit,MU] = polyfit(P1,P2,norder);
            %         if Sfit.normr < 5e3
            %
            %             if plot_curves == 1
            %                 P1 = sort(P1);
            %                 [Zfit,Dz] = polyval(Pfit,P1,Sfit,MU);
            %                 plot(P1,Zfit,'k-','LineWidth',3);
            %             end
            %         else
            %             warning(sprintf('Poor polyfit (normr = %10.0f). Skipping fault %s\n',Sfit.normr,S(i).name));
            %             %close gcf
            %         end
            %         % sort by P0 coordinate
            [P0,isort] = sort(P0,'descend');
            Ps1 = P1(isort);
            Ps2 = P2(isort);
            %         % sort by P1 coordinate
            %         [P1,isort] = sort(P1);
            %         P2 = P2(isort);
            %         P1 = smooth(P0,P1,min([25,ndata]));
            %         P2 = smooth(P0,P2,min([25,ndata]));
            %         P1 = smooth(P1,min([100,ndata]));
            %         P2 = smooth(P2,min([100,ndata]));
            %        P2 = smooth(P1,P2,min([50,floor(ndata/5)]));
            %         if plot_curves == 1
            %             plot(P1,P2,'k-','LineWidth',3);
            %         end
            %          idx = knnsearch([colvec(P1),colvec(P2)],[colvec(P1),colvec(P2)]);
            %          if plot_curves == 1
            %             plot(P1(idx(:,1)),P2(idx(:,1)),'k-','LineWidth',3);
            %         end
            Ps1=movmean(Ps1,min([10,floor(ndata/5)]));
            Ps2=movmean(Ps2,min([10,floor(ndata/5)]));
            if plot_curves == 1
                plot(Ps1,Ps2,'k-','LineWidth',3);
            end
        end
        %   draw intersections as dots
        if plot_points == 1
            plot(P1,P2,'w.','MarkerSize',1);
        end
    else
        warning(sprintf('bad fault: %s\n',S(i).faultname));
    end
end

% % return handle for current figure
% h = gcf;


return
end

