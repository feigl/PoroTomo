function figfilenames = plot_tomo_and_faults3(TOMO,FAULTS,title_str,SLICES,OPTIONS,funprint,rsearch, WELLS, elev_mean,BRADY2DGRID)
% make slices at SLICES of tomogram TOMO with faults in FAULTS
% 20170921 Kurt Feigl

% initialize
kfiles = 0;
figfilenames=[];

nx = size(TOMO.Xp)
ny = size(TOMO.Yp)
nz = size(TOMO.Zp)

tstart = tic;

% set faults
if exist('OPTIONS','var') == 1
    if isfield(OPTIONS,'vmin') == 1
        vmin=OPTIONS.vmin;
    else
        %vmin=nanmin(colvec(TOMO.V));
        vmin = nan;
    end
    if isfield(OPTIONS,'vmax') == 1
        vmax=OPTIONS.vmax;
    else
        %vmax=nanmax(colvec(TOMO.V));
        vmax = nan;
    end
    if isfield(OPTIONS,'cmap') == 1
        cmap=OPTIONS.cmap;
    else
        cmap=flipud(colormap('jet'));
    end
    if isfield(OPTIONS,'nmin') == 1
        nmin=OPTIONS.nmin;
    else
        nmin = 2;
    end
    if isfield(OPTIONS,'plot_points') == 1
        plot_points=OPTIONS.plot_points;
    else
        plot_points = 1;
    end
    if isfield(OPTIONS,'norder') == 1
        norder=OPTIONS.norder;
    else
        norder=2;
    end
    if isfield(OPTIONS,'plot_curves') == 1
        plot_curves=OPTIONS.plot_curves;
    else
        plot_curves = 1;
    end
    if isfield(OPTIONS,'interpolation_method') == 1
        interpolation_method=OPTIONS.interpolation_method;
    else
        interpolation_method = 'linear';
    end
    if isfield(OPTIONS,'extrapolation_method') == 1
        extrapolation_method=OPTIONS.extrapolation_method;
    else
        extrapolation_method = 'linear';
    end
end


% count number of slices
nxslices = numel(SLICES.Xp)
nyslices = numel(SLICES.Yp)
nzslices = numel(SLICES.Zp)



%% loop over normal direction
for knorm = [1,2,3]
    switch knorm
        case 1
            nslices=nxslices;
        case 2
            nslices=nyslices;
        case 3
            nslices=nzslices;
        otherwise
            error(sprintf('unknown knorm = %d\n',knorm));
    end
    %% loop over slices
    for kslice = 1:nslices
        % select sample points in bounds
        ikeep = 1:numel(TOMO.Xp);
        if knorm ~= 1
            ikeep = intersect(ikeep,find(TOMO.Xp <= nanmax(colvec(SLICES.Xp))));
            ikeep = intersect(ikeep,find(TOMO.Xp >= nanmin(colvec(SLICES.Xp))));
        end
        if knorm ~= 2
            ikeep = intersect(ikeep,find(TOMO.Yp <= nanmax(colvec(SLICES.Yp))));
            ikeep = intersect(ikeep,find(TOMO.Yp >= nanmin(colvec(SLICES.Yp))));
        end
        if knorm ~= 3
            ikeep = intersect(ikeep,find(TOMO.Zp <= nanmax(colvec(SLICES.Zp))));
            ikeep = intersect(ikeep,find(TOMO.Zp >= nanmin(colvec(SLICES.Zp))));
        end
        if numel(ikeep) < 10
            numel(ikeep)
            error('miscount 3');
        end
        
        switch knorm
            % build 3-D plaid grid using all 3 coordinates in slicing plane
            %             case 1
            %                 [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp(kslice),  SLICES.Yp,        SLICES.Zp);
            %             case 2
            %                 [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp,          SLICES.Yp(kslice),SLICES.Zp);
            %             case 3
            %                 [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp,SLICES.Yp,                  SLICES.Zp(kslice));
            
            % build 3-D plaid grid with original node points on slicing plane          
            case 1
                ikeep = intersect(ikeep,find(abs(TOMO.Xp - SLICES.Xp(kslice)) < rsearch));
                nkeep = numel(ikeep)
                [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp(kslice),unique(colvec(TOMO.Yp(ikeep))),unique(colvec(TOMO.Zp(ikeep))));
            case 2
                ikeep = intersect(ikeep,find(abs(TOMO.Yp - SLICES.Yp(kslice)) < rsearch));
                nkeep = numel(ikeep)
                [Q.x,Q.y,Q.z] = meshgrid(unique(colvec(TOMO.Xp(ikeep))),SLICES.Yp(kslice),unique(colvec(TOMO.Zp(ikeep))));
            case 3
                ikeep = intersect(ikeep,find(abs(TOMO.Zp - SLICES.Zp(kslice)) < rsearch));
                nkeep = numel(ikeep)
                [Q.x,Q.y,Q.z] = meshgrid(unique(colvec(TOMO.Xp(ikeep))),unique(colvec(TOMO.Yp(ikeep))),SLICES.Zp(kslice));
            otherwise
                error(sprintf('unknown knorm = %d\n',knorm));
        end
        
        % remove singletons to make 2-D meshes
        switch knorm
            case 1
                Q1 = squeeze(Q.y);
                Q2 = squeeze(Q.z);
            case 2
                Q1 = squeeze(Q.x);
                Q2 = squeeze(Q.z);
            case 3
                Q1 = squeeze(Q.x);
                Q2 = squeeze(Q.y);
            otherwise
                error(sprintf('unknown knorm = %d\n',knorm));
        end
 
        % interpolate values from 3D onto the 2-D mesh
        %image2d = griddata(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,Q.x,Q.y,Q.z,interpolation_method)
        mnpQx = size(Q.x)
        mnpQy = size(Q.y)
        mnpQz = size(Q.z)
        fprintf(1,'Starting 3-D scatteredInterpolant at t = %.001f seconds\n',toc(tstart));
        %extrapolation_method = 'none';
        %extrapolation_method = 'linear';
        Finterp2 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V ...
            ,interpolation_method,extrapolation_method);
        image2d = Finterp2(Q.x,Q.y,Q.z);
        
        % make a finer mesh with 10-meter spacing
        switch knorm
            case 1
                [M1,M2] = meshgrid([ceil(nanmin(colvec(SLICES.Yp))):10:floor(nanmax(colvec(SLICES.Yp)))]...
                    ,[ceil(nanmin(colvec(SLICES.Zp))):10:floor(nanmax(colvec(SLICES.Zp)))]);
            case 2
                [M1,M2] = meshgrid([ceil(nanmin(colvec(SLICES.Xp))):10:floor(nanmax(colvec(SLICES.Xp)))]...
                    ,[ceil(nanmin(colvec(SLICES.Zp))):10:floor(nanmax(colvec(SLICES.Zp)))]);
            case 3
                [M1,M2] = meshgrid([ceil(nanmin(colvec(SLICES.Xp))):10:floor(nanmax(colvec(SLICES.Xp)))]...
                    ,[ceil(nanmin(colvec(SLICES.Yp))):10:floor(nanmax(colvec(SLICES.Yp)))]);
            otherwise
                error(sprintf('unknown knorm = %d\n',knorm));
        end
        
       
        
        % interpolate onto 2-D mesh
        if numel(colvec(Q1)) > 10 && numel(colvec(Q2)) > 10
            %image2dinterp = interp2(Q1,Q2,squeeze(image2d),M1,M2);
            fprintf(1,'Starting 2-D scatteredInterpolant at t = %.001f seconds\n',toc(tstart));
            Finterpolant = scatteredInterpolant(colvec(Q1),colvec(Q2),colvec(squeeze(image2d))...
            ,interpolation_method,extrapolation_method);
            %whos Finterpolant
            image2dinterp = Finterpolant(colvec(M1),colvec(M2));
            if numel(image2dinterp) == numel(M1)
                image2dinterp = reshape(image2dinterp,size(M1));
                
                % plot the tomogram
                figure;
                hold on;
                
                % setting the velocity range for plotting
                %caxis([v1, v2]);
                if isfinite(vmin) == 1 && isfinite(vmax)==1
                    caxis([vmin,vmax]);
                    caxis manual;
                else
                    vmin1 = nanmin(colvec(image2dinterp));
                    vmax1 = nanmax(colvec(image2dinterp));
                    if isfinite(vmin1) == 1 && isfinite(vmax1)==1
                        caxis([vmin1,vmax1]);
                        caxis manual;
                    else
                        caxis auto;
                    end
                end
                
                %colormap(flipud(colormap('jet')));
                colormap(cmap);
                
                % plot in color
                pcolor(M1,M2,image2dinterp);
                shading interp;
                
                %% set the axes
                hold on;
                axis equal
                axis xy
                axis tight;
                %  same size as the slices
                switch knorm
                    case 1
                        axis([nanmin(SLICES.Yp), nanmax(SLICES.Yp), nanmin(SLICES.Zp), nanmax(SLICES.Zp)]);
                    case 2
                        axis([nanmin(SLICES.Xp), nanmax(SLICES.Xp), nanmin(SLICES.Zp), nanmax(SLICES.Zp)]);
                    case 3
                        axis([nanmin(SLICES.Xp), nanmax(SLICES.Xp), nanmin(SLICES.Yp), nanmax(SLICES.Yp)]);
                    otherwise
                        error(sprintf('unknown knorm = %d\n',knorm));
                end
                %             % same size as the tomogram
                %             switch knorm
                %                 case 1
                %                     axis([nanmin(TOMO.Yp), nanmax(TOMO.Yp), nanmin(TOMO.Zp), nanmax(TOMO.Zp)]);
                %                 case 2
                %                     axis([nanmin(TOMO.Xp), nanmax(TOMO.Xp), nanmin(TOMO.Zp), nanmax(TOMO.Zp)]);
                %                 case 3
                %                     axis([nanmin(TOMO.Xp), nanmax(TOMO.Xp), nanmin(TOMO.Yp), nanmax(TOMO.Yp)]);
                %                 otherwise
                %                     error(sprintf('unknown knorm = %d\n',knorm));
                %             end
                
                % clip outside the axes
                ax = gca;
                ax.Clipping='On';
                
                
                
                % draw contours
                %[c,h]=contour(M1,M2,image2dinterp,[v1:dv:v2],'w');
                
                %         S == structure containing fault coordinates
                %         nmin =5; % minimum number of intersecting points to qualify a fault for plot
                %         norder = 2;  % order of polynomial fit (norder = 1 is linear)
                %         plot_points = 1; % plot intersections as black dots
                %         plot_curves = 1; % plot polynomial as black lines
                switch knorm
                    case 1
                        constant_coordinate = SLICES.Xp(kslice);
                    case 2
                        constant_coordinate = SLICES.Yp(kslice);
                    case 3
                        constant_coordinate = SLICES.Zp(kslice);
                    otherwise
                        error(sprintf('unknown knorm = %d\n',knorm));
                end
                
                % sort according to Y coordinate
                ksor  = 2;
                
                % draw the faults
                %draw_faults(S,knorm,constant_coordinate,nmin,norder,plot_points,plot_curves,ksor);
                draw_faults2(FAULTS,knorm,constant_coordinate,nmin,norder,plot_points,plot_curves,ksor);
                
                % draw the wells
                switch knorm
                    case 1
                        for i=1:numel(WELLS.grdsurf.Yp)
                            if abs(WELLS.grdsurf.Xp(i)-SLICES.Xp(kslice)) < 10*rsearch
                                plot([WELLS.grdsurf.Yp(i),WELLS.openbtm.Yp(i)]...
                                    ,[WELLS.grdsurf.Zp(i),WELLS.openbtm.Zp(i)]...
                                    ,'k^:','LineWidth',1);
                            end
                        end
                    case 2
                        for i=1:numel(WELLS.grdsurf.Xp)
                            if abs(WELLS.grdsurf.Yp(i)-SLICES.Yp(kslice)) < 10*rsearch
                                plot([WELLS.grdsurf.Xp(i),WELLS.openbtm.Xp(i)]...
                                    ,[WELLS.grdsurf.Zp(i),WELLS.openbtm.Zp(i)]...
                                    ,'k^:','LineWidth',1);
                            end
                        end
                    case 3
                        plot(WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,'k^','MarkerSize',2,'MarkerFaceColor','k');
                    otherwise
                        error(sprintf('unknown knorm = %d\n',knorm));
                end
                
                %% draw the sample points
                if plot_points == 1
                    switch knorm
                        case 1
                            for i=1:numel(TOMO.Xp)
                                if abs(TOMO.Xp(i)-SLICES.Xp(kslice)) < 10*rsearch
                                    plot(TOMO.Yp(i),TOMO.Zp(i),'ok','MarkerSize',7);
                                end
                            end
                        case 2
                            for i=1:numel(TOMO.Yp)
                                if abs(TOMO.Yp(i)-SLICES.Yp(kslice)) < 10*rsearch
                                    plot(TOMO.Xp(i),TOMO.Zp(i),'ok','MarkerSize',7);
                                end
                            end
                        case 3
                            plot(TOMO.Xp,TOMO.Yp,'ok','MarkerSize',7);
                        otherwise
                            error(sprintf('unknown knorm = %d\n',knorm));
                    end
                end

                %% label axes
                switch knorm
                    case 1
                        % slicing plane is normal to X axis
                        title(sprintf('X_{PoroTomo} = %10.1f meters (%s)',constant_coordinate,strrep(title_str,'_',' ')));
                        xlabel('<-- SW       Y_{PoroTomo} [m]        NE -->');
                        ylabel('Z_{PoroTomo} [m]'); % (= elevation above WGS84 ellipsoid + 800 m'
                    case 2
                        % slicing plane is normal to Y axis
                        title(sprintf('Y_{PoroTomo} = %10.1f meters (%s)',constant_coordinate,strrep(title_str,'_',' ')));
                        xlabel('<-- NW       X_{PoroTomo} [m]        SE -->');
                        ylabel('Z_{PoroTomo} [m]'); % (= elevation above WGS84 ellipsoid + 800 m'
                    case 3
                        % slicing plane is normal to Z axis
                        title(sprintf('Z_{PoroTomo} = %10.1f meters (%s)',constant_coordinate,strrep(title_str,'_',' ')));
                        xlabel('<-- NW       X_{PoroTomo} [m]        SE -->');
                        ylabel('<-- SW       Y_{PoroTomo} [m]        NE -->');
                    otherwise
                        error(sprintf('unknown knorm %d\n',knorm));
                end
                
                %% set the axes to be the same size as the tomogram
                %       hold on;
                %         axis equal;
                %         axis xy
                %         axis tight;
                %         ax = gca;
                %         ax.Clipping='On';
                %        axis([nanmin(SLICES.Xp), nanmax(SLICES.Xp), nanmin(SLICES.Yp), nanmax(SLICES.Yp)]);
                
                
                % color bar
                h = colorbar;
                %xlabel(h,'Vp [m/s]');
                xlabel(h,OPTIONS.colorbarlabelstr);
                
                % save plot as PNG file with name
                switch knorm
                    case 1
                        fname_out = sprintf('%s_XnormX%05.0fm',title_str,SLICES.Xp(kslice));
                    case 2
                        fname_out = sprintf('%s_YnormY%05.0fm',title_str,SLICES.Yp(kslice));
                    case 3
                        fname_out = sprintf('%s_ZnormZ%05.0fm',title_str,SLICES.Zp(kslice));
                    otherwise
                        error(sprintf('unknown knorm %d\n',knorm));
                end
                feval(funprint,fname_out);
                kfiles = kfiles+1;
                figfilenames{kfiles} = fname_out;
            else
                size(image2dinterp)
                size(M1)
                warning('miscount 1');
                break;
            end
        else
            size(image2d)
            size(colvec(Q1))
            size(colvec(Q2))
            size(colvec(squeeze(image2d)))
            warning('miscount 2');
        end
    end
end
% return each file name as a row
figfilenames=figfilenames';

return
end

