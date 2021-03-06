function figfilenames = plot_tomo_and_faults5(TOMO,FAULTS,title_str,BOUNDS,OPTIONS,funprint...
    ,rsearch, WELLS, elev_mean, BRADY2DGRID, SLICES)
% make slices at BOUNDS of tomogram TOMO with faults in FAULTS
% 20170921 Kurt Feigl
% 20171203 Kurt Feigl - include topographic surface and depth
% 20180211 Kurt Feigl for GRL
% Begin forwarded message:
%
% From: Clifford Thurber <cthurber@wisc.edu>
% Subject: Re: Figures for SRL paper
% Date: February 10, 2018 at 10:26:32 AM CST
% To: Kurt Feigl <feigl@wisc.edu>
% Cc: Lesley Parker <lparker4@wisc.edu>
%
% If you can re-make them as eps and they aren?t awful like the pdfs, please also simplify the labeling:
%
% No plot title showing PoroTomo coordinate info
% No NE/SW/NW/SE labels
% No PoroTomo subscript on axis labels
% No depth on right axis of cross-sections
%
% Also, please no triangle at the bottom of wells on cross-sections
%
% I am going to use WGS84 for elevation for this paper, not PoroTomo Z,
% because I do not want to confuse SRL readers, so please add 800 to the Z
% values in the cross-sections.
% 20180224 Kurt Feigl
%  handle flattened cube case by changing SLICES to BOUNDS

% initialize
kfiles = 0;
figfilenames=[];

nx = size(TOMO.Xp)
ny = size(TOMO.Yp)
nz = size(TOMO.Zp)

tstart = tic;

% set faults
if exist('OPTIONS','var') == 1
    OPTIONS
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
    if isfield(OPTIONS,'draw_wells') == 1
        draw_wells=OPTIONS.draw_wells;
    else
        draw_wells = 1;
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
    if isfield(OPTIONS,'short_labels') == 1
        short_labels=OPTIONS.short_labels;
    else
        short_labels = 0;
    end
    if isfield(OPTIONS,'flattened_cube') == 1
        flattened_cube=OPTIONS.flattened_cube;
    else
        flattened_cube = 0;
    end
    if isfield(OPTIONS,'resolution_threshold') == 1
        resolution_threshold=OPTIONS.resolution_threshold;
    else
        resolution_threshold = 0.0;
    end
    if isfield(OPTIONS,'geologic_model') == 1
        geologic_model=OPTIONS.geologic_model;
    else
        geologic_model = '';
    end
    if isfield(OPTIONS,'draw_topo') == 1
        if isfield(TOMO,'topo') == 1
            draw_topo=OPTIONS.draw_topo;
        else
            warning('TOMO structure does not contain field named topo');
            draw_topo = 0;
        end
    else
        draw_topo = 0;
    end
    if isfield(OPTIONS,'contours') == 1
        contours=OPTIONS.contours;
    else
        contours = 0;
    end
end


%% count number of slices
nxslices = numel(SLICES.Xp)
nyslices = numel(SLICES.Yp)
nzslices = numel(SLICES.Zp)

%% clip out velocity values where resolution is poor
if resolution_threshold > 0.0 && isfield(TOMO,'R') == 1
    title_str = strcat(title_str,sprintf('_MRES_GT_%02d_percent',100*resolution_threshold))
    ilowres = find(TOMO.R < resolution_threshold);
    TOMO.V(ilowres) = NaN;
end

%% clip out velocity values where model voxel is above topographic surface
if resolution_threshold > 0.0 && isfield(TOMO,'R') == 1
    TOMO.topo=interp2(BRADY2DGRID.Xp,BRADY2DGRID.Yp,BRADY2DGRID.Zp,TOMO.Xp,TOMO.Yp);
    iinair = find(TOMO.Zp > TOMO.topo);
    TOMO.V(iinair) = NaN;
end

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
            ikeep = intersect(ikeep,find(TOMO.Xp <= nanmax(colvec(BOUNDS.Xp)) + rsearch));
            ikeep = intersect(ikeep,find(TOMO.Xp >= nanmin(colvec(BOUNDS.Xp)) - rsearch));
        end
        if knorm ~= 2
            ikeep = intersect(ikeep,find(TOMO.Yp <= nanmax(colvec(BOUNDS.Yp)) + rsearch));
            ikeep = intersect(ikeep,find(TOMO.Yp >= nanmin(colvec(BOUNDS.Yp)) - rsearch));
        end
        if knorm ~= 3
            ikeep = intersect(ikeep,find(TOMO.Zp <= nanmax(colvec(BOUNDS.Zp)) + rsearch));
            ikeep = intersect(ikeep,find(TOMO.Zp >= nanmin(colvec(BOUNDS.Zp)) - rsearch));
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
        mnV   = size(TOMO.V)
        fprintf(1,'Starting 3-D scatteredInterpolant at t = %.001f seconds\n',toc(tstart));
        %extrapolation_method = 'none';
        %extrapolation_method = 'linear';
        TOMO
        interpolation_method
        extrapolation_method
        numel(TOMO.V(isfinite(TOMO.V)))
        Finterp2 = scatteredInterpolant(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V ...
            ,interpolation_method,extrapolation_method);
        image2d = Finterp2(Q.x,Q.y,Q.z);
        
        % make a finer mesh with 10-meter spacing
        switch knorm
            case 1
                [M1,M2] = meshgrid([ceil(nanmin(colvec(BOUNDS.Yp))):10:floor(nanmax(colvec(BOUNDS.Yp)))]...
                    ,[ceil(nanmin(colvec(BOUNDS.Zp))):10:floor(nanmax(colvec(BOUNDS.Zp)))]);
            case 2
                [M1,M2] = meshgrid([ceil(nanmin(colvec(BOUNDS.Xp))):10:floor(nanmax(colvec(BOUNDS.Xp)))]...
                    ,[ceil(nanmin(colvec(BOUNDS.Zp))):10:floor(nanmax(colvec(BOUNDS.Zp)))]);
            case 3
                [M1,M2] = meshgrid([ceil(nanmin(colvec(BOUNDS.Xp))):10:floor(nanmax(colvec(BOUNDS.Xp)))]...
                    ,[ceil(nanmin(colvec(BOUNDS.Yp))):10:floor(nanmax(colvec(BOUNDS.Yp)))]);
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
                if flattened_cube == 1
                    switch knorm
                        case 3
                            ipanel = 1; %horizontal slice is in upper left
                            %'position',[left bottom width height])
                            panel_position = [0.1 0.5 0.4 0.4];
                        case 1
                            ipanel = 4; % Y-Z slice is in lower right
                            %'position',[left bottom width height])
                            panel_position = [0.5 0.0 0.4 0.4];
                        case 2
                            ipanel = 3; % X-Z slice is in lower left
                            panel_position = [0.1 0.0 0.4 0.4];
                        otherwise
                            error(sprintf('unknown knorm = %d\n',knorm));
                    end
                    
                    %subplot(2,2,ipanel);
                    subplot('position',panel_position);
                else
                    figure;
                end
                hold on;
                
                % setting the velocity range for plotting
                %caxis([v1, v2]);
                if isfinite(vmin) == 1 && isfinite(vmax)==1
                    caxis([vmin,vmax]);
                    caxis manual;
                else
                    vmin = nanmin(colvec(image2dinterp));
                    vmax = nanmax(colvec(image2dinterp));
                    if isfinite(vmin) == 1 && isfinite(vmax)==1
                        caxis([vmin,vmax]);
                        caxis manual;
                    else
                        caxis auto;
                    end
                end
                
                %colormap(flipud(colormap('jet')));
                colormap(cmap);
                
                % plot in color
                size(image2dinterp)
                iok = find(isfinite(image2dinterp)==1);
                nok = numel(iok)
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
                        axis([nanmin(BOUNDS.Yp), nanmax(BOUNDS.Yp), nanmin(BOUNDS.Zp), nanmax(BOUNDS.Zp)]);
                    case 2
                        axis([nanmin(BOUNDS.Xp), nanmax(BOUNDS.Xp), nanmin(BOUNDS.Zp), nanmax(BOUNDS.Zp)]);
                    case 3
                        axis([nanmin(BOUNDS.Xp), nanmax(BOUNDS.Xp), nanmin(BOUNDS.Yp), nanmax(BOUNDS.Yp)]);
                    otherwise
                        error(sprintf('unknown knorm = %d\n',knorm));
                end               
                % clip outside the axes
                ax = gca;
                ax.Clipping='On';
                              
                
                %% draw contours
                if contours == 1
                    %v_contours = linspace(vmin,vmax,3);
                    v_contours = nice_steps([vmin,vmax],5);
                    [c,h]=contour(M1,M2,image2dinterp,v_contours,'w');
                    if short_labels == 0
                       clabel(c,h,'LabelSpacing',8*72,'Color','White'); % spacing in points
                    %[c,h]=contour(M1,M2,image2dinterp,v_contours,'w','ShowText','on');
                    end
                end
                               
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
                            
                %% draw the faults              
                ksor  = 2; % sort according to Y coordinate
                draw_faults3(FAULTS,knorm,constant_coordinate,nmin,norder,plot_points,plot_curves,ksor);
                
                
                %% draw the topographic surface, i.e. the interface between rock and air
                if draw_topo == 1
                     switch knorm
                        case 1
                            for i=1:numel(TOMO.Xp)
                                if abs(TOMO.Xp(i)-SLICES.Xp(kslice)) < 1*rsearch
                                    plot(TOMO.Yp(i),TOMO.topo(i),'o-','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0.5]);
                                end
                            end
                        case 2
                            for i=1:numel(TOMO.Yp)
                                if abs(TOMO.Yp(i)-SLICES.Yp(kslice)) < 1*rsearch
                                    plot(TOMO.Xp(i),TOMO.topo(i),'o-','MarkerSize',7,'MarkerFaceColor',[0.5,0.5,0.5]);
                                end
                            end
                        case 3
                            fprintf('Not plotting topography on map view.\n');
                        otherwise
                            error(sprintf('unknown knorm = %d\n',knorm));
                    end
                end
                
                
                %% draw the wells
                if draw_wells == 1
                    nwells = numel(WELLS.WellName);
                    well_symbols = cell(nwells,1);
                    for i=1:nwells
                        if contains(WELLS.ActivityDuringPoroTomoTesting(i),'Inj') == 1
                            well_symbols{i} = 'kv';  % downward pointing triangle for injection well
                        elseif contains(WELLS.ActivityDuringPoroTomoTesting(i),'Pumping') == 1
                            well_symbols{i} = 'k^';  % upward pointing triangle for injection well
                        elseif contains(WELLS.ActivityDuringPoroTomoTesting(i),'Obs') == 1
                            well_symbols{i} = 'ko';  % open circle for observation well
                        else
                            well_symbols{i} = 'k+';  % 
                        end
                    end
                    switch knorm
                        case 1
                            for i=1:nwells
                                if abs(WELLS.grdsurf.Xp(i)-SLICES.Xp(kslice)) < rsearch/2.
                                    plot(WELLS.grdsurf.Yp(i),WELLS.grdsurf.Zp(i),well_symbols{i},'MarkerSize',7,'MarkerFaceColor','none'); % symbol
                                    plot([WELLS.grdsurf.Yp(i),WELLS.openbtm.Yp(i)]...
                                        ,[WELLS.grdsurf.Zp(i),min([WELLS.openbtm.Zp(i),WELLS.perfbtm.Zp(i)])]...
                                        ,'k-','LineWidth',1); % thin, solid line for borehole    
                                    plot([WELLS.perftop.Yp(i),WELLS.perfbtm.Yp(i)]...
                                        ,[WELLS.perftop.Zp(i),WELLS.perfbtm.Zp(i)]...
                                        ,'k:','LineWidth',2); % thin, solid line for borehole 
                                    plot([WELLS.opentop.Yp(i),WELLS.openbtm.Yp(i)]...
                                        ,[WELLS.opentop.Zp(i),WELLS.openbtm.Zp(i)]...
                                        ,'k:','LineWidth',4); % thick dashed line for open interval
                                end
                            end
                        case 2
                            for i=1:nwells
                                if abs(WELLS.grdsurf.Yp(i)-SLICES.Yp(kslice)) < rsearch/2.
                                    plot(WELLS.grdsurf.Xp(i),WELLS.grdsurf.Zp(i),well_symbols{i},'MarkerSize',7,'MarkerFaceColor','none'); % symbol
                                    plot([WELLS.grdsurf.Xp(i),WELLS.openbtm.Xp(i)]...
                                        ,[WELLS.grdsurf.Zp(i),min([WELLS.openbtm.Zp(i),WELLS.perfbtm.Zp(i)])]...
                                        ,'k-','LineWidth',1); % thin, solid line for borehole                                    
                                    plot([WELLS.perftop.Xp(i),WELLS.perfbtm.Xp(i)]...
                                        ,[WELLS.perftop.Zp(i),WELLS.perfbtm.Zp(i)]...
                                        ,'k:','LineWidth',2); % thin, dashed line for perforated interval
                                    plot([WELLS.opentop.Xp(i),WELLS.openbtm.Xp(i)]...
                                        ,[WELLS.opentop.Zp(i),WELLS.openbtm.Zp(i)]...
                                        ,'k:','LineWidth',4); % thick dashed line for open interval
                                end
                            end
                        case 3
                            for i=1:nwells
                                plot(WELLS.grdsurf.Xp(i),WELLS.grdsurf.Yp(i),well_symbols{i},'MarkerSize',5,'MarkerFaceColor','none');
                            end
                        otherwise
                            error(sprintf('unknown knorm = %d\n',knorm));
                    end
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
                
                %% title
                if flattened_cube ~= 1
                    switch knorm
                        case 1
                            title(sprintf('%s\nX_{PoroTomo} = %10.1f meters'...
                                ,strrep(title_str,'_',' ') ...
                                ,constant_coordinate));
                        case 2
                            title(sprintf('%s\nY_{PoroTomo} = %10.1f meters'...
                                ,strrep(title_str,'_',' ') ...
                                ,constant_coordinate));
                        case 3
                            title(sprintf('%s\nZ_{PoroTomo} = %7.0f m; depth = %7.0f m'...
                                ,strrep(title_str,'_',' ')...
                                ,constant_coordinate,elev_mean-constant_coordinate-800));
                        otherwise
                            error(sprintf('unknown knorm %d\n',knorm));
                    end
                end
                if short_labels == 1
                    switch knorm
                        case 1
                            title(sprintf('X = %10.1f m',constant_coordinate));
                        case 2
                            title(sprintf('Y = %10.1f m',constant_coordinate));
                            ax = gca;ax.Clipping='Off';
                            text(nanmin(BOUNDS.Xp),nanmax(BOUNDS.Zp),'SW','HorizontalAlignment','left','VerticalAlignment','Bottom');
                            text(nanmax(BOUNDS.Xp),nanmax(BOUNDS.Zp),'NE','HorizontalAlignment','right','VerticalAlignment','Bottom');
                            ax = gca;ax.Clipping='On';
                        case 3
                            title(sprintf('Z = %7.0f m',constant_coordinate+800));
                        otherwise
                            error(sprintf('unknown knorm %d\n',knorm));
                    end
                else
                    switch knorm
                        case 1
                            title(sprintf('%s\nX_{PoroTomo} = %10.1f meters',strrep(title_str,'_',' '),constant_coordinate));
                        case 2
                            title(sprintf('%s\nY_{PoroTomo} = %10.1f meters',strrep(title_str,'_',' '),constant_coordinate));
                        case 3
                            title(sprintf('%s\nZ_{PoroTomo} = %7.0f m; depth = %7.0f m'...
                                ,strrep(title_str,'_',' ') ...
                                ,constant_coordinate,elev_mean-constant_coordinate-800));
                        otherwise
                            error(sprintf('unknown knorm %d\n',knorm));
                    end
                end
                
                %% label axes
                switch knorm
                    case 1 % slicing plane is normal to X axis
                        if short_labels == 1
                            xlabel('Y [m]');
                            %ylabel('Z [m]'); % (= elevation above WGS84 ellipsoid - 800 m'
                            set(gca,'YTickLabel','');
                        else
                            xlabel('<-- SW       Y_{PoroTomo} [m]        NE -->');
                            ylabel('Z_{PoroTomo} [m]'); % (= elevation above WGS84 ellipsoid + 800 m'
                        end
                    case 2 % slicing plane is normal to Y axis
                        if short_labels == 1
                            xlabel('X [m]');
                            %ylabel('Z [m]'); % (= elevation above WGS84 ellipsoid - 800 m'
                            set(gca,'YTickLabel','');
                        else
                            xlabel(',<-- NW      X_{PoroTomo} [m]        SE -->');
                            ylabel('Z_{PoroTomo} [m]'); % (= elevation above WGS84 ellipsoid - 800 m'
                        end
                    case 3 % slicing plane is normal to Z axis
                        if short_labels == 1
                            xlabel('X [m]');
                            ylabel('Y [m]');
                        else
                            %  label with depth, too
                            xlabel('<-- NW       X_{PoroTomo} [m]        SE -->');
                            ylabel('<-- SW       Y_{PoroTomo} [m]        NE -->');
                        end
                    otherwise
                        error(sprintf('unknown knorm %d\n',knorm));
                end
                
                %% depth on second axis on right
                %                 switch knorm
                %                     case 1
                %                         %axis([nanmin(SLICES.Yp), nanmax(SLICES.Yp), nanmin(SLICES.Zp), nanmax(SLICES.Zp)]);
                %                         %                         y2 = elev_mean - y1 - 800;
                %                         yyaxis right
                %                         axis ij
                %                         plot([nanmin(SLICES.Yp), nanmax(SLICES.Yp)],elev_mean-[nanmin(SLICES.Zp), nanmax(SLICES.Zp)]-800,'w.');
                %                         ylabel('Depth [m]','color','k');
                %                         ax = gca;
                %                         ax.YAxis(1).Color='k';
                %                         ax.YAxis(2).Color='k';
                %                     case 2
                %                         %%axis([nanmin(SLICES.Xp), nanmax(SLICES.Xp), nanmin(SLICES.Zp), nanmax(SLICES.Zp)]);
                %                         yyaxis right
                %                         axis ij
                %                         plot([nanmin(SLICES.Xp), nanmax(SLICES.Xp)],elev_mean-[nanmin(SLICES.Zp), nanmax(SLICES.Zp)]-800,'w.');
                %                         ylabel('Depth [m]','color','k');
                %                         ax = gca;
                %                         ax.YAxis(1).Color='k';
                %                         ax.YAxis(2).Color='k';
                %                     case 3
                %                         %axis([nanmin(SLICES.Xp), nanmax(SLICES.Xp), nanmin(SLICES.Yp), nanmax(SLICES.Yp)]);
                %                     otherwise
                %                         error(sprintf('unknown knorm = %d\n',knorm));
                %                 end
                
                %% WGS elevation second axis on right
                if flattened_cube ~= 1
                    switch knorm
                        case 1
                            yyaxis right
                            axis xy                           
                            plot([nanmin(BOUNDS.Yp), nanmax(BOUNDS.Yp)],[nanmin(BOUNDS.Zp)+800, nanmax(BOUNDS.Zp)+800],'w.');
                            if short_labels == 1
                                ylabel('Z [m]','color','k');
                            else
                                ylabel('Elevation [m]','color','k');
                            end
                            ax = gca;
                            ax.YAxis(1).Color='k';
                            ax.YAxis(2).Color='k';
                        case 2
                            yyaxis right
                            axis xy
                            plot([nanmin(BOUNDS.Xp), nanmax(BOUNDS.Xp)],[nanmin(BOUNDS.Zp)+800, nanmax(BOUNDS.Zp)+800],'w.');
                            if short_labels == 1
                                ylabel('Z [m]','color','k');
                            else
                                ylabel('Elevation [m]','color','k');
                            end
                            ax = gca;
                            ax.YAxis(1).Color='k';
                            ax.YAxis(2).Color='k';
                        case 3
                            %axis([nanmin(SLICES.Xp), nanmax(SLICES.Xp), nanmin(SLICES.Yp), nanmax(SLICES.Yp)]);
                        otherwise
                            error(sprintf('unknown knorm = %d\n',knorm));
                    end
                end
                
                
                
                %% set the axes to be the same size as the tomogram
                %       hold on;
                %         axis equal;
                %         axis xy
                %         axis tight;
                %         ax = gca;
                %         ax.Clipping='On';
                %        axis([nanmin(BOUNDS.Xp), nanmax(BOUNDS.Xp), nanmin(BOUNDS.Yp), nanmax(BOUNDS.Yp)]);
                
                if flattened_cube == 1
                    subplot(2,2,2); % color bar in upper right
                end
                %% color bar
                h = colorbar;
                %xlabel(h,'Vp [m/s]');
                xlabel(h,OPTIONS.colorbarlabelstr);
                
%                 %% write contour levels on legend 
%                 if contours == 1 && short_labels == 0
%                     legend(sprintf('%4.1g\n',v_contours),'Location','bestoutside');
%                 end
                
                %% save plot as graphics file with name
                if flattened_cube == 1
                    fname_out = sprintf('%s_cubeX%05.0fmY%05.0fmZ%05.0fm',title_str...
                        ,SLICES.Xp(kslice),SLICES.Yp(kslice),SLICES.Zp(kslice));
                else
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
                end
                % append name of geologic model to file name
                fname_out = strcat(fname_out,'_',geologic_model);
                
                % call the appropriate printing function, which will append
                % file type (e.g., '.pdf')
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

