function figfilenames = plot_tomo_and_faults(TOMO,FAULTS,title_str,SLICES,OPTIONS,funprint)
% make slices at SLICES of tomogram TOMO with faults in FAULTS
% 20170921 Kurt Feigl

% initialize
kfiles = 0;

nx = size(TOMO.Xp)
ny = size(TOMO.Yp)
nz = size(TOMO.Zp)

% set faults
if exist('OPTIONS','var') == 1
    if isfield(OPTIONS,'vmin') == 1
        vmin=OPTIONS.vmin;
    else
        vmin=nanmin(colvec(TOMO.V));
    end
    if isfield(OPTIONS,'vmax') == 1
        vmax=OPTIONS.vmax;
    else
        vmax=nanmax(colvec(TOMO.V));
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
        % make a 2-D mesh
        switch knorm
            case 1
                [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp(kslice),  SLICES.Yp,        SLICES.Zp);
            case 2
                [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp,          SLICES.Yp(kslice),SLICES.Zp);
            case 3
                [Q.x,Q.y,Q.z] = meshgrid(SLICES.Xp,SLICES.Yp,                  SLICES.Zp(kslice));
            otherwise
                error(sprintf('unknown knorm = %d\n',knorm));
        end
        % interpolate values from 3D onto the 2-D mesh
        image2d = griddata(TOMO.Xp,TOMO.Yp,TOMO.Zp,TOMO.V,Q.x,Q.y,Q.z,interpolation_method);
        
        % make a finer mesh with 10-meter spacing
        switch knorm
            case 1
                [M1,M2] = meshgrid([unique(min(SLICES.Yp)):10:unique(max(SLICES.Yp))]...
                    ,[unique(min(SLICES.Zp)):10:unique(max(SLICES.Zp))]);
            case 2
                [M1,M2] = meshgrid([unique(min(SLICES.Xp)):10:unique(max(SLICES.Xp))]...
                    ,[unique(min(SLICES.Zp)):10:unique(max(SLICES.Zp))]);
            case 3
                [M1,M2] = meshgrid([unique(min(SLICES.Xp)):10:unique(max(SLICES.Xp))]...
                    ,[unique(min(SLICES.Yp)):10:unique(max(SLICES.Yp))]);
            otherwise
                error(sprintf('unknown knorm = %d\n',knorm));
        end
        
        % interpolate onto that mesh
        %         switch knorm
        %             case 1
        %                 image2dinterp = interp2(squeeze(Q.y),squeeze(Q.z),image2d,M1,M2);
        %             case 2
        %                 image2dinterp = interp2(squeeze(Q.x),squeeze(Q.z),image2d,M1,M2);
        %             case 3
        %                 image2dinterp = interp2(squeeze(Q.x),squeeze(Q.y),image2d,M1,M2);
        %             otherwise
        %                 error(sprintf('unknown knorm = %d\n',knorm));
        %         end
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
        
        if numel(colvec(Q1)) > 10 && numel(colvec(Q2)) > 10
            
            % interpolate onto 2-D mesh
            %image2dinterp = interp2(Q1,Q2,squeeze(image2d),M1,M2);
            Finterpolant = scatteredInterpolant(colvec(Q1),colvec(Q2),colvec(squeeze(image2d)));
            %whos Finterpolant
            image2dinterp = Finterpolant(colvec(M1),colvec(M2));
            image2dinterp = reshape(image2dinterp,size(M1));
            
            % plot the tomogram
            figure;
            hold on;
            
            % setting the velocity range for plotting
            %caxis([v1, v2]);
            caxis([vmin,vmax]);
            caxis manual;
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
            
            % label axes
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
                    fname_out = sprintf('%s_XnormX%05.0fm.png',title_str,SLICES.Xp(kslice));
                case 2
                    fname_out = sprintf('%s_YnormY%05.0fm.png',title_str,SLICES.Yp(kslice));
                case 3
                    fname_out = sprintf('%s_ZnormZ%05.0fm.png',title_str,SLICES.Zp(kslice));
                otherwise
                    error(sprintf('unknown knorm %d\n',knorm));
            end
            feval(funprint,fname_out);
            kfiles = kfiles+1;
            figfilenames{kfiles} = fname_out;
        end
    end   
end
% return each file name as a row
figfilenames=figfilenames';

return
end

