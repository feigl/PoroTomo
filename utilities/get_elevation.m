function elevation = get_elevation(xq,yq,grdfilename)
% extract elevation from a digital elevation model
% 20170925 Kurt Feigl

% check for existence
if exist(grdfilename,'file') ~= 2
    error(sprintf('Could not find GMT grid file named %s\n',grdfilename));
else    
    % read the grdfile in GMT format
    INFO = grdinfo3(grdfilename)
    xlab = INFO.xname;
    ylab = INFO.yname;
    zlab = INFO.zname;

    [xgrd,ygrd,zgrd] = grdread3(grdfilename);
    zgrd = double(zgrd);
    
    % find mean values
    xmean = nanmean(xgrd);
    ymean = nanmean(ygrd);
    zmean = nanmean(colvec(zgrd));
       
    %% if coordinates are in UTM meters
    if contains(INFO.xname,'meters') == 1 || contains(INFO.yname,'meters') == 1
        fprintf(1,'Coordinates are in meters\n');
    else
        [xutmmean,yutmmean,utmzone] = deg2utm(ymean,xmean);
    end
    fprintf(1,'Mean Coordinates: X = %.3f %s Y = %.3f %s Z = %.3f %s\n',xmean,xlab,ymean,ylab,zmean,zlab);
  
    
    % build a regular mesh
    [XGRD,YGRD] = meshgrid(xgrd,ygrd);
    
    if numel(XGRD) == numel(zgrd)
        ZGRD = reshape(zgrd,size(XGRD));
    else
        error('miscount');
    end
   
    
%     size(XGRD)
%     size(YGRD)
%     size(ZGRD)
    % build the 2-D interpolant function
    Finterpolant = scatteredInterpolant(colvec(XGRD),colvec(YGRD),colvec(ZGRD));
    Finterpolant.Method='linear';
    Finterpolant.ExtrapolationMethod ='none';
 
    
    % evaluate the interpolant function at the query points
    elevation = Finterpolant(colvec(xq),colvec(yq));
    return    
end

