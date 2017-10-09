function WELLS = read_wells(isite)
% return 3-coordinates of study area
% 20170925 Kurt Feigl

%% read the wells

switch isite
    case 'brady'
        % url = 'https://gdr.openei.org/submissions/828'
        % data_filename='Brady_Well_coordinates_Metadata.csv'
        % strcat(url,filesep,data_filename)
        % websave(data_filename,strcat(url,filesep,data_filename));
        %'Brady_Well_coordinates_Metadata.csv'
        dirname=strcat(getenv('HOME'),filesep,'PoroTomo',filesep,'metadata_txt_files')
        data_filename='Brady_Well_coordinates_Metadata.csv'
        table=readtable(strcat(dirname,filesep,data_filename));
        WELLS=table2struct(table,'ToScalar',true);
        % Transform coordinates into rotated PoroTomo frame (Xp,Yp,Zp)
        [WELLS.grdsurf.Xp,WELLS.grdsurf.Yp,WELLS.grdsurf.Zp] = utm2xyz_porotomo(WELLS.XCoordinate_UTMME_...
            ,WELLS.YCoordinate_UTMMN_...
            ,WELLS.LandSurfElev_mAboveWGS84Ellipsoid_);
        [WELLS.perftop.Xp,WELLS.perftop.Yp,WELLS.perftop.Zp] = utm2xyz_porotomo(WELLS.XCoordinate_UTMME_...
            ,WELLS.YCoordinate_UTMMN_...
            ,WELLS.PerfDepthTopElev_m_);
        [WELLS.perfbtm.Xp,WELLS.perfbtm.Yp,WELLS.perfbtm.Zp] = utm2xyz_porotomo(WELLS.XCoordinate_UTMME_...
            ,WELLS.YCoordinate_UTMMN_...
            ,WELLS.PerfDepthBtmElev_m_);
        [WELLS.opentop.Xp,WELLS.opentop.Yp,WELLS.opentop.Zp] = utm2xyz_porotomo(WELLS.XCoordinate_UTMME_...
            ,WELLS.YCoordinate_UTMMN_...
            ,WELLS.OpenHoleTopElev_m_);
        [WELLS.openbtm.Xp,WELLS.openbtm.Yp,WELLS.openbtm.Zp] = utm2xyz_porotomo(WELLS.XCoordinate_UTMME_...
            ,WELLS.YCoordinate_UTMMN_...
            ,WELLS.OpenHoleBtmEleve_m_);
    otherwise
        error(sprintf('Unknown isite = %d\n',isite));
end


return
end



