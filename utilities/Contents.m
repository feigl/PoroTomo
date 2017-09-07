%% UTILITIES
%  20170412 Elena Reinisch and Kurt Feigl 
%  
%
% Files
%   csv2struct            - function [DATA, DATA_fields] = csv2struct( csvname )
%   csv2struct2           - function [DATA, DATA_fields] = csv2struct( csvname )
%   csvimport             - reads the specified CSV file and stores the contents in a cell array or matrix
%   deg2utm               - project geographic coordinates (latitude, longitude in degrees) into cartographic coordinates (easting, northing in meters) using Universal Transverse Mercator (UTM) projection, assuming WGS84 ellipsoid.
%   fit_straight_line     - function [pest, psig, tfit, ymod, ymodl, ymodu, mse] = fit_straight_line(time,yobs,ysig)
%   fixticklabels         - function fixticklabels(xfmt,yfmt)
%   GDRcsv2struct         - [ DATA ] = GDRcsv2struct( gdr_csv_url )
%   get_das_utctime       - this version is show UTC time
%   PoroTomo_Get_Metadata - Get current version of metadata files from askja using anynomous ftp ftp://roftp.ssec.wisc.edu
%   printpdf              - write current graphics window to a PDF file
%   Read_PoroTomo_urls    - list of files available on GDR
%   rotate_2d             - rotate_2d: Code for rotating a set of points in a 2D plane
%   timetag_demo          - demonstrate time tags in Matlab
%   utm2deg               - inverse project  cartographic coordinates (UTM easting, northing in meters) to geographic coordinates (latitude, longitude in degrees) using inverse Universal Transverse Mercator (UTM) projection, assuming WGS84 ellipsoid.
%   utm2xy_porotomo       - function [x_rotated_in_meters, y_rotated_in_meters] = utm2xy_porotomo(UTM_easting_in_meters,UTM_northing_in_meters)
%   utm2xyz_porotomo      - function [x_rotated_in_meters, y_rotated_in_meters, z_rotated_in_meters] = utm2xyz_porotomo(UTM_easting_in_meters,UTM_northing_in_meters, UTM_height_in_meters)
%   xlsx2struct           - function [ DATA ] = xlsx2struct( xlsx_filename_list )
%   xy_porotomo2utm       - function [UTM_easting_in_meters,UTM_northing_in_meters] = xy_porotomo2utm(x_rotated_in_meters, y_rotated_in_meters)
%   xyz_porotomo2utm      - function [UTM_easting_in_meters,UTM_northing_in_meters, UTM_height_in_meters] = xyz_porotomo2utm(x_rotated_in_meters, y_rotated_in_meters, z_rotated_in_meters)
%   yymmdd2yr_mo_dy       - function [year,month,day] = yymmdd2yr_mo_dy(yymmdd)
