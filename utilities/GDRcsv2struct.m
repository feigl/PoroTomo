function [ METADATA ] = GDRcsv2struct( gdr_csv_url, data_type_short_name )
%[ DATA ] = GDRcsv2struct( gdr_csv_url )
%  Given a url to a csv file on GDR, download the file and save the
%  contents to a data structure
%
% INPUTS:
% gdr_csv_url - string (or vector of strings) containing the (full) url to the GDR submission. To
% get this url, just right click on the GDR file link and select "copy
% address"
% data_type_short_name - optional string for short name of data.  If not specified, GDRcsv2struct will try to select the appropriate short name
% based on the filename. 
%
% Current labels are:
% KEYWORD IN FILENAME      DATA SHORT NAME
%   (NOT Case Sens.)
% ------------------------------------------------------------
%  borehole + das --------> DASV
%  surface + das  --------> DASH
%  borehole --------------> BORE
%  nodal    --------------> NODE
%  reftek   --------------> REFT
%  vibroseis -------------> VIBRO
%  permaseis -------------> PERMS
%  insar    --------------> INSAR
%  pressure  -------------> PDAT
%  geol     --------------> GEOL
%  pump     --------------> PUMP
%  temp + dtsv -----------> DTSV
%  temp + dtsh -----------> DTSH
%  dtsv        -----------> DTSV
%  dtsh        -----------> DTSH
%  uav         -----------> UAV
%
% If GDRcsv2struct cannot find a match to any of the above, a prompt will appear that requires
% manual entry.
%
% OUTPUT:
% METADATA - data structure containing the full contents of the csv
% file(s).  If reading more than one url at a time, METADATA will have
% sub-strutures that contain the information for each csv file.
% Data is additionally saved to a matlab file with date of run
%
% Elena C. Reinisch, 20160827

 % check to see if data short name is supplied
   if nargin == 2
       isshortname = 1;
   elseif nargin == 1
       isshortname = 0;
   else
       error('wrong number of input variables');
   end

% build a structure
METADATA = struct;
class(gdr_csv_url)
size(gdr_csv_url)
gdr_csv_url = cellstr(gdr_csv_url)
for k = 1:numel(gdr_csv_url)
    % assign current url address to variable
   furl = gdr_csv_url{k}
   % pull name of csv file from url to save file as
   furl_split = regexp(furl, '/', 'split');
   csvname = furl_split{end}
   % check for % characters
   csv_split = regexp(csvname, '%', 'split')
   if numel(csv_split) > 1
    csvname = strcat(csv_split{1}, '.csv')
   end
       
   % assign shortname if needed
   if isshortname == 0
      if find(regexpi(csvname, 'borehole')) == 1
          if find(regexpi(csvname, 'das')) == 1
              data_short_name = 'DASV';
          else
              data_short_name = 'BORE';
          end
      elseif find(regexpi(csvname, 'das')) == 1
             data_short_name = 'DASH';
      elseif find(regexpi(csvname, 'nodal')) == 1
          data_short_name = 'NODE';
      elseif find(regexpi(csvname, 'reftek')) == 1
          data_short_name = 'REFT';
      elseif find(regexpi(csvname, 'vibroseis')) == 1
          data_short_name = 'VIBRO';
      elseif find(regexpi(csvname, 'permaseis')) == 1
          data_short_name = 'PERMS';
      elseif find(regexpi(csvname, 'insar')) == 1
          data_short_name = 'INSAR';
      elseif find(regexpi(csvname, 'pressure')) == 1
          data_short_name = 'PDAT';
      elseif find(regexpi(csvname, 'geol')) == 1
          data_short_name = 'GEOL';
      elseif find(regexpi(csvname, 'pump')) == 1
          data_short_name = 'PUMP';
      elseif find(regexpi(csvname, 'temp')) == 1
          if find(regexpi(csvname, 'dtsv')) == 1
            data_short_name = 'DTSV';  
          elseif find(regexpi(csvname, 'dtsh')) == 1
            data_short_name = 'DTSH'; 
          end
      elseif find(regexpi(csvname, 'dtsv')) == 1
          data_short_name = 'DTSV';    
      elseif find(regexpi(csvname, 'dtsh')) == 1
          data_short_name = 'DTSH';       
      elseif find(regexpi(csvname, 'uav')) == 1
          data_short_name = 'UAV';  
      else
          data_short_name = input('Manual entry of data short name needed:', 's')
      end
   else
       data_short_name = data_type_short_name{k};
   end
   data_short_name
   
   % download file to current directory under csvname
   urlwrite(furl, csvname);
   
   % read csv file and write to data structure
   METADATA.(data_short_name) = csv2struct(csvname);
    
end

save(strcat('METADATA_', date, '.mat'), 'METADATA')

end

