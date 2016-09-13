function [ METADATA ] = xlsx2struct(xlsx_filename_list, data_type_short_name)
%function [ DATA ] = xlsx2struct( xlsx_filename_list )
%  Given the filename for Excel files under Metdata folder on Askja:/data/PoroTomo, download the file using anynomous ftp and save the
%  contents to a data structure
%
% INPUTS:
% xlsx_filename_list - name of text file which contains name of xlsx files in column format
% data_type_short_name - optional string for short name of data.  If not specified, GDRcsv2struct will try to select the appropriate short name
% based on the filename. 
%
% Current labels are:
% KEYWORD IN FILENAME      DATA SHORT NAME
%   (NOT Case Sens.)
% ------------------------------------------------------------
%  borehole + das + dts --> DASV_DTS
%  surface + das + dts  --> DASH_DTS
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
% If xlsx2struct cannot find a match to any of the above, a prompt will appear that requires
% manual entry.
%
% OUTPUT:
% METADATA - data structure containing the full contents of the csv
% file(s).  If reading more than one url at a time, METADATA will have
% sub-strutures that contain the information for each csv file.
% Data is additionally saved to a matlab file with date of run
%   
% Elena C. Reinisch 20160913

% check to see if data short name is supplied
   if nargin == 1
       isshortname = 0;
   elseif nargin == 2
       isshortname = 1;
   else
       error('wrong number of input variables');
   end

% build a structure
METADATA = struct;



% read file names to text string
xlsx_filenames = textread(xlsx_filename_list, '%s');

% open ftp connection
  ftpsite = ftp('roftp.ssec.wisc.edu');
  
% loop through addresses and save all data to matlab structure METADATA as
% well as .mat file

for k = 1:numel(xlsx_filenames)
    xlsx_name = xlsx_filenames{k}
    xlsx_path = strcat('porotomo/PoroTomo/Metadata/', xlsx_name);
   
   % assign shortname if needed
   if isshortname == 0
%       if find(regexpi(xlsx_name, 'borehole')) == 1
%           if find(regexpi(xlsx_name, 'das')) == 1
%               data_short_name = 'DASV';
%           else
%               data_short_name = 'BORE';
%           end
      if find(regexpi(xlsx_name, 'nodal')) == 1
          data_short_name = 'NODE';
      elseif find(regexpi(xlsx_name, 'reftek')) == 1
          data_short_name = 'REFT';
      elseif find(regexpi(xlsx_name, 'vibroseis')) == 1
          data_short_name = 'VIBRO';
      elseif find(regexpi(xlsx_name, 'permaseis')) == 1
          data_short_name = 'PERMS';
      elseif find(regexpi(xlsx_name, 'insar')) == 1
          data_short_name = 'INSAR';
      elseif find(regexpi(xlsx_name, 'pressure')) == 1
          data_short_name = 'PDAT';
      elseif find(regexpi(xlsx_name, 'geol')) == 1
          data_short_name = 'GEOL';
      elseif find(regexpi(xlsx_name, 'pump')) == 1
          data_short_name = 'PUMP';
      elseif find(regexpi(xlsx_name, 'temp')) == 1
          if find(regexpi(xlsx_name, 'dtsv')) == 1
            data_short_name = 'DTSV';  
          elseif find(regexpi(xlsx_name, 'dtsh')) == 1
            data_short_name = 'DTSH'; 
          end
      elseif find(regexpi(xlsx_name, 'dtsv')) == 1
          data_short_name = 'DTSV'; 
      elseif find(regexpi(xlsx_name, 'surface')) == 1 & find(regexpi(xlsx_name, 'DAS')) == 1 & find(regexpi(xlsx_name, 'DTS')) == 1 
          data_short_name = 'DASH_DTS'; 
      elseif find(regexpi(xlsx_name, 'borehole')) == 1 & find(regexpi(xlsx_name, 'DAS')) == 1 & find(regexpi(xlsx_name, 'DTS')) == 1 
          data_short_name = 'DASV_DTS'; 
      elseif find(regexpi(xlsx_name, 'dtsh')) == 1
          data_short_name = 'DTSH';       
      elseif find(regexpi(xlsx_name, 'uav')) == 1
          data_short_name = 'UAV';  
      else
          data_short_name = input('Manual entry of data short name needed:', 's')
      end
   else
       data_short_name = data_type_short_name{k};
   end
   
   % pull file from ftp server 
     mget(ftpsite, xlsx_path);

    % load data from file to table
    data = readtable(xlsx_path);
    
    % convert data to structure 
    METADATA.(data_short_name) = table2struct(data,'ToScalar',true);

end

% close ftp connection
  close(ftpsite);
  
% save metadata to .mat file with date in name
 save(strcat('METADATA_', date, '.mat'), 'METADATA')
  
return

