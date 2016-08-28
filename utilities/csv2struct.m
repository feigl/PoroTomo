function DATA = csv2struct(csvname)
%function [DATA, DATA_fields] = csv2struct( csvname )
%   reads the header of a csv file and outputs the labels as strings.  To
%   be used with csvimport.m
%
% INPUT:
% csvname - string that contains name of csv data file (must be comma delimited)
% 
% OUTPUTS:
% DATA - data structure containing the data columns from the csv file
%              in field format, with labels from columns (modified if needed)
%
% NOTE: csvimport is sensitive to the format of the csv file. if a warning is thrown, 
% the data may not have been stored correctly in the DATA structure.
% Take special consideration with the warning that certain lines appear to
% have different numbers of columns than the rest.  This is indicative of
% commas with the "comment/notes" column of the csv file and need to be
% removed before the file can be read properly
%
% Elena C. Reinisch
% last update: 20160827 - streamline header extraction process and account
% for extra text lines

% OLD CODE
% get header information for csvimport
%vals = importdata(csvname); % pull text into cells by row (1 row per cell)

% check to see if data has been imported as data structure
% if strcmp(class(vals), 'struct') == 1
%     fld_list = fieldnames(vals);
%     % check to see if colheaders field exists
%     TF = strcmpi('colheaders', fld_list);
%     if isempty(find(TF)) == 0
%         column_labels = vals.colheaders;
%     else
%         text = vals.textdata;
%         column_labels = text(1,:);
%     end
% else %read as text file
%     column_labels = strsplit(char(vals(1)), ','); % separate header row into strings by commas
% end

% MOST RECENT
% get header information for csvimport
fid = fopen(csvname);
column_labels = fgetl(fid);
column_labels = regexp(column_labels, ',', 'split' )
% read data from file 
[data_matrix] = csvimport(csvname, 'columns', column_labels);

[nr,nc] = size(data_matrix)

% check to see that file was read correctly (i.e., no osx-unix
% conversions; if error, try removing <cntrl>M and replacing with new line
% characters
if nr <= 0 
    % write and execute command to replace Mac OS ^M with new line character; saves to
    % new file 
    command = sprintf('awk ''{ gsub("\\r", "\\n"); print $0;}'' %s > %s', csvname, strrep(csvname, '.csv', 'EOL.csv'));
    unix(command);
    
    % call new file for reading
    csvname = strrep(csvname, '.csv', 'EOL.csv');  
    
    %re-run csvimport on new file
    [data_matrix] = csvimport(csvname, 'columns', column_labels);
    
    % check if read; if still empty, throw error for further manual
    % inspection
    [nr,nc] = size(data_matrix)
    if nr <= 0 
        error(sprintf('file named %s appears to empty. Consider CRLF issues.\n',csvname));
    end
end
if strcmp(class(data_matrix), 'double') == 1
    data_matrix = num2cell(data_matrix);
end

% make sure column labels are appropriate field headers
% removes spaces, right parentheses/brackets, replaces left
% parentheses/brackets with underscore
column_labels = regexprep(column_labels, ' ', '');
column_labels = regexprep(column_labels, '(', '_');
column_labels = regexprep(column_labels, ')', '');
column_labels = regexprep(column_labels, '{', '_');
column_labels = regexprep(column_labels, '}', '');
column_labels = regexprep(column_labels, '[', '_');
column_labels = regexprep(column_labels, ']', '');
column_labels = regexprep(column_labels, '/', '');
column_labels = regexprep(column_labels, ':', '-');
column_labels = regexprep(column_labels, '\t', '');
column_labels = regexprep(column_labels, '', '_');
empty_label_ind = find(strcmp(column_labels, ''));
if isempty(empty_label_ind) == 0
    column_labels{empty_label_ind} = 'extra';
end

% store each dataset to structure
% if throws error, make sure that column labels do not contain special
% characters, etc.
 %if strcmp(class(column_labels), 'cell') == 1
 %    for i = 1:numel(column_labels)
%     test(i) = char(column_labels(i));
   
 %end

DATA = cell2struct(data_matrix, column_labels, 2);
%DATA_fields = fieldnames(DATA);
T=struct2table(DATA);
DATA=table2struct(T,'ToScalar',true);
end

