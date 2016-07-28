function [DATA, DATA_fields] = csv2struct(csvname)
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
% DATA_fields - cell array containing the names of the fields 


% get header information for csvimport
vals = importdata(csvname); % pull text into cells by row (1 row per cell)

% check to see if data has been imported as data structure
if strcmp(class(vals), 'struct') == 1
    text = vals.textdata;
    column_labels = text(1,:);
else %read as text file
    column_labels = strsplit(char(vals(1)), ','); % separate header row into strings by commas
% ncols_csv = numel(column_labels); % count number of columns in csv file
% nrows_csv_wo_header = numel(vals)-1; % count number of rows in csv file
end

% read data from file 
[data_matrix] = csvimport(csvname, 'columns', column_labels);

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

% store each dataset to structure
% if throws error, make sure that column labels do not contain special
% characters, etc.
DATA = cell2struct(data_matrix, column_labels, 2);
DATA_fields = fieldnames(DATA);
end

