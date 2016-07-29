function S = csv2struct2(csvname)
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

T=readtable(csvname,'ReadVariableNames',true,'HeaderLines',0,'TreatAsEmpty','');
S=table2struct(T,'ToScalar',true);
return


