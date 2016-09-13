%% Get current version of metadata files from askja using anynomous ftp ftp://roftp.ssec.wisc.edu
% Elena C. Reinisch
% 20160913

% Initialize
clear all; close all;

% Add path to PoroTomo software repository on askja
addpath(genpath('/data/PoroTomo/Software/PoroTomo'))

%% R

% set name of text file which contains name of xlsx files in column format 
% e.g., 
%      Borehole_DAS_DTS_UTM_coordinates.csv
%      Vibroseis_timingLog_PoroTomo_Mar2016_Phases_1234_v3.xlsx
% etc. 
% this can be made from command line by typing "ls *.xlsx > xlsx_list.txt"
% in Metadata subfolder
datfile = 'xlsx_list.txt'

% loop through addresses and save all data to matlab structure METADATA as
% well as .mat file
[ METADATA ] = xlsx2struct(datfile)



