%% Basic Script to read segy files
% Chelsea Lancelle
% 7/9/15

% clears all variables, closes all figures, clears the command window
clear all
close all
clc

%% Read segy files

% call to read the data from giving the filename
[Data,SegyTraceHeaders,SegyHeader]=ReadSegy('large shaker NEES_130910161319_01.sgy');

% the data will be stored in Data
% header information will be stored in SegyHeader
% individual trace headers will be stored in SegyTraceHeaders

% In the segy data, each column is a different trace.  Generally don't look
% at traces below number 189 (I think) because those channels were not
% buried in the ground.

%% Plot one DAS channel
% plots one channel (#250) as read from the segy file

% create time vector to plot against (1000 sps)
timevec = 0:0.001:30-0.001;
% plot the channel 250 time series over time (: plots all values in that
% direction)
plot(timevec,Data(:,250));
% create labels on the axes and a title
ylabel('DAS Units');
xlabel('Time (sec)');
title('Channel 250 Time Series');