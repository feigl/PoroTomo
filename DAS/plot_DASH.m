%% plot strain rate from summed, resampled DASH data
% 20180221 Kurt Feigl


%% initialize
close all
clear all
%% initialize
nf=0; % count figures
tbegin = tic;

%% set up path for Matlab
%addpath(genpath('/globus/PoroTomo/SOFTWARE'),'-begin');
%addpath('/Users/feigl/gipht/utils','-begin'); % needed for colortables
%addpath(genpath('/Users/feigl/PoroTomo'),'-begin');
addpath(genpath('/mnt/usr1/feigl/PoroTomo'),'-begin');
%addpath('/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_scripts/','-begin');
addpath('/mnt/usr1/dmiller/dmiller_utilities','-begin');

%% load header and data file Doug's way
datadir = '/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_files/';
filenam = 'DASHstrain'
hdrfile = strcat(datadir,filesep,filenam,'.dmih');
%load -mat DASHstrain.dmih
HDR = load(hdrfile,'-mat')
[nsamp,nchan] = size(HDR.io_size);

%ddk=masio_get_data('DASHstrain');
datfile = strcat(datadir,filesep,filenam);
ddk=masio_get_data(datfile);
% imagesc(ddk(1:12*60*2,:),[-1 1])
% caxis([-1 1]*1e3)

%% select time interval
% %[utime1,T_MAS] = mas_dstr2utime(tstr1,'dd-mm-yyyy HH:MM:SS.FFF');
% %time1_dstr = '2016-03-18 14:59:51.404';
 time1_dstr = '2016-03-18 03:00:00.000';
 [utime1,T_MAS] = mas_dstr2utime(time1_dstr,'yyyy-mm-dd HH:MM:SS.FFF');
 nseconds = utime1-HDR.utofs0
 nhalfminutes = floor(nseconds/30.)
 it1 = nhalfminutes;
 it2 = it1 + 120*9;  % 9 hours
 
 %% transpose and scale into nanostrain for display
Sdot = transpose(ddk(it1:it2,:)); % 1 row is now a channel
[nchannels,nsamples] = size(Sdot)
filter='nanostrain/second'

%% draw one trace
nf=nf+1;hh(nf)=figure;hold on
ichan = 2000;
trace1a = ddk(it1:it2,ichan);
trace1a = trace1a';
plot(trace1a,'k.-');
[tofs0_dstr,tofs0_dstruct] =  mas_utime2dstr(HDR.utofs0,'yyyy-mm-dd HH:MM:SS.FFF')
xlabel(sprintf('%s %s','half-minutes after ',time1_dstr));
ylabel('strain [units??]');
title(sprintf('channel %d from Doug',ichan));
delete(sprintf('%s_fig%02d.pdf',mfilename,nf));
printpdf(sprintf('%s_fig%02d.pdf',mfilename,nf));


