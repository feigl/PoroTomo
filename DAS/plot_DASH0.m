%% 

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

%% draw one trace
%[utime1,T_MAS] = mas_dstr2utime(tstr1,'dd-mm-yyyy HH:MM:SS.FFF');
time1_dstr = '2016-03-18 14:59:51.404';
[utime1,T_MAS] = mas_dstr2utime(time1_dstr,'yyyy-mm-dd HH:MM:SS.FFF');
nseconds = utime1-HDR.utofs0
nhalfminutes = floor(nseconds/30.)

% trace1 = ddk(1:12*60*2,5);
%trace1 = ddk(:,5);
trace1a = ddk(nhalfminutes:nhalfminutes+120*8,5);
trace1a = trace1a';

figure
plot(trace1a,'k.-');
[tofs0_dstr,tofs0_dstruct] =  mas_utime2dstr(HDR.utofs0,'yyyy-mm-dd HH:MM:SS.FFF')
xlabel(sprintf('%s %s','half-minutes after ',time1_dstr));
ylabel('strain [units??]');
title('channel 5 from Doug');


%% load metadata
% load ../dm_files/metadata
META = load('/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_files/metadata.mat','-mat')
% find labels for segments
META.itrn0(2) % index of first channel in segment 2
META.itrn0(72) % index of first channel in segment 72


%% plot original 
load sum_dash.mat
% imagesc(seconds(T_MATLAB-min(T_MATLAB)),1:nchan,log10(abs(sumdatH')));
% axis tight
% axis xy
% colormap(jet)
% colorbar
% xlabel(sprintf('%s %s','seconds after ',min(T_MATLAB)))
% ylabel('channel index');
% title(mfilename)
% printpng(sprintf('%s.png',mfilename));

figure
trace1b = sumdatH(1:120*8,5);
plot(trace1b,'r.-')
xlabel(sprintf('%s %s','half-minutes after ',T_MATLAB(1)));
ylabel('strain [units??]');
title('channel 5 from original');




