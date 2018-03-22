%% sum DASH from raw data
% roots
% porotomo.geology.wisc.edu:/mnt/t31/PoroTomo/ANALYSIS/feigl/DASH/20161112/DougExample20161113.m
% 
% 20160903 Doug Miller and Kurt Feigl
% 20160907 Kurt Feigl find use correct file name (21, not 51 seconds)
% 20161113 Kurt Feigl test paths 
% 201802128 Kurt Feigl

%% initialize
close all
clear all
nf=0; % count figures
tbegin = tic;
HOST = getenv('HOST');

%% set up paths
switch HOST
    case 'askja.ssec.wisc.edu'
        addpath(genpath('/globus/PoroTomo/SOFTWARE'),'-begin');
        datadir = '/s11/PoroTomo/DATA/DASH';
    case 'porotomo.geology.wisc.edu'
        addpath(genpath('/mnt/t31/PoroTomo/SOFTWARE'),'-begin');
        %addpath('/Users/feigl/gipht/utils','-begin'); % needed for colortables
        %addpath(genpath('/Users/feigl/PoroTomo'),'-begin');
        % addpath(genpath('/mnt/usr1/feigl/PoroTomo'),'-begin');
        %addpath('/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_scripts/','-begin');
        % addpath('/mnt/usr1/dmiller/dmiller_utilities','-begin');       
        datadir = '/mnt/t31/PoroTomo/DATA/DASH';
    otherwise
        warning(sprintf('Unknown HOST %s\n',HOST));
        addpath(genpath('/mnt/t31/PoroTomo/SOFTWARE'),'-begin');
        datadir = '/mnt/t31/PoroTomo/DATA/DASH';
end



%% start time of time series
% in Matlab datetime format
%tref_mdt=datetime(2016,3,21,07,37,00); % Hawthorne earthquake
tref_mdt=datetime(2016,3,18,0,0,0.0); % March 18th
tref_mdt.Format = 'yyyy/MM/dd hh:mm:ss.SSSSSSS';
tref_mdt.TimeZone = 'UTC';


%% set up to read files
fprintf(1,'At %.0f seconds after starting: Reading SEGY file...\n',toc(tbegin));
% %rsync -av gridftp1.ssec.wisc.edu:/globus/PoroTomo2/DASH/20160321/PoroTomo_iDAS16043_160321073721.sgy .
% %fn = '/data/PoroTomo/DATA/DASH/20160321/PoroTomo_iDAS16043_160321073721.sgy' % name of SEG-Y file
% fn = '/globus/PoroTomo2/DATA/DASH/20160321/PoroTomo_iDAS16043_160321073721.sgy' % name of SEG-Y file
%datadir = '/mnt/t31/PoroTomo/DATA/DASH';
yyyymmdd = '20160318'


%flist = ls(strcat(datadir,filesep,'*.sgy'));
%DIR = dir(strcat(datadir,filesep,'*.sgy'));
DIR = dir(strcat(datadir,filesep,yyyymmdd,filesep,'*.sgy'));
flist = {DIR.name};
flist = flist';
%%flist = flist(1:2)  % truncate list for testing
[nfiles,ndum] = size(flist)

%% read the header of the first file. 
kk=1; 
sgy_filename=strcat(datadir,filesep,yyyymmdd,filesep,char(flist{kk}))
HSEGY = GetSegyHeader(sgy_filename);
% get start time of time series
time_tag_mdt1 = get_das_utctime(HSEGY.TextualFileHeader);
time_tag_mdt1.Format = 'yyyy/MM/dd hh:mm:ss.SSSSSSS';
time_tag_mdt1.TimeZone = 'UTC';
time_tag_mdt1
% number of samples
nsamples=HSEGY.ns;
% number of channels
nchannels=HSEGY.DataTracePerEnsemble;
% sample interval in seconds
dt_das = HSEGY.dt*1e-6 


%% declare arrays 
strain_rate_summed_over30s_in_radians_per_second=zeros(nfiles,nchannels);
sample_standard_deviation_in_radians_per_second=zeros(nfiles,nchannels);
% time tag in unix time (seconds since epoch)
time_tag_uts=zeros(nfiles,1); % unix time of first sample in seconds
% time tag in Matlab datetime format 
%  NaT Not-a-Datetime.
%     NaT is the representation for Not-a-Datetime, a value that can be stored in
%     a datetime array to indicate an unknown or missing datetime value. DATETIME
%     creates a NaT automatically when reading text that cannot be parsed as a datetime,
%     or for elements in a datetime array where the Year, Month, 
%     Day, Hour, Minute, or Second properties are set to NaN.
%  
%     D = NaT with no inputs returns a scalar NaT datetime.
%     D = NaT(N) is an N-by-N matrix of NaTs.
%     D = NaT(M,N) or NaT([M,N]) is an M-by-N matrix of NaTs.
%     D = NaT(M,N,P,...) or NaT([M,N,P,...]) is an M-by-N-by-P-by-... array of NaTs.
time_tag_mdt = NaT(nfiles,1); % 
time_tag_mdt.TimeZone='UTC';
%time_tag_mdt.Format='yyyy/MM/dd hh:mm:ss.SSSSSSS';
time_tag_mdt.Format='yyyy/MM/dd HH:mm:ss.SSSSSSS'; % upper case H for 24-hour clock


tcalc0=tic;
fprintf(1,'At %.0f seconds after starting: Summing...\n',toc(tbegin));
for kk=1:nfiles
    sgy_filename=strcat(datadir,filesep,yyyymmdd,filesep,char(flist{kk}));
 
    % read header
    HSEGY1 = GetSegyHeader(sgy_filename);
    % get start time of time series
    time_tag_mdt1 = get_das_utctime(HSEGY1.TextualFileHeader);
   %time_tag_mdt1.Format = 'yyyy/MM/dd hh:mm:ss.SSSSSSS';
    time_tag_mdt1.Format = 'yyyy/MM/dd HH:mm:ss.SSSSSSS'; % 20180228 HH to displays hours out of 24, e.g., 3:00 PM as 13:00
    time_tag_mdt1.TimeZone = 'UTC';
    %ut1
    % number of samples
    nsamples1=HSEGY1.ns;
    % number of channels
    nchannels1=HSEGY1.DataTracePerEnsemble;
    % sample interval in seconds
    dt_das1 = HSEGY1.dt*1e-6;  %  value in SEGY header is in microseconds
    % time tag of middle of time interval spanned by file in Matlab datetime format
    time_tag1 = time_tag_mdt1 + seconds((nsamples-1)*dt_das/2.);
    time_tag_mdt(kk) = time_tag1;
    
    % load header and data file Doug's way
    [dd,fhdH,thdH]=masio_segy_get_data(sgy_filename);
    tstring = fhdH.texthead(9,[1:30]+35);
    
    % time tag of middle of time interval spanned by file in unix time
    time_tag_uts(kk)=mas_dstr2utime(tstring,'dd-mmm-yyyy HH:MM:SS.FFF')+(nsamples-1)*dt_das/2.;
    
    % difference between two time ways of tagging time
    time_diff_in_seconds = second(time_tag_mdt(1)) - mod(time_tag_uts(1),60.);

    %% check for sanity
    if nsamples1 == nsamples ...
            && nchannels1 == nchannels ...
            && abs(dt_das1-dt_das) <= eps ...
            && abs(time_diff_in_seconds) < dt_das1
         
        %% read data file 
        %[dd,fhdH,thdH]=masio_segy_get_data(fname); % get file header and time header
        [dd,fhdH]=masio_segy_get_data(sgy_filename); % get data only      

        %% sum over 30 seconds to calculate strain rate
        % Feigl and PoroTomo team [2017] Stanford Geothermal workshop says:
        % To understand the DAS recordings, we review a few of their
        % characteristics [Bakku, 2015; Daley et al., 2015]. DAS measures the
        % strain rate ?? in the fiber by averaging its elongation over a segment
        % of cable (called the ?gauge length?) during a temporal sampling interval.
        % The elongation represents the phase shift of the backscattered light
        % pulse. These data were written with dimensions of radians per
        % millisecond in the SEG-Y files. In this DAS system, one radian of phase
        % change corresponds to 116 nanometers of elongation. The wavelength of the
        % laser light is 1550 nanometer. The temporal sampling interval was set to
        % 1 millisecond and the spatial sampling length was set to 1 meter. The
        % spatial resolution of the DAS strain rate measurement equals the gauge
        % length of 10 m.
        % Assume that dd from SEG-Y files are in units of radians/millisecond
        srate1 = 1.0e3*sum(dd)/(nsamples*dt_das);                        % Units are radian/second.
        strain_rate_summed_over30s_in_radians_per_second(kk,:) = srate1; % Units are radian/second.
        std1=1.03*std(dd);                                               % units are radian/second.
        sample_standard_deviation_in_radians_per_second(kk,:) = std1;    % Units are radian/second.
        if mod(kk,1)==0
            fprintf(1,'file %5d of %5d %s done in %12.1f s rms(srate1) = %12.6E rms(std1) = %12.6E \n'...
                ,kk,nfiles,sgy_filename,toc(tbegin),rms(srate1),rms(std1));
        end

    else
        sgy_filename
        strain_rate_summed_over30s_in_radians_per_second(kk,:) = nan(size(dd)); 
        kk
        nsamp1
        nsamp
        nchan1
        nchan
        dt
        dt1
        time_diff_in_seconds
        warning('miscount 1');
    end
end

%% Save result using Matlab file format
outmat_filename = sprintf('DASHstrainrate%s.mat',yyyymmdd);
save(outmat_filename ...
    ,'flist'...
    ,'time_tag_mdt'...
    ,'time_tag_uts' ...
    , 'strain_rate_summed_over30s_in_radians_per_second'...
    , 'sample_standard_deviation_in_radians_per_second');

fprintf(1,'Saved strain rate to %s\n',outmat_filename);
telapsed=seconds(toc(tbegin));
telapsed.Format = 'hh:mm:ss';
fprintf(1,'Total Elapsed Time [hh:mm:ss] = %s\n',char(telapsed));
fprintf(1,'Average time per file %.3f seconds\n',toc(tcalc0)/nfiles);

%% Save result using Doug's file format
%% TODO Doug, please help us with your protocol
% dmid_filename = sprintf('DASHstrainrate%s',yyyymmdd);
% masio_appendto_header(dmid_filename,'flist', strcat(datadir, filesep, flist));
% masio_appendto_header(dmid_filename,'time_tag_uts',time_tag_uts);
% %masio_appendto_header(dmid_filename,'units',repmat({'radian/s'}, numel(time_tag_mdt), 1));
% masio_appendto_header(dmid_filename,'units','s');
% masio_appendto_header(dmid_filename,'units','radian/s');
% masio_appendto_header(dmid_filename,'units','radian/s');
% 
% %masio_write_data(dmid_filename,strain_rate_in_radians_per_second);
% masio_write_data(dmid_filename,[time_tag_uts ...
%     , strain_rate_in_radians_per_second...
%     , sample_standard_deviation_in_radians_per_second]);
% fprintf(1,'Saved strain rate to %s\n',dmid_filename);

%% Evaluate performance
telapsed=seconds(toc(tbegin));
[telapsed_HH, telapsed_mm, telapsed_ss] = hms(telapsed);
%telapsed.Format = 'HH:mm:ss'
fprintf(1,'Total Elapsed Time [HH:mm:ss] = %2i:%2i:%2.2f\n', telapsed_HH, telapsed_mm, telapsed_ss);
%fprintf(1,'Total Elapsed Time [seconds] = %s\n',char(telapsed));
fprintf(1,'Average time per file %.3f seconds\n',toc(tcalc0)/nfiles);




% %% plot a trace
% fprintf(1,'At %.0f seconds after starting: Plotting...\n',toc(tcalc0));
% nf=nf+1;HH(nf)=figure;
% plot(time_tag,strain_rate_in_radians_per_second(:,100)/1e-9,'k.-');
% xlabel(sprintf('seconds after %s',char(tref_mdt)),'Interpreter','none');
% xlabel('UTC time')
% %xtickformat('yyyy/MM/dd HH:mm:ss')
% xtickformat('HH:mm:ss')
% ylabel('strain rate [nanostrain/second]');
% title(sprintf('Channel 100 %s',dmid_filename),'Interpreter','none');
% printpdf(sprintf('%s_fig%02d.pdf',msgy_filename,nf));
% 

%% Verify
STR=load(outmat_filename)



