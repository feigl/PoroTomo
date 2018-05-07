%% initialize
close all
clear all
%% initialize
nf=0; % count figures
tbegin = tic;


%% set up path for Matlab
addpath(genpath('/mnt/t31/PoroTomo/SOFTWARE'),'-begin');
%addpath(genpath('/globus/PoroTomo/SOFTWARE'),'-begin');
%addpath('/Users/feigl/gipht/utils','-begin'); % needed for colortables
%addpath(genpath('/Users/feigl/PoroTomo'),'-begin');
%addpath(genpath('/mnt/usr1/feigl/PoroTomo'),'-begin');
%addpath('/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_scripts/','-begin');
%addpath('/mnt/usr1/dmiller/dmiller_utilities','-begin');

%% choose source of data
%isource = 1 % Kurt's integration from raw data
isource = 2 % Doug's integration of Xiangfang's resampled data -- fails
switch isource
    case 1
        
        %% file names
        datadir = '/mnt/t31/PoroTomo/DATA/DASH';
        yyyymmdd = '20160318'
        %flist = ls(strcat(datadir,filesep,'*.sgy'));
        %DIR = dir(strcat(datadir,filesep,'*.sgy'));
        DIR = dir(strcat(datadir,filesep,yyyymmdd,filesep,'*.sgy'));
        flist = {DIR.name};
        flist = flist';
        flist = flist(1:120)  % truncate list for testing
        [nfiles,ndum] = size(flist)
        
        %% select time intervals
        it1 =  3*120; % omit first 3 hours
        it2 = 12*120; % omit last 12 hours

        
        %% load result
        matfilename = sprintf('DASHstrainrate%s.mat',yyyymmdd)
        S = load(matfilename);
        
        % transpose and scale into nanostrain for display
        Sdot = transpose(S.strainrate30s)/1e-9; % 1 row is now a channel
        [nchannels,nsamples] = size(Sdot)
        
        % add processing step to documentation
        filter='nanostrain/second from raw data'
        unitlabel = 'nanonstrain/second';
    case 2
        %% load header and data file Doug's way
        
        datadir = '/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_files/';
        filenam = 'DASHstrain'
        matfilename = filenam;
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
        jt1 = nhalfminutes;
        jt2 = jt1 + 120*9;  % 9 hours
        
        dhour = 1/120; % halfminutes per hour
        S.time_tag = datetime(2016,03,18,3,0,0) + hours([0:dhour:9]);
        
        yyyymmdd = '20160318'
        
        %% transpose and scale into nanostrain for display
        Sdot = transpose(ddk(jt1:jt2,:)); % 1 row is now a channel
        %Sdot = Sdot*11.6/30/100/10; % 11.6 radians per guage length / 30 seconds / 100 samples per second / 10 m
        Sdot = Sdot/30; % strain over 30 seconds
        [nchannels,nsamples] = size(Sdot)
        filter='unknown units from Miller sum of resampled DASH'
        
        % set up for later time slice
        it1 = 1;
        it2 = jt2-jt1+1;
        
        unitlabel = 'unknown units';
    otherwise
        isource
        error('unknown isource')
end

%% number figures according to source
nf = isource*100;

%% load metadata
% load ../dm_files/metadata
META = load('/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_files/metadata.mat','-mat');
% find labels for segments
META.itrn0(2) % index of first channel in segment 2
META.itrn0(72) % index of first channel in segment 72
nsegments = numel(META.itrn0);


% %% remove mean over each row (time tag), returning nsample values
% Sdot = Sdot-mean(Sdot,1);
% filter = strcat(filter,':','mean1')
%
% %% remove mean over each column (channel), returning nchannel values
% Sdot = Sdot-mean(Sdot,2);
% filter = strcat(filter,':','mean2')

% %% filter by segment
% segmean1=nan(1,nsamples);
% SEG.mean     = nan(nsegments,nsamples);
% SEG.time_tag = S.time_tag;
% for i=1:nsegments
%     % take first channel in segment
%     i1=META.itrn0(i);
%     % take last channel in segment
%     if i==nsegments
%         i2=nchannels;
%     else
%         i2=META.itrn0(i+1);
%     end
%     % indices of segment
%     SEG.i1(i)=i1;
%     SEG.i2(i)=i2;
%     % go 5 channels inside corner
%     ii1 = i1+5;
%     ii2 = i2-5;
%     SEG.ii1(i)=ii1;
%     SEG.ii2(i)=ii2;
%     % take mean over only inside channels
%     segmean1=mean(Sdot(ii1:ii2,:),1);
%     SEG.mean(i,:) = segmean1;
%     % replace all channels with mean of inside channels
%     segblk = repmat(segmean1,i2-i1+1,1);
%     Sdot(i1:i2,:) = segblk;
% end
% filter = strcat(filter,':','segmean')
% matfilename2 = sprintf('DASHstrainrate%sSEG.mat',yyyymmdd);
% fprintf(1,'Saving %s\n',matfilename2);
% save(matfilename2,'SEG');


%% choose a channel in segment 50
% META.itrn0(50);ans =       5644
% META.itrn0(51);ans =       5800
ichannel = 5700

% anomolous channels exposed to air
% segment 42, channel 4700:
% META.itrn0(42); ans =         4673
% META.itrn0(43); ans =         4731
% segment 47, channel 5200:
% META.itrn0(47); ans =         5166
% META.itrn0(48); ans =         5214


%% histogram
nf=nf+1;hh(nf)=figure;
histogram(colvec(Sdot(:,it1:it2)));
set(gca,'YScale','log')
xlabel(sprintf('Strain rate  %s',unitlabel));
ylabel('N');
title(sprintf('%s\n%s',filter,matfilename),'Interpreter','none');
delete(sprintf('%s_fig%02d.pdf',mfilename,nf));
printpdf(sprintf('%s_fig%02d.pdf',mfilename,nf));



%% plot a trace
nf=nf+1;hh(nf)=figure;
plot(S.time_tag(it1:it2),Sdot(ichannel,it1:it2),'k.-');
%ylim([-20 0]); % clip
xlabel('UTC time')
%xtickformat('yyyy/MM/dd HH:mm:ss')
xtickformat('HH:mm:ss')
ylabel(sprintf('Strain rate  %s',unitlabel));
title(sprintf('%s\n%s',filter,matfilename),'Interpreter','none');
legend(sprintf('Channel %d',ichannel));
delete(sprintf('%s_fig%02d.pdf',mfilename,nf));
printpdf(sprintf('%s_fig%02d.pdf',mfilename,nf));

%% plot several traces
nf=nf+1;hh(nf)=figure;hold on
%ifavorites = [1000,2000,3000,4000,5000,6000,7000];
ifavorites = [4700,5200,5644+5:20:5800-5];
for i=1:numel(ifavorites)
    plot(S.time_tag(it1:it2),Sdot(ifavorites(i),it1:it2));
    legendlabel{i} = sprintf('channel %3d',ifavorites(i));
end
%ylim([-20 20]); % clip


% %% plot several segments
% %filter = strcat(filter,':','diff')
% nf=nf+1;hh(nf)=figure;hold on
% %ifavorites = [53:70];  % SE side?
% ifavorites = [51:60];  % not noisy
% for i=1:numel(ifavorites)
%     %undifferenced
%     plot(SEG.time_tag(it1:it2),SEG.mean(ifavorites(i),it1:it2));  
%     % difference with respect to first
%     %plot(SEG.time_tag(it1:it2),SEG.mean(ifavorites(i),it1:it2)-SEG.mean(ifavorites(1),it1:it2));
%     legendlabel{i} = sprintf('segment %3d',ifavorites(i));
% end
%ylim([-20 20]); % clip

xlabel('UTC time')
%xtickformat('yyyy/MM/dd HH:mm:ss')
%xtickformat('hh:mm:ss')
ylabel(sprintf('Strain rate  %s',unitlabel));
title(sprintf('%s\n%s',filter,matfilename),'Interpreter','none');
legend(legendlabel,'Location','bestoutside');
delete(sprintf('%s_fig%02d.pdf',mfilename,nf));
printpdf(sprintf('%s_fig%02d.pdf',mfilename,nf));

%% take log of strain rate
Sdotlog10 = log10(abs(Sdot));
filter = strcat(filter,':','log10(abs(Sdot)')


%% plot the whole array
nf=nf+1;hh(nf)=figure;
%subplot(2,1,1);
%imagesc(seconds(time_tag-min(time_tag)),1:nchannels,log10(abs(sumdatH')));

%imagesc(1:nsamples,1:nchannels,Sdot); 
imagesc(1:nsamples,1:nchannels,Sdotlog10); 
hold on;
xlabel('sample index');
ylabel('channel index');

% label with segment index on right side
% draw white line to expand margin
xp1 = (1.0                    )*nsamples; % stagger line
xp2 = (1.0 +            8*0.02)*nsamples; % stagger line
plot([xp1 xp2],max(META.itrn0)*[1 1],'w.-');
for i=1:numel(META.itrn0)
    xp1 = (1.0                     )*nsamples; % stagger line
    xp2 = (1.12 - (mod(i,4)+1)*0.02)*nsamples; % stagger line
    xt  = (1.12 - (mod(i,4)+1)*0.02)*nsamples; % stagger text position
    plot([xp1 xp2],META.itrn0(i)*[1 1],'k-');
    
    text(xt,META.itrn0(i),sprintf('%3d',i) ...
        ,'FontSize',9 ...
        ,'Margin',1 ...
        ,'HorizontalAlignment','left','VerticalAlignment','bottom');
    %    ,'EdgeColor','w','BackgroundColor','w'...

end

axis tight
axis xy
colormap(jet)
%caxis([-1,1])

%title(sprintf('log10(abs(strainrate/(1/s)) %s',matfilename),'Interpreter','none');
title(sprintf('%s\n%s',filter,matfilename),'Interpreter','none');

% label time on top
text(1,nchannels+1000,'UTC time','EdgeColor','w','BackgroundColor','w'...
    ,'HorizontalAlignment','left','VerticalAlignment','top')
for i=6*120:6*120:nsamples-1
     plot([i,i],[nchannels,nchannels+1000],'k-');

     ttick = S.time_tag(i);
    ttick.Format = 'HH:mm:ss';
    text(i,nchannels+1000,sprintf('%s',char(ttick)),'EdgeColor','w','BackgroundColor','w'...
    ,'HorizontalAlignment','center','VerticalAlignment','top');
end
colorbar;
delete(sprintf('%s_fig%02d.pdf',mfilename,nf));
printpdf(sprintf('%s_fig%02d.pdf',mfilename,nf));
%printpng(sprintf('%s_fig%02d.png',mfilename,nf));

