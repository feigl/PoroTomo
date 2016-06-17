%Basic Script for Brady Day 1 30-sec segy file
% clears all variables, closes all figures, clears the command window
clear all
close all
clc

%% Read segy files

[Data,SegyTraceHeaders,SegyHeader]=ReadSegy('BNL_IDAS__LINE_0_SID_0_FID_000030.sgy');
% the data will be stored in Data
% header information will be stored in SegyHeader
% individual trace headers will be stored in SegyTraceHeaders

% create time vector to plot against (1000 sps)
% Find length of Segy trace
nrows=size(Data,1);
timevec = 0:0.001:(nrows-1)*0.001;

% plot results

%% Plot one DAS channel
% plots one channel ch
ch=5000
DAS(:,ch)=11300*Data(:,ch);

figure(1)
plot(timevec,DAS(:,ch));
ylabel('DAS Strain Rate');
xlabel('Time (sec)');
title(['BNL_IDAS__LINE_0_SID_0_FID_000030 DAS Channel = ', num2str(ch)])

[ps,w] = spectrum(DAS(:,ch));
figure(2)
plot(w,ps);
ylabel('DAS Power Spectrum Amplitude');
xlabel('Frequency (Hz)');
title(['BNL_IDAS__LINE_0_SID_0_FID_000030 DAS Channel = ', num2str(ch)])

%{
%% Zoom into 10-second window between t1 and t2 seconds
t1=15;
t2=25;
figure(2)
plot(timevec(t1*1000:t2*1000),DAS(t1*1000:t2*1000,ch));
ylabel('DAS Strain Rate');
xlabel('Time (sec)');
title(['BNL_IDAS__LINE_0_SID_0_FID_000030 DAS Channel = ', num2str(ch)])
grid on
%}
