function t0=get_das_utctime(textual_header)
%this version is show UTC time
%for silixa ... 
%requirements: SegyMat. 
% 20160325 original Zeng Xiang-Fang UW-Madison
% 20160328 adapted and documented by Kurt Feigl and Chelsea Lancelle
% 20160509 Kurt Feigl with comments from Neal Lord
% 20180220 Kurt Feigl make less verbose

%segyf
%[data,trd,hd]=ReadSegy(segyf);

%trd(1)
%hd(1)

%% convert from EBCDIC to ASCII
atext = char(ebcdic2ascii(textual_header));

%% split into lines of 80 characters each
for i=1:ceil(numel(atext)/80)
    atextrows{i} = atext([((i-1)*80+1):(i-1)*80+80]);
    %fprintf(1,'%80s\n',char(atextrows{i}));
end



%  
% C01 Client: Uni. Wisconsin
% C02 Field: PoroTomo
% C03 Fibre: Trenched Surface Fibre
% C04 Data collected by Silixa iDAS, Distributed Fibre Optic Sensor
% C05 iDAS S/N: iDAS16043
% C06 SEGY Format: Rev.1, IEEE 32bit float, big endian
% C07 Field Recording Filename: BNL_IDAS__160318023921.tdms
% C08 Continuous acquisition data converted to SEGY
% C09 UTC Timestamp of first sample: 18-Mar-2016 02:39:21.404309911
% C10
% C11
% C12
% C13 Receiver positions are in true E, N, Elevation (m)
% C14 Number of Traces: 8721
% C15 Samples Per Trace: 30000
% C16 Sampling Interval (us): 1000
% C17 Record Length (sec): 30
% C18 Measurement Units: Depths = Metres, Coordinates = Metres
% C19
% C20 Trace amplitude is proportional to fibre strain-rate
% C21 For comparison with conventional geophones it is recommended to
% C22 time-integrate this data
% C23
% C24 Trace Header Byte Positions:
% C25 41-44: Receiver Elevation (m)
% C26 81-84: Receiver Easting (m)
% C27 85-88: Receiver Northing (m)
% C28 233-236: Samples per Trace
% C29 237-240: Fibre distance from beginning of trench (m)
% C30
% C31
% C32
% C33
% C34
% C35 Binary Header Byte Positions:
% C36 63-66: Samples per Trace
% C37 67-70: Number of Traces in file
% C38
% C39 Silixa Ltd, 230 Centennial Park, Elstree, UK, WD6 3SN
% C40 End Text Header
%% look at header 

% I looked at a DASH segy header and the Silixa time of the first sample is
% there to the nS. Below is the whole textual header. The portion
% pertaining to time is
%  
% C09 UTC Timestamp of first sample: 18-Mar-2016 02:39:21.404309911
%  
% I think you need to look at your code to make sure it is properly
% extracting the whole time time string. You are correct in that they seem
% to have added yet another version of the time in the header. I?ve had to
% deal with these formats from Silixa:
%  
% GPSTimeStamp: 10/09/2013 23:14:07.752 (UTC)
%  
% GPS Timestamp of first sample: 11/09/2013 23:08:49.888
%  
% GPS Time stamp: 2016/03/21 07:37:21.403999 (UTC)  çnew for PoroTomo for QC data in the field
%  
% UTC Timestamp of first sample: 18-Mar-2016 02:39:21.404309911     çnew for PoroTomo after field work
%  
% Fortunately they haven?t made a versions incompatible with the other versions.

atextrow09 = char(atextrows{9});

if numel(strfind(atextrow09,'GPSTimeStamp:')) > 0
    itype = 1;
elseif numel(strfind(atextrow09,'GPS Timestamp of first sample:')) > 0
    itype = 2;
elseif numel(strfind(atextrow09,'GPS Time stamp:')) > 0
    itype = 3;
elseif numel(strfind(atextrow09,'UTC Timestamp of first sample: ')) > 0
    itype = 4;
    ii = strfind(atextrow09,': ');
    timestamp = atextrow09(ii+2:ii+33);
    day = str2num(timestamp(1:2));
    monthstr = timestamp(4:6);
    switch monthstr
        case 'Jan'
            month = 1;
        case 'Feb'
            month = 2;
        case 'Mar'
            month = 3;
        case 'Apr'
            month = 4;
        case 'May'
            month = 5;
        case 'Jun'
            month = 6;
        case 'Jul'
            month = 7;
        case 'Aug'
            month = 8;
        case 'Sep'
            month = 9;
        case 'Oct'
            month = 10;
        case 'Nov'
            month = 11;
        case 'Dec'
            month = 12;
        otherwise
            error(sprintf('Unrecognized month %s',monthstr));
            month = 0;
    end
            
    year         = str2num(timestamp(8:11)); 
    hour         = str2num(timestamp(13:14)); 
    minute       = str2num(timestamp(16:17));  
    seconds      = str2num(timestamp(19:20));
    milliseconds = str2num(timestamp(21:32))*1.e3;
else
    itype = 0;
    timestamp = 'Unknown time stamp format'
end

%     %% older version used time stamp in trace header.
%     year = trd(1).YearDataRecorded;
%     doy = trd(1).DayOfYear;
%     day1 = datetime(year,1,1) + caldays(doy-1);
%     month = day1.Month;
%     %day2 = day(day1,'dayofmonth')
%     day = day1.Day;
%     h0=trd(1).HourOfDay;
%     m0=trd(1).MinuteOfHour;
%     s0=trd(1).SecondOfMinute;
% end

t0 = datetime(year,month,day,hour,minute,seconds,milliseconds...
    ,'TimeZone','UTC');
t0.Format='yyyy/MM/dd_HH:mm:ss.SSSSSSSSSS';


% % add time increment
% dt2=trd(1).dt * 1e-6; %unit is microsecond
% s1=(trd(1).ns-1)*dt2;

%% code time as class duration - 20160328 incorrect.
%times=duration(h0,m0,s0:dt2:s0+s1);

%% code time as absolute epoch (point in time) using Matlab day number
%times =datenum(year,month,day,h0,m0,s0:dt2:s0+s1);

%% code time as absolute epoch (point in time
% Matlab function datetime DOES support fractional seconds
% ss = s0:dt2:(s0+s1);
% 
% times = datetime(year*ones(size(ss)) ...
%     ,month*ones(size(ss)) ...
%     ,day*ones(size(ss))   ...
%     ,h0*ones(size(ss))    ...
%     ,m0*ones(size(ss))    ...
%     ,ss);
% times.Format='yyyy/MM/dd_hh:mm:ss.SSSSSSS';
% times.TimeZone = 'UTC';


return
