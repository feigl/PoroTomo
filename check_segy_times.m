% check times of SEG-Y files from Silixa
% 20160505 Kurt Feigl, with code from Xiangfang Zeng and Chelsea Lancelle

%dirname = '20160318';
%dirname = '20160319';
%dirname = '.'
%D = dir(strcat(dirname,'/','*.sgy'))
D = dir('*.sgy');

%% example
% regular file
%fn = '../Earthquake/BNL_IDAS_20160321_073621.sgy'
%fn = '20160319/PoroTomo_iDAS16043_160319000021.sgy'
% file that is too short
%fn = '/Volumes/PoroTomo2//PoroTomo/Task6_Deployment/Subtask6_2_Deploy_and_Operate_DAS_DTS/Silixa/Trenched_Fibre/TestSweeps/T85/SEGY/2016_03_16/BNL_IDAS_test sweeps_160310213743.sgy'

for ifile = 1:numel(D)
    if D(ifile).isdir == 0
        
        %fn = strcat(dirname,'/',D(ifile).name);
        fn = D(ifile).name;
        fprintf(1,'%s :',fn);
        
        H = GetSegyHeader(fn);
        t0 = get_das_utctime0(H.TextualFileHeader)
        %t0=times(1);   t0.Format = 'yyyy/MM/dd_hh:mm:ss.SSSSSSS';t0.TimeZone = 'UTC';
        
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


        t1=times(end); t1.Format = 'yyyy/MM/dd_hh:mm:ss.SSSSSSS';t1.TimeZone = 'UTC';
        n_epochs = numel(times);
        time_duration = seconds(t1 - t0);
        
        %% check that file is named properly
        % find index of period in file name
        idot = strfind(fn,'.sgy');
        idot = idot(1);
        
        nameisok = 0;
        if numel(idot) == 1
            
            % make structure to parse file name
            FN.name = fn;
            FN.ss   = str2num(fn(idot-2:idot-1));
            FN.mm   = str2num(fn(idot-4:idot-3));
            FN.hh   = str2num(fn(idot-6:idot-5));
            FN.DD   = str2num(fn(idot-8:idot-7));
            FN.MM   = str2num(fn(idot-10:idot-9));
            FN.YYYY = str2num(fn(idot-12:idot-11)) + 2000;
            
            if         int16(FN.ss) == int16(second(t0)) ...
                    && int16(FN.mm) == int16(minute(t0)) ...
                    && int16(FN.hh) == int16(hour(t0)) ...
                    && int16(FN.DD) == int16(day(t0)) ...
                    && int16(FN.MM) == int16(month(t0)) ...
                    && int16(FN.YYYY) == int16(year(t0))
                nameisok = 1;
            else
                nameisok = 0;
            end
        else           
            fprintf(1,'%s %s %12.6f seconds %5d epochs %s\n' ...
                , char(t0) ...
                , char(t1) ...
                , time_duration ...
                , n_epochs ...
                , fn ...
                );          
        end        
    end
    if nameisok == 1
        fprintf(1,' OK\n');
    else
         fprintf(1,' PROBLEM\n');
    end
end

