function [utime,dstruc] = mas_dstr2utime(dstr,fmt)

% function [utime,dstruc] = mas_dstr2utime(dstr,fmt)
% secperday=24*60*60;
% 
% if nargin<2
%     utime=mas_chop((datenum(dstr) - datenum('1-jan-1970'))*secperday,4);
% else
%     utime=mas_chop((datenum(dstr,fmt) - datenum('1-jan-1970'))*secperday,4);
% end
% if nargout>1
%     dstruc.utime=utime;
%     dstruc.year=str2num(mas_utime2dstr(utime,'yyyy'));
%     dstruc.month=str2num(mas_utime2dstr(utime,'mm'));
%     dstruc.day=str2num(mas_utime2dstr(utime,'dd'));
%     ut0=mas_dstr2utime(num2str(dstruc.year),'yyyy');
%     dstruc.dayofyear=floor((utime-ut0)/24/60/60)+1;
%     dstruc.datestr=mas_utime2dstr(utime);
% end

secperday=24*60*60;

if nargin<2
    utime=mas_chop((datenum(dstr) - datenum('1-jan-1970'))*secperday,4);
else
    utime=mas_chop((datenum(dstr,fmt) - datenum('1-jan-1970'))*secperday,4);
end

if nargout>1
    dstruc.utime=utime;
    dstruc.year=str2num(mas_utime2dstr(utime,'yyyy'));
    dstruc.month=str2num(mas_utime2dstr(utime,'mm'));
    dstruc.day=str2num(mas_utime2dstr(utime,'dd'));
    ut0=mas_dstr2utime(num2str(dstruc.year),'yyyy');
    dstruc.dayofyear=floor((utime-ut0)/24/60/60)+1;
    dstruc.datestr=mas_utime2dstr(utime);
end
