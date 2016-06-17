function [year,month,day] = yymmdd2yr_mo_dy(yymmdd)
%function [year,month,day] = yymmdd2yr_mo_dy(yymmdd)
% given date in 6-digit integer format, return 4-digit year, month, day
% example:
%
% yymmdd = 160311
% [year,month,day] = yymmdd2yr_mo_dy(yymmdd)
%
% year =
%     2016
% month =
%      3
% day =
%     11
%
% 2016/03/24 by Kurt Feigl

year  = nan(size(yymmdd));
month = nan(size(yymmdd));
day   = nan(size(yymmdd));

for i=1:numel(yymmdd)
    yy = floor(yymmdd(i)/10000);
    mm = floor((yymmdd(i) - 10000*yy)/100);
    dd = yymmdd(i) - 10000*yy - 100*mm;
    
    date = datevec(now);
    
    if yy > 80 && yy < 100 && abs(mod(yy,1)) < eps
        year = yy + 1900;
    elseif yy < 100 && yy <= date(1)-2000 && abs(mod(yy,1)) < eps
        year(i) = yy + 2000;
    else
        warning(sprintf('Suspicious year %f\n',yy));
        %year(i) = yy;
    end
    
    if mm >= 1 && mm <= 12 && (mm - abs(mod(mm,12))) < eps
        month(i) = mm;
    else
        warning(sprintf('Suspicious month %f\n',mm));
        %month(i) = floor(mm);
    end
    
    if dd >= 1 && dd <= 31 && (dd - abs(mod(dd,31))) < eps
        day(i) = dd;
    else
        warning(sprintf('Suspicious day %f\n',dd));
        %day(i) = floor(dd);
    end
end

return
