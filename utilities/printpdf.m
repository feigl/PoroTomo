function printpdf(filename,res)
% write current graphics window to a PostScript file

% if strcmp(getenv('NOPRINTpdf'),'TRUE') == 1
%     disp('Env variable NOPRINTpdf is set to TRUE');
%     return
% end

if nargin < 1
    filename = mfilename;
end
if nargin < 2
    res = 1200;
end
% print current figure to pdf file name with 1200 DPI and TIFF
t0=pwd;
t1=strrep(filename,'_','\_');
%t2=date;
t2 = datestr(now,31); %31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
tu=getenv('USER');
t3=sprintf('%s %s %s %s',t1,t0,t2,tu);
t4 = strrep(t3,'\','\\');
t5 = strrep(t4,'_','\_');
% Label the figure
% Does not work on Hengill
mycomputer = computer;
if strcmp(mycomputer, 'GLNXA64')==0
    subplot('Position',[0 0 10 0.05],'Units','Centimeters');
    axis off
    text(1,0,t5...
        ,'Units','Centimeters'...
        ,'VerticalAlignment','Bottom'...
        ,'HorizontalAlignment','Left'...
        ,'Clipping','off'...
        ,'FontName','Courier','FontSize',9 ...
        ,'Rotation',0);
end

h=gcf;
old_unit = get(h,'Units');
set(h,'Units','inches','PaperPositionMode','auto','PaperUnits','inches')
set(h,'PaperSize',get(h,'position')*[0 0;0 0;1 0;0 1])
print(h,filename,'-dpdf',sprintf('-r%04d',res));
fprintf('%s created\n',filename)
set(h,'Units',old_unit)
return
end



