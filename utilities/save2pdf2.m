%SAVE2PDF Saves a figure as a properly cropped pdf
%
%   save2pdf(pdfFileName,handle,dpi)
%
%   - pdfFileName: Destination to write the pdf to.
%   - handle:  (optional) Handle of the figure to write to a pdf.  If
%              omitted, the current figure is used.  Note that handles
%              are typically the figure number.
%   - dpi: (optional) Integer value of dots per inch (DPI).  Sets
%          resolution of output pdf.  Note that 150 dpi is the Matlab
%          default and this function's default, but 600 dpi is typical for
%          production-quality.
%
%   Saves figure as a pdf with margins cropped to match the figure size.

%   (c) Gabe Hoffmann, gabe.hoffmann@gmail.com
%   Written 8/30/2007
%   Revised 9/22/2007
%   Revised 1/14/2007
%   Revised 2018/05/27 Kurt Feigl

function save2pdf2(pdfFileName,handle,dpi)

% Verify correct number of arguments
%error(nargchk(0,3,nargin));
narginchk(0,3);

% If no handle is provided, use the current figure as default
if nargin<1
    [fileName,pathName] = uiputfile('*.pdf','Save to PDF file:');
    if fileName == 0; return; end
    pdfFileName = [pathName,fileName];
end
if nargin<2 || ishandle(handle)==0  
    handle = gcf;
end
if nargin<3
    dpi = 150;
end

% Backup previous settings
prePaperType = get(handle,'PaperType');
prePaperUnits = get(handle,'PaperUnits');
preUnits = get(handle,'Units');
prePaperPosition = get(handle,'PaperPosition');
prePaperSize = get(handle,'PaperSize');

% Make changing paper type possible
set(handle,'PaperType','<custom>');

% Set units to all be the same
set(handle,'PaperUnits','inches');
set(handle,'Units','inches');

% Set the page size and position to match the figure's dimensions
PaperPosition = get(handle,'PaperPosition')

% get the size of the white space
Position = get(handle,'Position')
OuterPosition = get(handle,'OuterPosition')

t1=sprintf('%s',strrep(pdfFileName,'_','\_'));
t2 = datestr(now,31); %31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17
tu=getenv('USER');
t3=sprintf('%s %s %s',t1,t2,tu);
t4 = strrep(t3,'\','\\');
t5 = strrep(t4,'_','\_');

dirname = sprintf('%s',strrep(pwd,'_','\_'));

%% write text string horizontally at lower left
subplot('Position',[0. 0. OuterPosition(3)-1 0.02],'Units','Inches');
text(1.,0.,t5 ...
            ,'Units','inches'...
            ,'VerticalAlignment','Bottom'...
            ,'HorizontalAlignment','Left'...
            ,'Clipping','off'...
            ,'FontName','Courier','FontSize',7 ...
            ,'Rotation',0);
axis off


%% print the name of the working directory vertically on the left hand side
subplot('Position',[0. 0. 0.04 OuterPosition(4)],'Units','Inches');
text(0.,1.,dirname ...
            ,'Units','inches'...
            ,'VerticalAlignment','Bottom'...
            ,'HorizontalAlignment','Left'...
            ,'Clipping','off'...
            ,'FontName','Courier','FontSize',7 ...
            ,'Rotation',90);
axis off

set(handle,'PaperPosition',[0,0,OuterPosition(3),OuterPosition(4)]);
set(handle,'PaperSize',[OuterPosition(3),OuterPosition(4)]);

% Save the pdf (this is the same method used by "saveas")
print(handle,'-dpdf',pdfFileName,sprintf('-r%d',dpi))

% Restore the previous settings
set(handle,'PaperType',prePaperType);
set(handle,'PaperUnits',prePaperUnits);
set(handle,'Units',preUnits);
set(handle,'PaperPosition',prePaperPosition);
set(handle,'PaperSize',prePaperSize);
return
end
