clear all
close all
histogram(randn([1000,1]))

handle = gcf;
pdfFileName='test.pdf'
res=1200;
% Set the page size and position to match the figure's dimensions
PaperPosition = get(handle,'PaperPosition');

% get the size of the white space
Position = get(handle,'Position');
OuterPosition = get(handle,'OuterPosition');

% Mx = max([1,OuterPosition(3)-Position(3)]);
% My = max([1,OuterPosition(4)-Position(4)]);
Mx = max([Position(3),PaperPosition(3)-Position(3)]);
My = max([1,PaperPosition(4)-Position(4)]);
% Mx = PaperPosition(3)-Position(3)
% My = PaperPosition(4)-Position(4)
%MarginPosition = [OuterPosition(1)-Mx, OuterPosition(2)-My, Mx, My]
%MarginPosition = [OuterPosition(1)-Mx, OuterPosition(2)-My, Mx, My]
MarginPosition = [PaperPosition(1), PaperPosition(2)+PaperPosition(4)-My, Mx, My]
%MarginPosition = [4, 4, Mx, My]
%set(handle,'NextPlot','add');
%subplot('position',MarginPosition,'Units','inches');
% text(MarginPosition(1),MarginPosition(2),pdfFileName ...
%             ,'Units','inches'...
%             ,'VerticalAlignment','Bottom'...
%             ,'HorizontalAlignment','Left'...
%             ,'Clipping','off'...
%             ,'FontName','Courier','FontSize',9 ...
%             ,'Rotation',0);

%subplot('Position',[0 0 8 0.05],'Units','inches');
subplot('Position',MarginPosition,'Units','inches');
axis off
text(0,0,pdfFileName...
    ,'Units','inches'...
    ,'VerticalAlignment','Bottom'...
    ,'HorizontalAlignment','Left'...
    ,'Clipping','off'...
    ,'FontName','Courier','FontSize',9 ...
    ,'Rotation',0);
%print(handle,pdfFileName,'-dpdf',sprintf('-r%04d',res));


%save2pdf('test.pdf',gcf,1200);