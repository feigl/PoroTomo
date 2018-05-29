% LABELFIG
% label a figure with metadata
%   Revised 2018/05/27 Kurt Feigl

function labelfig(hlabel,vlabel)

% Verify correct number of arguments
%error(nargchk(0,3,nargin));
narginchk(0,2);

%% set up horizontal label
if nargin < 1
   hlabel=sprintf('%s %s',datestr(now,31),getenv('USER'));
end
hlabel=sprintf('%s',strrep(strrep(hlabel,'\','\\'),'_','\_'));
   
%% set up horizontal label
if nargin < 2
    vlabel = pwd;
end
vlabel=sprintf('%s',strrep(strrep(vlabel,'\','\\'),'_','\_'));

%% write text string horizontally at lower left
subplot('Position',[0. 0. 5 0.02],'Units','Inches');
text(1.,0.,hlabel ...
            ,'Units','inches'...
            ,'VerticalAlignment','Bottom'...
            ,'HorizontalAlignment','Left'...
            ,'Clipping','off'...
            ,'FontName','Courier','FontSize',7 ...
            ,'Rotation',0);
axis off


%% print the label vertically on the left hand side
subplot('Position',[0. 0. 0.04 5],'Units','Inches');
text(0.,1.,vlabel ...
            ,'Units','inches'...
            ,'VerticalAlignment','Bottom'...
            ,'HorizontalAlignment','Left'...
            ,'Clipping','off'...
            ,'FontName','Courier','FontSize',7 ...
            ,'Rotation',90);
axis off

return
end
