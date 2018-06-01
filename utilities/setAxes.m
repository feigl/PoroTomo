function setAxes(varargin)
%SETAXES automatically sets the position and size of the axes (subplots) in
%the given figure so that white spaces are removed
%
%Optional Inputs
%   pad: A vector of size 1,2,3,4. The padding method follows that of HTML
%       If size 1, then each subplot is padded with white space by the 
%       given amount
%       If size 2, then for each subplot, the top and bottom are padded by 
%       pad[0], and the right and left are padded by pad[1]
%       If size 3, then for each subplot, the top is padded by pad[0], the
%       right and left are padded by pad[1], and the bottom is padded by 
%       pad[2]
%       If size 4, then for each subplot, the top is padded by pad[0], the 
%       right by pad[1], the bottom by pad[2], the left by pad[3]
%   fh: figure handle. By default, == gcf
%
%Examples
%   setAxes
%       set axes positions for current figure (gcf)
%   setAxes(0.1)
%       each subplot in current figure has boundary white spaces of 0.1 (in
%       normalized units)
%   setAxes([0.05,0,0,0])
%       each subplot in current figure has top boundary (white space) of 
%       0.05 and 0 in the right, bottom and left boundary
%   setAxes([0.05,0,0,0],fh)
%       set axes for subplots in the figure with figure handle fh
%
%Custom Functions
%   cell2var
%   figsize
%
%   last updated 09/22/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%
%% optional inputs
input = {[0.05,0,0,0], gcf};
input_N = length(input);
for jj = 1:input_N
    if nargin >= jj && ~isempty(varargin{jj})
        input{jj} = varargin{jj};
    end
end
[pad, fh] = cell2var(input);

%% parameter preparation
switch length(pad)
    case 1
        pad = [pad, pad, pad, pad]; % top right bottom left
    case 2
        pad = [pad(1),pad(2),pad(1),pad(2)]; % top right bottom left
    case 3
        pad = [pad(1),pad(2),pad(3),pad(2)]; % top right bottom left
    case 4
        % do nothing
    otherwise
        error('padding size is invalid');
end

%% main
num = figsize(fh);
ax = findobj(fh,'type','axes');
axN = length(ax);

axX = 1/num(2)-pad(2)-pad(4); % X direction
axY = 1/num(1)-pad(1)-pad(3); % Y direction


for ii = 1:num(1)
    posY = pad(3)+1-ii/num(1);
    for jj = 1:num(2)
        posX = pad(4)+(jj-1)/num(2);
        
        kk = (ii-1)*num(2)+jj;
        set(ax(axN+1-kk),'units','normalized',...
                'position',[posX,posY,axX,axY],...
                'XTick','','YTick','');
    end
end
end

% ----------------------- custom functions --------------------------------
%% figsize
function num = figsize(varargin)
%FIGSIZE returns the number of figures in the given figure handle (usually
%gcf)
%   num: 1-by-2 vector specifying the number of rows and columns of
%   subplots in the given figure
%
%Optional Inputs
%   fh: figure handle. By default, == gcf
%
%   last updated 09/22/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%
%% optional input
fh = gcf;
if nargin >= 1 && ~isempty(varargin{1})
    fh = varargin{1};
end

%%
ax = findobj(fh,'type','axes');
pos = get(ax,'position');
if iscell(pos)
    pos = cell2mat(pos);
    colNum = numel(unique(pos(:,1))); % same X positions
    rowNum = numel(unique(pos(:,2))); % same Y positions
else
    colNum = 1;
    rowNum = 1;
end

num = [rowNum, colNum];

end

%% cell2var
function [varargout] = cell2var(cellA, varargin)
%CELL2VAR stands for cell to variable. It converts a cell of length N to N
%specified variables
%   cellA: 1D cell
%
%Optional Inputs
%   num: a # specifying which cell to convert to variable
%
%   last updated 06/05/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%

%% possible case
if nargin >= 2 && ~isempty(varargin{1})
    varargout{1} = cellA{varargin{1}};
    return;
end
%% main
if nargout <= length(cellA)
   for kk = 1:nargout
       varargout{kk} = cellA{kk};
   end
else
    error('number of outputs exceeds length of input cell');
end

end
