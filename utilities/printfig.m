function printeps(figfilename,ghandle)
% write current graphics window to a EPS file

if nargin < 1
    figfilename = mfilename;
end

if numel(strfind(figfilename,'.fig')) <= 0
    figfilename = sprintf('%s.fig',figfilename);
end

%   savefig(H,FILENAME,'compact') saves the figures identified by the graphics 
%    handle array H to a MATLAB figure file called FILENAME. This MATLAB figure 
%    file can be opened only in R2014b or later version of MATLAB. Using the 
%    'compact' option reduces the size of the MATLAB figure file and the 
%    time required to create the file.


if exist('ghandle','var') == 1
    savefig(ghandle,figfilename,'compact'); %
else
    savefig(gcf,figfilename,'compact');
end

return


