function [data,fhd,thd]=masio_segy_get_data(file,trace_indices,sample_indices,endian)

%function [data,fhd,thd]=masio_segy_get_data(file,trace_indices,sample_indices,endian)
%
% if nargout>1 fhd = masio_segy_get_fhead(file);
%


if nargin==0 
    [file,indir] = uigetfile('*.sgy;*.segy','Select segy file');
    file=[indir,file];
end

if nargin<4 endian='b';end

fhd=masio_segy_get_fhead(file,endian);

file=fhd.filename;
ntrace=fhd.ntrace;
ltrace=fhd.ltrace;

ltbuff=ltrace*4+240;

if (nargin<3 || isempty(sample_indices) ) sample_indices=[1:ltrace];end
if (nargin<2 || isempty(trace_indices)) trace_indices=[1:ntrace];end


fid=fopen(file,'r',endian);


ltrout=length(sample_indices);
ntrout=length(trace_indices);
data=zeros(ltrout,ntrout);
trace=zeros(ltrout,1);
%efficiently read a simple block:
if  (ltrout==(sample_indices(end)-sample_indices(1)+1))
    for itr=1:ntrout
    dstart=3600+240+ltbuff*(trace_indices(itr)-1)+ 4*(sample_indices(1)-1);   
    fseek(fid,dstart,'bof'); 
    if fhd.formatcode==5
 [trace,count]=fread(fid,ltrout,'real*4');
    elseif fhd.formatcode==2
   [trace,count]=fread(fid,ltrout,'int32');    
    elseif fhd.formatcode==1
   [utrace,count]=fread(fid,ltrout,'uint');
   trace=ibm2ieee(utrace);
    else
        ['format code ' num2str(fhd.formatcode) ' not supported']
        break
    end
        
    try
     data(:,itr)=trace;
    catch
    
    [num2str(itr) '  requested: '  num2str(ltrout) ';  read:  ' num2str(count)]
    break
    end 
    end

else
    for itr=1:ntrout
    for jsamp=1:ltrout
    dstart=3600+240+ltbuff*(trace_indices(itr)-1)+ 4*(sample_indices(jsamp)-1);
    fseek(fid,dstart,'bof');
    data(jsamp,itr)=fread(fid,1,'real*4');
end
end
end

fclose(fid);

if nargout>2 thd=masio_segy_get_thead(file,endian);end



