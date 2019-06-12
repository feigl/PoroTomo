function [thd,fhd]=masio_segy_get_thead(file,endian)

% function [thd,fhd]=masio_segy_get_thead(file,endian)
% 
% if nargin<2 endian='b';end

if nargin<2 endian='b';end

hdbytes=...
    [1:4:28 ...
     29:2:36 ...
     37:4:68 ...
     69:2:72 ...
     73:4:88 ...
     89:2:180 ...
     181:4:200 ...
     201:2:204 ...
     205:4:208 ...
     209:2:240 ...
     ]';

dbytes=diff([hdbytes; 241]);
blist=[hdbytes [dbytes]];
bstart=hdbytes([0; find(diff([dbytes]))]+1);
b2=find(dbytes==2);
b4=find(dbytes==4);

if nargin==0 
    [file,indir] = uigetfile('*.sgy;*.segy','Select segy file');
    file=[indir,file];
end
fhd=masio_segy_get_fhead(file,endian);

file=fhd.filename;
ntrace=fhd.ntrace;
ltrace=fhd.ltrace;

trbuff=zeros(ntrace,95);
lskip2=(ltrace*4+240-2);
lskip4=(ltrace*4+240-4);

fid=fopen(file,'r',endian);

for kk=1:95
fseek(fid,3600+hdbytes(kk)-1,'bof');
if dbytes(kk)==2
    trbuff(:,kk)=fread(fid,ntrace,'int16',lskip2);
else
    trbuff(:,kk)=fread(fid,ntrace,'int32',lskip4);
end
end

fclose(fid);

tmx=max(trbuff);
tmn=min(trbuff);
fg=find(abs(tmx)+abs(tmn)>0);
f1=find(tmx(fg)==tmn(fg));
f2=find(tmx(fg)>tmn(fg));

thd.ucodes=hdbytes(fg(f1));
thd.uvals=tmx(fg(f1));
thd.bcodes=hdbytes(fg(f2));
thd.bvals=trbuff(:,fg(f2));

for kk=1:length(thd.ucodes)
    eval(['thd.b' num2str(thd.ucodes(kk)) '= thd.uvals(' num2str(kk) ');']);
end
for kk=1:length(thd.bcodes)
    eval(['thd.b' num2str(thd.bcodes(kk)) '= thd.bvals(:,' num2str(kk) ');']);
end

    


