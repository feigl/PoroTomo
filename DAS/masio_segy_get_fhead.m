function fhd=masio_segy_get_fhead(file,endian)

%function fhd=masio_segy_get_fhead(file,endian)
%
%  returns a struct with fields 
% fhd.texthead
% fhd.ltrace  (byte 3223 or 3593)
% fhd.ntrace  (byte 3213 or 3597)
% fhd.srate   (1e6/byte 3217; = sample rate in Hz)
% fhd.list    (all nonzero entries)
% fhd.filename 
%
% if nargin<2 endian='b'

if nargin==0 
    [file,indir] = uigetfile('*.sgy;*.segy','Select segy file');
    file=[indir,file];
end

if nargin<2 endian='b';end

fid=fopen(file,'r',endian);

GENHDR	= fread (fid,3200,'uchar');
if strcmp(char(ebcdic2ascii(GENHDR(1))),'C')
    ebchdr = char(ebcdic2ascii(GENHDR));% for rev0
else
    ebchdr = char(GENHDR);% for rev1
end

fhd.texthead=reshape(ebchdr,80,40)';

hblock4a=3201:4:3212; %rev-1 standard
hblock2=3213:2:3520;  %
hblock4b=3521:4:3600; %need 4-byte integers for large ntrace
lb4a=length(hblock4a); % that's a lower case "L'
lb2=length(hblock2);
lb4b=length(hblock4b);
hdbytes=[hblock4a hblock2 hblock4b]';

hd=0*hdbytes;
hd(1:lb4a)=fread(fid,lb4a,'uint32');
hd(lb4a+[1:lb2])=fread(fid,lb2,'uint16');
hd(lb4a+lb2+[1:lb4b])=fread(fid,lb4b,'uint32');


fhd.ltrace=hd(8);
fhd.ntrace=hd(4);
fhd.srate=1e6/hd(6);
fhd.formatcode=hd(10);

%If 4-byte ltrace and ntrace are in the unassigned block, use them:
if hd(end) fhd.ntrace=hd(end); end
if hd(end-1) fhd.ltrace=hd(end-1); end

%list of filled header slots:
ff=find(hd);
fhd.list=[hdbytes(ff) hd(ff)];

% check whether file length, ntrace and ltrace are consistent
fseek(fid,0,'eof');
fhd.flength=ftell(fid);

predicted_length=3600+(fhd.ltrace*4+240)*fhd.ntrace;
if fhd.flength~=predicted_length %test whether ntrace or ltrace can be replaced
    dsize=fhd.flength-3600; %bytes of data with headers
    % we should have dsize = (ltrace*4 + 240)*ntrace = (ltrace+60)*4*ntrace
    if rem(dsize/4,fhd.ntrace) ==0
        fhd.ltrace=dsize/4/fhd.ntrace-60;
    elseif rem(dsize/4,fhd.ltrace+6) ==0
        fhd.ntrace=dsize/4/(fhd.ltrace+60);
    else
       error(sprintf('Bad Size: file length %u does not match 3600+(240+%u)*%u',fhd.length,fhd.ltrace,fhd.ntrace));
    end
end

fclose(fid);

fhd.filename=file;
