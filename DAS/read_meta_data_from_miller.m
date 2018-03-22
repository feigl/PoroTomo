%% load metadata
% load ../dm_files/metadata
META = load('/mnt/t31/PoroTomo/ANALYSIS/dmiller/dm_files/metadata.mat','-mat');
% find labels for segments
META.itrn0(2) % index of first channel in segment 2
META.itrn0(72) % index of first channel in segment 72
nsegments = numel(META.itrn0);


%% try to build look up table of channel numbers
%              jjt: [1×2195 double]
%              jjl: [1×2211 double]
%              jjv: [1×2218 double]
% 
nchanmin = min([numel(META.jjl) numel(META.jjt) numel(META.jjt)])
fprintf(1,'ii, META.jjl(ii), META.jjt(ii), META.jjl(ii)\n');
for ii=1:nchanmin
    fprintf(1,'%5d %5d %5d %5d\n',ii,META.jjl(ii),META.jjt(ii),META.jjl(ii));
end
%  ii, META.jjl(ii), META.jjt(ii), META.jjl(ii)
%    ii   jjl   jjt   jjl
%     1     4     7     4
%     2     5     8     5
%     3     6     9     6
%     4    13    16    13
%     5    14    17    14
%     6    15    18    15
%     7    22    25    22
%     8    23    26    23
%     9    24    27    24
%    10    33    28    33
%    11    34    36    34
% ....
%  2187  6556  6607  6556
%  2188  6563  6611  6563
%  2189  6564  6612  6564
%  2190  6565  6613  6565
%  2191  6572  6614  6572
%  2192  6573  6615  6573
%  2193  6574  6628  6574
%  2194  6581  6629  6581
%  2195  6582  6630  6582