function out = mas_chop(in,precision)  
                                       
% function out = mas_chop(in,precision)
%                                      
if nargin<2 precision=15; end;         
ee=10^precision;ei=1./ee;              
out=ei*round(ee*in);                   
